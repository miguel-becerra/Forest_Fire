/*
 * main.cpp
 *
 *  Created on: 20-06-2015
 *      Author: samsung
 */

#include <sqlite3.h>
#include <string>

#include "sitios.h"
#include "radiation.h"

using namespace std;

int sitio::contador = 0;

int Nx=300;
int Ny=300;
int mass=0;
int FPRINT=0; //si se imprime el fuel de cada sitio
int DIMS=0;
int xi0,yj0; // semillas
double dzf=5.0;
double zk0; // semillas
double BD=0.0;
double dxmalla,slope;
double dT=0.5; //paso temporal
double PP=0.1;
double l,l2;
double Tmean=20.0; // vida media
double walpha=0.00001; // coef de degradacion
double Cth=0.9; // degradation/burning threshold
double dCth=-0.0005; // degradation/burning threshold of the derivative
double front,quem,activ;
ofstream out1;
unsigned long limit=numeric_limits<unsigned long>::max();
double PI=acos(-1.0);
double dpi=2.0*acos(-1.0);

int nPP=100;
double Pini=0.0, Pfin=1.0;

int ENDB=0; // end the calculations when reaching the border
int REAB=0; // flag that indicates when the border is reached

// physical distances! in [m]
double pdx=3.3;
double pdy=3.3;
double pdz=1.0;

//vegetation types
double vch=3.0;
double vcd=1.5;
double vfh=4.2;
double vfd=1.5;
//             c h  c d  f h  f d
vegetation veg(vch, vcd, vfh, vfd, Tmean, walpha, Cth);


int test_radiation(){
  double dx=0.0;
  double dy=0.0;
  double dz=0.0;

  ofstream aux;
  aux.open("aux.dat");
  //init_vfactor();

  double dist=0.0;
  for(dx=-10;dx<=10;dx+=0.1){
    dist=sqrt(dx*dx+dy*dy);
    if(dist>0){
      for(dz=-20.0; dz<=20.0;dz+=0.5){
	aux<<dx<<" "
	   <<dz<<" "
	  //hay truco en radiation.cpp para poner la llama mas arriba?
	   <<veg.radiation(dx,dy,dz,0.0,0.0)<<" "
	   <<veg.radiation(dx,dy,dz,3.0,0.0)<<" "
	   <<veg.radiation(dx,dy,dz,6.0,0.0)<<" "
	  //<<veg.FF(dist,dz)<<" "
	   <<endl;
      }
      aux<<endl;
    }
  }

  aux.close();
  return 0;
}

bool isSeed(int ii, int jj, vector<double>& sX, vector<double>& sY){
  int sI,sJ;
  for(int i=0; i<sX.size(); i++){
    sI=round(Nx*sX[i]);
    sJ=round(Ny*sY[i]);
    if(ii==sI && jj==sJ){
      return true;
    }
  }
  return false;
}
bool isTime(int contador, vector<double>& times){
  double secs=contador*dT;
  for(int i=0; i<times.size(); i++){
    if( abs(secs-times[i])<dT ){
      return true;
    }
  }
  return false;
}

string parseOption(sqlite3 * db, string env_var){
  char * valor;
  string sql_exp_base = "insert into input (key, param) values ('%s', '%s')";
  int len = (int)sql_exp_base.size() + 20;
  char * sql_exp = new char[len];
  valor = getenv (env_var.c_str());
  if (valor==NULL){
    cerr<<"error "<<" "<<env_var<<endl;
    return "";
  }
  printf(sql_exp_base.c_str(), env_var.c_str(), valor);
  cout<<endl;
  sprintf(sql_exp, sql_exp_base.c_str(), env_var.c_str(), valor);
  sqlite3_exec(db, sql_exp,0, 0, 0);
  return valor;
}

int main(){

  //Hacemos las base de datos para el output
  sqlite3 *db;
  sqlite3_stmt *res;
  char *zErrMsg = 0;
  int rc;
  rc = sqlite3_open("job.db", &db);
  if( rc ){
    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
    sqlite3_close(db);
    exit(1);
  }
  sqlite3_exec(db,"PRAGMA synchronous=OFF",0,0,0);

  //Benchmarks
  clock_t start,finish;
  start = clock();

  //resto
  double x0,y0;
  int kw; //en que k escribir los archivos
  int kmax; //numero de ejecuciones

  //recogemos las opciones desde el entorno
  string landfile;
  string vegfile;
  string windfile;
  bool boolveg=false;
  bool boolwind=false;

  out1.open("output.dat");
  rc = sqlite3_exec(db,
		    "create table input (\
                     key     text,\
                     param   text);",
		    0, 0, 0);

  l=atof(parseOption(db,"l").c_str());
  l2=l*l;

  x0=atof(parseOption(db,"x0").c_str());
  y0=atof(parseOption(db,"y0").c_str());

  //en lo siguiente hay ke tener cuidado con el punto 0.0
  //strtod devuelve 0.0 cuando no hay nada mas en la string
  vector <double> seedsX,seedsY;
  string seedaux;
  double value;
  char * pEnd;
  seedaux=parseOption(db,"sx");
  value=strtod(seedaux.c_str(),&pEnd);
  while(value!=0.0){
    seedsX.push_back(value);
    value=strtod(pEnd,&pEnd);
  }
  seedaux=parseOption(db,"sy");
  value=strtod(seedaux.c_str(),&pEnd);
  while(value!=0.0){
    seedsY.push_back(value);
    value=strtod(pEnd,&pEnd);
  }
  if(seedsX.size()!=seedsY.size()){
    cerr<<"Problema fatal con las semillas: sx,sy"<<endl;
    return 5;
  }

  //pTimes: vector con los tiempos a los ke se escribe el estado
  // del boske. en segundos!
  vector <double> pTimes;
  seedaux=parseOption(db,"pTimes");
  value=strtod(seedaux.c_str(),&pEnd);
  while(value!=0.0){
    pTimes.push_back(value);
    value=strtod(pEnd,&pEnd);
  }

  kw=atof(parseOption(db,"kw").c_str());
  kmax=atof(parseOption(db,"kmax").c_str());

  landfile=parseOption(db,"landf");

  vegfile =parseOption(db,"vegf");
  if (vegfile==""){
    cerr<<"error vegf"<<endl;
    boolveg=false;
  }else{
    boolveg=true;
  }

  windfile=parseOption(db,"windf");
  if (windfile==""){
    cerr<<"error windf"<<endl;
    boolwind=false;
  }else{
    boolwind=true;
  }
  //
  Pini=atof(parseOption(db,"Pini").c_str());
  Pfin=atof(parseOption(db,"Pfin").c_str());
  nPP =atof(parseOption(db,"nPP" ).c_str());
  DIMS=atoi(parseOption(db,"DIMS").c_str());
  ENDB=atof(parseOption(db,"ENDB").c_str());
  //
  vch =atof(parseOption(db,"vCh").c_str());
  vcd =atof(parseOption(db,"vCd").c_str());
  vfh =atof(parseOption(db,"vFh").c_str());
  vfd =atof(parseOption(db,"vFd").c_str());

  pdx =atof(parseOption(db,"pdx").c_str());
  pdy=pdx;
  //
  Tmean =atof(parseOption(db,"Tmean" ).c_str());
  walpha=atof(parseOption(db,"alpha" ).c_str());
  Cth   =atof(parseOption(db,"Cth"   ).c_str());
  dCth  =atof(parseOption(db,"dCth"  ).c_str());
  FPRINT=atof(parseOption(db,"FPRINT").c_str());

  cerr<<l<<endl;
  cerr<<x0<<endl;
  cerr<<y0<<endl;
  cerr<<kw<<endl;
  cerr<<landfile<<endl;
  cerr<<vegfile<<endl;


  //-------------------------------------------------------------------

  // mas adelante, tratar de cambiar esto por una base de datos.
  // hay que cambiarlo en landf también
  double alpha,lambda;
  cerr<<"-Leemos el landscape de "<<landfile<<endl;
  ifstream bland(landfile.c_str(), ios::binary);
  bland.read((char *)(&Nx), sizeof(Nx));
  bland.read((char *)(&Ny), sizeof(Ny));
  bland.read((char *)(&dxmalla), sizeof(dxmalla));
  bland.read((char *)(&slope), sizeof(slope));
  bland.read((char *)(&lambda), sizeof(lambda));
  bland.read((char *)(&alpha), sizeof(alpha));
  bland.read((char *)(&BD), sizeof(BD));

  int NN=Nx*Ny;
  double malla[Nx][Ny];//para guardar la altura
  double x,y,uk,vk,lk,lak;

  bland.read((char *)(&malla), sizeof(malla));
  bland.close();
  cerr<<"   ( leida una matriz de "<<Nx<<"x"<<Ny<<" )"<<endl;
  cerr<<"   ( lambda,alpha,dim= "
      <<lambda<<" , "
      <<alpha<<" "
      <<BD<<" )"
      <<endl;
  imprime(BD,"dim-superficie");
  //-------------------------------------------------------------------
  int wNx=0;
  int wNy=0;
  double wdx;
  if(boolwind){
    cerr<<"  -Leemos el viento de "<<windfile<<endl;
    ifstream bwind(windfile.c_str(), ios::binary);
    bwind.read((char *)(&wNx), sizeof(wNx));
    bwind.read((char *)(&wNy), sizeof(wNy));
    bwind.read((char *)(&wdx), sizeof(wdx));
    double wxmalla[wNx][wNy];
    double wymalla[wNx][wNy];
    bwind.read((char *)(&wxmalla), sizeof(wxmalla));
    bwind.read((char *)(&wymalla), sizeof(wymalla));
    bwind.close();
    cerr<<"   (viento: leida dos matrices de "<<wNx<<"x"<<wNy<<" )"<<endl;
    if(wNx!=Nx || wNy!=Ny || wdx!=dxmalla){
      cerr<<"Dimenssion mismatch!!"<<endl;
      cerr<<"problem with the wind file"<<endl;
      return 2;
    }
  }
  //else{
    double wxmalla[wNx][wNy];
    double wymalla[wNx][wNy];
    wxmalla[0][0]=0.0;
    wymalla[0][0]=0.0;
    //}
  //-------------------------------------------------------------------

  //guardamos los valores leidos en la vegetacion
  veg.set(vch, vcd, vfh, vfd, Tmean, walpha, Cth);
  veg.Tmean=Tmean;
  veg.alpha=walpha;
  veg.Cth=Cth;
  // Initialization of the flux factors
  veg.init_vfactor();
  test_radiation();

  int seed=time(NULL);
  srand(seed);
  //VSLStreamStatePtr stream;
  double *maux;
  maux=new double[NN];

  double dPP=(Pfin-Pini)/(nPP-1);

  int contador=0;
  vector <sitio> forest;
  vector <double> errordFF,errordQ,errordA;
  vector <double> massFF,massQ,massA,rhoA,rhoFF;
  vector <double>dens      (nPP,0.0);
  vector <double>time      (nPP,0.0);
  vector <double>pPP       (nPP,0.0);
  vector <double>massQT    (nPP,0.0);
  vector <double>rhoQ      (nPP,0.0);
  vector <double>meanDMFF  (nPP,0.0);
  vector <double>meanDMQ   (nPP,0.0);
  vector <double>meanDMA   (nPP,0.0);
  vector <double>sigmaDMFF (nPP,0.0);
  vector <double>sigmaDMQ  (nPP,0.0);
  vector <double>sigmaDMA  (nPP,0.0);
  vector <double>Pci       (kmax,0.0); // umbral de percolacion para cada realizacion
  vector <vector<double> > Ps;  // survival probability
  vector <vector<double> > Chi; // mean cluster size
  vector <vector<double> > R2;  // mean cluster radius
  vector <vector<double> > R2z; // mean cluster radius with Z
  vector <vector<double> > Tdin;// holds the survival runs

  int **diccionario;//para cambiar de posicion a orden
  diccionario = new int*[Nx];
  for(int i=0; i<Nx; i++){
    diccionario[i] = new int[Ny];
  }
  for(int i=0; i<Nx; i++){
    for(int j=0; j<Ny; j++){
      diccionario[i][j]=-1;
    }
  }

  vector <double> paux;
  for(int iPP=0; iPP<nPP; ++iPP){
    pPP[iPP]=Pini+iPP*dPP;
    Ps.push_back(paux);
    Chi.push_back(paux);
    R2.push_back(paux);
    R2z.push_back(paux);
    Tdin.push_back(paux);
  }

  set <int> ffo; // fire front order
  set <int>::iterator iffo,iffoe; // fire front order iterator
  vector <sitio>::iterator tree1,tree2,forestend;

  cerr<<endl;

  double delta=0.0;
  int tot=0;
  double **mallaMap;
  mallaMap=new double*[Nx];
  for(int i=0; i<Nx; ++i) mallaMap[i]=new double[Ny];

  //calculo del promedio de diferencia de altura entre sitios adyacentes
  for(int i=1; i<Nx;++i){
    for(int j=1; j<Ny;++j){
      mallaMap[i][j]=malla[i][j];
      delta+=fabs(malla[i][j]-malla[i-1][j]);
      delta+=fabs(malla[i][j]-malla[i][j-1]);
      tot+=2;
    }
  }
  cerr<<"<dh> ~ "<<delta/tot<<endl;
  //generamos un mapa vectorial
  vectorMap(mallaMap);

  //===========================================
  //             Ciclo principal
  //===========================================
  seed+=895;
  //    vslNewStream( &stream, VSL_BRNG_MRG32K3A, seed );

  //dejamos listos los char* para escribir en la base de datos
  rc = sqlite3_exec(db,
		    "create table Map (\
                           x      real,		   \
                           y      real,		   \
                           z      real)",
		    0, 0, 0);
  rc = sqlite3_exec(db,
		    "create table FireFront (\
                           PP     real,		   \
                           t      real,		   \
                           x      real,		   \
                           y      real,		   \
                           z      real,		   \
                           state  real)",
		    0, 0, 0);
  rc = sqlite3_exec(db,
		    "create table ForestState (\
                           PP     real,		     \
                           t      real,		     \
                           x      real,		     \
                           y      real,		     \
                           z      real,		     \
                           state  real)",
		    0, 0, 0);

  int printMap=1;
  char sql_exp_base_map[] = "insert into Map (x,y,z) values (%f, %f, %f)";
  char *sql_exp_map = new char[(int)(string(sql_exp_base_map).length())*3];

  char sql_exp_base[] = "insert into %s (PP,t,x,y,z,state) values (%f, %f, %f, %f, %f, %f)";
  char *sql_exp = new char[(int)(string(sql_exp_base).length())*3];

  int Nxfile=Nx; //despues modificamos el numero de puntos totales,
  int Nyfile=Ny; // pero es importante saber cuantos puntos se leyeron del archivo
  for(int kn=0; kn<kmax; ++kn){
    cerr<<"\r";
    cerr<<"iteracion: "
        <<kn+1<<"/"
        <<kmax<<"         "
        <<endl;

    //se genera un arreglo de NN numeros aleatorios entre 0 y 1
    //vdRngUniform(VSL_METHOD_DUNIFORM_STD_ACCURATE,stream,NN,maux,0.0,1.0);
    mrandom(NN,maux);

    for(int iPP=0; iPP<nPP; ++iPP){
      PP=Pini+iPP*dPP;
      cerr<<"\r  densidad "<<iPP+1<<" de "<<nPP<< " ( "<<PP<<" )";
      cerr.flush();


	//Dos tipos de vegetacion:
	// -random, dependiendo de P
	// -leida de un archivo.
	//  Si se lee de un archivo, ademas se especifican los "tipos"
	//  de vegetacion, a si ke se cambian los valores de los miembros de
	//  la clase sitio.
	contador=0;
	forest.clear();
	for(int i=0; i<Nx; i++){
	  for(int j=0; j<Ny; j++){
	    diccionario[i][j]=-1;
	  }
	}
	if(boolveg){
	  //Generamos las mallas usando la distribucion del archivo de vegetacion
	  //Como Yann tiene el archivo separado en cuadrados de 50metros, tenemos
	  //ke interpolar para llenar esos cuadrados.
	  //Adicionalmente, las posiciones las voy a llenar con probabilidad
	  //correspondiente a PP. Aunke parece ke ellos no usan esa probabilidad
	  //creo ke va a kedar bonito, asi hay una distribucion menos homogenea.
	  //Seguramente tendre ke forzar ke PP sea un numero alto de todas formas.
	  cerr<<"  -Leemos la vegetacion de "<<vegfile<<endl;
	  ifstream bveg(vegfile.c_str(), ios::binary);
	  int vNx,vNy;
	  double vdx;
	  bveg.read((char *)(&vNx), sizeof(vNx));
	  bveg.read((char *)(&vNy), sizeof(vNy));
	  bveg.read((char *)(&vdx), sizeof(vdx));
	  double vmalla[vNx][vNy];//guardamos la vegetacion
	  bveg.read((char *)(&vmalla), sizeof(vmalla));
	  bveg.close();
	  cerr<<"   (vegetacion: leida una matriz de "<<vNx<<"x"<<vNy<<" )"<<endl;
	  if(vNx!=Nxfile || vNy!=Nyfile || vdx!=dxmalla){
	    cerr<<"Dimenssion mismatch!!"<<endl;
	    cerr<<"problem with the vegetation file"<<endl;
	    return 2;
	  }
	  //tenemos ke borrar diccionario y reacondicionarlo a sus nuevas dimensiones.
	  //lo mas seguro es reescribir los valores de Nx y Ny, para no romper nada mas
	  //abajo. Con lo unico ke habria ke tener cuidado despues es de no llamar a
	  //malla[][] con los nuevo valores, pero deberia ser facil, creo ke no uso
	  //mucho las matrices, si no ke solo el vector forest.
	  for(int i=0; i<Nx; i++){
	    delete[] diccionario[i];
	  }
	  delete[] diccionario;

	  double xpt,ypt;
	  xpt=(vNx-1)*vdx; //largos totales
	  ypt=(vNy-1)*vdx;
	  Nx=floor(xpt/pdx) +1;//nueva cantidad de puntos
	  Ny=floor(ypt/pdy) +1;//...
	  diccionario = new int*[Nx];
	  for(int i=0; i<Nx; i++){
	    diccionario[i] = new int[Ny];
	  }
	  for(int i=0; i<Nx; i++){
	    for(int j=0; j<Ny; j++){
	      diccionario[i][j]=-1;
	    }
	  }

	  contador=0;
	  //srand(); esta inicializado mas arriba
	  double urmax=1.0/((double)RAND_MAX);
	  int ip,jp;
	  double distBx,distBy,indiceX, indiceY;
	  for(int i=0; i<Nx; i++){
	    distBx=modf(i*pdx/vdx,&indiceX);//encontramos el borde de la izkierda
	    ip=(int)indiceX;
	    for(int j=0; j<Ny; j++){
	      distBy=modf(j*pdy/vdx,&indiceY);//encontramos el borde de abajo
	      jp=(int)indiceY;
	      if( (urmax*((double)rand())<=PP) || isSeed(i,j,seedsX,seedsY) ){
		switch(int(vmalla[ip][jp])){//eskina de la izkierda abajo
		case 0:
		  break;
		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
		case 6:
		  //agregamos el sitio con una altura y viento interpolados
		  forest.push_back(sitio(i,j,
					 malla[ip][jp]*(1.0-distBx)*(1.0-distBy) +
					 malla[ip+1][jp]*distBx*(1.0-distBy) +
					 malla[ip][jp+1]*(1.0-distBx)*distBy +
					 malla[ip+1][jp+1]*distBx*distBy
					 ));
		  if(boolwind){
		    forest[contador].wx=
		      wxmalla[ip][jp]*(1.0-distBx)*(1.0-distBy) +
		      wxmalla[ip+1][jp]*distBx*(1.0-distBy) +
		      wxmalla[ip][jp+1]*(1.0-distBx)*distBy +
		      wxmalla[ip+1][jp+1]*distBx*distBy;
		    forest[contador].wy=
		      wymalla[ip][jp]*(1.0-distBx)*(1.0-distBy) +
		      wymalla[ip+1][jp]*distBx*(1.0-distBy) +
		      wymalla[ip][jp+1]*(1.0-distBx)*distBy +
		      wymalla[ip+1][jp+1]*distBx*distBy;
		  }
		  diccionario[i][j]=contador;
		  //forest[contador].??=??;
		  ++contador;
		  break;
		default:
		  break;
		}
	      }
	    }
	  }
	}else{ // !boolveg
	  // Generamos la malla con la distribucion
	  // de sitios dependiendo de P.
	  cerr<<endl<<"Se genera la vegetacion de forma aleatoria"<<endl;
	  cerr.flush();
	  for(int i=0; i<Nx; ++i){
	    for(int j=0; j<Ny; ++j){
	      cerr<<"sitio "<<i<<" "<<j<<endl;
	      cerr.flush();
	      forest.push_back(sitio(i,j,malla[i][j]));
	      if(maux[i*Ny+j]<=PP || isSeed(i,j,seedsX,seedsY) ){
		forest.push_back(sitio(i,j,malla[i][j]));
		diccionario[i][j]=contador;
		++contador;
	      }
	    }
	  }

	  // guardamos la altura de la semilla
	  //nota: no se para ke... O_O
	  xi0=round(Nx*x0);
	  yj0=round(Ny*y0);
	  if(0<=x0 && x0<=1 && 0<=y0 && y0<=1){
	    zk0=forest[diccionario[xi0][yj0]].z;
	  }else{
	    xi0=0.0;
	    yj0=0.0;
	    zk0=0.0;
	  }
	}

	//limpiamos las conexiones de los sitios
	for(tree1=forest.begin(); tree1!=forest.end(); ++tree1){
	  int i=tree1->orden;
	  forest[i].linksTo.clear();
	  forest[i].linked=false;
	  //forest[i].findNeigh(forest,diccionario);
	}

	/* Lo hacemos mas abajo, junto antes de propaga()r y solo para los ffo
	//encontramos las conexiones de los sitios
	for(tree1=forest.begin(); tree1!=forest.end(); ++tree1){
	  int i=tree1->orden;
	  //forest[i].linksTo.clear();
	  forest[i].findNeigh(forest,diccionario);
	}
	*/

	double xxx,xxxx;
	count(forest,front,quem,activ,0.01,Cth,xxx,xxxx);
	double massT=activ;
	cerr<<"-El boske final tiene "<<massT<<" sitios"
	    <<endl;
	cerr<<"-El boske final tiene "<<forest.size()<<" sitios"
	    <<endl;

	// Escribimos la CI
	// Si x,y e [0,1] entonces iniciamos como punto
	// si no, como linea
	mass=0;
	ffo.clear();
	forestend=forest.end();
	int semilla=0;
	int added=0;
	// puntos
	for(tree1=forest.begin(); tree1<forestend; ++tree1){
	  if( isSeed(tree1->x,tree1->y,seedsX,seedsY) ){
	    tree1->estado=Cth*0.95;
	    ffo.insert(tree1->orden);
	    added++;
	    semilla=tree1->orden;
	    cerr<<"semilla "
		<<tree1->x*pdx<<", "
		<<tree1->y*pdy<<", "
		<<tree1->z<<" "
		<<endl;
	  }
	}
	// si es ke las condiciones iniciales
	//son con los bordes prendidos
	if(x0+y0 != 0.0){
	  if(x0<0 && 0<=y0 && y0<=1){
	    for(tree1=forest.begin(); tree1<forestend; ++tree1){
	      if( tree1->x==0 ){
		tree1->estado=Cth*0.95;
		ffo.insert(tree1->orden);
		added++;
	      }
	    }
	  }else if(x0>1 && 0<=y0 && y0<=1){
	    for(tree1=forest.begin(); tree1<forestend; ++tree1){
	      if( tree1->x==(Nx-1) ){
		tree1->estado=Cth*0.95;
		ffo.insert(tree1->orden);
		added++;
	      }
	    }
	  }else if(y0<0 && 0<=x0 && x0<=1){
	    for(tree1=forest.begin(); tree1<forestend; ++tree1){
	      if( tree1->y==0 ){
		tree1->estado=Cth*0.95;
		ffo.insert(tree1->orden);
		added++;
	      }
	    }
	  }else if(y0>1 && 0<=x0 && x0<=1){
	    for(tree1=forest.begin(); tree1<forestend; ++tree1){
	      if( tree1->y==(Ny-1) ){
		tree1->estado=Cth*0.95;
		ffo.insert(tree1->orden);
		added++;
	      }
	    }
	  }else{
	    cerr<<"condicion inicial invalida"<<endl;
	    return 3;
	  }
	}
	// si no se anhadio, devolvemos error
	if(added==0){
	  cerr<<"CI error"<<endl;
	  return 5;
	}


	forest[0].contador=0;

	// buscamos los puntos en los bordes
	//hay ke mejorar esto, podria no haber ningun arbol justo en el borde.
	//deberia buscar los puntos mas "extremos"
	for(tree1=forest.begin(); tree1<forestend; ++tree1){
	  if( tree1->x==0      ||
	      tree1->x==(Nx-1) ||
	      tree1->y==0      ||
	      tree1->y==(Ny-1)  )
	    forest[tree1->orden].borde=1;
	}

	// the real-deal
	errordFF.clear();
	errordQ.clear();
	errordA.clear();
	double error=0.0,r2=0.0,r2z=0.0;
	massFF.clear();
	massQ.clear();
	massA.clear();
	rhoA.clear();
	rhoFF.clear();
	contador=0;
	REAB=0; // flag that indicates when the border is reached

	while(ffo.size()>0 &&
	      (ENDB*REAB)==0 &&
	      //    contador*dT<=8400){
	      contador*dT<=18000){
	  //	    contador*dT<24601){

	  if(contador%500==0) cerr<<"paso: "<<contador<<" tiempo: "<<contador*dT
				  <<" #(ffo): "<<ffo.size()<<endl;

	  // mientras haya fuego, propagamos...
	  // pero primero,
	  // si se cumplen las condiciones escribimos el avance
	  // en el output
	  // luego propagamos el fuego
	  //cerr<<contador*dT<<" "<<pTimes[0]<<"\n";
	  if(isTime(contador, pTimes)
	     || contador==0){
	    //|| contador%5000==0){
	    //lo pongo aca por mientras, para ke no se demore tanto
	    //sacamos los datos solo para algunos momentos
	    cerr<<"+ "
		<<contador*dT<<"s, "
		<<ffo.size()<<" sitios quemandose, en "
		<<(double(clock())-double(start))/CLOCKS_PER_SEC/60.0<<" m"
		<<endl;

	    msr(forest,quem,r2,r2z);
	    if( kn==kw && abs(PP)>0.0 ){
	      cerr<<"  Tiempo de escritura alcanzado:"<<contador*dT<<endl;

	      //escribimos el estado de los sitios quemados y quemandose,
	      // y encontramos el fire front
	      int i0,i1,i2,i3,i4,i5,i6,i7,i8;
	      int j0,j1,j2,j3,j4,j5,j6,j7,j8;
	      double estado=0.0;
	      for(tree1=forest.begin(); tree1<forestend; ++tree1){

		if (printMap==1){
		  sprintf(sql_exp_map, sql_exp_base_map,
			  tree1->x*pdx,tree1->y*pdy,tree1->z*pdz);
		  sqlite3_exec(db, sql_exp_map, 0, 0, 0);
		}

		estado=tree1->estado;
		if(estado < 1){
		  //escritura del sitio a la base de datos
		  sprintf(sql_exp, sql_exp_base, "ForestState", PP, contador*dT,
			  tree1->x*pdx,tree1->y*pdy,tree1->z*pdz,estado);
		  sqlite3_exec(db, sql_exp,0, 0, 0);
		  //ahora buscamos el firefront buscando sitios
		  // quemandose con vecinos intactos. Lo hacemos asi
		  // en vez de usando el vector linksTo por que
		  // queremos los vecinos inmediatamente adyacentes,
		  // mientras que linksTo podria tener contribuciones
		  // mas lejanas.
		  i0=tree1->x;
		  j0=tree1->y;
		  if(i0>0 && j0>0 && i0<Nx-1 && j0<Ny-1){
		    i1=diccionario[i0+1][j0]; i5=diccionario[i0+1][j0+1];
		    i2=diccionario[i0][j0+1]; i6=diccionario[i0+1][j0-1];
		    i3=diccionario[i0-1][j0]; i7=diccionario[i0-1][j0+1];
		    i4=diccionario[i0][j0-1]; i8=diccionario[i0-1][j0-1];
		    if ((i1>=0 && forest[i1].estado==1) || (i5>=0 && forest[i5].estado==1) ||
			(i2>=0 && forest[i2].estado==1) || (i6>=0 && forest[i6].estado==1) ||
			(i3>=0 && forest[i3].estado==1) || (i7>=0 && forest[i7].estado==1) ||
			(i4>=0 && forest[i4].estado==1) || (i8>=0 && forest[i8].estado==1) ){
		      sprintf(sql_exp, sql_exp_base, "FireFront", PP, contador*dT,
			      tree1->x*pdx,tree1->y*pdy,tree1->z*pdz,estado);
		      sqlite3_exec(db, sql_exp,0, 0, 0);
		    }
		  }
		}
	      }

	      printMap=0; //imprimimos el mapa una sola vez


	      //imprime(forest,diccionario,"estadoP");
	      //imprime(forest,"estado", *malla);
	      //imprime(forest,diccionario,"estado");
	    }
	  }

	  if(FPRINT==1 && PP>0){
	    imprime(forest,contador*dT,"fuel");
	  }

	  iffoe=ffo.end();
	  for(iffo=ffo.begin(); iffo!=iffoe; ++iffo){
	    if(! forest[*iffo].linked){
	      forest[*iffo].findNeigh(forest,diccionario);
	    }
	  }
	  propaga(forest,ffo);

	  //                          count(forest,bso,front,quem,activ,0.1,0.9,r2,r2z);
	  //                          massFF.push_back(front);
	  //                          massQ.push_back(quem);
	  //                          massA.push_back(activ);
	  //quem=bso.size();

	  /*
	  //LO SACO SOLO PARA TESTEAR!! HAY KE PONERLO DE VUELTA PARA TENER VALORES!!
	  msr(forest,quem,r2,r2z);
	  pushadd(Chi,iPP,contador,quem);
	  pushadd(R2,iPP,contador,r2);
	  pushadd(R2z,iPP,contador,r2z);
	  pushadd(Tdin,iPP,contador,1.0);
	  */
	  contador++;
	}

	// we add to the survival probability counter
	Ps[iPP].push_back(contador-1);
	// time[iPP]=contador*dT;

	//imprime(massFF,massQ,massA,"Masas");
	//imprime(massFF,massT,"Mfrente");
	//imprime(massQ,massT,"Mquemados");
	//imprime(massA,massT,"Mactivos");

	// Revisamos si percola (llega a todos los bordes)
	// para eso primero necesitamos seguir
	// el proceso, por si paro por ENDB
	//while( ffo.size()>0 ){
	//        propaga(forest,ffo);
	//}
	/*                        int b1=0;
				  int b2=0;
				  int b3=0;
				  int b4=0;
				  for(tree1=forest.begin(); tree1<forestend; ++tree1){
				if(tree1->estado<1.0){
				if( tree1->x==0      )b1++;
				if( tree1->x==(Nx-1) )b2++;
				if( tree1->y==0      )b3++;
				if( tree1->y==(Ny-1) )b4++;
				}
				}
				// si percola anhadimos 1
				// al contador para esta densidad
				if( (b1+b2+b3+b4)>0 )dens[iPP]+=1.0;
	*/
	if( REAB==1 ){
	  dens[iPP]+=1.0/kmax;
	  //imprime(forest,"estado", *malla);
	}

	// sumamos la masa kemada y la
	// densidad=(masa kemada)/(masa total)
	// a los contadores para esta densidad
	//                        massQT[iPP]+=massQ[massQ.size()-1];
	//                        rhoQ[iPP]+=massQ[massQ.size()-1]/massT;
    }// for iPP

    // encontramos el umbral de percolacion Pc
    // para cada realizacion (Pci)
    double dermax=(dens[1]-dens[0])/(pPP[1]-pPP[0]);
    double der=0;
    for(int i=1; i<nPP; i++){
      der=(dens[i]-dens[i-1])/(pPP[i]-pPP[i-1]);
      if(abs(der)>=abs(dermax)){
	dermax=der;
	Pci[kn]=pPP[i];
      }
    }

  }// k

  /*
  //we take the mean of the curves Chi, R2, R2z

  for(int iPP=0; iPP<nPP; ++iPP){
    int imax=Chi[iPP].size();
    double aaa=1.0;
    for(int i=0; i<imax; i++){
      aaa=1.0/Tdin[iPP][i];
      Chi[iPP][i]=Chi[iPP][i]*aaa;
      R2[iPP][i]=R2[iPP][i]*aaa;
      R2z[iPP][i]=R2z[iPP][i]*aaa;
    }
  }

  for(int iPP=0; iPP<nPP; ++iPP){
    massQT[iPP]  =Chi[iPP].back();
    rhoQ[iPP]    /=kmax;

    meanDMFF[iPP]/=kmax;
    meanDMQ[iPP] /=kmax;
    meanDMA[iPP] /=kmax;
    sigmaDMFF[iPP]/=kmax;
    sigmaDMQ[iPP] /=kmax;
    sigmaDMA[iPP] /=kmax;
    sigmaDMFF[iPP] =sqrt(sigmaDMFF[iPP]-meanDMFF[iPP]*meanDMFF[iPP]);
    sigmaDMQ[iPP]  =sqrt(sigmaDMQ[iPP]-meanDMQ[iPP]*meanDMQ[iPP]);
    sigmaDMA[iPP]  =sqrt(sigmaDMA[iPP]-meanDMA[iPP]*meanDMA[iPP]);
  }

  imprime(pPP,massQT,"massQT");
  imprime(pPP,dens,"percolacion");
  imprime(Pci,walpha,Tmean,"Pci");
  // imprime(pPP,time,"time");

  // Calculamos y escribimos la probabilidad
  // de supervivencia dependiente del tiempo
  double pt=0.0;
  int Tmax=0;
  for(int ii=0; ii<nPP;ii++){
    if(Chi[ii].size()>Tmax)Tmax=Chi[ii].size();
  }
  //int incr=max(1.0,(double)Tmax*0.01);
  int incr=10000/200;
  for(int ii=0; ii<nPP;ii++){
    Tmax=Chi[ii].size();
    for(int iT=0; iT<Tmax; iT+=incr){
      pt=survive(Ps,ii,iT);
      out1<<"<dinamica> "
          <<pPP[ii]<<" "
          <<iT*dT<<" "
          <<pt<<" "
          <<Chi[ii][iT]<<" "
          <<R2[ii][iT]<<" "
          <<R2z[ii][iT]<<" "
          <<"</dinamica>"<<endl;
    }
    out1<<"<dinamica> </dinamica>"<<endl;
  }
  */

  finish = clock();

  cerr<<endl;
  cerr<<(double(finish)-double(start))/CLOCKS_PER_SEC<<" s = ";
  cerr<<(double(finish)-double(start))/CLOCKS_PER_SEC/60.0<<" m";
  cerr<<endl;

  sqlite3_close(db);
  return 0;
}

void pushadd(vector<vector<double> >&Chi,int iPP,int contador,double ratio){
  try{
    Chi[iPP].at(contador)+=ratio;
  }catch (...){
    if(Chi[iPP].size()==contador){
      Chi[iPP].push_back(ratio);
    }else{
      cerr<<"Err FATAL:pushadd "<<contador<<" "<<Chi[iPP].size()<<endl;
    }
  }
}

double survive(vector<vector<double> >& Ps,int iPP,int tiempo){
  int total=Ps[iPP].size();
  int cont=0;
  for(int i=0; i<total; i++){
    if(tiempo<=Ps[iPP][i])cont++;
  }
  return ((double)cont)/((double)total);
}

void propaga(vector<sitio>&forest,set<int>&ffo){
  // 0     = quemado
  // (0,1) = qumandose
  // 1     = activo

  // creamos un vector para guardar los cambios antes de pasarlos
  // a forest. de esta forma nos aseguramos de una propagacion "a un tiempo"
  set <int> nuevo;
  set <int>::iterator orden,ffoend;
  //vector <sitio>::iterator tree1;
  vector <int>::iterator tree1,linkend;

  // iteramos sobre todos los sitios quemandose,
  // identificamos a sus vecinos
  // iteramos sobre ellos y cambiamos sus estados (si estan activos)
  // Los anhadimos al nuevo vector de FireFront
  //
  // por mientras definimos 0.01 como la cota minima
  // para la cual hay fuego
  ffoend=ffo.end();
  for(orden=ffo.begin(); orden!=ffoend; ++orden){
    int i1=*orden;
    linkend=forest[i1].linksTo.end();
    for(tree1=forest[i1].linksTo.begin();
        tree1<linkend;
        ++tree1){
      if( forest[*tree1].estado>0.01 ){
        nuevo.insert(*tree1);
        forest[*tree1].burnedBy(forest[i1]);
      }
    }
  }

  // insertamos los nuevos sitios
  // lo hacemos despues del loop principal
  // para evitar loops infinitos
  for(orden=nuevo.begin(); orden!=nuevo.end(); ++orden){
    ffo.insert(*orden);
  }

  // avanzamos un paso temporal para cada sitio que participa.
  // usamos "nuevo" para guardar los lugares que deberiamos borrar
  // (pues no se keman).
  nuevo.clear();
  for(orden=ffo.begin(); orden!=ffo.end(); ++orden){
    forest[*orden].avanza();
    //    if(forest[*orden].umbral <= forest[*orden].estado ||
    if(forest[*orden].isBurning < 1  ||
       forest[*orden].estado<=0.01){
      nuevo.insert(*orden);
    }
  }

  // eliminamos los sitios "muertos"
  // lo hacemos despues de avanzar pues no sabemos con
  // certeza que pasa con el puntero al eliminar un sitio
  for(orden=nuevo.begin(); orden!=nuevo.end(); ++orden){
    ffo.erase(*orden);
  }

  return;}
