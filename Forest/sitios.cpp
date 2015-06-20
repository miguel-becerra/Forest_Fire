/*
 * sitios.cpp
 *
 *  Created on: 20-06-2015
 *      Author: samsung
 */

#include "sitios.h"
#include "radiation.h"

using namespace std;

// definidos en main.cpp
extern int Nx;
extern int Ny;
extern int xi0;
extern int yj0;
extern double dzf;
extern double zk0;
extern double BD;
extern double PP;
extern double l;
extern double dT;
extern double Tmean; // vida media
extern double walpha; // coef de degradacion "termica"
extern double Cth; // degradation/burning threshold
extern double dCth; // degradation/burning threshold for the derivative
extern double dxmalla, slope;
extern ofstream out1;
extern int REAB; // flag that indicates when the border is reached

extern double pdx;
extern double pdy;
extern double pdz;

extern vegetation veg;

sitio::sitio(int posx, int posy, double altura){
  x=posx;
  y=posy;
  z=altura;
  orden=contador;
  contador++;
  estado=1.0;
  estado0=1.0;
  linked=false;
  lx=l; //rangos maximos
  ly=l;

  // vientos
  wx=0.0;
  wy=0.0;

  brate=1.0/Tmean*dT; // u_i*dt para facilitar solucion numerica

  f0=0.0; // radiacion recibida el paso anterior
  f1=0.0; // radiacion recibida este paso
  umbral=Cth; // C_th
  dtumbral=dCth*dT; // para facilitar el avance en .avanza

  isBurning=0;

  borde=0; // esta en el borde?
  //se usa para saber cuando detener la iteracion
}

bool sitio::pointsTo(sitio b){
  int dx=abs(b.x-x);
  int dy=abs(b.y-y);
  double lx2=lx*lx;
  // primera seleccion de elementos
  // nos quedamos solo con los elementos
  // dentro del rango deseado
  //    if( dx>lx || dy>ly)
  //        return false;

  // si el elemento es el mismo
  // no lo consideramos
  if( dx==0 && dy==0 )
    return false;

  // si hay l pasos de distancia y b esta
  // mas alto, entonces tambien es vecino
  //    if(dx+dy==lx && z<=b.z)
  //        return true;

  // es vecino si hay 'l' pasos de distancia
  if( (dx*dx+dy*dy) > lx2 ){
    return false;
  }else{
    return true;
  }


  /*
  // seleccion como elipse
  double dx2=(x-b.x)*(x-b.x);
  double dy2=(y-b.y)*(y-b.y);
  double lx2=lx*lx;
  double ly2=ly*ly;
  if(dx2/lx2 + dy2/ly2 <= 1) return true;
  */

  // en cualquier otro caso, no es vecino
  return false;
}

void sitio::addNeigh(sitio b){
  linksTo.push_back(b.orden);
}

void sitio::addedTo(sitio b){
  b.linkedBy.push_back((*this).orden);
}

/*
  void sitio::findNeigh(vector <sitio>&forest){
  static vector <sitio>::iterator tree1,forestend;
  forestend=forest.end();
  for(tree1=forest.begin(); tree1!=forestend; ++tree1){
  if( (*this).pointsTo(*tree1) ){
  //(*this).addNeigh(*tree1);
  //(*tree2).addedTo(*this);
  (*this).linksTo.push_back(tree1->orden);
  }
  }
  }
*/
void sitio::findNeigh(vector <sitio>&forest, int **diccionario){
  int dx=0;
  int dy=0;
  int finx=min((int)ceil(x+lx),Nx-1);
  int finy=min((int)ceil(y+ly),Ny-1);
  int inix=max((int)floor(x-lx),0);
  int iniy=max((int)floor(y-ly),0);
  double dist=0;
  for(int i=inix; i<=finx; i++){
    dx=abs(x-i);
    for(int j=iniy; j<=finy; j++){
      dy=abs(y-j);
      dist=sqrt(dx*dx+dy*dy);
      if( dist>0 && // no son los 2 iguales a cero
	  dist<=lx && // el N de pasos es a lo mas lx
	  diccionario[i][j]>=0){
	(*this).linksTo.push_back(diccionario[i][j]);
      }
    }
  }
  (*this).linked=true;
}

void sitio::burnedBy(sitio b){
  int    dx=(x-b.x);// desde el sitio ke emite al ke recibe
  int    dy=(y-b.y);
  double dz=(z-b.z);
  double cap=12000.0;   // cap calorica [joule/kilo]
  double mass=17.5;     // masa del arbol [kilo]
  double r=0.35;        // % de energia como radiacion [1]
  // walpha esta en los extern al comienzo del archivo
  double ww=r*walpha*mass*cap;
  double W=veg.radiation(dx*pdx, dy*pdy, dz*pdz, wx, wy); // Factor de forma (geometrico)

  f1+=(b.estado-b.estado0)*ww*W;
}

void sitio::avanza(){
  //We use predictor method Adam-Bashforth 3
  //to solve the differential equation
  //for the fuel
  double estadoA;
  isBurning=(int)Theta(dtumbral-estado+estado0);
  estadoA = estado
    + 0.5*(
	   +3.0*estado*f1
	   -3.0*isBurning*estado*(1.0-estado)*brate
	   -estado0*f0
	   +isBurning*estado0*(1.0-estado0)*brate
	   );
  estado0	= estado;
  estado	= estadoA;
  if(estado<0)estado=0.0;
  f0		= f1;
  f1		= 0.0;

  // after the advance in time and if we are in the border
  // we mark the border as reached
  if(borde==1)REAB=1;
  }
// FIN METODOS DE CLASE SITIO

double Theta(double x){
  /*	if(x>=0){
	return 1.0;
	}else{
	return 0.0;
	}
  */
  return (x>=0?1.0:0.0);
}

double Sign(double x){
  return (x>=0?1.0:-1.0);
}

void count(vector <sitio>&forest, double &ff, double &q,
	   double &a, double li, double lf, double &r2, double &r2z){
  double dx1=0.0,dy1=0.0,dz1=0.0;
  ff =0.0;
  q  =0.0;
  a  =0.0;
  r2 =0.0;
  r2z=0.0;
  vector <sitio>::iterator tree1,forestend;
  forestend=forest.end();
  for(tree1=forest.begin(); tree1<forestend; ++tree1){
    //if(0<=tree1->estado && tree1->estado<=li){
    //we include firefront and burned places
		//to make the curve smooth
    if(0<=tree1->estado && tree1->estado<=lf){
      q++;
      // medimos r2 y r2z
      // usamos el 3 correspondiente a los 3m de radiation.cpp
      dx1=pdx*(tree1->x-xi0);
      dx1*=dx1;
      dy1=pdy*(tree1->y-yj0);
      dy1*=dy1;
      dz1=pdz*(tree1->z-zk0);
      dz1*=dz1;
      dx1+=dy1;
      r2+=dx1;
      r2z+=dx1 + dz1;
    }else if(li<tree1->estado && tree1->estado<=lf){
      ff++;
    }else if(lf<tree1->estado && tree1->estado<=1.0){
      a++;
    }
  }
  double uq=1.0/max(0.1,q);
  r2*=uq;
  r2z*=uq;

  return;
}

void msr(vector <sitio>&forest,
	 double &quem,
	 double &r2,
	 double &r2z){
  double dx1=0.0,dy1=0.0,dz1=0.0,zk1=0.0,q=0.0;
  int xi1=0,yj1=0;
  r2 =0.0;
  r2z=0.0;
  quem=0.0;
  vector <sitio>::iterator tree1,fend;
  fend=forest.end();
  for(tree1=forest.begin(); tree1!=fend; ++tree1){
    if(tree1->estado<tree1->umbral){
      // medimos r2 y r2z
      // usamos el 3 correspondiente a los 3m de radiation.cpp
      xi1=tree1->x;
      yj1=tree1->y;
      zk1=tree1->z;
      dx1=pdx*(xi1-xi0);
      dx1*=dx1;
      dy1=pdy*(yj1-yj0);
      dy1*=dy1;
      dz1=pdz*(zk1-zk0);
      dz1*=dz1;
      dx1+=dy1;
      r2+=dx1;
      r2z+=dx1 + dz1;
      q++;
    }
  }
  double uq=1.0/max(0.1,q);
  r2*=uq;
  r2z*=uq;
  quem=q;

  return;
}
//---------------
void imprime(vector <sitio>&forest,int **diccionario, string label){
  if(label=="estado" || label=="estadoP"){
    //    set <int>::iterator iffo,iffoend;
    //iffoend=ffo.end();

    out1<<"<"<<label<<"> "
	<<"#nt "
	<<"</"<<label<<"> "<<endl;
    /*
    for(iffo=ffo.begin(); iffo!=iffoend; ++iffo){
      out1<<"<"<<label<<"> "
	  <<PP<<" "
	  <<forest[*iffo].x*pdx<<" "
	  <<forest[*iffo].y*pdy<<" "
	  <<forest[*iffo].z*pdz<<" "
	  <<forest[*iffo].estado<<" "
	  <<"</"<<label<<">"<<endl;
    }
    */

    int nn=0;
    int i1,i2,i3,i4,i0;
    for(int i=2; i<Nx-2; i++){
      for(int j=2; j<Ny-2; j++){
	i0=diccionario[i][j];
	if(i0>=0 && forest[i0].estado<1){
	  i1=diccionario[i+1][j];
	  i2=diccionario[i][j+1];
	  i3=diccionario[i-1][j];
	  i4=diccionario[i][j-1];
	  if ((i1>=0 && forest[i1].estado==1) ||
	      (i2>=0 && forest[i2].estado==1) ||
	      (i3>=0 && forest[i3].estado==1) ||
	      (i4>=0 && forest[i4].estado==1) ){
	    out1<<"<"<<label<<"> "
		<<PP<<" "
		<<forest[i0].x*pdx<<" "
		<<forest[i0].y*pdy<<" "
		<<forest[i0].z*pdz<<" "
		<<forest[i0].estado<<" "
		<<"</"<<label<<">"<<endl;
	  }
	}
      }
    }
    cerr<<"----"<<nn<<endl;
    out1<<"<"<<label<<"> </"<<label<<">"<<endl;
    out1<<"<"<<label<<"> </"<<label<<">"<<endl;
  }
}
//---------------
void imprime(vector <sitio>&forest,string label, double *malla){
  if(label=="estado"){
        double mallabin[Nx][Ny];
        for(int i=0;i<Nx;i++)
	  for(int j=0;j<Ny;j++)
	    mallabin[i][j]=0.0;

        vector <sitio>::iterator tree1,forestend;
	forestend=forest.end();
        for(tree1=forest.begin(); tree1<forestend; ++tree1){
	  if(tree1->estado==1.0){
	    mallabin[tree1->x][tree1->y]=1.0;
	  }else{
	    mallabin[tree1->x][tree1->y]=-tree1->estado;
	  }
        }

        out1<<"<"<<label<<"> "
            <<"#nt "
            <<"</"<<label<<"> "<<endl;
        for(int i=0;i<Nx;i++){
	  for(int j=0;j<Ny;j++){
            out1<<"<"<<label<<"> "
                <<PP<<" "
                <<i*dxmalla*cos(slope)<<" "
                <<j*dxmalla<<" "
		<<mallabin[i][j]<<" "
		<<malla[i*Ny+j]<<" "
                <<"</"<<label<<">"<<endl;
	  }out1<<"<"<<label<<"> </"<<label<<">"<<endl;
        }out1<<"<"<<label<<"> </"<<label<<">"<<endl
             <<"<"<<label<<"> </"<<label<<">"<<endl;
  }
}

void imprime(vector <sitio>&forest,double T, string label){
  out1<<"<"<<label<<"> "
		<<T<<" ";
  for(int i=0; i<forest.size(); i++){
    out1<<forest[i].estado<<" ";
  }
  out1<<"</"<<label<<">"<<endl;
}

void imprime(vector <double>&Pci, double al, double tau, string label){
  out1<<"<"<<label<<"> "
      <<BD<<" "
      <<al<<" "
      <<tau<<" ";
  for(int i=0; i<Pci.size(); i++){
    out1<<Pci[i]<<" ";
  }
  out1<<"</"<<label<<">"<<endl;

  out1<<"<"<<label<<"> ";
  out1<<"</"<<label<<">"<<endl;
  out1<<"<"<<label<<"> ";
  out1<<"</"<<label<<">"<<endl;

  double s=0.0;
  double s2=0.0;
  for(int i=0; i<Pci.size(); i++){
    s+=Pci[i];
    s2+=Pci[i]*Pci[i];
  }
  s=s/Pci.size();
  s2=s2/Pci.size();

  out1<<"<"<<label<<"> "
      <<BD<<" "
      <<s<<" "
      <<sqrt(s2-s*s)<<" ";
  out1<<"</"<<label<<">"<<endl;

}

void imprime(vector <double>&dim,string label){
  for(int i=0; i<dim.size(); i+=1){
    out1<<"<"<<label<<"> "
	<<PP<<" "
	<<i<<" "
	<<dim[i]<<" "
	<<"</"<<label<<">"<<endl;
  }
  out1<<"<"<<label<<"> </"<<label<<">"<<endl;
  out1<<"<"<<label<<"> </"<<label<<">"<<endl;
}

void imprime(vector <double>&dim,double aa,string label){
  for(int i=0; i<dim.size(); i+=1){
    out1<<"<"<<label<<"> "
	<<PP<<" "
	<<i<<" "
	<<dim[i]<<" "
	<<aa<<" "
	<<"</"<<label<<">"<<endl;
  }
  out1<<"<"<<label<<"> </"<<label<<">"<<endl;
  out1<<"<"<<label<<"> </"<<label<<">"<<endl;
}

void imprime(vector <double>&m1,
	     vector <double>&m2,
	     vector <double>&m3,
	     string label){
  for(int i=0; i<m1.size(); i+=1){
    out1<<"<"<<label<<"> "
	<<i*dT<<" "
	<<m1[i]<<" "
	<<m2[i]<<" "
	<<m3[i]<<" "
	<<"</"<<label<<">"<<endl;
  }
  out1<<"<"<<label<<"> </"<<label<<">"<<endl;
   out1<<"<"<<label<<"> </"<<label<<">"<<endl;
}

void imprime(vector <double>&dim, vector <double>&error,string label){
  for(int i=0; i<dim.size(); i+=1){
    out1<<"<"<<label<<"> "
	<<PP<<" "
	<<i<<" "
	<<dim[i]<<" "
	<<error[i]<<" "
	<<"</"<<label<<">"<<endl;
  }
  out1<<"<"<<label<<"> </"<<label<<">"<<endl;
  out1<<"<"<<label<<"> </"<<label<<">"<<endl;
}

void imprime( double a,string label){
  out1<<"<"<<label<<"> "
      <<a
      <<" </"<<label<<">"<<endl;
}

void imprime( double a,double b, double c, double d,string label){
  out1<<"<"<<label<<"> "
      <<a<<" "
      <<b<<" "
      <<c<<" "
      <<d<<" "
      <<" </"<<label<<">"<<endl;
}

void imprime( string label){
  out1<<"<"<<label<<"> "
      <<" </"<<label<<">"<<endl;
}

void imprime( double dT, vector<double>& pPP,
	      vector<vector<double> >&Ps, string label){
  int nT,nPP=Ps.size();
  for(int i=0; i<nPP;i++){
    nT=Ps[i].size();
    double s=0.0;
    double s2=0.0;
    for(int j=0; j<nT; j++){
      s+=Ps[i][j];
      s2+=Ps[i][j]*Ps[i][j];
    }
    s=s/nT;
    s2=s2/nT;
    out1<<"<"<<label<<"> "
	<<BD<<" "
	<<pPP[i]<<" "
	<<s<<" "
	<<sqrt(s2-s*s)<<" "
	<<" </"<<label<<">"<<endl;
  }
  out1<<"<"<<label<<"> </"<<label<<">"<<endl;
  out1<<"<"<<label<<"> </"<<label<<">"<<endl;
}

void vectorMap(double **mallaMap){
  //Normaliza la altura de la superficie
  int rango=l;
  ofstream out2;
  out2.open("vector.dat");
  double gx[Nx][Ny];
  double gy[Nx][Ny];
  double gmax=0;

  for(int i=1;i<Nx-1;i++){
    for(int j=1;j<Ny-1;j++){
      gx[i][j]=mallaMap[i+1][j]-mallaMap[i-1][j];
      gy[i][j]=mallaMap[i][j+1]-mallaMap[i][j-1];
      if(fabs(gx[i][j])>=gmax)gmax=fabs(gx[i][j]);
      if(fabs(gy[i][j])>=gmax)gmax=fabs(gy[i][j]);
    }
  }
  for(int i=1;i<Nx-1;i++){
    for(int j=1;j<Ny-1;j++){
      out2<<i<<" "
	  <<j<<" "
	  <<gx[i][j]/gmax<<" "
	  <<gy[i][j]/gmax<<" "
	  <<endl;
    }
  }
  out2.close();
}

void mass_analysis(vector <double>&mass,double BD, double PP,string label){
  int nmass=mass.size();
  double total=0,total2=0, max=0;
  for(int i=0; i<nmass; i++){
    if(mass[i]>0){ //get rid of nan's
      total+=mass[i];
      total2+=mass[i]*mass[i];
      if(mass[i]>max)max=mass[i];
    }
  }

  out1<<"<"<<label<<"> "
      <<BD<<" "
      <<PP<<" "
      <<total<<" "
      <<total2<<" "
      <<max<<" "
      <<nmass<<" "
      <<"</"<<label<<">"<<endl;
}

void mrandom(int NN, double* maux){
  /*	sin GSL
   */
 	static double abajo=1.0/(double)RAND_MAX;
	for(int i=0;i<NN;i++){
	maux[i]=((double)rand())*abajo;
	}
  /*	con GSL

  const gsl_rng_type * T;
  gsl_rng * r;
  //	gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r,time(NULL));

  cerr<<"random "<<gsl_rng_uniform(r)<<endl;
  for(int i=0;i<NN;i++){
    maux[i]=gsl_rng_uniform(r);
  }
  gsl_rng_free (r);
  */
}
