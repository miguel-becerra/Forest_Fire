#include "sitios.h"
#include "radiation.h"

using namespace std;

// Cambio de estrategia con respecto a la ultima version:
// -Antes:
//  . Se calcula la radiaiacion recibida por un un elemento pekenho, por parte
//    de 2 mitades de la llama, la superior e inferior.
// -Ahora:
//  . La idea es incorporar el viento, por lo ke hay ke hacer calculos con la llama
//    inclinada. Para lograr esto, la aproximacion ke usamos es la siguiente:
//    -> Guardamos en un arreglo la cantidad de radiacion recibida por un arbol
//       (cilindro) completo proveniente de un disco (cilindro pekenho) a una distancia
//       d y altura z
//    -> Sumamos sobre todos los discos puestos formando una diagonal dada por el viento
//       segun "Numerical study of wind effects on the characteristics of fames from
//              non-propagating vegetation fires.
//              F. Nmira, J.L. Consalvi, P. Boulet, B. Porterie".
//
// Otro cambio tecnico es el de pasar de namespace a una clase, asi despues 
// podria cambiarse facilmente a distintos tipos de vegetacion.
        

extern double l,l2;
int ll=ceil(l)+1;

extern double pdx;
extern double pdy;
extern double pdz;

vegetation::vegetation (double ch, double cd, 
			double fh, double fd, 
			double Tmeam, double alpha, double cth){
  isStarted=0;

  canopiH=ch;
  canopiD=cd;
  flameH=fh;
  flameD=fd;
  Dzf=flameH/(100.0-1.0);  //altura del disco ke emite
  Dzc=canopiH/(100.0-1.0); //altura del disco ke recibe

  Tmean=20;
  alpha=0.000005;
  Cth=0.9;

  nDd=100; //solo valores default, ver init_vfactor()
  Dd =0.1;
  nDz=100;
  Dz=0.1;
}

void vegetation::set (double ch, double cd, 
		      double fh, double fd, 
		      double Tmeam, double alpha, double cth){
  isStarted=0;

  canopiH=ch;
  canopiD=cd;
  flameH=fh;
  flameD=fd;
  Dzf=flameH/(100.0-1.0);  //altura del disco ke emite
  Dzc=canopiH/(100.0-1.0); //altura del disco ke recibe

  Tmean=20;
  alpha=0.000005;
  Cth=0.9;

  nDd=100; //solo valores default, ver init_vfactor()
  Dd =0.1;
  nDz=100;
  Dz=0.1;
}

void vegetation::init_vfactor(){
  if(isStarted){
    //factor.clear();
    for(int i=0; i<nDd; i++){
      delete [] factor2[i];
    }
    delete [] factor2;
  }
  isStarted=1;
  
  //numero de segmentos en los ke se divide la distancia al foco.
  // cada arista de la malla se divide en ~100 puntos.
  // calculamos para 2*(maximo aristas) para dar cuenta de los lejanos ke estan en
  // diagonal. Mas alla de esa distancia lo hacemos 0.
  nDd =(2*l)*150+1; // l viene de un extern mas arriba!
  Dd  =(2.0*l)*pdx/nDd;

  //numero de segmentos en los ke dividimos la diferencia de altura.
  // calculamos la influencia de cada disco en un cilindro completo,
  // para hasta 4 alturas completas, 2 para arriba y 2 para abajo
  nDz =int(4*canopiH)*150+1;
  Dz  =(4.0*canopiH)/nDz;

  factor2 = new double*[nDd];
  for(int i=0; i<nDd; i++){
    //factor.push_back(vector<double>());
    factor2[i]=new double [nDz];
  }
  
  double vf=0.0;
  double dist=0.0;
  for(int i=0; i<nDd; i++){
    dist=i*Dd; //phys distance
    for(int k=0; k<nDz; k++){
      vf=FF(dist,-2.0*canopiH+k*Dz);// partimos de 2 alturas hacia abajo
      //factor[i].push_back(vf);
      factor2[i][k]=vf;
    }
  }
}

double vegetation::FF(double dist, double z){
  //da el flujo sobre un cilindro (dado por la vegetacion) proveniente
  //de un disco de altura dada, a una distancia dist, y una diferencia de altura z
  //tomado de:
  //SFPE - Handbook of - Fire Protection Engineering
  //p 3-276
  if(dist<flameD*0.5){
    return 1.0;
  }
  extern double PI;
  double L=dist;    //distance from center of cylinder to target
  double D=flameD;  //cylinder diameter
  //double H1=dzf;    //cylinder height
  
  double S=2.0*L/D;
  double uSP=1.0/S/PI;
  double B=(1+S*S)*0.5/S;
  //double h=2.0*H1/D;
  //double A=(h*h+S*S+1.0)*0.5/S;
  double hp,Ap,hm,Am;
  double f=0.0;
  //integramos sobre el cilindro ke recibe
  //la idea es, para cada disco ke recibe, calcular la radiacion de un cilindro
  //de altura H+dzf y restar la ke recibe de un cilindro de altura H.
  //parece idiota, pero no se como hacerlo de otra forma
  for(double z0=0.0; z0<=canopiH; z0+=Dzc){
    hp=(abs(z+z0)+Dzf)*2.0/D; 
    Ap=(hp*hp+S*S+1.0)*0.5/S;
    hm=abs(z+z0)*2.0/D;
    Am=(hm*hm+S*S+1.0)*0.5/S;
    // esta geometria funciona si pensamos en ke Dzf->0, 
    //habria ke pensarla un poco mejor para incorporar bien el disco en todas
    //las configuraciones

    f+=uSP*(
	    (atan(hp/sqrt(S*S-1.0))
	    -hp*atan(sqrt((S-1.0)/(S+1.0)))
	    +Ap*hp/sqrt(Ap*Ap-1.0)
	     *atan(sqrt((Ap+1.0)*(S-1.0)/(Ap-1.0)/(S+1.0))))
	    - 
	    (atan(hm/sqrt(S*S-1.0))
	    -hm*atan(sqrt((S-1.0)/(S+1.0)))
	    +Am*hm/sqrt(Am*Am-1.0)
	     *atan(sqrt((Am+1.0)*(S-1.0)/(Am-1.0)/(S+1.0))))
	    );
  }
  double arc=2.0*dist*(PI-2.0*acos(D*0.25/dist));
  return f*arc*Dzc;
}

double vegetation::radiation(double dx, double dy, double dz, double wx, double wy){
  //radiation received by a section* at distance (dx,dy,dz)
  // todos los d? deben apuntar desde la base del sitio kemandose 
  // a la base del sitio por kemar!!!

  //artificially make the fire start towards the top of the canopy
  dz=dz-canopiH*0.7;
  
  //calculamos los angulos de los vectores necesarios: posicion del sitio kemado,
  // direccion del viento e inclinacion de la llama
  double rw  =sqrt(wx*wx + wy*wy);
  double thw =(wx==0.0&&wy==0.0?0.0:atan2(wy,wx));//tiene cuidado con los cuadrantes

  static double uu=1.0/9.8/flameH;
  double Frn =rw*rw*uu; //Froude number=U_{\infty}/g/Hf
  double beta=atan(1.22*sqrt(Frn)); 
  //see F.A. Albini, A model for the window flame from a line fire,  
  //Combustion and Flame 43 (1981) 155-174

  double tb=tan(beta);
  double flx=tb*cos(thw);//(no) rotamos el vector (se expande)
  double fly=tb*sin(thw);
  //double flz=cos(beta);
  //notar ke la llama no "rota", si no ke es como inclinar una torre de monedas

  static double uDd=1.0/Dd;
  static double uDz=1.0/Dz;
  double hx=0.0;
  double hy=0.0;
  int idistXY=0;
  int idistZ=0;
  double sumr=0.0;
  for(double hz=0.0; hz<=flameH; hz+=Dzf){
    hx=hz*flx;
    hy=hz*fly;
    //notar ke la llama no "rota", si no ke es como inclinar una torre de monedas
    idistXY=round(sqrt((dx-hx)*(dx-hx) + (dy-hy)*(dy-hy)) * uDd);
    idistZ =round((2.0*canopiH + (dz-hz)) * uDz); 
    //             ^toma en cuenta el offset de los indices dados en ini_vfactor.
    if(idistXY<nDd && idistZ>=0 && idistZ<nDz){
      //sumr+=FF(idistXY*Dd, idistZ*Dz);
      //      sumr+=factor[idistXY][idistZ];
      sumr+=factor2[idistXY][idistZ];
    }
  }
  
  return sumr;
}
