#ifndef SITIOS_H
#define SITIOS_H

#include <iostream>
#include <fstream>
#include <cstdio> // fwrite()
#include <math.h>
#include <string>
#include <vector>
#include <set>
#include <complex>
#include <ctime>
#include <time.h>
#include <iomanip>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_rng.h>

//#include "mkl.h"
//#include "mkl_types.h"

using namespace std;

extern int Nx;
extern int Ny;

//exportamos las distancias fisicas
extern double pdx;
extern double pdy;
extern double pdz;

class sitio;
void propaga(vector<sitio>&,set<int>&);
double boxD(vector<sitio>&, double, double,double&);

void pushadd(vector<vector<double> >&,int,int,double);
double survive(vector<vector<double> >&, int, int);
//cuidado!! box2D cambia el vector entregado!!
double box2D(vector<double>&, double&);
void count(vector<sitio>&,double&,double&,double&,double,double,
		double&,double&);
void msr(vector <sitio>&,double &,double &, double &);

//double count(vector<sitio>&, double, double);
void mass_analysis(vector<double>&, double, double,string);

void imprime(string);
void imprime(vector <sitio>&, int **, string);
void imprime(vector <sitio>&, string, double *);
void imprime(vector <double>&, string);        
void imprime(vector <double>&, double, string);        
void imprime(vector <double>&, vector <double>&, string);        
void imprime( double, string);                
void imprime( double, vector<double>&, vector<vector<double> >&, string);
void imprime( double, double, double, double, string);                
void imprime(vector <sitio>&, double, string);
void imprime(vector <double>&, double, double, string);

void vectorMap(double **);

double Theta( double);                
double Sign( double);                

// in radiation.cpp
//void distm(double dxmN, double dymN,double dzmN);
//void init_vfactor();
//double FF(double dist, double z);
//double radiation(int dx, int dy, double dz);
//

void mrandom(int, double*);

class sitio{
 public:
  //generamos la clase sitio.
  //x = posicion en x
  //y = posicion en y
  //z = altura del punto
  //orden = numeracion del sitio
  //estado = estado del sitio:  	0 inactivo 
  //				1 activo (opcion por defecto)
  //				(0,1) qumandose
  //lx , ly = definen la zona de influencia. ver pointsTo
  //linksTo = vector con los sitios apuntados por this.
  //linkedBy =  vector con los sitios ke apuntan a this.
  
  int x, y, orden, borde, isBurning;
  static int contador;
  double estado,umbral,dtumbral;
  double z, lx, ly;
  double estado0,f0,f1,brate;
  bool linked;
  vector <int> linksTo,linkedBy;
  
  double wx,wy;
  
  // sitio x, y, z
  sitio (int,int,double);
  bool pointsTo(sitio b);
  void addNeigh(sitio b);
  void addedTo(sitio b);
  void burnedBy(sitio b);
  void avanza();
  void findNeigh(vector<sitio>&, int **);
};


#endif
