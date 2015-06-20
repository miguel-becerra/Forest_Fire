/*
 * radiation.h
 *
 *  Created on: 20-06-2015
 *      Author: samsung
 */

#ifndef RADIATION_H_
#define RADIATION_H_

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

//distancias fisicas
extern double pdx;
extern double pdy;
extern double pdz;

class vegetation;

class vegetation{
 public:
  int isStarted;

  // variables fisicas
  double canopiH;
  double canopiD;
  double flameH;
  double flameD;
  double Tmean;
  double alpha;
  double Cth;

  //variables de discretizacion
  double Dzf;
  double Dzc;

  int nDd;  //numero de puntos para interpolar la distancia
  double Dd;
  int nDz;  //numero de puntos parainterpolar la altura
  double Dz;

  //matriz para interpolar
  vector<vector<double> > factor;
  double **factor2;

  // metodos
  //          c hei   c diam  f hei   f dia   Tmean   alpha   cth
  vegetation (double, double, double, double, double, double, double);
  void init_vfactor();
  void set(double, double, double, double, double, double, double);
  double FF(double, double);
  //               dx      dy      dz      wx      wy
  double radiation(double, double, double, double, double);
};

#endif /* RADIATION_H_ */
