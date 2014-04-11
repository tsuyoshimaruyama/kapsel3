#ifndef PERIODIC_SPLINE_H
#define PERIODIC_SPLINE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string>
#include <cstring>

#include "alloc.h"
#include "macro.h"

using std::string;
/*!
  \file periodic_spline.h
  \author J. Molina
  \date 2014/04/03
  \version 1.0
  \brief Spline interpolation for 1D periodic systems using equally
  spaced sample points
 */

extern int     splN;       //number of interpolation points
extern int     splDim;     //interpolation order + 1 (= 4)
extern double  splDx;      //grid spacing for interpolating points
extern double *splQ;       //solution vector
extern double *splAii;     //diagonal of coefficient matrix
extern double *splAin;     //last column of coefficient matrix
extern double **splC;      //spline coefficienst a,b,c,d

void splineInit(const int &n, const double &dx);
void splineFree();
void splineReset(const double *fx);
double fspline(const double &x);
#endif
