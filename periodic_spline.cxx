#include "periodic_spline.h"

int     splN;
int     splDim;
double  splDx;
double *splQ;
double *splAii;
double *splAin;
double **splC;

void splineInit(const int &n, const double &dx){
  assert(n >= 4);
  splDim= 4;
  splDx = dx;
  splN  = n;

  splQ     = alloc_1d_double(splN);
  splAii   = alloc_1d_double(splN);
  splAin   = alloc_1d_double(splN);

  splC = (double**)malloc(4*sizeof(double*));
  for(int d = 0; d < splDim; d++){
    splC[d] = alloc_1d_double(splN);
  }
}
void splineFree(){
  for(int d = 0; d < splDim; d++){
    free_1d_double(splC[d]);
  }
  free(splC);

  free_1d_double(splQ);
  free_1d_double(splAii);
  free_1d_double(splAin);
  
  splN = splDim = 0;
  splDx = 0.0;
}

void splineReset(const double *fx){
  
  //initialize
  {  
    int inext;
    int iprev;
    double *a = splC[0];

    for(int i = 0; i < splN; i++){
      inext = MOD(i+1, splN);
      iprev = MOD(i-1, splN);
      splQ[i]   = 3.0*(fx[inext] - 2.0*fx[i] + fx[iprev])/SQ(splDx);
      splAii[i] = 4.0;
      splAin[i] = 0.0;

      //a_i coefficient
      a[i] = fx[i];
    }

    splAin[0] = splAin[splN-2] = 1.0;
    splAin[splN-1] = 4.0;
  }

  //reduce to upper triangular form
  //only non-zero elements are on diagonal and last column
  {
    int i;
    int i0;
    double alpha;
    //lower part (exluding last row)
    for(i = 1; i < splN - 1; i++){
      i0 = i-1;
      alpha = - 1.0 / splAii[i0];
      splAii[i] += alpha;
      splAin[i] += alpha*splAin[i0];
      splQ[i]   += alpha*splQ[i0];
    }
    
    //upper part
    for(i = splN-3; i >= 0; i--){
      i0 = i+1;
      alpha = -1.0 / splAii[i0];
      splAin[i] += alpha*splAin[i0];
      splQ[i]   += alpha*splQ[i0];
    }
    
    //last row (two bottom edges)
    i = splN - 1;
    i0= 0;
    alpha = - 1.0 / splAii[i0];
    splAii[i] += alpha*splAin[i0];
    splAin[i] += alpha*splAin[i0];
    splQ[i]   += alpha*splQ[i0];
    i0= splN - 2;
    alpha = - 1.0 / splAii[i0];
    splAii[i] += alpha*splAin[i0];
    splAin[i] += alpha*splAin[i0];
    splQ[i]   += alpha*splQ[i0];
  }

  //back substitute for c_i coefficient
  {
    double *c = splC[2];
    double cn = splQ[splN - 1] / splAii[splN - 1];

    c[splN - 1] = cn;
    for(int i = splN - 2; i >= 0; i--){
      c[i] = (splQ[i] - splAin[i]*cn) / splAii[i];
    }
  }

  //b_i & d_i coefficient
  {
    int inext;
    int iprev;

    double *b = splC[1];
    double *d = splC[3];
    const double *a = splC[0];
    const double *c = splC[2];

    for(int i = 0; i < splN; i++){
      inext = MOD(i+1, splN);
      iprev = MOD(i-1, splN);
      d[i] = (c[inext] - c[i])/(3.0*splDx);
      b[i] = (a[inext] - a[i])/splDx - c[i]*splDx - d[i]*SQ(splDx);
    }
  }
}
double fspline(const double &x){
  assert(x >= 0.0 && x < splN * splDx);
  int i = x / splDx;
  double fx    = 0;
  double dmy   = 1;
  double delta = (x - (double)i*splDx);
  for(int d = 0; d < splDim; d++){
    fx += splC[d][i]*dmy;
    dmy*= delta;
  }
  return fx;
}
