#include "periodic_spline.h"

const int splineDim = 4;   //spline coefficients
const int splineDim_= 7;   //spline coefficients + working memory


void splineFree(splineSystem* &spl){
  free_1d_double(spl->a);
  free_1d_double(spl->b);
  free_1d_double(spl->c);
  free_1d_double(spl->d);
  free_1d_double(spl->Q);
  free_1d_double(spl->Aii);
  free_1d_double(spl->Ain);
  delete spl;
  spl = NULL;
}
void splineInit(splineSystem* &spl, const int &n, const double &dx){
  assert(n >= splineDim);
  spl = new splineSystem;
  spl->n   = n;
  spl->dx  = dx;
  spl->a   = alloc_1d_double(n);
  spl->b   = alloc_1d_double(n);
  spl->c   = alloc_1d_double(n);
  spl->d   = alloc_1d_double(n);
  spl->Q   = alloc_1d_double(n);
  spl->Aii = alloc_1d_double(n);
  spl->Ain = alloc_1d_double(n);
  for(int i = 0; i < n; i++){
    spl->a[i] = 0.0;
    spl->b[i] = spl->c[i] = spl->d[i] = 0.0;
  }
}
void splineCompute(splineSystem* spl, const double* fx){
  int n    = spl->n;
  double dx   = spl->dx;
  double* a   = spl->a;
  double* b   = spl->b;
  double* c   = spl->c;
  double* d   = spl->d;
  double* Q   = spl->Q;
  double* Aii = spl->Aii;
  double* Ain = spl->Ain;

  //initialize
  {  
    int inext;
    int iprev;
    for(int i = 0; i < n; i++){
      inext = MOD(i+1, n);
      iprev = MOD(i-1, n);
      Q[i]   = 3.0*(fx[inext] - 2.0*fx[i] + fx[iprev])/SQ(dx);
      Aii[i] = 4.0;
      Ain[i] = 0.0;

      //a_i coefficient
      a[i] = fx[i];
    }

    Ain[0] = Ain[n-2] = 1.0;
    Ain[n-1] = 4.0;
  }

  //reduce to upper triangular form
  //only non-zero elements are on diagonal (Aii) and last column (Ain)
  {
    int i;
    int i0;
    double alpha;
    //lower part (exluding last row)
    for(i = 1; i < n - 1; i++){
      i0 = i-1;
      alpha = - 1.0 / Aii[i0];
      Aii[i] += alpha;
      Ain[i] += alpha*Ain[i0];
      Q[i]   += alpha*Q[i0];
    }
    
    //upper part
    for(i = n-3; i >= 0; i--){
      i0 = i+1;
      alpha = -1.0 / Aii[i0];
      Ain[i] += alpha*Ain[i0];
      Q[i]   += alpha*Q[i0];
    }
    
    //last row (two bottom edges)
    i  = n - 1;
    i0 = 0;
    alpha = - 1.0 / Aii[i0];
    Aii[i] += alpha*Ain[i0];
    Ain[i] += alpha*Ain[i0];
    Q[i]   += alpha*Q[i0];

    i0 = n - 2;
    alpha = - 1.0 / Aii[i0];
    Aii[i] += alpha*Ain[i0];
    Ain[i] += alpha*Ain[i0];
    Q[i]   += alpha*Q[i0];
  }

  //back substitute for c_i coefficient
  {
    double cn = Q[n - 1] / Aii[n - 1];

    c[n - 1] = cn;
    for(int i = n - 2; i >= 0; i--){
      c[i] = (Q[i] - Ain[i]*cn) / Aii[i];
    }
  }

  //b_i & d_i coefficient
  {
    int inext;
    int iprev;
    for(int i = 0; i < n; i++){
      inext = MOD(i+1, n);
      iprev = MOD(i-1, n);
      d[i] = (c[inext] - c[i])/(3.0*dx);
      b[i] = (a[inext] - a[i])/dx - c[i]*dx - d[i]*SQ(dx);
    }
  }
}
