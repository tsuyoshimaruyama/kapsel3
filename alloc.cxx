/*!
  \file alloc.c
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Memory allocation routines
 */

#include "alloc.h"

template<typename T>  T*  alloc_1d(int n1){
  T* p = (T*) std::malloc(n1*sizeof(T));
  alloc_error_check(p);
  return p;
}
template<typename T>  void free_1d(T* p){
  std::free(p);
}

template<typename T> T** initview_2d(int n1, int n2,  T * const p){
  T** pp = (T **) alloc_1d<T*>(n1);
  alloc_error_check(pp);
  pp[0] = p;
  for(int j = 1; j < n1; j++) pp[j] = pp[j-1] + n2;
  return pp;
}
template<typename T> void freeview_2d(T **pp){
  free_1d<T*>(pp);
}

template<typename T> T** alloc_2d(int n1, int n2){
  T*  p  = (T *)  alloc_1d<T>(n1 * n2);
  alloc_error_check(p);
  T** pp = initview_2d<T>(n1, n2, p);
  return pp;
}
template<typename T> void free_2d(T **pp){
  free_1d<T>(pp[0]);
  freeview_2d<T>(pp);
}

template<typename T> T*** initview_3d(int n1, int n2, int n3, T* const p){
  T ***ppp = (T ***) alloc_1d<T**>(n1);
  alloc_error_check(ppp);
  T **pp   = (T **)  alloc_1d<T*>(n1 * n2);
  alloc_error_check(pp);

  ppp[0] = pp;
  for(int j = 1; j < n1; j++) ppp[j] = ppp[j-1] + n2;

  pp[0]  = p;
  for(int j = 1; j < n1*n2; j++) pp[j] = pp[j-1] + n3;
  return ppp;
}
template<typename T> void freeview_3d(T ***ppp){
  free_1d<T*>(ppp[0]);
  free_1d<T**>(ppp);
}

template<typename T> T ***alloc_3d(int n1, int n2, int n3){
  T*     p = (T *) alloc_1d<T>(n1 * n2 * n3);
  alloc_error_check(p);
  T*** ppp = initview_3d<T>(n1, n2, n3, p);
  return ppp;
}
template<typename T> void free_3d(T ***ppp){
  free_1d<T>(ppp[0][0]);
  freeview_3d<T>(ppp);
}

int*       alloc_1d_int(int n1){return alloc_1d<int>(n1);}
int**      alloc_2d_int(int n1, int n2){return alloc_2d<int>(n1, n2);}
int***     alloc_3d_int(int n1, int n2, int n3){return alloc_3d<int>(n1, n2, n3);}
int**      initview_2d_int(int n1, int n2, int* const i){return initview_2d<int>(n1, n2, i);}
int***     initview_3d_int(int n1, int n2, int n3, int* const i){return initview_3d<int>(n1, n2, n3, i);}

void       free_1d_int(int *i){return free_1d<int>(i);}
void       free_2d_int(int **ii){return free_2d<int>(ii);}
void       free_3d_int(int ***iii){return free_3d<int>(iii);}
void       freeview_2d_int(int **ii){return freeview_2d<int>(ii);}
void       freeview_3d_int(int ***iii){return freeview_3d<int>(iii);}


double*    alloc_1d_double(int n1){return alloc_1d<double>(n1);}
double**   alloc_2d_double(int n1, int n2){return alloc_2d<double>(n1, n2);}
double***  alloc_3d_double(int n1, int n2, int n3){return alloc_3d<double>(n1, n2, n3);}
double**   initview_2d_double(int n1, int n2, double* const d){return initview_2d<double>(n1, n2, d);}
double***  initview_3d_double(int n1, int n2, int n3, double* const d){return initview_3d<double>(n1, n2, n3, d);}

void       free_1d_double(double *d){return free_1d<double>(d);}
void       free_2d_double(double **dd){return free_2d<double>(dd);}
void       free_3d_double(double ***ddd){return free_3d<double>(ddd);}
void       freeview_2d_double(double **dd){return freeview_2d<double>(dd);}
void       freeview_3d_double(double ***ddd){return freeview_3d<double>(ddd);}





