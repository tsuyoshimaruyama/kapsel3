/*!
  \file alloc.c
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Memory allocation routines
 */

#include "alloc.h"

template<typename T>  T*  alloc_1d(int n1){
  T* i = (T*) std::malloc(n1*sizeof(T));
  alloc_error_check(i);
  return i;
}
template int*    alloc_1d(int n1);
template double* alloc_1d(int n1);

template<typename T>  void free_1d(T *i){
  std::free(i);
}
template void free_1d(int *i);
template void free_1d(double *i);


template<typename T>
T** initview_2d(int n1, int n2,  T* i){
  T** ii = (T **) alloc_1d<T*>(n1);
  alloc_error_check(ii);
  ii[0] = i;
  for(int j = 1; j < n1; j++) ii[j] = ii[j-1] + n2;
  return ii;
}
template int**    initview_2d(int n1, int n2, int* i);
template double** initview_2d(int n1, int n2, double* i);

template<typename T>
void freeview_2d(T **ii){
  free(ii);
}
template void freeview_2d(int **ii);
template void freeview_2d(double **ii);

template<typename T>
T** alloc_2d(int n1, int n2){
  T*  i  = (T *)  alloc_1d<T>(n1 * n2);
  alloc_error_check(i);
  T** ii = initview_2d(n1, n2, i);
  ii[0] = i;
  for(int j = 1; j < n1; j++) ii[j] = ii[j-1] + n2;
  return ii;
}
template int**    alloc_2d(int n1, int n2);
template double** alloc_2d(int n1, int n2);

template<typename T>
void free_2d(T **ii){
    free(ii[0]);
    freeview_2d(ii);
}
template void free_2d(int **ii);
template void free_2d(double **ii);

template<typename T>
T*** initview_3d(int n1, int n2, int n3, T* i){
  T ***iii = (T ***) alloc_1d<T**>(n1);
  alloc_error_check(iii);
  T **ii   = (T **)  alloc_1d<T*>(n1 * n2);
  alloc_error_check(ii);

  iii[0] = ii;
  for(int j = 1; j < n1; j++) iii[j] = iii[j-1] + n2;

  ii[0]  = i;
  for(int j = 1; j < n1*n2; j++) ii[j] = ii[j-1] + n3;
  return iii;
}
template int***    initview_3d(int n1, int n2, int n3, int* i);
template double*** initview_3d(int n1, int n2, int n3, double* i);

template<typename T>
void freeview_3d(T ***iii){
  free(iii[0]);
  free(iii);
}
template void freeview_3d(int ***i);
template void freeview_2d(double ***i);

template<typename T>
T ***alloc_3d(int n1, int n2, int n3){
  T*     i = (T *) alloc_1d<T>(n1 * n2 * n3);
  alloc_error_check(i);
  T*** iii = initview_3d(n1, n2, n3, i);
  return iii;
}
template int***    alloc_3d(int n1, int n2, int n3);
template double*** alloc_3d(int n1, int n2, int n3);

template<typename T>
void free_3d(T ***iii){
    free(iii[0][0]);
    freeview_3d(iii);
}
template void free_3d(int ***iii);
template void free_3d(double ***iii);


