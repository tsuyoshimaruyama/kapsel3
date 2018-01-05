/*!
  \file alloc.h
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Memory allocation routines (header file)
 */

#ifndef ALLOC_H
#define ALLOC_H

#include <cstdlib>
#include <cstdio>
#include <mm_malloc.h>

#define alloc_error_check(p) {			\
    if ((p) == NULL) { \
        fprintf(stderr, "Allocation Failure!\n"); \
        exit(1); \
    } \
}

/*!
  \brief Allocate memory for 1d array
 */
template<typename T> T* alloc_1d(int n1);

/*!
  \brief Free memory for 1d array
 */
template<typename T> void free_1d(T *p);

/*!
  \brief Create 2d view of 1d array referenced by i
  \details Allocates array ii of size n1, each element of ii give view of rows of i
  ii[I][J] = i[I*n2 + J]
 */
template<typename T> T** initview_2d(int n1, int n2,  T * const p);

/*!
  \brief Destroy 2d view
 */
template<typename T> void freeview_2d(T **pp);
/*!
  \brief Allocate memory for 2d array
  \details Allocates linear array of size n1*n2, position of element
  (I,J) is im = I*n2 + J
 */
template<typename T>  T** alloc_2d(int n1, int n2);
/*!
  \brief Free memory for 2d array
 */
template<typename T> void free_2d(T **ii);

/*!
  \brief Create 3d view of 1d array referenced by i
  \details Allocates array iii of size n1, each iii[.] points to 2D array ii such that
  iii[I][J][K] = i[I*n2*n3 + J*n3 + K]
 */
template<typename T>  T*** initview_3d(int n1, int n2, int n3, T * const p);
/*!
  \brief Destroy 3d view
 */
template<typename T> void  freeview_3d(T ***ppp);
/*!
  \brief Allocate memory for 3d array
  \details Allocates linear array of size n1*n2*n3, position of element
  (I,J,K) is im = I*n2*n3 + J*n3 + K
 */
template<typename T>  T*** alloc_3d(int n1, int n2, int n3);
/*!
  \brief Free memory for 3d array
 */
template<typename T> void  free_3d(T ***ii);

int*      alloc_1d_int(int n1);
int**     alloc_2d_int(int n1, int n2);
int***    alloc_3d_int(int n1, int n2, int n3);
int**     initview_2d_int(int n1, int n2, int * const i);
int***    initview_3d_int(int n1, int n2, int n3, int * const i);

void      free_1d_int(int *i);
void      free_2d_int(int **ii);
void      free_3d_int(int ***iii);
void      freeview_2d_int(int **ii);
void      freeview_3d_int(int ***iii);

double*   alloc_1d_double(int n1);
double**  alloc_2d_double(int n1, int n2);
double*** alloc_3d_double(int n1, int n2, int n3);
double**  initview_2d_double(int n1, int n2, double * const d);
double*** initview_3d_double(int n1, int n2, int n3, double * const d);

void      free_1d_double(double *d);
void      free_2d_double(double **dd);
void      free_3d_double(double ***dddd);
void      freeview_2d_double(double **dd);
void      freeview_3d_double(double ***dddd);


#endif
