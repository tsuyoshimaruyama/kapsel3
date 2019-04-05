#ifndef __FFT_BASE__
#define __FFT_BASE__
#ifdef _MPI
#include <mpi.h>
#endif
#include "variable.h"
#include "alloc.h"
#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <emmintrin.h> //sseëgÇ›çûÇ›ä÷êîåƒÇ—èoÇµÇÃÇΩÇﬂ

#ifdef _MPI
#include "operate_mpi.h"
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#if (defined (_OPENMP) || defined(_MPI))
#include "mkl.h"
#include <mkl_dfti.h>
#endif

#ifdef _FFT_IMKL
#include <mkl_dfti.h>
#elif _FFT_FFTW
#include <fftw3.h>
#endif

extern void rdft3d(int n1, int n2, int n3, int sign, double ***a, double *t, int *ip, double *w);
extern void rdft3dsort(int, int, int, int, double ***);
void Init_fft_omp_mpi_imkl(void);
void Dfti_finalize(void);
void A2a_k_1D(double *a);
void A_k2a_1D(double *a);
void A2a_k_nD(double **a, int dim);
void A_k2a_nD(double **a, int dim);
#endif
