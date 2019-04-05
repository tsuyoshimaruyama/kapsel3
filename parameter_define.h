/*!
  \file parameter_define.h
  \brief Define global system parameters for FFT routines
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */
#ifndef PARAMETER_DEFINE_H_
#define PARAMETER_DEFINE_H_

#define DIM 3
#define QDIM ((DIM*(DIM+1))/2-1)
#define SPACE 2 //REAL or SPECTRUM SPACE
/////// for SW_FFT
#define Ooura 0
#define RFFTW 1
#define MPI_RFFTW 2
#define IMKL_FFT 3

#if !defined (NDEBUG)
/* for Mersennetwister Parameter (MT_ARRAY = (521/32) + 1) */
#define MTW 32
#define MTP 521
#define MT_ARRAY (MTP/MTW + 2)
#endif
#define Cell_length 16
/* division number*/
#define DIV_N 2
#ifdef _MPI
#define ID(a,b) idtbl[(offset_x) + (a)][(offset_y) + (b)]
#define LJ_CUTOFF_REFID(a,b) lj_ref[(a)][(b)]
#define SEKIBUN_REFID(a,b) sekibun_ref[(a)][(b)]
#define TAG(a,b) tags[(1 + (a))][(1 + (b))]
#endif
#define REALMODE_ARRAYINDEX(i,j,k) (((i) * NPs[REAL][1] * NZ_) + ((j) * NZ_) + (k))
#define REALMODE_ARRAYINDEX_MESH(p,i,j,k) (((p) * NPs[REAL][0] * NPs[REAL][1] * NZ_) + ((i) * NPs[REAL][1] * NZ_) + ((j) * NZ_) + (k))
#define SPECTRUMMODE_ARRAYINDEX(i,j,k) (((i) * NPs[SPECTRUM][1] * NPs[SPECTRUM][2]) + ((j) * NPs[SPECTRUM][2]) + (k))
/* Boolean */
#define ANS_TRUE 1
#define ANS_FALSE 0

#endif
