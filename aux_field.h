#ifndef AUX_FIELD_H
#define AUX_FIELD_H

#include "input.h"

/* change: spatial decomposition for MPI */

enum FIELD_SPACE{R_SPACE, K_SPACE};

inline void Copy_v1_Primitive(double *ucp, double const *u, const FIELD_SPACE &flag){

    int im;

    if(flag == R_SPACE){
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, ucp, u) private(im)
#endif
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    im = REALMODE_ARRAYINDEX(i, j, k);
                    ucp[im] = u[im];
                }
            }
        }
    }else if(flag == K_SPACE){
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, ucp, u) private(im)
#endif
        for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
            for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
                for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
                    im = SPECTRUMMODE_ARRAYINDEX(i, j, k);
                    ucp[im] = u[im];
                }
            }
        }
    }
}

inline void Copy_v2_Primitive(double **ucp, double const* const* u, const FIELD_SPACE &flag){

    int im;

	if(flag == R_SPACE){

#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, ucp, u) private(im)
#endif
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    im = REALMODE_ARRAYINDEX(i, j, k);
                    ucp[0][im] = u[0][im];
                    ucp[1][im] = u[1][im];
                }
            }
        }
    }else if(flag == K_SPACE){

#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, ucp, u) private(im)
#endif
        for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
            for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
                for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
                    im = SPECTRUMMODE_ARRAYINDEX(i, j, k);
                    ucp[0][im] = u[0][im];
                    ucp[1][im] = u[1][im];
                }
            }
        }
    }
}

inline void Copy_v3_Primitive(double **ucp, double const* const* u, const FIELD_SPACE &flag){

    int im;

	if(flag == R_SPACE){

#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, ucp, u) private(im)
#endif
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    im = REALMODE_ARRAYINDEX(i, j, k);
                    ucp[0][im] = u[0][im];
                    ucp[1][im] = u[1][im];
                    ucp[2][im] = u[2][im];
                }
            }
        }
    }else if(flag == K_SPACE){

#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, ucp, u) private(im)
#endif
        for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
            for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
                for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
                    im = SPECTRUMMODE_ARRAYINDEX(i, j, k);
                    ucp[0][im] = u[0][im];
                    ucp[1][im] = u[1][im];
                    ucp[2][im] = u[2][im];
                }
            }
        }
    }
}

inline void Copy_v5_Primitive(double **ucp, double const* const* u, const FIELD_SPACE &flag){

    int im;

	if(flag == R_SPACE){

#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, ucp, u) private(im)
#endif
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    im = REALMODE_ARRAYINDEX(i, j, k);
                    ucp[0][im] = u[0][im];
                    ucp[1][im] = u[1][im];
                    ucp[2][im] = u[2][im];
                    ucp[3][im] = u[3][im];
                    ucp[4][im] = u[4][im];
                }
            }
        }
    }else if(flag == K_SPACE){

#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, ucp, u) private(im)
#endif
        for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
            for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
                for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
                    im = SPECTRUMMODE_ARRAYINDEX(i, j, k);
                    ucp[0][im] = u[0][im];
                    ucp[1][im] = u[1][im];
                    ucp[2][im] = u[2][im];
                    ucp[3][im] = u[3][im];
                    ucp[4][im] = u[4][im];
                }
            }
        }
    }
}

inline void Copy_v1(double *ucp, double const* u){
    Copy_v1_Primitive(ucp, u, R_SPACE);
}
inline void Copy_v1_k(double *ucp, double const* u){
    Copy_v1_Primitive(ucp, u, K_SPACE);
}

inline void Copy_v2(double **ucp, double const* const* u){
    Copy_v1_Primitive(ucp[0], u[0], R_SPACE);
    Copy_v1_Primitive(ucp[1], u[1], R_SPACE);
}
inline void Copy_v2_k(double **ucp, double const* const* u){
    Copy_v1_Primitive(ucp[0], u[0], K_SPACE);
    Copy_v1_Primitive(ucp[1], u[1], K_SPACE);
}

inline void Copy_v3(double **ucp, double const* const* u){
    Copy_v1_Primitive(ucp[0], u[0], R_SPACE);
    Copy_v1_Primitive(ucp[1], u[1], R_SPACE);
    Copy_v1_Primitive(ucp[2], u[2], R_SPACE);
}
inline void Copy_v3_k(double **ucp, double const* const* u){
    Copy_v1_Primitive(ucp[0], u[0], K_SPACE);
    Copy_v1_Primitive(ucp[1], u[1], K_SPACE);
    Copy_v1_Primitive(ucp[2], u[2], K_SPACE);
}

inline void Copy_v5(double **ucp, double const* const* u){
    Copy_v1_Primitive(ucp[0], u[0], R_SPACE);
    Copy_v1_Primitive(ucp[1], u[1], R_SPACE);
    Copy_v1_Primitive(ucp[2], u[2], R_SPACE);
    Copy_v1_Primitive(ucp[3], u[3], R_SPACE);
    Copy_v1_Primitive(ucp[4], u[4], R_SPACE);
}
inline void Copy_v5_k(double **ucp, double const* const* u){
    Copy_v1_Primitive(ucp[0], u[0], K_SPACE);
    Copy_v1_Primitive(ucp[1], u[1], K_SPACE);
    Copy_v1_Primitive(ucp[2], u[2], K_SPACE);
    Copy_v1_Primitive(ucp[3], u[3], K_SPACE);
    Copy_v1_Primitive(ucp[4], u[4], K_SPACE);
}

inline void Swap_mem(double* &unew, double* &uold){
    double* temp = unew;
    unew = uold;
    uold = temp;
}

inline void Swap_mem(double** &unew, double ** &uold){
    double** temp = unew;
    unew = uold;
    uold = temp;
}

#endif
