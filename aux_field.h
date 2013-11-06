#ifndef AUX_FIELD_H
#define AUX_FIELD_H

#include "input.h"


enum FIELD_SPACE{R_SPACE, K_SPACE};

inline void Copy_v1_Primitive(double *ucp, double const *u, 
                              const FIELD_SPACE &flag){

  int im;
  int NZMAX = (flag == R_SPACE ? NZ : NZ_);
#pragma omp parallel for schedule(dynamic, 1) private(im)
  for(int i = 0; i < NX; i++){
    for(int j = 0; j < NY; j++){
      for(int k = 0; k < NZMAX; k++){
        im = (i * NY * NZ_) + (j * NZ_) + k;
        ucp[im] = u[im];
      }
    }
  }
}

inline void Copy_v2_Primitive(double **ucp, double const* const* u, 
                              const FIELD_SPACE &flag){

  int im;
  int NZMAX = (flag == R_SPACE ? NZ : NZ_);
#pragma omp parallel for schedule(dynamic, 1) private(im)
  for(int i = 0; i < NX; i++){
    for(int j = 0; j < NY; j++){
      for(int k = 0; k < NZMAX; k++){
        im = (i * NY * NZ_) + (j * NZ_) + k;
        ucp[0][im] = u[0][im];
        ucp[1][im] = u[1][im];
      }
    }
  }
}

inline void Copy_v3_Primitive(double **ucp, double const* const* u, 
                              const FIELD_SPACE &flag){

  int im;
  int NZMAX = (flag == R_SPACE ? NZ : NZ_);
#pragma omp parallel for schedule(dynamic, 1) private(im)
  for(int i = 0; i < NX; i++){
    for(int j = 0; j < NY; j++){
      for(int k = 0; k < NZMAX; k++){
        im = (i * NY * NZ_) + (j * NZ_) + k;
        ucp[0][im] = u[0][im];
        ucp[1][im] = u[1][im];
        ucp[2][im] = u[2][im];
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
  Copy_v2_Primitive(ucp, u, R_SPACE);
}
inline void Copy_v2_k(double **ucp, double const* const* u){
  Copy_v2_Primitive(ucp, u, K_SPACE);
}

inline void Copy_v3(double **ucp, double const* const* u){
  Copy_v3_Primitive(ucp, u, R_SPACE);
}
inline void Copy_v3_k(double **ucp, double const* const* u){
  Copy_v3_Primitive(ucp, u, K_SPACE);
}

#endif
