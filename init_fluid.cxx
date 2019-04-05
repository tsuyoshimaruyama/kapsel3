/*!
  \file init_fluid.cxx
  \brief Initialize fluid velocity fields
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */
#include "init_fluid.h"

void Init_zeta_k(double **zeta, double *uk_dc){
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(NPs, NZ_, u)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                int im = REALMODE_ARRAYINDEX(i, j, k);
                u[0][im] = 0.0e0;
                u[1][im] = 0.0e0;
                u[2][im] = 0.0e0;
            }
        }
    }
    U2u_k (u);
    Solenoidal_uk (u);
    U_k2zeta_k (u, zeta, uk_dc);

}


