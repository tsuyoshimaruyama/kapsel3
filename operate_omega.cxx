/*!
  \file operate_omega.cxx
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief High-level routines to handle basic operations on velocity field
  \see \ref page_design_fsolver section of manual for further details
 */

#include "operate_omega.h"

const int DIM2 = 2*DIM;

void U2advection_k(double **u, double **advection){
    double *u_advection[] = { u[0], u[1], u[2], advection[0], advection[1] };
    double us[3][DIM];
    double ks[3][DIM];
    int im;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, u, advection) private(us, im)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
                us[0][0] = u[0][im];
                us[0][1] = u[1][im];
                us[0][2] = u[2][im];
                u[0][im] = us[0][1] * us[0][2];
                u[1][im] = us[0][2] * us[0][0];
                u[2][im] = us[0][0] * us[0][1];
                advection[0][im] = (us[0][1] - us[0][2]) * (us[0][1] + us[0][2]);
                advection[1][im] = (us[0][2] - us[0][0]) * (us[0][2] + us[0][0]);
            }
        }
    }
    U2u_k (u_advection, (DIM + 2) );
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, KX_int, KY_int, KZ_int, WAVE_X, WAVE_Y, WAVE_Z, u, advection) private(ks, us, im)
#endif
    for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
        for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
            for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
                im = SPECTRUMMODE_ARRAYINDEX(i, j, k);
                ks[0][0] = KX_int[im];
                ks[0][1] = KY_int[im];
                ks[0][2] = KZ_int[im];
                ks[1][0] = ks[0][0] * WAVE_X;
                ks[1][1] = ks[0][1] * WAVE_Y;
                ks[1][2] = ks[0][2] * WAVE_Z;
                us[1][0] = u[0][im];
                us[1][1] = u[1][im];
                us[1][2] = u[2][im];
                us[2][0] = advection[0][im];
                us[2][1] = advection[1][im];
                us[2][2] = -us[2][1] - us[2][0];
                ks[2][0] = ks[1][0] * ks[1][1];
                ks[2][1] = ks[1][1] * ks[1][2];
                ks[2][2] = ks[1][2] * ks[1][0];
                u[0][im] = - (ks[2][1] * us[2][0] + (ks[1][2] - ks[1][1])
                              * (ks[1][2] + ks[1][1]) * us[1][0] + ks[2][2]
                              * us[1][2] - ks[2][0] * us[1][1]);
                u[1][im] = - (ks[2][2] * us[2][1] + (ks[1][0] - ks[1][2])
                              * (ks[1][0] + ks[1][2]) * us[1][1] + ks[2][0]
                              * us[1][0] - ks[2][1] * us[1][2]);
                u[2][im] = - (ks[2][0] * us[2][2] + (ks[1][1] - ks[1][0])
                              * (ks[1][1] + ks[1][0]) * us[1][2] + ks[2][1]
                              * us[1][1] - ks[2][2] * us[1][0]);
                if (ks[0][0] != 0) {
                    advection[0][im] = u[1][im];
                    advection[1][im] = u[2][im];
                } else if (ks[0][1] != 0) {
                    advection[0][im] = u[2][im];
                    advection[1][im] = u[0][im];
                } else {
                    advection[0][im] = u[0][im];
                    advection[1][im] = u[1][im];
                }
            }
        }
    }
    if (procid == 0) {
        assert ( (advection[0][0] == 0.0) && (advection[1][0] == 0.0) );
    }
}

void Zeta_k2advection_k(double **zeta, double uk_dc[DIM], double **advection){
    for (int d = 0; d < (DIM - 1); d++) {
        Truncate_two_third_rule_ooura(zeta[d]);
    }
    Zeta_k2u_k (zeta, uk_dc, u);
    U_k2u (u);
    U2advection_k (u, advection);
}

void Zeta_k2advection_k_OBL(double **zeta, double uk_dc[DIM], double **advection){
    //Truncate_vector_two_third_rule(zeta, DIM-1);
    //Zeta_k2u(zeta, uk_dc, u);
    for(int d=0; d<(DIM-1); d++){
        Truncate_two_third_rule_ooura(zeta[d]);
    }

    Zeta_k2u_cpuky(zeta, uk_dc, u, work_v1);//contra
    Zeta_k2omega_OBL(zeta, work_v3);//contra

    int im;
    double u1;
    double u2;
    double u3;
    double w1;
    double w2;
    double w3;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, u, work_v3) private(u1,u2,u3,w1,w2,w3,im)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);

                u1 = u[0][im];
                u2 = u[1][im];
                u3 = u[2][im];

                w1 = work_v3[0][im];
                w2 = work_v3[1][im];
                w3 = work_v3[2][im];

                u[0][im] = u2*w3 - u3*w2;//co
                u[1][im] = u3*w1 - u1*w3;//co
                u[2][im] = u1*w2 - u2*w1;//co
            }
        }
    }

    for(int d=0;d<DIM;d++){
        A2a_k(u[d]);
    }

#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, u, work_v1, Shear_rate_eff, degree_oblique) private(im)
#endif
    for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
        for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
            for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
                im = SPECTRUMMODE_ARRAYINDEX(i, j, k);
                u[0][im] -= 2.*Shear_rate_eff*work_v1[im];//co
                u[1][im] -= 2.*Shear_rate_eff*degree_oblique*work_v1[im];//co
            }
        }
    }
    
    U_k2rotation_k(u);//contra
    Omega_k2zeta_k_OBL(u, advection);//contra
}

void Add_zeta_viscous_term(double ** zeta, double **f, const Index_range &ijk_range){
    int im;
    Index_range renge;

    if (Range_check (&ijk_range, &renge) ) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (renge, NPs, f, NU, K2, zeta) private(im)
#endif
        for (int i = renge.istart; i <= renge.iend; i++) {
            for (int j = renge.jstart; j <= renge.jend; j++) {
                for (int k = renge.kstart; k <= renge.kend; k++) {
                    im = SPECTRUMMODE_ARRAYINDEX(i, j, k);
                    f[0][im] += - (NU * K2[im] * zeta[0][im]);
                    f[1][im] += - (NU * K2[im] * zeta[1][im]);
                }
            }
        }
    }
}

void Solenoidal_uk(double **u){
    double ks[DIM];
    double dmy;
    int im;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, KX_int, KY_int, KZ_int, WAVE_X, WAVE_Y, WAVE_Z, IK2, u) private(ks, dmy, im)
#endif
    for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
        for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
            for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
                im = SPECTRUMMODE_ARRAYINDEX(i ,j, k);
                ks[0] = KX_int[im] * WAVE_X;
                ks[1] = KY_int[im] * WAVE_Y;
                ks[2] = KZ_int[im] * WAVE_Z;
                dmy = IK2[im] * (u[0][im] * ks[0] + u[1][im] * ks[1] + u[2][im] * ks[2]);
                u[0][im] -= dmy * ks[0];
                u[1][im] -= dmy * ks[1];
                u[2][im] -= dmy * ks[2];
            }
        }
    }
}

void Solenoidal_uk_OBL(double **u){
    double kx;
    double ky;
    double kz;
    double kx_contra;
    double ky_contra;
    double dmy;
    int im;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, WAVE_X, WAVE_Y, WAVE_Z, KX_int, KY_int, KZ_int, IK2, u, degree_oblique) private(kx,ky,kz,kx_contra,ky_contra,dmy,im) 
#endif
    for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
        for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
            for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
                im = SPECTRUMMODE_ARRAYINDEX(i, j, k);

                kx =KX_int[im] * WAVE_X;
                ky =KY_int[im] * WAVE_Y;
                kz =KZ_int[im] * WAVE_Z;

                kx_contra = (1. + degree_oblique*degree_oblique)*KX_int[im]*WAVE_X - degree_oblique*KY_int[im]*WAVE_Y;
                ky_contra = -degree_oblique*KX_int[im]*WAVE_X + KY_int[im]*WAVE_Y;
                dmy = IK2[im] * (u[0][im] * kx + u[1][im] * ky + u[2][im] * kz);

                u[0][im] -= dmy * kx_contra;
                u[1][im] -= dmy * ky_contra;
                u[2][im] -= dmy * kz;
            }
        }
    }
}

