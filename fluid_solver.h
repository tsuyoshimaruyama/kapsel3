/*!
  \file fluid_solver.h
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Solver for Navier-Stokes equations (header file)
 */
#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#ifdef _MPI
#include <mpi.h>
#endif
#include <math.h>
#include "variable.h"
#include "input.h"
#include "f_particle.h"
#include "operate_omega.h"
#include "solute_rhs.h"
#include "fluct.h"
#include "memory_model.h"

extern double *Pressure;
extern double **Shear_force;
extern double **Shear_force_k;
extern double **f_ns0;
extern double **f_ns1;


//Documentation for inline functions defined in fluid_solver.cxx

/*!
  \fn void Field_solver_Euler(const int &dim, double **zeta_k, const CTime &jikan, double **rhs, const Index_range &ijk_range)
  \brief solver euler
 */

/*!
  \fn void Rhs_NS(double **zeta, double uk_dc[DIM], double **rhs, const Index_range *ijk_range, const int &n_ijk_range, Particle *p)
  \brief rhs of ns
 */

/*!
  \fn void Rhs_solvent(double **u_solvent, double **rhs_solvent, double rhs_uk_dc[DIM])
  \brief rhs
*/

/*!
  \fn void Rhs_NS_solute(Particle *p, double **zeta_k, double uk_dc[DIM], double **u, double **concentration_k, double **rhs_ns, double **rhs_solute, const Index_range *ijk_range, const int &n_ijk_range, double **solute_flux, double **grad_potential, double **surface_normal, double rhs_uk_dc[DIM])
  \brief rhs + ns
 */

/*!
  \fn Add_constant_field_k(double **grad_potential_k, double e_ext[DIM], const CTime &jikan)
  \brief field
 */

/*!
  \fn Rhs_NS_Nernst_Planck(Particle *p, double **zeta, double uk_dc[DIM], double **u, double **concentration_k, double **rhs_ns, double **rhs_solute, const Index_range *ijk_range, const int &n_ijk_range, double **solute_flux, double **grad_potential, double **surface_normal, double rhs_uk_dc[DIM], const CTime &jikan)
  \brief ns + np
 */

// End of inline documentation

/*!
  \brief Allocate workspace variables 
 */

//function is integrated to Mem_alloc_var()
void Mem_alloc_NS_solver(void);

/*!
  \brief Solve Navier-Stokes equation to update reduced vorticity field
  \details \f[
  \ft{\vec{\zeta}} \longrightarrow \ft{\vec{\zeta}} + \left(e^{-\nu (2\pi k)^2 h} - 1\right)\left[\ft{\vec{\zeta}} + \frac{\ft{\vec{\Omega}}^*}{\nu(2\pi k)^2}\right]
  \f]
  Analytic solution is valid in the absence of solute terms, with no shear.
  \param[in,out] zeta reduced vorticity field
  \param[in] jikan time data
  \param[in] uk_dc zero-wavenumber Fourier transform of the velocity field
  \param[in] ijk_range field iterator parameters for update
  \param[in] n_ijk_range field iterator parameters for update
  \param[in] p particle data (unused)
  \see \ref page_design_fsolver section of manual for further details.
  \todo Specify the difference between ijk_range and n_ijk_range
 */
void NS_solver_slavedEuler(double **zeta, const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);

/*! 
  \brief Solve Navier-Stokes equation under zig-zag shear flow to update
  reduced vorticity field
  \param[in,out] zeta reduced vorticity field (reciprocal space)
  \param[in] jikan time data
  \param[in] uk_dc zero-wavenumber Fourier transform of the velocity
  field
  \param[in] ijk_range field iterator parameters for update
  \param[in] p particle data (unused)
 */
void NS_solver_slavedEuler_Shear_PBC(double **zeta, const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p, double **force);

/*!
  \brief Solve Navier-Stokes equations with Lees-Edwards PBC under
  shear flow to update reduced vorticity field
  \details \f[
  \ft{\zeta}^\alpha \longrightarrow \ft{\zeta}^{\alpha} + 
  \left(e^{-\nu(2\pi k)^2 h} - 1\right)\left[\ft{\zeta}^\alpha +
  \frac{\ft{\Omega}^{\alpha*}}{\nu (2\pi k)^2}\right]
  \f]
  \param[in,out] zeta contravariant reduced vorticity field
  (reciprocal space)
  \param[in] jikan time data
  \param[in] uk_dc zero-wavenumber Fourier transform of the velocity
  field (contravariant)
  \param[in] ijk_range field iterator parameters for update
  \param[in] n_ijk_range field iterator parameters for update
  \param[in] p particle data (unused)
  \see \ref page_design_fsolverOBL section of manual for further details
 */
void NS_solver_slavedEuler_Shear_OBL(double **zeta, const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p, double **force);

//Type of Rheology 
void Mean_shear_sustaining_yforce_PBC(double **u, double **force,const CTime &jikan);

void Mean_shear_sustaining_kforce_PBC(double **u, double **force,const CTime &jikan);

//void Mean_shear_sustaining_kforce_PBC(Value u[DIM], Value force[DIM],const CTime &jikan);
void Mean_shear_sustaining_force_PBC_OBL(double **u);

void Ion_diffusion_solver_Euler(double **zeta
				,const CTime &jikan
				,double uk_dc[DIM]
				,double **concentration_k
				,Particle *p
				,const Index_range *ijk_range
				,const int &n_ijk_range
				);
/*!
  \brief Solves Navier-Stokes and Poisson-Nernst Plank equation to update the fluid velocity and the electrolyte concentration
  \details Calls Rhs_NS_Nernst_Plank to compute the various terms appearing in the right hand side of the relevant equations, then proceeds to call Field_solver_Euler to update the reduced vorticity field \f$\vec{\zeta}\f$ and the solute concentration \f$C_{\alpha}\f$.
  \param[in,out] zeta reduced vorticity field (reciprocal space)
  \param[in] jikan time data
  \param[in,out] uk_dc zero-wavenumber Fourier transform of velocity field
  \param[in,out] concentration_k solute concentration field (reciprocal space)
  \param[in] p particle data
  \param[in] ijk_range field iterator data
  \param[in] n_ijk_range field iterator data
  \note The solute source term appearing in the rhs of the Navier-Stokes equation has not been included. This section of the code should probably be removed for clarity...
  \see \ref page_design_ssolver section of manual for furter details.
 */
void NSsolute_solver_Euler(double **zeta
			   ,const CTime &jikan
			   ,double uk_dc[DIM]
			   ,double **concentration_k
			   ,Particle *p
			   ,const Index_range *ijk_range
			   ,const int &n_ijk_range
			   );

inline void Calc_Reynolds_shear_stress(double const* const* u, double &reynolds_shear_stress){

    int im;
	static const double ivolume = Ivolume * POW3(DX);

    double stress_yx = 0.0;
    double stress_yx_local = 0.0;

//#ifdef _OPENMP
//#pragma omp parallel for default(none) reduction(+:stress_yx_local) shared (NPs, NY, NZ_, DX, Shear_rate_eff, u, reynolds_shear_stress) private (im)
//#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            double delta_y = (double)((j + PREV_NPs[REAL][1]) - NY/2)*DX;
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
                //only consider fluctuating velocities
#ifdef _MPI
                stress_yx_local += (u[0][im] - Shear_rate_eff * delta_y)*(u[1][im]); 
#else
				stress_yx += (u[0][im] - Shear_rate_eff * delta_y)*(u[1][im]);
#endif

            }
        }
    }
#ifdef _MPI
    MPI_Allreduce(&stress_yx_local, &stress_yx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    reynolds_shear_stress = -ivolume*RHO*stress_yx;
}

inline void Calc_shear_stress(const CTime &jikan, Particle *p, double *phi //,double **force
                                 ,double **Shear_force, double stress[DIM][DIM]){
    static double factor = Ivolume / jikan.dt_fluid;
    static const double ivolume = Ivolume * POW3 (DX);
    const double dmy = ivolume / jikan.dt_fluid; //Select this
    double mean_force_x = 0.0;
    double stress_yx = 0.0;
    int im;
    double dmy_rhop;
#ifdef _MPI
    double *dmy_mpi;
#endif
    Reset_phi (phi);
    Make_rho_field (phi, p);
#ifdef _MPI
    dmy_mpi = alloc_1d_double (procs);
#endif
    if (Shear_AC) {
        Inertia_stress = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:Inertia_stress) shared (NPs, NZ_, PREV_NPs, phi, ucp, DX, rhop) private(im, dmy_rhop)
#endif
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    im = REALMODE_ARRAYINDEX(i, j, k);
                    dmy_rhop = phi[im] * ucp[0][im];
                    Inertia_stress += ((double)(j + PREV_NPs[REAL][1]) * DX) * (dmy_rhop - rhop[im]);
                    rhop[im] = dmy_rhop;
                }
            }
        }
        Inertia_stress *= factor;
    }
#ifdef _MPI
    ierr = MPI_Allgather (&Inertia_stress, 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_ERRORCHECK (ierr);
    Inertia_stress = 0.0;
    for (int n = 0; n < procs; n++) {
        Inertia_stress += dmy_mpi[n];
    }
#endif
    mean_force_x = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:mean_force_x) shared (NPs, NZ_, Shear_force, phi) private(im)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
                Shear_force[0][im] *= phi[im];
                mean_force_x += Shear_force[0][im];
            }
        }
    }
#ifdef _MPI
    ierr = MPI_Allgather (&mean_force_x, 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_ERRORCHECK (ierr);
    mean_force_x = 0.0;
    for (int n = 0; n < procs; n++) {
        mean_force_x += dmy_mpi[n];
    }
#endif
    mean_force_x *= ivolume;
    stress_yx = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:stress_yx) shared (NPs, NZ_, PREV_NPs, DX, Shear_force, mean_force_x) private(im)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
                stress_yx += ((double)(j + PREV_NPs[REAL][1]) * DX) * (Shear_force[0][im] - mean_force_x);
            }
        }
    }
#ifdef _MPI
    ierr = MPI_Allgather (&stress_yx, 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_ERRORCHECK (ierr);
    stress_yx = 0.0;
    for (int n = 0; n < procs; n++) {
        stress_yx += dmy_mpi[n];
    }
#endif
    stress[1][0] = stress_yx * dmy;
#ifdef _MPI
    free_1d_double (dmy_mpi);
#endif
}

inline void Calc_hydro_stress(const CTime &jikan, const Particle *p, double *phi, double *force, double stress[DIM][DIM]){

    int im;
    double stress_yx = 0.0;
    double stress_yx_local = 0.0;

    static const double ivolume = Ivolume * POW3(DX);
//#ifdef _OPENMP
//#pragma omp parallel for default(none) reduction(+:stress_yx_local) shared (NPs, NZ_, jikan, p, phi, force, stress) private(im)
//#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
            /*
                int jj = (int)fmod((double)targetJ + j + NY, NY); 

                double y = (double)j*DX;
                double fx = force[i][jj][k] - mean_force_x;
                stress_yx += y * fx;
             */
#ifdef _MPI
                stress_yx_local += force[im];
#else
				stress_yx += force[im];
#endif
            }
        }
    }
#ifdef _MPI
    MPI_Allreduce(&stress_yx_local, &stress_yx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    static const double dmy = (double)ivolume/jikan.dt_fluid; //Select this 
    stress[1][0]=-stress_yx*dmy;
}

/*!
  \brief Reset phi field when oblique degree is equal to one
 */
inline void Reset_phi_OBL(double *phi, double const* work_v1){
    int im, im_obl;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NX, NY, NZ_, phi, work_v1, PREV_NPs) private(im, im_obl)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            double sign = j - NPs[REAL][1]/2;
            if (!(sign == 0)) {
                sign = sign/fabs(sign);
            }
            int i_oblique = (int)(sign*((j + PREV_NPs[REAL][1]) - NY/2))*sign;
            i_oblique     = (int) fmod(i + i_oblique + 4.*NX, NX);
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
                im_obl = REALMODE_ARRAYINDEX((i_oblique), j, k);

                //Warning: reset grid points AND oblique basis vectors
                phi[im_obl] = work_v1[im];
            }
        }
    }
}

/*!
  \brief Reset velocity field when oblique degree is equal to one
  \details Reset 1-> 0 (flag == 0), otherwise reset 0->1
 */

inline void Reset_U_OBL(double **u, double const* const* ucp, const int &flg){
    const double delta = (flg == 0 ? 1.0 : -1.0);
    int im, im_obl;
	int num_mesh_transfer = NX * NY * NZ_ / yprocs;

    for(int i = 0; i < DIM; i++){
		for (int j = 0; j < num_mesh_transfer; j++) {
			work_v3[i][j] = 0.0;
        }
	}

    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {

            double sign = (j + PREV_NPs[REAL][1]) - NY/2;
            if (!(sign == 0)) {
                sign = sign/fabs(sign);
            }
            int i_oblique = (int)(sign*((j + PREV_NPs[REAL][1]) - NY/2))*sign;
            //i_oblique      = (int) fmod(i + delta*i_oblique + 4.*NX, NX);
            int i_real = i + PREV_NPs[REAL][0];
            i_oblique = (int) fmod(i_real + delta*i_oblique + 4.*NX, NX);
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
                im_obl = REALMODE_ARRAYINDEX((i_oblique), j, k);

#ifdef _MPI
                work_v3[0][im_obl] = ucp[0][im] + delta*ucp[1][im];
                work_v3[1][im_obl] = ucp[1][im];
                work_v3[2][im_obl] = ucp[2][im];
#else
				//Warning: reset grid points AND oblique basis vectors
				u[0][im_obl] = ucp[0][im] + delta*ucp[1][im];
				u[1][im_obl] = ucp[1][im];
				u[2][im_obl] = ucp[2][im];
#endif
            }
        }
    }

#ifdef _MPI
	MPI_Allreduce(work_v3[0], work_v4[0], num_mesh_transfer, MPI_DOUBLE, MPI_SUM, OWN_X_COMM);
	MPI_Allreduce(work_v3[1], work_v4[1], num_mesh_transfer, MPI_DOUBLE, MPI_SUM, OWN_X_COMM);
	MPI_Allreduce(work_v3[2], work_v4[2], num_mesh_transfer, MPI_DOUBLE, MPI_SUM, OWN_X_COMM);


	int offset_mesh = NX * NY * NZ_ / xprocs / yprocs * xid;
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
				im_obl = im + offset_mesh;
                u[0][im] = work_v4[0][im_obl];
                u[1][im] = work_v4[1][im_obl];
                u[2][im] = work_v4[2][im_obl];
			}
		}
	}
#endif
}


/*!
  \brief Update magnitude of k vectors 
 */

inline void Update_K2_OBL(void){
    int im;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, WAVE_X, WAVE_Y, WAVE_Z, KX_int, KY_int, KZ_int, K2, IK2, degree_oblique) private(im)
#endif
    for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
        for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
            for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
                im = SPECTRUMMODE_ARRAYINDEX(i, j, k);
                K2[im] = SQ(WAVE_X * KX_int[im]) +
                         SQ(WAVE_Y * KY_int[im] -
                            WAVE_X * KX_int[im] * degree_oblique) +
                         SQ(WAVE_Z * KZ_int[im]);

                if(K2[im] > 0.0){
                    IK2[im] = 1.0/K2[im];
                }else{
                    IK2[im] = 0.0;
                }

            }
        }
    }
}

inline void Update_Obl_Coord(double **u, const double &delta_gamma){
    int im;
//#pragma omp parallel for shared (NPs, NZ_, u, delta_gamma) private(im)
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
                u[0][im] -= delta_gamma * u[1][im];
            }
        }
    }
}

#endif


