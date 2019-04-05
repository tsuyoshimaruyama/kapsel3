/*!
  \file solute_rhs.cxx
  \brief Routines to compute the terms appearing in the right hand side of the solute advection diffusion equation
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */

#include "solute_rhs.h"

//////

void Solute_impermeability(Particle *p, double **solute_flux_x, double **surface_normal){
    double xp[DIM];
    int x_int[DIM];
    double residue[DIM];
    int r_mesh[DIM];
    double r[DIM];
    double normal[DIM];
    double norm2;
    double dmy;
    int im;
    int dmy_mesh[DIM];

//<COMM 要検討>#ifdef _OPENMP
//<COMM 要検討>#pragma omp parallel for default(none) \
//<COMM 要検討>    private(xp, x_int, residue, dmy_mesh, r_mesh, r, im, norm2, normal, dmy) \
//<COMM 要検討>    shared(Local_Particle_Number, p, NP_domain, Sekibun_cell, Ns, DX, NPs, NZ_, \
//<COMM 要検討>           surface_normal, solute_flux_x)
//<COMM 要検討>#endif
    for (int n = 0; n < Local_Particle_Number; n++) {
        for (int d = 0; d < DIM; d++) {
            xp[d] = p[n].x[d];
        }
        int sw_in_cell = Particle_cell (xp, DX, x_int, residue); // {1, 0} が返ってくる
        sw_in_cell = 1;
        for (int mesh = 0; mesh < NP_domain; mesh++) {
            Relative_coord (Sekibun_cell[mesh], x_int, residue, sw_in_cell, Ns, DX, r_mesh, r);
            if (Range_coord (r_mesh, dmy_mesh) ) {
                im = REALMODE_ARRAYINDEX(dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
                norm2 = 0.0;
                for (int d = 0; d < DIM; d++) {
                    normal[d] = surface_normal[d][im];
                    norm2 += SQ (normal[d]);
                }
                if (norm2 > 0.0) {
                    dmy = 0.0;
                    for (int d = 0; d < DIM; d++) {
                        dmy += solute_flux_x[d][im] * normal[d];
                    }
                    dmy /= norm2;
                    for (int d = 0; d < DIM; d++) {
                        solute_flux_x[d][im] -= normal[d] * dmy;
                    }
                }
            }
        }
    }
}

void Rescale_solute(double *rescale_factor, double *total_solute, double **conc_k, Particle *p
            ,double *phi // working memory
            ,double *conc_x // working memory
		    ){
    int im;

    Reset_phi (phi);
    Make_phi_particle (phi, p);
    for (int n = 0; n < N_spec; n++) {
        A_k2a_out (conc_k[n], conc_x);
        rescale_factor[n] = total_solute[n] / Count_single_solute (conc_x, phi);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(NPs, conc_k, n, rescale_factor) private(im)
#endif
        for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
            for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
                for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
                    im = SPECTRUMMODE_ARRAYINDEX(i, j, k);
                    conc_k[n][im] *= rescale_factor[n];
                }
            }
        }
    }
}

double Count_single_solute(double *conc_x, double *phi_p){
    // set phi(x) before calling
    double dmy = 0.0;
#ifdef _MPI
    double *dmy_mpi;
#endif
    int im;
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:dmy) shared(NPs, NZ_, phi_p, conc_x) private(im)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
                dmy += (1.0 - phi_p[im]) * conc_x[im];
            }
        }
    }
#ifdef _MPI
    dmy_mpi = alloc_1d_double (procs);

    ierr = MPI_Allgather (&dmy, 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_ERRORCHECK (ierr);
    dmy = 0.0;
    for (int n = 0; n < procs; n++) {
        dmy += dmy_mpi[n];
    }
    free_1d_double (dmy_mpi);
#endif
    return dmy * DX3;
}

void Solute_solver_rhs_nonlinear_x_single(double **grad_potential, double *concentration_x, double **solute_flux
                    ,double &valency_e, double &onsager_coeff){
    int im;
    double dmy_grad_pot[DIM];
    double dmy_conc;
    double dmy_interaction;
//#ifdef _OPENMP
//#pragma omp parallel for default(none) shared(NPs, NZ_, grad_potential, concentration_x, valency_e, onsager_coeff, solute_flux, NS_source, omega_rhs, IRHO) private(im, dmy_grad_pot, dmy_conc, dmy_interaction)
//#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
                for (int d = 0; d < DIM; d++) {
                    dmy_grad_pot[d] = grad_potential[d][im];
                }
                dmy_conc = concentration_x[im];
                dmy_interaction = valency_e * onsager_coeff * dmy_conc;
                for (int d = 0; d < DIM; d++) {
                    solute_flux[d][im] += dmy_interaction * dmy_grad_pot[d];
                }
            }
        }
    }
} 
void Add_advection_flux(double **solute_flux, double **u_solvent, double *concentration_x){
    int im;
    double dmy_conc;
//#ifdef _OPENMP
//#pragma omp parallel for default(none) shared(NPs, NZ_, concentration_x, solute_flux, u_solvent) private(im, dmy_conc)
//#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
                dmy_conc = concentration_x[im];
                for (int d = 0; d < DIM; d++) {
                    solute_flux[d][im] += dmy_conc * -u_solvent[d][im];
                }
            }
        }
    }
}

void Diffusion_flux_single(double **diff_flux_x, double *conc_k, double &onsager_coeff, double *dmy_value // working memory
){
    int im;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(NPs, dmy_value, kBT, onsager_coeff, conc_k) private(im)
#endif
    for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
        for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
            for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
                im = SPECTRUMMODE_ARRAYINDEX(i, j, k);
                dmy_value[im] = kBT * onsager_coeff * conc_k[im];
            }
        }
    }
    A_k2da_k (dmy_value, diff_flux_x);
    U_k2u (diff_flux_x);
}

