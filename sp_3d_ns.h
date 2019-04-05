/*!
  \file sp_3d_ns.h
  \author Y. Nakayama
  \date 2006/11/30
  \version 1.9
  \brief Main program file (header)
  \todo documentation
 */
#ifndef SP_3D_NS_H
#define SP_3D_NS_H

#ifdef _MPI
#include <mpi.h>
#endif
#define  NDEBUG
#include <assert.h> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <math.h>
#include <float.h>
#include <time.h>
#include "macro.h"
#include "input.h"
#include "variable.h"
#include "aux_field.h"
#include "make_phi.h"
#include "output.h"
#include "resume.h"
#include "md_force.h"
#include "fluid_solver.h"
#include "particle_solver.h"
#include "f_particle.h"
#include "init_fluid.h"
#include "init_particle.h"
#include "operate_omega.h"
#include "operate_electrolyte.h"
#include "operate_surface.h"
#include "memory_model.h"
#ifdef _MPI
#include "operate_mpi.h"
#include "operate_mpi_particle.h"
#endif

#include "dc.h"
#include "mt19937ar.h"

inline void Electrolyte_free_energy(const Count_SW &OPERATION, FILE *fout, Particle *p, double **Concentration_rhs1, const CTime &jikan){
    static const char *labels[]={"", "nu", "radius", "xi", "M", "I", "kBT", "kBT/M"};
    static char line_label[1<<10];
    static const char *labels_free_energy_nosalt[]={"", "total", "ideal_gas", "electrostatic", "liquid_charge", "total_counterion"};
    static const char *labels_free_energy_salt[]={"", "total", "ideal_gas", "electrostatic", "liquid_charge", "total_positive_ion", "total_negative_ion"};
    //static char line_label_free_energy[1<<10];

    if(OPERATION == INIT){
        sprintf(line_label,"#");
        for(int d=1;d<sizeof(labels)/sizeof(char *);d++){
	        sprintf(line_label,"%s%d:%s ",line_label, d, labels[d]);
        }
    }else if(OPERATION == SHOW || OPERATION == MEAN){
        if(OPERATION == MEAN){
            fprintf_single(fout,"%s\n",line_label);
            fprintf_single(fout, "%g %g %g %g %g %g %g\n", NU, A*DX, XI*DX, MASS[p[0].spec], MOI[p[0].spec], kBT, kBT*IMASS[p[0].spec]);
        }

        double free_energy[3];
        Calc_free_energy_PB(Concentration_rhs1, p, free_energy, up[0], up[1], up[2], jikan);
        double ion_density = 0.; 
        double *n_solute = new double[N_spec];
        Count_solute_each(n_solute, Concentration_rhs1, p, phi, up[0]);
        for(int n=0;n < N_spec;n++){
            ion_density += n_solute[n] * Valency_e[n];
        }
        char line0[1<<10];
        char line1[1<<10];
        int d;
        int dstart;

        if(OPERATION == SHOW){
            dstart = 1;
            sprintf(line0,"#%d:%s ",dstart,"time");
            sprintf(line1, "%d ",jikan.ts);
        } else {
            dstart = 0;
            sprintf(line0,"#");
            sprintf(line1, "");
        }

        if(N_spec == 1){
            for(d=1;d<sizeof(labels_free_energy_nosalt)/sizeof(char *);d++){
                sprintf(line0,"%s%d:%s ", line0, dstart+d, labels_free_energy_nosalt[d]);
            }
            sprintf(line1, "%s%.15g %.15g %.15g %.15g %.15g", line1, free_energy[0], free_energy[1], free_energy[2], ion_density, n_solute[0]);
            sprintf(line0,"%s\n",line0);
            sprintf(line1,"%s\n",line1);
        }else{
            for(d=1;d<sizeof(labels_free_energy_salt)/sizeof(char *);d++){
                sprintf(line0,"%s%d:%s ", line0, dstart+d, labels_free_energy_salt[d]);
            }
        }
	    sprintf(line1, "%s%.15g %.15g %.15g %.15g %.15g %.15g", line1
               ,free_energy[0], free_energy[1], free_energy[2], ion_density, n_solute[0], n_solute[1]);
        sprintf(line0,"%s\n",line0);
        sprintf(line1,"%s\n",line1);

        fprintf_single(fout, "%s%s", line0, line1);
        delete [] n_solute;

	} else {
        fprintf_single(stderr, "invalid OPERATION in Electrolyte_free_energy().\n"); 
        exit_job(EXIT_FAILURE);
    }
}

inline void Calc_shear_rate_eff(){
    static const double ivolume = Ivolume * POW3(DX);
    double s_rate_eff = 0.0;
    double s_rate_eff_local = 0.0;
	int im, im_next;

    //get above xz-slab
#ifdef _MPI
    MPI_Status stat;

    int send_pid, recv_pid, recv_tag;
    int mesh_slab = NPs[REAL][0] * NPs[REAL][2];

    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int k = 0; k < NPs[REAL][2]; k++) {
            im = REALMODE_ARRAYINDEX(i, 0, k);
            im_next =  i * NPs[REAL][2] + k;
            work_v4[0][im_next] = ucp[0][im];
		}
	}
    if(yid == 0) recv_pid = yprocs * (xid + 1) - 1;
	else         recv_pid = procid - 1;
	if(yid == (yprocs-1)) send_pid = yprocs * xid;
	else              send_pid = procid + 1;

    ierr = MPI_Sendrecv(work_v4[0], mesh_slab, MPI_DOUBLE, recv_pid, procid, work_v3[0], mesh_slab, MPI_DOUBLE, send_pid, send_pid, MPI_COMM_WORLD, &stat);
    MPI_ERRORCHECK (ierr);

	if(yid == (yprocs-1)){
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]-1; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    im = REALMODE_ARRAYINDEX(i, j, k);
                    im_next = REALMODE_ARRAYINDEX(i, (j+1), k);
                    s_rate_eff_local += (ucp[0][im_next] - ucp[0][im])/DX;
                }
            }
        }
	} else {
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
					if( j == (NPs[REAL][1] - 1) ){
                        im = REALMODE_ARRAYINDEX(i, j, k);
                        im_next =  i * NPs[REAL][2] + k;
                        s_rate_eff_local += (work_v3[0][im_next] - ucp[0][im])/DX;
					} else {
                        im = REALMODE_ARRAYINDEX(i, j, k);
                        im_next = REALMODE_ARRAYINDEX(i, (j+1), k);
                        s_rate_eff_local += (ucp[0][im_next] - ucp[0][im])/DX;
					}
                }
            }
        }
	}

    MPI_Allreduce(&s_rate_eff_local, &s_rate_eff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_ERRORCHECK (ierr);
#else
	{
#pragma omp parallel for reduction(+:s_rate_eff) 

		for (int i = 0; i < NX; i++) {
			for (int j = 0; j < NY - 1; j++) {
				for (int k = 0; k < NZ; k++) {

					s_rate_eff +=
						(ucp[0][(i*NY*NZ_) + ((j + 1)*NZ_) + k] - ucp[0][(i*NY*NZ_) + (j*NZ_) + k]) / DX;

				}
			}
		}
	}
#endif
    s_rate_eff *= ivolume*NY/(NY - 1);
    Shear_rate_eff = s_rate_eff;
}

inline double Update_strain (double &shear_strain_realized, const CTime &jikan, double **zeta, double uk_dc[DIM], double **u) {
    static const double hivolume = Ivolume * POW3 (DX) * 2.0;
    static const int dmy_ny0 = NY / 4;
    static const int dmy_ny1 = 3 * NY / 4;
    static int update_strain_1st = 0;
    static int ny0, ny1;
    double srate_eff = 0.0;
#ifdef _MPI
    double *dmy_mpi;
#endif
    if (!update_strain_1st) {
        if (dmy_ny1 < PREV_NPs[REAL][1] || NEXT_NPs[REAL][1] <= dmy_ny0) {
            ny0 = -1;
            ny1 = -2;
        } else {
            if (PREV_NPs[REAL][1] <= dmy_ny0 && dmy_ny0 < NEXT_NPs[REAL][1]) {
                ny0 = (int)lrint(dmy_ny0) - PREV_NPs[REAL][1];
            } else {
                ny0 = 0;
            }
            if (PREV_NPs[REAL][1] <= dmy_ny1 && dmy_ny1 < NEXT_NPs[REAL][1]) {
                ny1 = (int)lrint(dmy_ny1) - PREV_NPs[REAL][1];
            } else {
                ny1 = NPs[REAL][1];
            }
        }
        update_strain_1st++;
    }
    Zeta_k2u_k (zeta, uk_dc, u);
    A_k2dya_k (u[0], u[1]);
    A_k2a (u[1]);
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:srate_eff) shared(NPs, NZ_, u, ny0, ny1)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = ny0; j < ny1; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                srate_eff += u[1][REALMODE_ARRAYINDEX(i, j, k)];
            }
        }
	}
#ifdef _MPI
    dmy_mpi = alloc_1d_double (procs);
    ierr = MPI_Allgather (&srate_eff, 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_ERRORCHECK (ierr);
    srate_eff = 0.0;
    for (int n = 0; n < procs; n++) {
        srate_eff += dmy_mpi[n];
    }
    free_1d_double (dmy_mpi);
#endif
    srate_eff *= - (hivolume);
    shear_strain_realized += srate_eff * jikan.dt_fluid;
    return srate_eff;
}

inline void Mean_shear_stress(const Count_SW &OPERATION, FILE *fout, Particle *p, const CTime &jikan, const double &srate_eff){
    static const char *labels_zz_dc[]={"", "time", "shear_rate_temporal", "shear_strain_temporal", "shear_stress_temporal", "viscosity"};
    static const char *labels_zz_ac[]={"", "time", "shear_rate_temporal", "shear_strain_temporal"
                 ,"shear_stress_temporal", "shear_inertia_stress_temporal", "apparent_shear_stress"};
    static const char *labels_le[]={"", "time", "shear_rate", "degree_oblique", "shear_strain_temporal", "lj_dev_stress_temporal"
                 , "shear_stress_temporal_old", "shear_stress_temporal_new", "reynolds_stress", "viscosity"};

    static char line_label[1<<10];

    if(OPERATION == INIT){
        sprintf(line_label,"#");
        if(SW_EQ==Shear_Navier_Stokes_Lees_Edwards){
	        for(int d=1; d<sizeof(labels_le)/sizeof(char *); d++)
                sprintf(line_label, "%s%d:%s ", line_label, d, labels_le[d]);
        }else if(SW_EQ==Shear_Navier_Stokes){
            if(!Shear_AC){
                for(int d=1; d<sizeof(labels_zz_dc)/sizeof(char *); d++)
                    sprintf(line_label, "%s%d:%s ", line_label, d, labels_zz_dc[d]);
            }else{
                for(int d=1; d<sizeof(labels_zz_dc)/sizeof(char *); d++)
                    sprintf(line_label, "%s%d:%s ", line_label, d, labels_zz_ac[d]);
            }
        }else{
            fprintf_single(stderr, "Error: Incorrect Shear calculation\n");
            exit_job(EXIT_FAILURE);
        }
        fprintf_single(fout,"%s\n",line_label);
	}else if(OPERATION == SHOW){
        double stress[DIM][DIM] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
        double hydro_stress[DIM][DIM] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
        double hydro_stress_new[DIM][DIM] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};

        double strain_output = Shear_strain_realized;

        if (SW_EQ == Shear_Navier_Stokes_Lees_Edwards) {
            Calc_hydro_stress(jikan, p, phi, Hydro_force, hydro_stress);
            Calc_hydro_stress(jikan, p, phi, Hydro_force_new, hydro_stress_new);
            double dev_stress = (SW_PT == rigid ? rigid_dev_shear_stress_lj : dev_shear_stress_lj);
            double dev_stress_rot = (SW_PT == rigid ? rigid_dev_shear_stress_rot : dev_shear_stress_rot);

            fprintf_single(fout, "%16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n"
                      , jikan.time, srate_eff, degree_oblique, strain_output, dev_stress
                      , hydro_stress[1][0], hydro_stress_new[1][0], Inertia_stress
                      ,(hydro_stress_new[1][0] + Inertia_stress + dev_stress)/srate_eff + ETA);
        } else if(SW_EQ == Shear_Navier_Stokes){ 
            if(!Shear_AC){
                Calc_shear_stress(jikan, p, phi, Shear_force, stress);
                fprintf_single(fout, "%16.8g %16.8g %16.8g %16.8g %16.8g\n"
                      , jikan.time, srate_eff, strain_output, -stress[1][0], -stress[1][0]/srate_eff);
            }else{
                Calc_shear_stress(jikan, p, phi, Shear_force, stress);
                fprintf_single(fout, "%16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n"
                      , jikan.time, srate_eff, strain_output, -stress[1][0], Inertia_stress, -stress[1][0]+Inertia_stress);
            }
		}
	}else {
        fprintf_single(stderr, "invalid OPERATION in Mean_shear_stress().\n"); 
        exit_job(EXIT_FAILURE);
	}
}

#endif

