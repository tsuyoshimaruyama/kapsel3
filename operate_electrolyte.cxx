/*!
  \file operate_electrolyte.cxx
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Routines to compute the charge distributions and forces
 */
#include "operate_electrolyte.h"

inline void Add_external_electric_field_x(double *potential, const CTime &jikan){
    double external[DIM];
    double time = jikan.time;
    int im;
    for (int d = 0; d < DIM; d++) {
        external[d] = E_ext[d];
        if (AC) {
            external[d] = sin (Angular_Frequency * time);
        }
    }
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, PREV_NPs, NZ_, external, potential, DX) private(im)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
                potential[im] -= ( (external[0] * (double) (i + PREV_NPs[REAL][0]) +
                                    external[1] * (double) (j + PREV_NPs[REAL][1]) +
                                    external[2] * (double) (k + PREV_NPs[REAL][2])) * DX);
            }
        }
    }
}
void Conc_k2charge_field(Particle *p, double **conc_k, double *charge_density
			 ,double *phi_p // working memory
			 ,double *dmy_value // working memory
){
    int im;
    Reset_phi (phi_p);
    Reset_phi (charge_density);
    Make_phi_qq_particle (phi_p, charge_density, p);
    for (int n = 0; n < N_spec; n++) {
        A_k2a_out (conc_k[n], dmy_value);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, charge_density, Valency_e, n, dmy_value, phi_p) private(im)
#endif
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    im = REALMODE_ARRAYINDEX(i, j, k);
                    charge_density[im] += Valency_e[n] * dmy_value[im] * (1.0 - phi_p[im]);
                }
            }
        }
    }
}

void Charge_field_k2Coulomb_potential_k_PBC(double *potential){
    const double iDielectric_cst = 1.0 / Dielectric_cst;
    int im;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, potential, IK2, iDielectric_cst) private(im)
#endif
    for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
        for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
            for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
                im = SPECTRUMMODE_ARRAYINDEX(i, j, k);
                potential[im] *= (IK2[im] * iDielectric_cst);
            }
        }
    }
}

void Make_Coulomb_force_x_on_fluid(double **force, Particle *p, double **conc_k
				    ,double *charge_density // working memory
				    ,double *potential // working memory
				    ,const CTime &jikan){
    int im;
    double electric_field;
    double external[DIM];
    double time = jikan.time;
    Conc_k2charge_field (p, conc_k, charge_density, force[0], force[1]);
    A2a_k_out (charge_density, potential);
    Charge_field_k2Coulomb_potential_k_PBC (potential);
    A_k2da_k (potential, force);
    U_k2u (force);
    if (Fixed_particle) {
        Reset_phi (phi);
        Reset_phi (charge_density);
        Make_phi_qq_fixed_particle (phi, charge_density, p);
        for (int n = 0; n < N_spec; n++) {
            A_k2a_out (conc_k[n], potential);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, charge_density, Valency_e, n, potential, phi) private(im)
#endif
            for (int i = 0; i < NPs[REAL][0]; i++) {
                for (int j = 0; j < NPs[REAL][1]; j++) {
                    for (int k = 0; k < NPs[REAL][2]; k++) {
                        im = REALMODE_ARRAYINDEX(i, j, k);
                        charge_density[im]
                          += Valency_e[n] * potential[im] * (1.0 - phi[im]);
                    }
                }
            }
        }
    }
    for (int d = 0; d < DIM; d++) {
        external[d] = E_ext[d];
        if (AC) {
            external[d] = sin (Angular_Frequency * time);
        }
    }
#ifdef _OPENMP
#pragma omp parallel default(none) shared(NPs, NZ_, External_field, force, charge_density, external)
#endif
    for (int d = 0; d < DIM; d++) {
#ifdef _OPENMP
#pragma omp for nowait private(im, electric_field)
#endif
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    im = REALMODE_ARRAYINDEX(i, j, k);
                    electric_field = -force[d][im];
                    if (External_field) {
                        electric_field += external[d];
                    }
                    force[d][im] = charge_density[im] * electric_field;
                }
            }
        }
    }
}

void Make_phi_qq_particle(double *phi, double *surface, Particle *p){
    double abs_total_surface_charge = 0.0;
    int im;
    int r_mesh[DIM];
    int dmy_mesh[DIM];
    int x_int[DIM];
    double dmy = 0.0;
    double dmy_surface_charge;
    double r[DIM];
    double rescale;
    double residue[DIM];
    double x[DIM];
    double xp[DIM];
#ifdef _MPI
    double *dmy_mpi;
#endif

#ifdef _OPENMP
    Reset_mesh(tmp_buffer1);
    Reset_mesh(tmp_buffer2);
#else
    tmp_buffer1 = phi;
    tmp_buffer2 = surface;
#endif
#ifdef _OPENMP
#pragma omp parallel default(none) shared(Local_Particle_Number, Sekibun_cell, NP_domain, THREADNUM, Ns, NPs, NZ_, L, DX, RADIUS, p, Surface_charge_e, abs_total_surface_charge, phi, surface, tmp_buffer1, tmp_buffer2, XI)
    {
#pragma omp for reduction(+:abs_total_surface_charge) private(r, r_mesh, residue, x, xp, x_int, dmy, dmy_mesh, dmy_surface_charge, im)
#endif
        for (int n = 0; n < Local_Particle_Number; n++) {
            dmy_surface_charge = Surface_charge_e[p[n].spec];
            for (int d = 0; d < DIM; d++) {
                xp[d] = p[n].x[d];
                assert ( (p[n].x[d] >= 0) || (p[n].x[d] < L[d]) );
            }
            int sw_in_cell = Particle_cell (xp, DX, x_int, residue); // {1, 0} が返ってくる
            sw_in_cell = 1;
            for (int mesh = 0; mesh < NP_domain; mesh++) {
                Relative_coord (Sekibun_cell[mesh], x_int, residue, sw_in_cell, Ns, DX, r_mesh, r);
                if (Range_coord (r_mesh, dmy_mesh) ) {
#ifdef _OPENMP
                    im = REALMODE_ARRAYINDEX_MESH(omp_get_thread_num(), dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
#else
                    im = REALMODE_ARRAYINDEX(dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
#endif
                    for (int d = 0; d < DIM; d++) {
                        x[d] = r_mesh[d] * DX;
                    }
                    dmy = Distance (x, xp);
                    tmp_buffer1[im] += Phi (dmy);
                    tmp_buffer2[im] += dmy_surface_charge * DPhi_compact_sin (dmy);
                }
            }
            abs_total_surface_charge += fabs (dmy_surface_charge);
        }
#ifdef _OPENMP
#pragma omp for
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int t = 0; t < THREADNUM; t++) {
                for (int j = 0; j < NPs[REAL][1]; j++) {
                    for (int k = 0; k < NPs[REAL][2]; k++) {
                        phi[REALMODE_ARRAYINDEX(i, j, k)] += tmp_buffer1[REALMODE_ARRAYINDEX_MESH(t, i, j, k)];
                        surface[REALMODE_ARRAYINDEX(i, j, k)] += tmp_buffer2[REALMODE_ARRAYINDEX_MESH(t, i, j, k)];
                    }
                }
            }
        }
    }
#endif
    // volume of surface section is normalized to unity
    dmy = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:dmy) shared(NPs, NZ_, surface)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                dmy += fabs (surface[REALMODE_ARRAYINDEX(i, j, k)]);
            }
        }
    }
#ifdef _MPI
    dmy_mpi = alloc_1d_double (procs);

    ierr = MPI_Allgather (&dmy, 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE,
                          MPI_COMM_WORLD);

    MPI_ERRORCHECK (ierr);
    dmy = 0.0;
    for (int n = 0; n < procs; n++) {
        dmy += dmy_mpi[n];
    }
    free_1d_double (dmy_mpi);
#endif
    rescale = abs_total_surface_charge / (dmy * DX3);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, surface, rescale) private(im)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                surface[REALMODE_ARRAYINDEX(i, j, k)] *= rescale;
            }
        }
    }
}

void Make_phi_qq_fixed_particle(double *phi, double *surface, Particle *p){
    double abs_total_surface_charge = 0.0;
    int im;
    int r_mesh[DIM];
    int dmy_mesh[DIM];
    int x_int[DIM];
    double dmy;
    double dmy_surface_charge;
    double r[DIM];
    double rescale;
    double residue[DIM];
    double x[DIM];
    double xp[DIM];
#ifdef _MPI
    double *dmy_mpi;
#endif

#ifdef _OPENMP
    Reset_mesh(tmp_buffer1);
    Reset_mesh(tmp_buffer2);
#else
    tmp_buffer1 = phi;
    tmp_buffer2 = surface;
#endif
#ifdef _OPENMP
#pragma omp parallel default(none) shared(Local_Particle_Number, Sekibun_cell, NP_domain, THREADNUM, Ns, NPs, NZ_, L, DX, RADIUS, p, Surface_charge_e, abs_total_surface_charge, phi, surface, tmp_buffer1, tmp_buffer2, XI)
    {
#pragma omp for reduction(+:abs_total_surface_charge) private(r, r_mesh, residue, x, xp, x_int, dmy, dmy_mesh, dmy_surface_charge, im)
#endif
        for (int n = 0; n < Local_Particle_Number; n++) {
            dmy_surface_charge = Surface_charge_e[p[n].spec];
            for (int d = 0; d < DIM; d++) {
                xp[d] = p[n].x[d];
                assert ( (p[n].x[d] >= 0) || (p[n].x[d] < L[d]) );
            }
            int sw_in_cell = Particle_cell (xp, DX, x_int, residue); // {1, 0} が返ってくる
            sw_in_cell = 1;
			for (int mesh = 0; mesh < NP_domain; mesh++) {
                Relative_coord (Sekibun_cell[mesh], x_int, residue, sw_in_cell, Ns, DX, r_mesh, r);
                if (Range_coord (r_mesh, dmy_mesh) ) {
#ifdef _OPENMP
                    im = REALMODE_ARRAYINDEX_MESH(omp_get_thread_num(), dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
#else
                    im = REALMODE_ARRAYINDEX(dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
#endif
                    for (int d = 0; d < DIM; d++) {
                        x[d] = r_mesh[d] * DX;
                    }
                    dmy = Distance (x, xp);
                    tmp_buffer1[im] += Phi (dmy);
                    tmp_buffer2[im] += dmy_surface_charge * DPhi_compact_sin (dmy);
                }
            }
            abs_total_surface_charge += fabs (dmy_surface_charge);
        }
#ifdef _OPENMP
#pragma omp for
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int t = 0; t < THREADNUM; t++) {
                for (int j = 0; j < NPs[REAL][1]; j++) {
                    for (int k = 0; k < NPs[REAL][2]; k++) {
                        phi[REALMODE_ARRAYINDEX(i, j, k)] += tmp_buffer1[REALMODE_ARRAYINDEX_MESH(t, i, j, k)];
                        surface[REALMODE_ARRAYINDEX(i, j, k)] += tmp_buffer2[REALMODE_ARRAYINDEX_MESH(t, i, j, k)];
                    }
                }
            }
        }
    }
#endif
    dmy = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:dmy) shared (NPs, NZ_, surface, phi) private(im)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
                dmy += fabs (surface[im]) * phi[im];
            }
        }
    }
#ifdef _MPI
    dmy_mpi = alloc_1d_double (procs);

    ierr = MPI_Allgather (&dmy, 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE,
                          MPI_COMM_WORLD);

    MPI_ERRORCHECK (ierr);
    dmy = 0.0;
    for (int n = 0; n < procs; n++) {
        dmy += dmy_mpi[n];
    }
    free_1d_double (dmy_mpi);
#endif
    rescale = abs_total_surface_charge * 0.5 / (dmy * DX3);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(NPs, NX, NY, NZ_, surface, rescale)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                surface[REALMODE_ARRAYINDEX(i, j, k)] *= rescale;
            }
        }
    }
}
void Calc_free_energy_PB(double **conc_k, Particle *p, double *free_energy
			  ,double *phi // working memory
			  ,double *charge_density // working memory
			  ,double *dmy_value // working memory
			  ,const CTime &jikan){

    double free_energy_ideal = 0.0;
    double free_energy_electrostatic = 0.0;
    int im;
    double dmy_conc;
#ifdef _MPI
    double *dmy_mpi;
#endif

    Reset_phi (phi);
    Make_phi_particle (phi, p);
    for (int n = 0; n < N_spec; n++) {
        A_k2a_out (conc_k[n], dmy_value);
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:free_energy_ideal) shared (NPs, NZ_, dmy_value, phi) private(im, dmy_conc)
#endif
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    im = REALMODE_ARRAYINDEX(i, j, k);
                    dmy_conc = dmy_value[im];
                    free_energy_ideal += (1.0 - phi[im]) * dmy_conc * (log (dmy_conc) - 1.0);
                }
            }
        }
    }
#ifdef _MPI
    dmy_mpi = alloc_1d_double (procs);

    ierr = MPI_Allgather (&free_energy_ideal, 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_ERRORCHECK (ierr);
    free_energy_ideal = 0.0;
    for (int n = 0; n < procs; n++) {
        free_energy_ideal += dmy_mpi[n];
    }
#endif
    free_energy_ideal *= (kBT * DX3);
    Conc_k2charge_field (p, conc_k, charge_density, phi, dmy_value);
    A2a_k_out (charge_density, dmy_value);
    Charge_field_k2Coulomb_potential_k_PBC (dmy_value);
    A_k2a (dmy_value);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(NPs, NZ_, dmy_value)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                dmy_value[REALMODE_ARRAYINDEX(i, j, k)] *= 0.5;
            }
        }
    }
    if (External_field) {
        Add_external_electric_field_x (dmy_value, jikan);
    }
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:free_energy_electrostatic) shared (NPs, NZ_, charge_density, dmy_value) private(im)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
                free_energy_electrostatic += charge_density[im] * dmy_value[im];
            }
        }
    }
#ifdef _MPI

    ierr = MPI_Allgather (&free_energy_electrostatic, 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_ERRORCHECK (ierr);
    free_energy_electrostatic = 0.0;
    for (int n = 0; n < procs; n++) {
        free_energy_electrostatic += dmy_mpi[n];
    }
    free_1d_double (dmy_mpi);
#endif
    free_energy_electrostatic *= DX3;
    free_energy[1] = free_energy_ideal;
    free_energy[2] = free_energy_electrostatic;
    free_energy[0] = free_energy[1] + free_energy[2];
}

inline void Set_steadystate_ion_density(double **Concentration, Particle *p, CTime &jikan){
    double free_energy[3];
    double free_energy_previous;
    double dmy_uk_dc[DIM];
    double dum;
    int STEP = 10000;
    const Index_range ijk_range[] = {
        {0, TRN_X - 1, 0, TRN_Y - 1, 0, 2 * TRN_Z - 1},
        {0, TRN_X - 1, NY - TRN_Y + 1, NY - 1, 0, 2 * TRN_Z - 1},
        {NX - TRN_X + 1, NX - 1, 0, TRN_Y - 1, 0, 2 * TRN_Z - 1},
        {NX - TRN_X + 1, NX - 1, NY - TRN_Y + 1, NY - 1, 0, 2 * TRN_Z - 1}
    };
    const int n_ijk_range = sizeof (ijk_range) / sizeof (Index_range);
    Calc_free_energy_PB (Concentration, p, free_energy, up[0], up[1], up[2], jikan);
    free_energy_previous = free_energy[0];
    for (int m = 1; m < STEP; m++) {
        for (int d = 0; d < DIM; d++) {
            Reset_phi (u[d]);
        }
        U_k2zeta_k (u, up, dmy_uk_dc);
        Ion_diffusion_solver_Euler (up, jikan, dmy_uk_dc, Concentration, p,
                                    ijk_range, n_ijk_range);
        Calc_free_energy_PB (Concentration, p, free_energy, up[0], up[1], up[2],
                             jikan);
        dum = fabs (free_energy[0] - free_energy_previous);
        if (dum < TOL) {
            break;
        }
        free_energy_previous = free_energy[0];
    }
}

inline void Set_uniform_ion_charge_density_nosalt(double *Concentration, double *Total_solute, Particle *p){
    double Counterion_density;
    double volume_phi = 0.0;
#ifdef _MPI
    double *dmy_mpi;
#endif
    if (N_spec != 1) {          //N_spec = 1 にかぎることに注意
        fprintf_single (stderr, "invalid (number of ion spencies) = %d\n", N_spec);
        exit_job (EXIT_FAILURE);
    }
    Reset_phi (phi);
    Make_phi_particle (phi, p);
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:volume_phi) shared(NPs, NZ_, phi)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                volume_phi += (1.0 - phi[REALMODE_ARRAYINDEX(i, j, k)]);
            }
        }
    }
#ifdef _MPI
    dmy_mpi = alloc_1d_double (procs);

    ierr = MPI_Allgather (&volume_phi, 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_ERRORCHECK (ierr);
    volume_phi = 0.0;
    for (int n = 0; n < procs; n++) {
        volume_phi += dmy_mpi[n];
    }
    free_1d_double (dmy_mpi);
#endif
    Counterion_density = Counterion_number / (volume_phi * DX3);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(NPs, NZ_, Concentration, Counterion_density)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                Concentration[REALMODE_ARRAYINDEX(i, j, k)]
                  = Counterion_density;
            }
        }
    }
    Total_solute[0] = Count_single_solute (Concentration, phi);
    A2a_k (Concentration);
}

inline void Set_Poisson_Boltzmann_ion_charge_density_nosalt(double **Concentration, double *Total_solute, Particle *p){
    int STEP = 100000;
    int im;
    double alp = 0.01; //0.001;//0.1;
    double dmy1 = Valency[0] * Elementary_charge / kBT;
    double *e_potential = up[2];
    double rho_0 = 0.0;
    double dmy;
    double sum = 0.0;
    double error = 0.0;
    double *rescale_factor;
#ifdef _MPI
    double *dmy_mpi;
    double *dmy_mpilocal;
    double *dmy_mpiall;
#endif

    if (N_spec != 1) {          //N_spec = 1 にかぎることに注意
        fprintf (stderr, "invalid (number of ion spencies) = %d\n", N_spec);
        exit_job (EXIT_FAILURE);
    }
    Reset_phi (phi);
    Make_phi_particle (phi, p);
#ifdef _MPI
    dmy_mpi = alloc_1d_double (procs);
    dmy_mpilocal = alloc_1d_double (2);
    dmy_mpiall = alloc_1d_double (2 * procs);
#endif
    for (int m = 1; m < STEP; m++) {
        Conc_k2charge_field (p, Concentration, e_potential, up[0], up[1]);
        A2a_k (e_potential);
        Charge_field_k2Coulomb_potential_k_PBC (e_potential);
        A_k2a (e_potential);
        rho_0 = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:rho_0) shared (NPs, NZ_, phi, dmy1, e_potential) private(im)
#endif
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    im = REALMODE_ARRAYINDEX(i, j, k);
                    rho_0 += (1.0 - phi[im]) * exp (-dmy1 * e_potential[im]);
                }
            }
        }
#ifdef _MPI

        ierr = MPI_Allgather (&rho_0, 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE,
                              MPI_COMM_WORLD);

        MPI_ERRORCHECK (ierr);
        rho_0 = 0.0;
        for (int n = 0; n < procs; n++) {
            rho_0 += dmy_mpi[n];
        }
#endif
        rho_0 = Counterion_number / (rho_0 * DX3);
        A_k2a (Concentration[0]);
        sum = error = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:sum) reduction(+:error) shared (NPs, NZ_, rho_0, dmy1, e_potential, Concentration, alp) private(im, dmy)
#endif
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    im = REALMODE_ARRAYINDEX(i, j, k);
                    dmy = rho_0 * exp (-dmy1 * e_potential[im]);
                    Concentration[0][im]
                      = (1.0 - alp) * Concentration[0][im] + alp * dmy;
                    sum += fabs (dmy);
                    error += fabs (Concentration[0][im] - dmy);
                }
            }
        }
        A2a_k (Concentration[0]);
        rescale_factor = new double[N_spec];
        Rescale_solute (rescale_factor, Total_solute, Concentration, p, up[0],
                        up[1]);
        delete[]rescale_factor;
#ifdef _MPI
        dmy_mpilocal[0] = sum;
        dmy_mpilocal[1] = error;

        ierr = MPI_Allgather (dmy_mpilocal, 2, MPI_DOUBLE, dmy_mpiall, 2, MPI_DOUBLE, MPI_COMM_WORLD);

        MPI_ERRORCHECK (ierr);
        sum = error = 0.0;
        for (int n = 0; n < procs; n++) {
            sum += dmy_mpiall[2 * n];
            error += dmy_mpiall[2 * n + 1];
        }
#endif
        if (error / sum < TOL) {
            break;
        }
    }
#ifdef _MPI
    free_1d_double (dmy_mpi);
    free_1d_double (dmy_mpilocal);
    free_1d_double (dmy_mpiall);
#endif
    A_k2a_out (Concentration[0], up[0]);
    Total_solute[0] = Count_single_solute (up[0], phi);
}
inline void Set_uniform_ion_charge_density_salt(double **Concentration, double *Total_solute, Particle *p){
    int im;
    double volume_phi = 0.0;
    double dmy_surface_charge;
    double dmy_particle_number;
    double positive_ion_number = 0.0;
    double negative_ion_number = 0.0;
    double positive_ion_density;
    double negative_ion_density;
    double Rho_inf = kBT * Dielectric_cst / SQ (Elementary_charge * Debye_length);
    double Rho_inf_positive_ion = Rho_inf / (Valency[0] * (Valency[0] - Valency[1]) );
    double Rho_inf_negative_ion = -Rho_inf / (Valency[1] * (Valency[0] - Valency[1]) );
#ifdef _MPI
    double *dmy_mpi;
#endif

    if (N_spec != 2) { //N_spec = 2 にかぎることに注意
        fprintf_single (stderr, "invalid (number of ion spencies) = %d\n", N_spec);
        exit_job (EXIT_FAILURE);
    }
    Reset_phi (phi);
    Make_phi_particle (phi, p);
    for (int i = 0; i < Component_Number; i++) {
        dmy_surface_charge = Surface_charge[i];
        dmy_particle_number = (double) Particle_Numbers[i];
        if (dmy_surface_charge < 0.0) {
            positive_ion_number
              += fabs (dmy_surface_charge * dmy_particle_number / Valency[0]);
        } else if (dmy_surface_charge > 0.0) {
            negative_ion_number
              += fabs (dmy_surface_charge * dmy_particle_number / Valency[1]);
        }
    }
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:volume_phi) shared(NPs, NZ_, phi)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                volume_phi += (1.0 - phi[REALMODE_ARRAYINDEX(i, j, k)]);
            }
        }
    }
#ifdef _MPI
    dmy_mpi = alloc_1d_double (procs);

    ierr = MPI_Allgather (&volume_phi, 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_ERRORCHECK (ierr);
    volume_phi = 0.0;
    for (int n = 0; n < procs; n++) {
        volume_phi += dmy_mpi[n];
    }
    free_1d_double (dmy_mpi);
#endif
    positive_ion_density = positive_ion_number / (volume_phi * DX3);
    negative_ion_density = negative_ion_number / (volume_phi * DX3);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, Concentration, positive_ion_density, Rho_inf_positive_ion, negative_ion_density, Rho_inf_negative_ion) private(im)
#endif
    for (int i = 0; i < NPs[REAL][0]; i++) {
        for (int j = 0; j < NPs[REAL][1]; j++) {
            for (int k = 0; k < NPs[REAL][2]; k++) {
                im = REALMODE_ARRAYINDEX(i, j, k);
                Concentration[0][im] = positive_ion_density + Rho_inf_positive_ion;
                Concentration[1][im] = negative_ion_density + Rho_inf_negative_ion;
            }
        }
    }
    Total_solute[0] = Count_single_solute (Concentration[0], phi);
    Total_solute[1] = Count_single_solute (Concentration[1], phi);
    U2u_k (Concentration, N_spec);
    fprintf_single (stderr, "# density of positive and negative ions %g %g\n",
                    positive_ion_density + Rho_inf_positive_ion,
                    negative_ion_density + Rho_inf_negative_ion);
}

inline void Set_Poisson_Boltzmann_ion_charge_density_salt(double **Concentration, double *Total_solute, Particle *p){
    int STEP = 100000;
    int im;
    double *e_potential0 = up[0];
    double *e_potential1 = up[1];
    double *dmy_value0 = f_particle[0];
    double *dmy_value1 = f_particle[1];
    double alp = 0.01; //0.001;//0.1;
    double a, b, c, nu;
    double dmy;
    double dmy_pot;
    double dmy_phi;
    double sum;
    double error;
    double dmy0 = Elementary_charge * Valency[0] / kBT;
    double dmy1 = Elementary_charge * Valency[1] / kBT;
    double dmy_rho_positive_ion = 0.0;
    double dmy_rho_negative_ion = 0.0;
    double Rho_inf = kBT * Dielectric_cst / SQ (Elementary_charge * Debye_length);
    double Rho_inf_positive_ion = Rho_inf / (Valency[0] * (Valency[0] - Valency[1]) );
    double Rho_inf_negative_ion = -Rho_inf / (Valency[1] * (Valency[0] - Valency[1]) );
#ifdef _MPI
    double *dmy_mpiall;
    double *dmy_mpilocal;
#endif

    if (N_spec != 2) { //N_spec = 2 にかぎることに注意
        fprintf (stderr, "invalid (number of ion spencies) = %d\n", N_spec);
        exit_job (EXIT_FAILURE);
    }
    //対称z:z電解質溶液 Valency_positive_ion=z, Valency_negative_ion=-z の場合にかぎる
    if (Valency[0] != -Valency[1]) {
        fprintf (stderr, "invalid (Valency_positive_ion, Valency_negative_ion) = %g %g\n",
                 Valency[0], Valency[1]);
        fprintf (stderr, "set Valency_positive_ion = - Valency_negative_ion\n");
        exit_job (EXIT_FAILURE);
    }
    //仕込みのイオン濃度(バルクでカウンターイオン, 共イオンは
    //濃度Rho_inf_positive_ion, Rho_inf_negative_ionで平衡であると仮定)
    fprintf_single (stderr, "# bulk density of positive and negative ions  %g %g\n",
                    Rho_inf_positive_ion, Rho_inf_negative_ion);
    Reset_phi (phi);
    Make_phi_particle (phi, p);
    Reset_phi (e_potential0);
#ifdef _MPI
    dmy_mpilocal = alloc_1d_double (2);
    dmy_mpiall = alloc_1d_double (2 * procs);
#endif
    for (int m = 1; m < STEP; m++) {
        dmy_rho_positive_ion = 0.0;
        dmy_rho_negative_ion = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:dmy_rho_positive_ion) reduction(+:dmy_rho_negative_ion) shared (NPs, NZ_, phi, dmy0, dmy1, e_potential0) private(im, dmy_phi)
#endif
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    im = REALMODE_ARRAYINDEX(i, j, k);
                    dmy_phi = phi[im];
                    dmy_rho_positive_ion
                      += (1.0 - dmy_phi) * exp (-dmy0 * e_potential0[im]);
                    dmy_rho_negative_ion
                      += (1.0 - dmy_phi) * exp (-dmy1 * e_potential0[im]);
                }
            }
        }
#ifdef _MPI
        dmy_mpilocal[0] = dmy_rho_positive_ion;
        dmy_mpilocal[1] = dmy_rho_negative_ion;

        ierr = MPI_Allgather (dmy_mpilocal, 2, MPI_DOUBLE, dmy_mpiall, 2, MPI_DOUBLE, MPI_COMM_WORLD);

        MPI_ERRORCHECK (ierr);
        dmy_rho_positive_ion = dmy_rho_negative_ion = 0.0;
        for (int n = 0; n < procs; n++) {
            dmy_rho_positive_ion += dmy_mpiall[2 * n];
            dmy_rho_negative_ion += dmy_mpiall[2 * n + 1];
        }
#endif
        dmy_rho_positive_ion *= DX3;
        dmy_rho_negative_ion *= DX3;
        a = Valency[0] * Rho_inf_positive_ion * dmy_rho_positive_ion;
        b = -Surface_ion_number;
        c = Valency[1] * Rho_inf_negative_ion * dmy_rho_negative_ion;
        nu = (-b + sqrt (b * b - 4.0 * a * c) ) / (2.0 * a);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, e_potential0, e_potential1, dmy0, dmy1, nu, Rho_inf_positive_ion, Rho_inf_negative_ion, Concentration) private(im, dmy_pot)
#endif
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    im = REALMODE_ARRAYINDEX(i, j, k);
                    dmy_pot = e_potential0[im];
                    Concentration[0][im]
                      = Rho_inf_positive_ion * exp (-dmy0 * dmy_pot) * nu;
                    Concentration[1][im]
                      = Rho_inf_negative_ion * exp (-dmy1 * dmy_pot) / nu;
                    e_potential1[im] = dmy_pot;
                }
            }
        }
        Reset_phi (dmy_value0);
        Reset_phi (dmy_value1);
        Make_phi_qq_particle (dmy_value0, dmy_value1, p);
        for (int n = 0; n < N_spec; n++) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (NPs, NZ_, dmy_value1, Valency_e, n, Concentration, dmy_value0) private(im)
#endif
            for (int i = 0; i < NPs[REAL][0]; i++) {
                for (int j = 0; j < NPs[REAL][1]; j++) {
                    for (int k = 0; k < NPs[REAL][2]; k++) {
                        im = REALMODE_ARRAYINDEX(i, j, k);
                        dmy_value1[im]
                          += Valency_e[n] * Concentration[n][im] * (1.0 - dmy_value0[im]);
                    }
                }
            }
        }
        A2a_k (dmy_value1);
        Charge_field_k2Coulomb_potential_k_PBC (dmy_value1);
        A_k2a (dmy_value1);
        sum = error = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:sum) reduction(+:error) shared (NPs, NZ_, e_potential0, dmy_value1, alp) private(im, dmy)
#endif
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    im = REALMODE_ARRAYINDEX(i, j, k);
                    dmy = alp * (dmy_value1[im] - e_potential0[im]);
                    e_potential0[im] += dmy;
                    sum += fabs (e_potential0[im]);
                    error += fabs (dmy);
                }
            }
        }
#ifdef _MPI
        dmy_mpilocal[0] = sum;
        dmy_mpilocal[1] = error;

        ierr = MPI_Allgather (dmy_mpilocal, 2, MPI_DOUBLE, dmy_mpiall, 2, MPI_DOUBLE, MPI_COMM_WORLD);

        MPI_ERRORCHECK (ierr);
        sum = error = 0.0;
        for (int n = 0; n < procs; n++) {
            sum += dmy_mpiall[2 * n];
            error += dmy_mpiall[2 * n + 1];
        }
#endif
        if (error / sum < TOL) {
            break;
        }
    }
#ifdef _MPI
    free_1d_double (dmy_mpilocal);
    free_1d_double (dmy_mpiall);
#endif
    Total_solute[0] = Count_single_solute (Concentration[0], phi);
    Total_solute[1] = Count_single_solute (Concentration[1], phi);
    U2u_k (Concentration, N_spec);
}

void Init_rho_ion(double **Concentration, Particle *p, CTime &jikan){
    double dmy0;
    double dmy1;
    double volume_phi = 0.0;
    double Counterion_density;
#ifdef _MPI
    double *dmy_mpi;
#endif

    //コロイド表面電荷は0以外に設定する
    for (int i = 0; i < Component_Number; i++) {
        if (Surface_charge[i] == 0.0) {
            fprintf (stderr, "invalid (Surface_charge) = %g\n",
                     Surface_charge[i]);
            exit_job (EXIT_FAILURE);
        }
        Surface_charge_e[i] = Surface_charge[i] * Elementary_charge;
    }
    Bjerrum_length = SQ (Elementary_charge) / (PI4 * kBT * Dielectric_cst);
    Surface_ion_number = 0.0;
    for (int i = 0; i < Component_Number; i++) {
        Surface_ion_number -= Surface_charge[i] * Particle_Numbers[i];
    }
    if (N_spec == 1) {
        Valency[0] = Valency_counterion;
        Onsager_coeff[0] = Onsager_coeff_counterion;
        Counterion_number = Surface_ion_number / Valency_counterion;
    } else if (N_spec == 2) {
        Valency[0] = Valency_positive_ion;
        Valency[1] = Valency_negative_ion;
        Onsager_coeff[0] = Onsager_coeff_positive_ion;
        Onsager_coeff[1] = Onsager_coeff_negative_ion;
    }
    for (int n = 0; n < N_spec; n++) {
        Valency_e[n] = Valency[n] * Elementary_charge;
    }
    if (N_spec == 1) {
        if (Component_Number == 1) {
            //イオン価数はコロイド表面電荷とは異符号にする
            if (Valency_counterion * Surface_charge[0] >= 0.0) {
                fprintf (stderr, "invalid (Valency_counterion) = %g\n",
                         Valency_counterion);
                exit_job (EXIT_FAILURE);
            }
        } else {
            dmy0 = Surface_charge[0];
            for (int i = 1; i < Component_Number; i++) {
                dmy1 = Surface_charge[i];
                //saltfreeのとき正負コロイド混合系には未対応
                if (dmy0 * dmy1 <= 0.0) {
                    fprintf (stderr, "invalid (number of ion spencies) = %d\n",
                             N_spec);
                    fprintf (stderr, "select salt in Add_salt\n");
                    exit_job (EXIT_FAILURE);
                }
            }
        }
        Reset_phi (phi);
        Make_phi_particle (phi, p);
#ifdef _OPENMP
#pragma omp parallel for default(none) reduction(+:volume_phi) shared(NPs, NZ_, phi)
#endif
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    volume_phi += (1.0 - phi[REALMODE_ARRAYINDEX(i, j, k)]);
                }
            }
        }
#ifdef _MPI
        dmy_mpi = alloc_1d_double (procs);

        ierr = MPI_Allgather (&volume_phi, 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);

        MPI_ERRORCHECK (ierr);
        volume_phi = 0.0;
        for (int n = 0; n < procs; n++) {
            volume_phi += dmy_mpi[n];
        }
        free_1d_double (dmy_mpi);
#endif
        Counterion_density = Counterion_number / (volume_phi * DX3);
        fprintf_single (stderr, "############################ initial state for counterions\n");
        fprintf_single (stderr, "# counterion_density %g\n", Counterion_density);
        fprintf_single (stderr, "# Debye length (1/(4 pi Bjerrum_length rho_c)^{1/2}) %g\n",
                        1.0 / sqrt (PI4 * SQ (Valency[0]) * Bjerrum_length * Counterion_density) );
        fprintf_single (stderr, "# radius of particle / Debye length  %g\n",
                        RADIUS * sqrt (PI4 * SQ (Valency[0]) * Bjerrum_length * Counterion_density) );
        if (Poisson_Boltzmann) {
            Set_uniform_ion_charge_density_nosalt (Concentration[0],
                                                   Total_solute, p);
            Set_Poisson_Boltzmann_ion_charge_density_nosalt (Concentration,
                                                             Total_solute, p);
            fprintf_single (stderr, "# initialized by Poisson-Boltzmann distribution\n");
            if (External_field) {
                Set_steadystate_ion_density (Concentration, p, jikan);
                fprintf_single (stderr, "# initialized by Poisson-Boltzmann distribution under external field\n");
            }
        } else {
            Set_uniform_ion_charge_density_nosalt (Concentration[0],
                                                   Total_solute, p);
            fprintf_single (stderr, "# initialized by uniform distribution\n");
        }
        fprintf_single (stderr, "############################\n");
    } else {
        if (Valency_positive_ion <= 0.0) {  //正イオン価数 > 0
            fprintf (stderr, "invalid (Valency_positive_ion) = %g\n",
                     Valency_positive_ion);
            exit_job (EXIT_FAILURE);
        }
        if (Valency_negative_ion >= 0.0) {  //負イオン価数 < 0
            fprintf (stderr, "invalid (Valency_negative_ion) = %g\n",
                     Valency_negative_ion);
            exit_job (EXIT_FAILURE);
        }
        fprintf_single (stderr, "############################ initial state for positive and negative ions\n");
        fprintf_single (stderr, "# radius of particle / Debye length  %g\n",
                        RADIUS / Debye_length);
        if (Poisson_Boltzmann) {
            if (Particle_Number == 1) {
                Set_Poisson_Boltzmann_ion_charge_density_salt (Concentration,
                                                               Total_solute, p);
            } else {
                Set_uniform_ion_charge_density_salt (Concentration,
                                                     Total_solute, p);
                Set_steadystate_ion_density (Concentration, p, jikan);
            }
            if (External_field != 0 && Particle_Number != 1) {
                fprintf_single (stderr, "# initialized by Poisson-Boltzmann distribution under external field\n");
            } else {
                fprintf_single (stderr, "# initialized by Poisson-Boltzmann distribution\n");
            }
        } else {
            Set_uniform_ion_charge_density_salt (Concentration, Total_solute, p);
            fprintf_single (stderr, "# initialized by uniform distribution\n");
        }
        fprintf_single (stderr, "############################\n");
    }
    Count_solute_each (Total_solute, Concentration, p, phi, up[0]);

}