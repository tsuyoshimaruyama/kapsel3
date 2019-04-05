/*!
  \file sp_3d_ns.cxx
  \author Y. Nakayama
  \date 2006/11/14
  \version 1.5
  \brief Main program file
 */
#include "sp_3d_ns.h"
#ifdef TIME_MEASURE
wall_timer evolM1_timer, evolM2_timer, evolM3_timer, evolM4_timer, evolM5_timer;
double evolM1_time, evolM2_time, evolM3_time, evolM4_time, evolM5_time;
wall_timer evolM6_timer, evolM7_timer, evolM8_timer, evolM9_timer,
evolM10_timer;
double evolM6_time, evolM7_time, evolM8_time, evolM9_time, evolM10_time;
wall_timer evolM11_timer, evolM12_timer, evolM13_timer, evolM14_timer,
evolM15_timer;
double evolM11_time, evolM12_time, evolM13_time, evolM14_time, evolM15_time;
wall_timer evolM16_timer, evolM17_timer, evolM18_timer, evolM19_timer,
evolM20_timer;
double evolM16_time, evolM17_time, evolM18_time, evolM19_time, evolM20_time;
#endif

void(*Time_evolution)(double **zeta, double uk_dc[DIM], double **f, Particle *p, CTime &jikan);

void Time_evolution_noparticle(double **zeta, double uk_dc[DIM], double **f, Particle *p, CTime &jikan) {
	const Index_range* ijk_range = ijk_range_two_third_filter;
	//const Index_range ijk_range[] = {
	//    {0, TRN_X - 1, 0, TRN_Y - 1, 0, 2 * TRN_Z - 1},
	//    {0, TRN_X - 1, NY - TRN_Y + 1, NY - 1, 0, 2 * TRN_Z - 1},
	//    {NX - TRN_X + 1, NX - 1, 0, TRN_Y - 1, 0, 2 * TRN_Z - 1},
	//    {NX - TRN_X + 1, NX - 1, NY - TRN_Y + 1, NY - 1, 0, 2 * TRN_Z - 1}};
	const int n_ijk_range = n_ijk_range_two_third_filter;
	// const int n_ijk_range = sizeof (ijk_range) / sizeof (Index_range);
	if (SW_EQ == Navier_Stokes) {
		NS_solver_slavedEuler(zeta, jikan, uk_dc, ijk_range, n_ijk_range, p);
	} else if (SW_EQ == Shear_Navier_Stokes) {
		NS_solver_slavedEuler_Shear_PBC(zeta, jikan, uk_dc, ijk_range, n_ijk_range, p, Shear_force);
	} else if (SW_EQ == Electrolyte) {
		NSsolute_solver_Euler(zeta, jikan, uk_dc, Concentration, p, ijk_range, n_ijk_range);
	}
}

void Time_evolution_hydro(double **zeta, double uk_dc[DIM], double **f, Particle *p, CTime &jikan) {
	int im;

#ifdef TIME_MEASURE
	evolM1_timer.start();
#endif
	// Update of Fluid Velocity Field
	Time_evolution_noparticle(zeta, uk_dc, f, p, jikan);
#ifdef TIME_MEASURE
	evolM1_time += evolM1_timer.stop();
#endif   

	if (Particle_Number >= 0) {
		if (FIX_CELL) { // time-dependent average pressure gradient
			for (int d = 0; d < DIM; d++) {
				if (FIX_CELLxyz[d]) {
					uk_dc[d] = 0.0;
				}
			}
		}

#ifdef TIME_MEASURE
		evolM2_timer.start();
#endif	
		for (int d = 0; d < (DIM - 1); d++) {
			Truncate_two_third_rule_ooura(zeta[d]);
		}
		Zeta_k2u(zeta, uk_dc, u);
#ifdef TIME_MEASURE
		evolM2_time += evolM2_timer.stop();
#endif

#ifdef TIME_MEASURE
		evolM3_timer.start();
#endif
#ifdef TIME_MEASURE
		evolM3_time += evolM3_timer.stop();
#endif

#ifdef TIME_MEASURE
		evolM4_timer.start();
#endif
#ifdef TIME_MEASURE
		evolM4_time += evolM4_timer.stop();
#endif

#ifdef TIME_MEASURE
		evolM5_timer.start();
#endif
#ifdef TIME_MEASURE
		evolM5_time += evolM5_timer.stop();
#endif
#ifdef TIME_MEASURE
		evolM6_timer.start();
#endif
#ifdef TIME_MEASURE
		evolM6_time += evolM6_timer.stop();
#endif

#ifdef TIME_MEASURE
		evolM7_timer.start();
#endif
		if (!Fixed_particle) {
			if (jikan.ts == 0) {
				MD_solver_position_Euler(p, jikan);
			} else {
				MD_solver_position_AB2(p, jikan);
			}
		}

#ifdef TIME_MEASURE
		evolM7_time += evolM7_timer.stop();
#endif

#ifdef TIME_MEASURE
		evolM8_timer.start();
#endif
		Reset_phi(phi);
		Reset_phi(phi_sum);
#ifdef TIME_MEASURE
		evolM8_time += evolM8_timer.stop();
#endif

#ifdef TIME_MEASURE
		evolM9_timer.start();
#endif
		Make_phi_particle_sum(phi, phi_sum, p);
#ifdef TIME_MEASURE
		evolM9_time += evolM9_timer.stop();
#endif

		if (SW_EQ == Electrolyte) {
			double * rescale_factor = new double[N_spec];
			Rescale_solute(rescale_factor, Total_solute, Concentration, p, phi, up[0]);
			delete[] rescale_factor;

			Reset_u(up);
			Calc_f_hydro_correct_precision(p, phi_sum, u, jikan);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(Local_Particle_Number, p, Particle_Number)
#endif
			for (int n = 0; n < Particle_Number; n++) {
				for (int d = 0; d < DIM; d++) {
					p[n].f_hydro1[d] = p[n].f_hydro[d];
					p[n].torque_hydro1[d] = p[n].torque_hydro[d];
				}
			}
			Make_Coulomb_force_x_on_fluid(f, p, Concentration, up[0], up[1], jikan);

			double dmy = jikan.dt_fluid * IRHO;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(NPs, NZ_, u, f, dmy) private(im)
#endif
			for (int i = 0; i < NPs[REAL][0]; i++) {
				for (int j = 0; j < NPs[REAL][1]; j++) {
					for (int k = 0; k < NPs[REAL][2]; k++) {
						im = REALMODE_ARRAYINDEX(i, j, k);
						u[0][im] += (f[0][im] * dmy);
						u[1][im] += (f[1][im] * dmy);
						u[2][im] += (f[2][im] * dmy);
					}
				}
			}
			Solenoidal_u(u);
		}
#ifdef TIME_MEASURE
		evolM10_timer.start();
#endif	
		// Calculation of hydrodynamic force
		Reset_u(up);


		Calc_f_hydro_correct_precision(p, phi_sum, u, jikan); //hydrodynamic force
#ifdef TIME_MEASURE
		evolM10_time += evolM10_timer.stop();
#endif
#ifdef TIME_MEASURE
		evolM11_timer.start();
#endif
#ifdef TIME_MEASURE
		evolM11_time += evolM11_timer.stop();
#endif

		if (!SW_JANUS_SLIP) {
#ifdef TIME_MEASURE
			evolM12_timer.start();
#endif
			if (!Fixed_particle) {
				if (jikan.ts == 0) {
					MD_solver_velocity_Euler(p, jikan);
				} else {
					MD_solver_velocity_AB2_hydro(p, jikan);
				}
			}
#ifdef TIME_MEASURE
			evolM12_time += evolM12_timer.stop();
#endif
		} else { // Self-Consistent slip force

			int slip_converge = 0;
			int slip_iter = 0;
			Make_particle_momentum_factor(u, p);
			Update_slip_particle_velocity(p, slip_iter); // initial particle velocity for slip profile
			while (!slip_converge) {
				Reset_u(up);
				Make_force_u_slip_particle(up, u, p, jikan);
				Solenoidal_u(up);
				Calc_f_slip_correct_precision(p, up, jikan); //slip force
				Add_f_particle(up, u); //up += u

				// Update particle velocity
				if (!Fixed_particle) {
					if (slip_iter == 0) {
						MD_solver_velocity_slip_iter(p, jikan, start_iter);
					} else {
						MD_solver_velocity_slip_iter(p, jikan, new_iter);
					}
				}

				if (PINNING) {
					Pinning(p);
				}
				slip_iter++;

				if (Slip_particle_convergence(p) < MAX_SLIP_TOL || slip_iter == MAX_SLIP_ITER) {
					slip_converge = 1;
					MD_solver_velocity_slip_iter(p, jikan, end_iter);
				} else {
					Update_slip_particle_velocity(p, slip_iter); // use new particle velocity for new slip profile
					MD_solver_velocity_slip_iter(p, jikan, reset_iter);
				}
			}//slip_convergence

			if (slip_iter == MAX_SLIP_ITER) {
				fprintf(stderr, "#Warning: increase MAX_SLIP_ITER (%d)\n", jikan.ts);
			}

			Swap_mem(u, up);
		} // slip 
#ifdef TIME_MEASURE
		evolM13_timer.start();
#endif
		if (kBT > 0. && SW_EQ != Electrolyte) {
			Add_random_force_thermostat(p, jikan);
		}
#ifdef TIME_MEASURE
		evolM13_time += evolM13_timer.stop();
#endif

		if (PINNING) {
			Pinning(p);
		}
		//Reset_phi_u(phi, up);
		//Make_phi_u_particle(phi, up, p);
#ifdef TIME_MEASURE
		evolM14_timer.start();
#endif
		Reset_u(up);
#ifdef TIME_MEASURE
		evolM14_time += evolM14_timer.stop();
#endif
#ifdef TIME_MEASURE
		evolM15_timer.start();
#endif
		Make_u_particle_sum(up, phi_sum, p);
#ifdef TIME_MEASURE
		evolM15_time += evolM15_timer.stop();
#endif
#ifdef TIME_MEASURE
		evolM16_timer.start();
#endif
		Make_f_particle_dt_sole(f, u, up, phi);
#ifdef TIME_MEASURE
		evolM16_time += evolM16_timer.stop();
#endif

#ifdef TIME_MEASURE
		evolM17_timer.start();
#endif
#ifdef TIME_MEASURE
		evolM17_time += evolM17_timer.stop();
#endif

#ifdef TIME_MEASURE
		evolM18_timer.start();
#endif
		Add_f_particle(u, f);
#ifdef TIME_MEASURE
		evolM18_time += evolM18_timer.stop();
#endif
#ifdef TIME_MEASURE
		evolM19_timer.start();
#endif
#ifdef TIME_MEASURE
		evolM19_time += evolM19_timer.stop();
#endif

		if (Shear_AC) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(NPs, NZ_, u, ucp) private(im)
#endif  
			for (int i = 0; i < NPs[REAL][0]; i++) {
				for (int j = 0; j < NPs[REAL][1]; j++) {
					for (int k = 0; k < NPs[REAL][2]; k++) {
						im = REALMODE_ARRAYINDEX(i, j, k);
						ucp[0][im] = u[0][im];
					}
				}
			}
		}

#ifdef TIME_MEASURE
		evolM20_timer.start();
#endif	
		U2zeta_k(zeta, uk_dc, u);
#ifdef TIME_MEASURE
		evolM20_time += evolM20_timer.stop();
#endif

	}

}

void Time_evolution_hydro_OBL(double **zeta, double uk_dc[DIM], double **f, Particle *p, CTime &jikan) {
	//AC
	Angular_Frequency = PI2*Shear_frequency;
	double ifreq = 1 / Angular_Frequency;
	double shear_amp = Shear_rate * ifreq;
	const Index_range* ijk_range = ijk_range_two_third_filter;
	const int n_ijk_range = n_ijk_range_two_third_filter;

	NS_solver_slavedEuler_Shear_OBL(zeta, jikan, uk_dc, ijk_range, n_ijk_range, p, Shear_force);

	if (Particle_Number >= 0) {
		if (FIX_CELL) { // time-dependent average pressure gradient
			for (int d = 0; d < DIM; d++) {
				if (FIX_CELLxyz[d]) {
					uk_dc[d] = 0.0;
				}
			}
		}

		Zeta_k2u_k_OBL(zeta, uk_dc, u);
		U_k2u(u);

		if (!Shear_AC) {
			Shear_rate_eff = Shear_rate;
		} else {
			Shear_rate_eff = Shear_rate * cos(Angular_Frequency*jikan.dt_fluid*jikan.ts);
		}
		degree_oblique += Shear_rate_eff*jikan.dt_fluid;

		if (degree_oblique >= 1.) {
			int flag = 0;
			Reset_U_OBL(ucp, u, flag);
			Swap_mem(u, ucp);
			degree_oblique -= 1.;
		} else if (degree_oblique < 0.) {
			int flag = 1;
			Reset_U_OBL(ucp, u, flag);
			Swap_mem(u, ucp);
			degree_oblique += 1.;
		}

		Update_K2_OBL();

		Copy_v3(ucp, u);

		Transform_obl_u(ucp, oblique2cartesian);


		// u   -> velocity field in oblique coordinates
		// ucp -> velocity fild in cartesian coordinates
		// Calc_shear_rate_eff();

		//AC
		if (!Shear_AC) Calc_shear_rate_eff();
		//End Deformation

		if (!Fixed_particle) {
			if (jikan.ts == 0) {
				MD_solver_position_Euler_OBL(p, jikan);
			} else {
				MD_solver_position_AB2_OBL(p, jikan);
			}
		}

		// Calculation of hydrodynamic force
		Reset_phi(phi);
		Reset_phi(phi_sum);
		//Reset_u(up);

		Make_phi_particle_sum_OBL(phi, phi_sum, p);

		Calc_f_hydro_correct_precision_OBL(p, phi_sum, ucp, jikan);

		Calc_Reynolds_shear_stress(ucp, Inertia_stress);


	}//end if(Particle_Number)

	if (!Fixed_particle) {// Update of Particle Velocity
		if (jikan.ts == 0) {
			MD_solver_velocity_Euler_OBL(p, jikan);
		} else {
			MD_solver_velocity_AB2_hydro_OBL(p, jikan);
		}
	}

	if (kBT > 0. && SW_EQ != Electrolyte) {
		Add_random_force_thermostat(p, jikan);
	}

	if (PINNING) {
		Pinning(p);
	}

	Reset_u(up);

	Make_u_particle_sum_OBL(up, phi_sum, p);

	Make_f_particle_dt_nonsole(f, ucp, up, phi);

	Transform_obl_u(f, cartesian2oblique);

	Add_f_particle(u, f);

	U2u_k(u);
	Solenoidal_uk_OBL(u);

	U_k2zeta_k_OBL(u, zeta, uk_dc);

}

inline void init_threads() {
#ifdef _OPENMP
	{
		int nthreads, tid, procs, maxt, inpar, dynamic, nested;
#pragma omp parallel private(nthreads, tid)
		{
			tid = omp_get_thread_num();
			if (tid == 0) {
				cerr << "#" << endl;
				cerr << "# OMP RUNTIME : " << endl;
				cerr << "# Number of processors        : " << omp_get_thread_num()
					<< endl;
				cerr << "# Number of threads           : " << omp_get_num_threads()
					<< endl;
				cerr << "# Max OMP threads             : " << omp_get_max_threads()
					<< endl;
#ifdef _FFT_IMKL
				cerr << "# Max MKL threads             : " << mkl_get_max_threads()
					<< endl;
#endif
				cerr << "# Dynamic thread enabled?     : " << omp_get_dynamic() << endl;
				cerr << "# Nested parallelism enabled? : " << omp_get_nested() << endl;
				cerr << "#" << endl;
			}
		}
	}
#endif
}

int main(int argc, char *argv[]) {
	double uk_dc[DIM];
	double **zeta;
	Particle *particles;

#if defined (_MPI)
	int flg = 0;
	int iprovided;
	MPI_Status status;

#ifdef _OPENMP
#pragma omp parallel
	{
#pragma omp master
		{
			ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &iprovided);
		}
#pragma omp barrier
	}
#else
	ierr = MPI_Init(&argc, &argv);
#endif
	ierr = MPI_Barrier(MPI_COMM_WORLD);
	MPI_ERRORCHECK(ierr);
	Set_MPI_initalize();
	flg = (procid == 0) ? 1 : 0;
#else
	procs = 1;
	procid = 0;

#endif
#ifndef NDEBUG 
	if (procid == root) {
		cerr << "###########################" << endl;
		cerr << "#  " << endl;
		cerr << "# built by debug mode " << endl;
		cerr << "# (NDEBUG is not defined in sp_3d_ns.h)" << endl;
		cerr << "#  " << endl;
		cerr << "###########################" << endl;
	}
#endif
	init_threads();

	wall_timer global_timer, block_timer;

	global_timer.start();
	block_timer.start();

	if (argc > 0) {
#ifdef _MPI
		if (procid != 0) {
			ierr = MPI_Recv(&flg, 1, MPI_INT, procid - 1, TAG(procid, procid - 1), MPI_COMM_WORLD, &status);
			MPI_ERRORCHECK(ierr);
		}
#endif
		file_get(argc, argv);
		Gourmet_file_io(In_udf, Out_udf, Sum_udf, Def_udf, Ctrl_udf, Res_udf);
#ifdef _MPI
		if (procid != procs - 1) {
			ierr = MPI_Send(&flg, 1, MPI_INT, procid + 1, TAG(procid + 1, procid), MPI_COMM_WORLD);
			MPI_ERRORCHECK(ierr);
		}
#endif
	}

#ifdef TIME_MEASURE
	evolM1_time = 0.0;
	evolM2_time = 0.0;
	evolM3_time = 0.0;
	evolM4_time = 0.0;
	evolM5_time = 0.0;
	evolM6_time = 0.0;
	evolM7_time = 0.0;
	evolM8_time = 0.0;
	evolM9_time = 0.0;
	evolM10_time = 0.0;
	evolM11_time = 0.0;
	evolM12_time = 0.0;
	evolM13_time = 0.0;
	evolM14_time = 0.0;
	evolM15_time = 0.0;
	evolM16_time = 0.0;
	evolM17_time = 0.0;
	evolM18_time = 0.0;
	evolM19_time = 0.0;
	evolM20_time = 0.0;
#endif

	// Main time evolution type 
	if (SW_EQ == Shear_Navier_Stokes_Lees_Edwards) {
		fprintf_single(stdout, "#Evolution type Shear_Navier_Stokes_Lees_Edwards\n");
		Time_evolution = Time_evolution_hydro_OBL;
	} else {
		if (SW_EQ == Navier_Stokes) {
			fprintf_single(stdout, "#Evolution type Navier_Stokes\n");
		} else if (SW_EQ == Shear_Navier_Stokes) {
			fprintf_single(stdout, "#Evolution type Shear_Navier_Stokes\n");
		} else if (SW_EQ == Electrolyte) {
			fprintf_single(stdout, "#Evolution type Electrolyte\n");
		}
		Time_evolution = Time_evolution_hydro;
	}
	zeta = (double **)malloc(sizeof(double *) * (DIM - 1));
	alloc_error_check(zeta);
	Range_division();
	Mem_alloc_var(zeta);

	MT_seed(GIVEN_SEED, 0);
	//MT_seed(RANDOM_SEED,0);
	Init_fft();
	Init_Transform_obl();

	static CTime jikan = { 0, 0.0, DT, DT*0.5, DT, DT*0.5 };

	particles = new Particle[Particle_Number];
	//particles = (Particle *)malloc(sizeof(Particle) * Particle_Number);

	if (Particle_Number > 0) {
#ifdef _MPI
		Build_Particle_Datatypes();
#endif

		MT_Init_Particle(particles);

		Init_Particle(particles);

		//if(SW_PT == chain){
		if ((SW_PT == chain) && !(DISTRIBUTION == user_specify)) {
			Init_Chain(particles);
		} else if ((SW_PT == rigid) && !(DISTRIBUTION == user_specify)) {
			Init_Rigid(particles);
		}
	}

	Init_output(particles);

	Init_zeta_k(zeta, uk_dc);

	Reset_phi_u(phi, up);
	Reset_phi(phi_sum);

	Make_phi_particle_sum(phi, phi_sum, particles);
	Make_u_particle_sum(up, phi_sum, particles);

	Zeta_k2u(zeta, uk_dc, u);
	Make_f_particle_dt_sole(f_particle, u, up, phi);
	Add_f_particle(u, f_particle);
	U2zeta_k(zeta, uk_dc, u);

	for (int d = 0; d < DIM; d++) {
		uk_dc[d] = .0;
	}

	if (SW_EQ == Electrolyte) {
		Init_rho_ion(Concentration, particles, jikan);
	}

	//  return EXIT_SUCCESS;

	Show_parameter(particles);
	Show_output_parameter();

	if ((SW_EQ == Shear_Navier_Stokes) || (SW_EQ == Shear_Navier_Stokes_Lees_Edwards)) {
		Mean_shear_stress(INIT, stderr, particles, jikan, Shear_rate_eff);
		degree_oblique = 0.0;
		p_info = (Particle_Info **)malloc(sizeof(Particle_Info *) * Particle_Number);
		for (int i = 0; i < Particle_Number; i++) {
			p_info[i] = (Particle_Info *)malloc(sizeof(Particle_Info) * NP_domain);
		}
	} else if (SW_EQ == Electrolyte) {
		Electrolyte_free_energy(INIT, stderr, particles, Concentration, jikan);
	}

	if (procid == root) {
		fprintf_single(stderr, "# Initialization time (s): %12.3f\n", block_timer.stop());
		block_timer.start();
	}

	//////////////////////////////////////////////////////////
	// Main Loop //
	int resumed_ts = 0;
	if (RESUMED) resumed_ts = last_ts;
	for (jikan.ts = resumed_ts; jikan.ts <= MSTEP; jikan.ts++) {

		int resumed_and_1st_loop = 0;

		if (RESUMED && (jikan.ts == resumed_ts)) {
			resumed_and_1st_loop = 1;
		}
		if (jikan.ts % GTS == 0) {

			if (!resumed_and_1st_loop) {
				if (SW_OUTFORMAT != OUT_NONE) {// Output field & particle data
					Output_open_frame();
					if (SW_EQ != Electrolyte) {
						Output_field_data(zeta, uk_dc, particles, jikan);
					} else if (SW_EQ == Electrolyte) {
						Output_charge_field_data(zeta, uk_dc, Concentration, particles, jikan);
					}
					Output_particle_data(particles, jikan);
					Output_close_frame();
				}
				if (SW_UDF) {// Output_UDF
					Output_udf(ufout, zeta, uk_dc, particles, jikan);
				}
				if (SW_EQ == Electrolyte) {
					Electrolyte_free_energy(SHOW, stderr, particles, Concentration, jikan);
				}
				if (jikan.ts != resumed_ts) {
					double block_time = block_timer.stop();
					double global_time = global_timer.stop();
					fprintf_single(stderr,
						"# Step: %9d/%9d\t  Block time (m): %8.3f\t Global time (m): %8.3f/%8.3f\n",
						jikan.ts,
						MSTEP,
						block_time / 60.0,
						global_time / 60.0,
						((double)(MSTEP - jikan.ts + 1)) / ((double)GTS) * block_time / 60.0);
					fflush(stderr);
					block_timer.start();
				}
			}
		}

		Time_evolution(zeta, uk_dc, f_particle, particles, jikan);

		jikan.time += jikan.dt_fluid;

		if (SW_EQ == Shear_Navier_Stokes) {
			Shear_rate_eff = Update_strain(Shear_strain_realized, jikan, zeta, uk_dc, u);
			Mean_shear_stress(SHOW, stderr, particles, jikan, Shear_rate_eff);
		} else if (SW_EQ == Shear_Navier_Stokes_Lees_Edwards) {
			Shear_strain_realized += Shear_rate_eff * jikan.dt_fluid;
			Mean_shear_stress(SHOW, stderr, particles, jikan, Shear_rate_eff);
		}

		if (jikan.ts == MSTEP) {
			Save_Restart_udf(zeta, uk_dc, particles, jikan, Concentration);
		}
		if (resumed_and_1st_loop) {
			Force_restore_parameters(zeta, uk_dc, particles, jikan, Concentration);
			delete ufin;
			fprintf_single(stderr, "############################ Parameters are restored.\n");
		}
	} //end of main loop

#ifdef TIME_MEASURE
	fprintf_single(stderr, "# evolM1_time(s)=%5.3f\n", evolM1_time);
	fprintf_single(stderr, "# evolM2_time(s)=%5.3f\n", evolM2_time);
	fprintf_single(stderr, "# evolM3_time(s)=%5.3f\n", evolM3_time);
	fprintf_single(stderr, "# evolM4_time(s)=%5.3f\n", evolM4_time);
	fprintf_single(stderr, "# evolM5_time(s)=%5.3f\n", evolM5_time);
	fprintf_single(stderr, "# evolM6_time(s)=%5.3f\n", evolM6_time);
	fprintf_single(stderr, "# evolM7_time(s)=%5.3f\n", evolM7_time);
	fprintf_single(stderr, "# evolM8_time(s)=%5.3f\n", evolM8_time);
	fprintf_single(stderr, "# evolM9_time(s)=%5.3f\n", evolM9_time);
	fprintf_single(stderr, "# evolM10_time(s)=%5.3f\n", evolM10_time);
	fprintf_single(stderr, "# evolM11_time(s)=%5.3f\n", evolM11_time);
	fprintf_single(stderr, "# evolM12_time(s)=%5.3f\n", evolM12_time);
	fprintf_single(stderr, "# evolM13_time(s)=%5.3f\n", evolM13_time);
	fprintf_single(stderr, "# evolM14_time(s)=%5.3f\n", evolM14_time);
	fprintf_single(stderr, "# evolM15_time(s)=%5.3f\n", evolM15_time);
	fprintf_single(stderr, "# evolM16_time(s)=%5.3f\n", evolM16_time);
	fprintf_single(stderr, "# evolM17_time(s)=%5.3f\n", evolM17_time);
	fprintf_single(stderr, "# evolM18_time(s)=%5.3f\n", evolM18_time);
	fprintf_single(stderr, "# evolM19_time(s)=%5.3f\n", evolM19_time);
	fprintf_single(stderr, "# evolM20_time(s)=%5.3f\n", evolM20_time);
	fprintf_single(stderr, "\n");
	double particle_time =
		evolM2_time + evolM3_time + evolM4_time + evolM5_time + evolM6_time +
		evolM7_time + evolM8_time + evolM9_time + evolM10_time + evolM11_time;
	fprintf_single(stderr, "# evolM particle solving time(s)=%5.3f\n", particle_time);
	double fluid_time = evolM1_time + evolM12_time + evolM13_time +
		evolM14_time + evolM15_time + evolM16_time +
		evolM17_time + evolM18_time + evolM19_time +
		evolM20_time;
	fprintf_single(stderr, "# evolM fluid solving time(s)=%5.3f\n", fluid_time);
	double loop_time =
		evolM1_time + evolM2_time + evolM3_time + evolM4_time + evolM5_time +
		evolM6_time + evolM7_time + evolM8_time + evolM9_time + evolM10_time +
		evolM11_time + evolM12_time + evolM13_time + evolM14_time +
		evolM15_time + evolM16_time + evolM17_time + evolM18_time +
		evolM19_time + evolM20_time;
	fprintf_single(stderr, "# evolM loop time(s)=%5.3f\n", loop_time);
	fprintf_single(stderr, "\n");
#endif

	if (SW_UDF) {
		if (procid == root) {
			ufout->write();
		}
	}
	delete ufout;
	fprintf_single(stderr, "#%s end.\n", Out_udf);

	//Always write restart file
	if (procid == root) {
		ufres->write();
	}
	delete ufres;
	fprintf_single(stderr, "#%s end.\n", Res_udf);

	Free_output();

	double global_time = global_timer.stop();
	fprintf_single(stderr, "#Simulation has ended!\n");
	fprintf_single(stderr, "#Total Running Time (s): %10.2f\n", global_time);
	fprintf_single(stderr, "#                   (m): %10.2f\n", global_time / 60.0);
	fprintf_single(stderr, "#                   (h): %10.2f\n", global_time / 3600.0);

	global_time /= (double)(MSTEP - resumed_ts + 1);
	fprintf_single(stderr, "#Average Step Time  (s): %10.2f\n", global_time);
	fprintf_single(stderr, "#                   (m): %10.2f\n", global_time / 60.0);
	fprintf_single(stderr, "#                   (h): %10.2f\n", global_time / 3600.0);
#if ((defined _OPENMP) || (defined _MPI))	
	Dfti_finalize();
#endif
	Free_Transform_obl();

	if (SW_EQ == Shear_Navier_Stokes_Lees_Edwards) {
		for (int i = 0; i < Particle_Number; i++) {
			free(p_info[i]);
		}
		free(p_info);
	}

	Mem_free_var(zeta);
	Free_fft();
	delete[] particles;
#ifdef _MPI
	Set_MPI_finalize();
#endif

	return EXIT_SUCCESS;
}

