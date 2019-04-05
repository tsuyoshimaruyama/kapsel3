/*!
  \file make_phi.cxx
  \brief Routines to compute smooth profile and grid particle properties
  \author T. Iwashita
  \date 2006/10/20
  \version 1.3
 */

#include "make_phi.h"

void(*Angular2v)(const double *omega, const double *r, double *v);

/////////////
void Make_surface_normal(double **surface_normal, const Particle *p) {
	int im;
	int dmy_mesh[DIM];
	int r_mesh[DIM];
	int x_int[DIM];
	double dmy;
	double dmy_r;
	double ir;
	double r[DIM];
	double residue[DIM];
	int sw_in_cell;
	double xp[DIM];
	for (int d = 0; d < DIM; d++) {
		Reset_phi(surface_normal[d]);
	}
#ifdef _OPENMP
	Reset_mesh(tmp_buffer_dim[0]);
	Reset_mesh(tmp_buffer_dim[1]);
	Reset_mesh(tmp_buffer_dim[2]);
#else
	tmp_buffer_dim = surface_normal;
#endif
#ifdef _OPENMP
	//#pragma omp parallel default(none) shared(Local_Particle_Number, Sekibun_cell, NP_domain, THREADNUM, Ns, NPs, NZ_, L, DX, HXI, RADIUS, p, surface_normal, tmp_buffer_dim, sw_in_cell)
	{
		//#pragma omp for private(r, ir, r_mesh, residue, xp, x_int, dmy_mesh, dmy_r, dmy, im)
#endif
		for (int n = 0; n < Local_Particle_Number; n++) {
			for (int d = 0; d < DIM; d++) {
				xp[d] = p[n].x[d];
				assert((xp[d] >= 0) && (xp[d] < L[d]));
			}
			sw_in_cell = Particle_cell(xp, DX, x_int, residue); // {1, 0} が返ってくる
			sw_in_cell = 1;
			for (int mesh = 0; mesh < NP_domain; mesh++) {
				Relative_coord(Sekibun_cell[mesh], x_int, residue, sw_in_cell, Ns, DX, r_mesh, r);
				if (Range_coord(r_mesh, dmy_mesh)) {
#ifdef _OPENMP
					im = REALMODE_ARRAYINDEX_MESH(omp_get_thread_num(), dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
#else
					im = REALMODE_ARRAYINDEX(dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
#endif
					dmy_r = sqrt(SQ(r[0]) + SQ(r[1]) + SQ(r[2]));
					dmy = ABS(dmy_r - RADIUS);
					if (dmy < HXI) {
						ir = 1.0 / dmy_r;
						for (int d = 0; d < DIM; d++) {
							tmp_buffer_dim[d][im] += r[d] * ir;
						}
					}
				}
			}
		}
#ifdef _OPENMP
		//#pragma omp barrier
		//#pragma omp for
		for (int i = 0; i < NPs[REAL][0]; i++) {
			for (int t = 0; t < THREADNUM; t++) {
				for (int d = 0; d < DIM; d++) {
					for (int j = 0; j < NPs[REAL][1]; j++) {
						for (int k = 0; k < NPs[REAL][2]; k++) {
							surface_normal[d][REALMODE_ARRAYINDEX(i, j, k)] += tmp_buffer_dim[d][REALMODE_ARRAYINDEX_MESH(t, i, j, k)];
						}
					}
				}
			}
		}
	}
#endif
}

/////////////

inline void Make_rho_field_primitive(double *phi, Particle *p, const double &dx,
	const int &np_domain, int **sekibun_cell, int Nlattice[DIM]) {
	int im;
	int r_mesh[DIM];
	int dmy_mesh[DIM];
	int x_int[DIM];
	double drho;
	double r[DIM];
	double residue[DIM];
	double x[DIM];
	double xp[DIM];
#ifdef _OPENMP
	Reset_mesh(tmp_buffer1);
#else
	tmp_buffer1 = phi;
#endif
#ifdef _OPENMP
	//#pragma omp parallel default(none) shared(Local_Particle_Number, sekibun_cell, np_domain, THREADNUM, Ns, NPs, NZ_, L, dx, RHO, RHO_particle, RADIUS, p, phi, tmp_buffer1)
	{
		//#pragma omp for private(drho, r, r_mesh, residue, x, xp, x_int, dmy_mesh, im)
#endif
		for (int n = 0; n < Local_Particle_Number; n++) {
			drho = RHO_particle[p[n].spec] - RHO;
			for (int d = 0; d < DIM; d++) {
				xp[d] = p[n].x[d];
			}
			int sw_in_cell = Particle_cell(xp, dx, x_int, residue); // {1, 0} が返ってくる
			sw_in_cell = 1;
			for (int mesh = 0; mesh < np_domain; mesh++) {
				Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Ns, dx, r_mesh, r);
				if (Range_coord(r_mesh, dmy_mesh)) {
#ifdef _OPENMP
					im = REALMODE_ARRAYINDEX_MESH(omp_get_thread_num(), dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
#else
					im = REALMODE_ARRAYINDEX(dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
#endif
					for (int d = 0; d < DIM; d++) {
						x[d] = r_mesh[d] * dx;
					}
					tmp_buffer1[im] += Phi(Distance(x, xp)) * drho;
				}
			}
		}
#ifdef _OPENMP
		//#pragma omp barrier
		//#pragma omp for
		for (int i = 0; i < NPs[REAL][0]; i++) {
			for (int t = 0; t < THREADNUM; t++) {
				for (int j = 0; j < NPs[REAL][1]; j++) {
					for (int k = 0; k < NPs[REAL][2]; k++) {
						phi[REALMODE_ARRAYINDEX(i, j, k)] += tmp_buffer1[REALMODE_ARRAYINDEX_MESH(t, i, j, k)];
					}
				}
			}
		}
#endif
#ifdef _OPENMP
		//#pragma omp barrier
		//#pragma omp for
#endif
		for (int i = 0; i < NPs[REAL][0]; i++) {
			for (int j = 0; j < NPs[REAL][1]; j++) {
				for (int k = 0; k < NPs[REAL][2]; k++) {
					phi[REALMODE_ARRAYINDEX(i, j, k)] += RHO;
				}
			}
		}
#ifdef _OPENMP
	}
#endif
}

void Make_rho_field(double *phi, Particle *p) {
	int *nlattice;
	nlattice = Ns;
	Make_rho_field_primitive(phi, p, DX, NP_domain, Sekibun_cell, nlattice);
}

//
inline double janus_geometry(const Particle &p, const double normal[DIM]) {
	double body_normal[DIM];
	double cos_theta;

	rigid_body_rotation(body_normal, normal, p.q, SPACE2BODY);
	if (janus_axis[p.spec] == x_axis) {
		cos_theta = body_normal[0];
	} else if (janus_axis[p.spec] == y_axis) {
		cos_theta = body_normal[1];
	} else if (janus_axis[p.spec] == z_axis) {
		cos_theta = body_normal[2];
	} else if (janus_axis[p.spec] == no_axis) {
		cos_theta = 1.0;
	} else {
		fprintf_single(stderr, "Error: %d not a janus particle\n", p.spec);
		exit_job(EXIT_FAILURE);
	}

	return ((cos_theta >= 0.0) ? 1.0 : -1.0);
}

inline void Make_phi_u_primitive(double *phi, double **up, Particle *p, const int &SW_UP, const double &dx
	, const int &np_domain, int **sekibun_cell, const int Nlattice[DIM], const double radius = RADIUS) {

	double xp[DIM], vp[DIM], omega_p[DIM];
	int x_int[DIM];
	double residue[DIM];
	int sw_in_cell;
	int r_mesh[DIM];
	int dmy_mesh[DIM];
	double r[DIM];
	double x[DIM];
	double dmy;
	double dmy_phi;
	double v_rot[DIM];
	int im;
#ifdef _MPI
#ifdef _OPENMP
	Reset_mesh(tmp_buffer1);
	Reset_mesh(tmp_buffer_dim[0]);
	Reset_mesh(tmp_buffer_dim[1]);
	Reset_mesh(tmp_buffer_dim[2]);
#else
	tmp_buffer1 = phi;
	tmp_buffer_dim = up;
#endif
#endif
	for (int n = 0; n < Local_Particle_Number; n++) {
		for (int d = 0; d < DIM; d++) {
			xp[d] = p[n].x[d];
			vp[d] = p[n].v[d];
			omega_p[d] = p[n].omega[d];
		}
		sw_in_cell = Particle_cell(xp, dx, x_int, residue);// {1,0} が返ってくる
		sw_in_cell = 1;
		for (int mesh = 0; mesh < NP_domain; mesh++) {
			Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
#ifdef _MPI
			if (Range_coord(r_mesh, dmy_mesh)) {
#endif
				//dmy = 0.;
				for (int d = 0; d < DIM; d++) {
					x[d] = r_mesh[d] * dx;
					//dmy += SQ(r[d]);
				}
				dmy = Distance(x, xp);
				//dmy = sqrt(dmy);
				dmy_phi = Phi(dmy, radius);
#ifdef _MPI
#ifdef _OPENMP
				im = REALMODE_ARRAYINDEX_MESH(omp_get_thread_num(), dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
#else
				im = REALMODE_ARRAYINDEX(dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
#endif
				tmp_buffer1[im] += dmy_phi;
#else
				im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
#ifdef _OPENMP
#pragma omp atomic
#endif
				phi[im] += dmy_phi;
#endif

				if (SW_UP) {
					Angular2v(omega_p, r, v_rot);
#ifdef _MPI
					tmp_buffer_dim[0][im] += ((vp[0] + v_rot[0]) * dmy_phi);
					tmp_buffer_dim[1][im] += ((vp[1] + v_rot[1]) * dmy_phi);
					tmp_buffer_dim[2][im] += ((vp[2] + v_rot[2]) * dmy_phi);
#else
#ifdef _OPENMP
#pragma omp atomic
#endif
					up[0][im] += ((vp[0] + v_rot[0]) * dmy_phi);
#ifdef _OPENMP
#pragma omp atomic
#endif
					up[1][im] += ((vp[1] + v_rot[1]) * dmy_phi);
#ifdef _OPENMP
#pragma omp atomic
#endif
					up[2][im] += ((vp[2] + v_rot[2]) * dmy_phi);
#endif

				}
#ifdef _MPI
			}
#endif
		}
	}
#ifdef _MPI
	//#ifdef _OPENMP
	//#pragma omp barrier
	//#pragma omp for
	for (int i = 0; i < NPs[REAL][0]; i++) {
		for (int t = 0; t < THREADNUM; t++) {
			for (int j = 0; j < NPs[REAL][1]; j++) {
				for (int k = 0; k < NPs[REAL][2]; k++) {
					phi[REALMODE_ARRAYINDEX(i, j, k)] += tmp_buffer1[REALMODE_ARRAYINDEX_MESH(t, i, j, k)];
				}
			}
		}
	}
	//#pragma omp barrier
	//#pragma omp for
	for (int i = 0; i < NPs[REAL][0]; i++) {
		for (int t = 0; t < THREADNUM; t++) {
			for (int d = 0; d < DIM; d++) {
				for (int j = 0; j < NPs[REAL][1]; j++) {
					for (int k = 0; k < NPs[REAL][2]; k++) {
						up[d][REALMODE_ARRAYINDEX(i, j, k)] += tmp_buffer_dim[d][REALMODE_ARRAYINDEX_MESH(t, i, j, k)];
					}
				}
			}
		}
	}
	//#endif
#endif

		// koba code //
	if (SW_UP) {
		double idmy_phi;
#ifdef _OPENMP
#pragma omp for private(im, idmy_phi)
#endif
		for (int i = 0; i < NPs[REAL][0]; i++) {
			for (int j = 0; j < NPs[REAL][1]; j++) {
				for (int k = 0; k < NPs[REAL][2]; k++) {
					im = REALMODE_ARRAYINDEX(i, j, k);
					if (phi[im] > 1.0) {
						idmy_phi = 1.0 / phi[im];
						up[0][im] *= idmy_phi;
						up[1][im] *= idmy_phi;
						up[2][im] *= idmy_phi;
						phi[im] = 1.0;
					}
				}
			}
		}
	}
}

inline void Make_phi_particle_sum_primitive(double *phi, double *phi_sum, Particle *p,
	const double &dx, const int &np_domain, int **sekibun_cell, const int Nlattice[DIM], const double radius) {
	int im;
#ifdef _MPI
#ifdef _OPENMP
	Reset_mesh(tmp_buffer1);
#else
	tmp_buffer1 = phi;
#endif
#endif
#ifdef _OPENMP
#pragma omp parallel for private(im)
#endif
	for (int n = 0; n < Local_Particle_Number; n++) {
		double xp[DIM];
		for (int d = 0; d < DIM; d++) xp[d] = p[n].x[d];

		int x_int[DIM];
		double residue[DIM];
		int sw_in_cell = Particle_cell(xp, dx, x_int, residue);
		sw_in_cell = 1;

		int r_mesh[DIM];
		double dmy, dmy_phi;
		double r[DIM], x[DIM];
#ifdef _MPI
		int dmy_mesh[DIM];
#endif
		for (int mesh = 0; mesh < NP_domain; mesh++) {
			Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
#ifdef _MPI
			if (Range_coord(r_mesh, dmy_mesh)) {
#endif
				for (int d = 0; d < DIM; d++) x[d] = r_mesh[d] * DX;
#ifdef _MPI
#ifdef _OPENMP
				im = REALMODE_ARRAYINDEX_MESH(omp_get_thread_num(), dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
#else
				im = REALMODE_ARRAYINDEX(dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
#endif
#endif
				dmy = Distance(x, xp);
				dmy_phi = Phi(dmy, radius);
#ifdef _MPI
#ifdef _OPENMP
#pragma omp atomic
#endif
				tmp_buffer1[im] += dmy_phi;
#else
#ifdef _OPENMP
#pragma omp atomic
#endif
				phi_sum[(r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2]] += dmy_phi;
#endif
#ifdef _MPI
			}
#endif
		}
	}

#ifdef _MPI
#ifdef _OPENMP
#pragma omp barrier
#pragma omp for
#endif
	for (int i = 0; i < NPs[REAL][0]; i++) {
		for (int t = 0; t < THREADNUM; t++) {
			for (int j = 0; j < NPs[REAL][1]; j++) {
				for (int k = 0; k < NPs[REAL][2]; k++) {
					phi_sum[REALMODE_ARRAYINDEX(i, j, k)] += tmp_buffer1[REALMODE_ARRAYINDEX_MESH(t, i, j, k)];
				}
			}
		}
	}
#endif
	for (int i = 0; i < NPs[REAL][0]; i++) {
		for (int j = 0; j < NPs[REAL][1]; j++) {
			for (int k = 0; k < NPs[REAL][2]; k++) {
				im = REALMODE_ARRAYINDEX(i, j, k);
				phi[im] = MIN(phi_sum[im], 1.0);
			}
		}
	}
}

inline void Make_phi_particle_sum_primitive_OBL(double *phi, double *phi_sum, Particle *p, const double &dx,
	const int &np_domain, int **sekibun_cell, const int Nlattice[DIM], const double radius) {
	int im;

	for (int n = 0; n < Particle_Number; n++) {
		double xp[DIM];
		for (int d = 0; d < DIM; d++) xp[d] = p[n].x[d];

		int x_int[DIM];
		double residue[DIM];
		int sw_in_cell = Particle_cell(xp, dx, x_int, residue);
		sw_in_cell = 1;

		/*
		double mesh_width = (NY / yprocs) * DX;
		double max_radius = radius + HXI + 2.0 * DX;
		double coord_high = xp[1] + max_radius;
		int pbc_flip_high = 1;
		if (coord_high > (double)(NY*DX)) {
			coord_high -= (double)(NY*DX);
			pbc_flip_high = 0;
		}
		int high_yid = (int)(coord_high / mesh_width);
		double coord_low = xp[1] - max_radius;
		int pbc_flip_low = 1;
		if (coord_low < 0.0) {
			coord_low += (double)(NY*DX);
			pbc_flip_low = 0;
		}
		int low_yid = (int)(coord_low / mesh_width);

		if (pbc_flip_high == 1 && pbc_flip_low == 1) {
			if (yid >= pbc_flip_low && yid <= pbc_flip_high) {
				;//blanck
			} else {
				continue;
			}
		} else if (pbc_flip_high == 0 && pbc_flip_low == 1) {
			if (yid <= pbc_flip_low && yid >= pbc_flip_high) {
				;//blanck
			} else {
				continue;
			}
		} else if (pbc_flip_high == 1 && pbc_flip_low == 0) {
			if (yid <= pbc_flip_low && yid >= pbc_flip_high) {
				;//blanck
			} else {
				continue;
			}
		}*/

		int sign;
		int r_mesh[DIM];
		int dmy_mesh[DIM];
		double dmy, dmy_phi;
		double r[DIM], x[DIM];
		for (int mesh = 0; mesh < NP_domain; mesh++) {
			sign = Relative_coord_check_stepover_Y(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
			if (Range_coord(r_mesh, dmy_mesh)) {
				for (int d = 0; d < DIM; d++) x[d] = r_mesh[d] * dx;

				im = REALMODE_ARRAYINDEX(dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);

				dmy = Distance_OBL(x, xp);
				dmy_phi = Phi(dmy, radius);
				phi_sum[im] += dmy_phi;
			}
		}
	}

	for (int i = 0; i < NPs[REAL][0]; i++) {
		for (int j = 0; j < NPs[REAL][1]; j++) {
			for (int k = 0; k < NPs[REAL][2]; k++) {
				im = REALMODE_ARRAYINDEX(i, j, k);
				phi[im] = MIN(phi_sum[im], 1.0);
			}
		}
	}
}

inline void Make_u_particle_sum_primitive(double **up,
	double const* phi_sum, Particle *p, const double &dx, const int &np_domain,
	int const* const* sekibun_cell, const int Nlattice[DIM], const double radius) {
#ifdef _MPI
#ifdef _OPENMP
	Reset_mesh(tmp_buffer_dim[0]);
	Reset_mesh(tmp_buffer_dim[1]);
	Reset_mesh(tmp_buffer_dim[2]);
#else
	tmp_buffer_dim = up;
#endif
#endif
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int n = 0; n < Local_Particle_Number; n++) {
		double xp[DIM], vp[DIM], omega_p[DIM];
		for (int d = 0; d < DIM; d++) {
			xp[d] = p[n].x[d];
			vp[d] = p[n].v[d];
			omega_p[d] = p[n].omega[d];
		}

		int im, sw_in_cell;
		int x_int[DIM], r_mesh[DIM];
		int dmy_mesh[DIM];
		double residue[DIM], r[DIM], x[DIM], v_rot[DIM];
		double dmy, dmy_phi;
		sw_in_cell = Particle_cell(xp, dx, x_int, residue);
		sw_in_cell = 1;
		for (int mesh = 0; mesh < NP_domain; mesh++) {
			Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
#ifdef _MPI
			if (Range_coord(r_mesh, dmy_mesh)) {
#endif
				for (int d = 0; d < DIM; d++) x[d] = r_mesh[d] * dx;
#ifdef _MPI
#ifdef _OPENMP
				im = REALMODE_ARRAYINDEX_MESH(omp_get_thread_num(), dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
#else
				im = REALMODE_ARRAYINDEX(dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
#endif
#else
				im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
#endif
				dmy = Distance(x, xp);
				dmy_phi = Phi(dmy, radius) / MAX(phi_sum[im], 1.0);

				Angular2v(omega_p, r, v_rot);
#ifdef _MPI
				tmp_buffer_dim[0][im] += ((vp[0] + v_rot[0]) * dmy_phi);
				tmp_buffer_dim[1][im] += ((vp[1] + v_rot[1]) * dmy_phi);
				tmp_buffer_dim[2][im] += ((vp[2] + v_rot[2]) * dmy_phi);
#else
#ifdef _OPENMP
#pragma omp atomic
#endif
				up[0][im] += ((vp[0] + v_rot[0]) * dmy_phi);
#ifdef _OPENMP
#pragma omp atomic
#endif
				up[1][im] += ((vp[1] + v_rot[1]) * dmy_phi);
#ifdef _OPENMP
#pragma omp atomic
#endif
				up[2][im] += ((vp[2] + v_rot[2]) * dmy_phi);
#endif
#ifdef _MPI
			}
#endif
		}
	}

#ifdef _MPI
#ifdef _OPENMP
	//#pragma omp barrier
	//#pragma omp for
	for (int i = 0; i < NPs[REAL][0]; i++) {
		for (int t = 0; t < THREADNUM; t++) {
			for (int d = 0; d < DIM; d++) {
				for (int j = 0; j < NPs[REAL][1]; j++) {
					for (int k = 0; k < NPs[REAL][2]; k++) {
						up[d][REALMODE_ARRAYINDEX(i, j, k)] += tmp_buffer_dim[d][REALMODE_ARRAYINDEX_MESH(t, i, j, k)];
					}
				}
			}
		}
	}
#endif
#endif
}

inline void Make_u_particle_sum_primitive_OBL(double **up, double const* phi_sum, Particle* p,
	const double &dx, const int &np_domain, int const* const* sekibun_cell, const int Nlattice[DIM], const double radius) {

	for (int n = 0; n < Particle_Number; n++) {
		double xp[DIM], vp[DIM], omega_p[DIM];
		for (int d = 0; d < DIM; d++) {
			xp[d] = p[n].x[d];
			vp[d] = p[n].v[d];
			omega_p[d] = p[n].omega[d];
		}

		int x_int[DIM];
		double residue[DIM];
		int sw_in_cell = Particle_cell(xp, dx, x_int, residue);
		sw_in_cell = 1;

		/*
		double mesh_width = (NY / yprocs) * DX;
		double max_radius = radius + HXI + 2.0 * DX;
		double coord_high = xp[1] + max_radius;
		int pbc_flip_high = 1;
		if (coord_high > (double)(NY*DX)) {
			coord_high -= (double)(NY*DX);
			pbc_flip_high = 0;
		}
		int high_yid = (int)(coord_high / mesh_width);
		double coord_low = xp[1] - max_radius;
		int pbc_flip_low = 1;
		if (coord_low < 0.0) {
			coord_low += (double)(NY*DX);
			pbc_flip_low = 0;
		}
		int low_yid = (int)(coord_low / mesh_width);

		if (pbc_flip_high == 1 && pbc_flip_low == 1) {
			if (yid >= pbc_flip_low && yid <= pbc_flip_high) {
				;//blanck
			} else {
				continue;
			}
		} else if (pbc_flip_high == 0 && pbc_flip_low == 1) {
			if (yid <= pbc_flip_low && yid >= pbc_flip_high) {
				;//blanck
			} else {
				continue;
			}
		} else if (pbc_flip_high == 1 && pbc_flip_low == 0) {
			if (yid <= pbc_flip_low && yid >= pbc_flip_high) {
				;//blanck
			} else {
				continue;
			}
		}
		*/
		int im, sign;
		int r_mesh[DIM];
		int dmy_mesh[DIM];
		double dmy, dmy_phi;
		double r[DIM], x[DIM], v_rot[DIM];
		for (int mesh = 0; mesh < NP_domain; mesh++) {
			sign = Relative_coord_check_stepover_Y(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);

			if (Range_coord(r_mesh, dmy_mesh)) {
				for (int d = 0; d < DIM; d++) x[d] = r_mesh[d] * dx;

				im = REALMODE_ARRAYINDEX(dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);

				dmy = Distance_OBL(x, xp);
				dmy_phi = Phi(dmy, radius) / MAX(phi_sum[im], 1.0);

				Angular2v(omega_p, r, v_rot);

				up[0][im] += ((vp[0] - sign*Shear_rate_eff*L_particle[1] + v_rot[0]) * dmy_phi);
				up[1][im] += (vp[1] + v_rot[1]) * dmy_phi;
				up[2][im] += (vp[2] + v_rot[2]) * dmy_phi;
			}
		}
	}
}

inline void Make_phi_u_primitive_OBL(double *phi, double **up, Particle *p, const int &SW_UP, const double &dx
	, const int &np_domain, int **sekibun_cell, const int Nlattice[DIM], const double radius = RADIUS) {
	double xp[DIM], vp[DIM], omega_p[DIM];
	int x_int[DIM];
	double residue[DIM];
	int sw_in_cell;
	int r_mesh[DIM];
	int dmy_mesh[DIM];
	double r[DIM];
	double x[DIM];
	double dmy;
	double dmy_phi;
	double v_rot[DIM];
	int sign;
	int im;
#ifdef _OPENMP
	Reset_mesh(tmp_buffer1);
	Reset_mesh(tmp_buffer_dim[0]);
	Reset_mesh(tmp_buffer_dim[1]);
	Reset_mesh(tmp_buffer_dim[2]);
#else
	tmp_buffer1 = phi;
	tmp_buffer_dim = up;
#endif
	//#ifdef _OPENMP
	//#pragma omp parallel default(none) shared(Local_Particle_Number, sekibun_cell, np_domain, THREADNUM, Ns, NPs, NZ_, dx, RADIUS, Angular2v, p, phi, up, tmp_buffer1, tmp_buffer_dim, sw_in_cell, sign, Nlattice, dmy, SW_UP, Shear_rate_eff, L_particle)
	{
		//#pragma omp for private(r, r_mesh, residue, x, xp, x_int, vp, v_rot, omega_p, dmy_mesh, dmy_phi, im)
		//#endif
				//for (int n = 0; n < Local_Particle_Number; n++) {
		for (int n = 0; n < Particle_Number; n++) {
			for (int d = 0; d < DIM; d++) {
				xp[d] = p[n].x[d];
				vp[d] = p[n].v[d];
				omega_p[d] = p[n].omega[d];
			}
			sw_in_cell = Particle_cell(xp, dx, x_int, residue);// {1,0} が返ってくる
			sw_in_cell = 1;
			for (int mesh = 0; mesh < np_domain; mesh++) {
				sign = Relative_coord_check_stepover_Y(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
				if (Range_coord(r_mesh, dmy_mesh)) {
					dmy = 0.;
					for (int d = 0; d < DIM; d++) {
						x[d] = r_mesh[d] * dx;
						dmy += SQ(r[d]);
					}
					dmy = sqrt(dmy);

					dmy_phi = Phi(dmy, radius);
					//#ifdef _OPENMP
					//                    im = REALMODE_ARRAYINDEX_MESH(omp_get_thread_num(), dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
					//#else
					im = REALMODE_ARRAYINDEX(dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
					//#endif
					tmp_buffer1[im] += dmy_phi;

					if (SW_UP) {
						Angular2v(omega_p, r, v_rot);
						tmp_buffer_dim[0][im] += ((vp[0] - sign*Shear_rate_eff*L_particle[1] + v_rot[0])*dmy_phi);
						tmp_buffer_dim[1][im] += (vp[1] + v_rot[1])*dmy_phi;
						tmp_buffer_dim[2][im] += (vp[2] + v_rot[2])*dmy_phi;
					}
				}
			}
		}
		//#ifdef _OPENMP
		//#pragma omp barrier
		//#pragma omp for
		for (int i = 0; i < NPs[REAL][0]; i++) {
			for (int t = 0; t < THREADNUM; t++) {
				for (int j = 0; j < NPs[REAL][1]; j++) {
					for (int k = 0; k < NPs[REAL][2]; k++) {
						phi[REALMODE_ARRAYINDEX(i, j, k)] += tmp_buffer1[REALMODE_ARRAYINDEX_MESH(t, i, j, k)];
					}
				}
			}
		}
		//#pragma omp barrier
		//#pragma omp for
		for (int i = 0; i < NPs[REAL][0]; i++) {
			for (int t = 0; t < THREADNUM; t++) {
				for (int d = 0; d < DIM; d++) {
					for (int j = 0; j < NPs[REAL][1]; j++) {
						for (int k = 0; k < NPs[REAL][2]; k++) {
							up[d][REALMODE_ARRAYINDEX(i, j, k)] += tmp_buffer_dim[d][REALMODE_ARRAYINDEX_MESH(t, i, j, k)];
						}
					}
				}
			}
		}
		//#endif
				// koba code //
		if (SW_UP) {
			double idmy_phi;

			//#ifdef _OPENMP
			//#pragma omp for private(im, idmy_phi)
			//#endif
			for (int i = 0; i < NPs[REAL][0]; i++) {
				for (int j = 0; j < NPs[REAL][1]; j++) {
					for (int k = 0; k < NPs[REAL][2]; k++) {
						im = REALMODE_ARRAYINDEX(i, j, k);
						if (phi[im] > 1.0) {
							idmy_phi = 1.0 / phi[im];
							up[0][im] *= idmy_phi;
							up[1][im] *= idmy_phi;
							up[2][im] *= idmy_phi;
							phi[im] = 1.0;
						}
					}
				}
			}
		}
	}
}

void Make_phi_particle(double *phi, Particle *p, const double radius) {
	const int SW_UP = 0;
	double **dmy_up = work_v3;
	int *nlattice;
	nlattice = Ns;
	Make_phi_u_primitive(phi, dmy_up, p, SW_UP, DX, NP_domain, Sekibun_cell, nlattice, radius);
}

void Make_phi_particle_sum(double *phi, double* phi_sum, Particle *p, const double radius) {
	int *nlattice;
	nlattice = Ns;
	Make_phi_particle_sum_primitive(phi, phi_sum, p, DX, NP_domain, Sekibun_cell, nlattice, radius);
}

void Make_phi_particle_sum_OBL(double *phi, double* phi_sum, Particle *p, const double radius) {
	int *nlattice;
	nlattice = Ns;
	Make_phi_particle_sum_primitive_OBL(phi, phi_sum, p, DX, NP_domain, Sekibun_cell, nlattice, radius);
}

void Make_u_particle_sum(double **up, double const* phi_sum, Particle *p, const double radius) {
	int *nlattice;
	nlattice = Ns;
	Make_u_particle_sum_primitive(up, phi_sum, p, DX, NP_domain, Sekibun_cell, nlattice, radius);
}

void Make_u_particle_sum_OBL(double **up, double const* phi_sum, Particle *p, const double radius) {
	int *nlattice;
	nlattice = Ns;
	Make_u_particle_sum_primitive_OBL(up, phi_sum, p, DX, NP_domain, Sekibun_cell, nlattice, radius);
}

void Make_phi_u_particle(double *phi, double **up, Particle *p) {
	const int SW_UP = 1;
	int *nlattice;
	nlattice = Ns;
	Make_phi_u_primitive(phi, up, p, SW_UP, DX, NP_domain, Sekibun_cell, nlattice);
}

void Make_phi_particle_OBL(double *phi, Particle *p, const double radius) {
	const int SW_UP = 0;
	double **dmy_up = work_v3;
	int *nlattice;
	nlattice = Ns;
	Make_phi_u_primitive_OBL(phi, dmy_up, p, SW_UP, DX, NP_domain, Sekibun_cell, nlattice, radius);
}
void Make_phi_u_particle_OBL(double *phi, double **up, Particle *p) {
	const int SW_UP = 1;
	int *nlattice;
	nlattice = Ns;
	Make_phi_u_primitive_OBL(phi, up, p, SW_UP, DX, NP_domain, Sekibun_cell, nlattice);
}

//not used
void Make_phi_u_advection(double *phi, double **up, Particle *p) {
	// map only V_p, excepting \Omega_p to the field up
	int *nlattice;
	if (SW_EQ == Shear_Navier_Stokes || SW_EQ == Shear_Navier_Stokes_Lees_Edwards) {
		nlattice = Ns_shear;
	} else {
		nlattice = Ns;
	}

	double xp[DIM], vp[DIM];
	int x_int[DIM];
	double residue[DIM];
	int sw_in_cell;
	int r_mesh[DIM];
	double r[DIM];
	double x[DIM];
	double dmy;
	double dmy_phi;
#pragma omp parallel for private(xp,vp,x_int,residue,sw_in_cell,r_mesh,r,x,dmy,dmy_phi)
	for (int n = 0; n < Particle_Number; n++) {
		for (int d = 0; d < DIM; d++) {
			xp[d] = p[n].x[d];
			vp[d] = p[n].v[d];
			{
				assert(p[n].x[d] >= 0);
				assert(p[n].x[d] < L[d]);
			}
		}

		sw_in_cell
			= Particle_cell(xp, DX, x_int, residue);// {1,0} が返ってくる
		sw_in_cell = 1;
		for (int mesh = 0; mesh < NP_domain; mesh++) {
			Relative_coord(Sekibun_cell[mesh], x_int, residue, sw_in_cell, nlattice, DX, r_mesh, r);
			for (int d = 0; d < DIM; d++) {
				x[d] = r_mesh[d] * DX;
			}
			dmy = Distance(x, xp);
			dmy_phi = Phi(dmy);
			int im = (r_mesh[0] * NY*NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
#pragma omp atomic
			phi[im] += dmy_phi;
#pragma omp atomic
			up[0][im] += (vp[0] * dmy_phi);
#pragma omp atomic
			up[1][im] += (vp[1] * dmy_phi);
#pragma omp atomic
			up[2][im] += (vp[2] * dmy_phi);
		}
	}
}

void Make_phi_rigid_mass(const double *phi_sum, Particle* p) {
	const double dx = DX;
	const double dx3 = DX3;
	const int np_domain = NP_domain;
	int const* const* sekibun_cell = Sekibun_cell;
	int const* nlattice = Ns;

#pragma omp parallel for 
	for (int rigidID = 0; rigidID < Rigid_Number; rigidID++) {
		double dmy, dmy_phi, dmy_mass;
		int x_int[DIM], r_mesh[DIM];
		double dmy_com[DIM], residue[DIM], r[DIM], x[DIM];
		dmy_mass = dmy_com[0] = dmy_com[1] = dmy_com[2] = 0.0;

		for (int n = Rigid_Particle_Cumul[rigidID]; n < Rigid_Particle_Cumul[rigidID + 1]; n++) {
			double xp[DIM];
			for (int d = 0; d < DIM; d++) xp[d] = p[n].x[d];

			int sw_in_cell = Particle_cell(xp, dx, x_int, residue);
			sw_in_cell = 1;
			for (int mesh = 0; mesh < np_domain; mesh++) {
				Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, nlattice, dx, r_mesh, r);
				for (int d = 0; d < DIM; d++) x[d] = r_mesh[d] * dx;
				int im = (r_mesh[0] * NY*NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];

				dmy = Distance(x, xp);
				dmy_phi = Phi(dmy) / MAX(phi_sum[im], 1.0);

				dmy_mass += dmy_phi;
				for (int d = 0; d < DIM; d++) {
					dmy_com[d] += dmy_phi*(xp[d] + r[d]);
				}
			}
		}


		Rigid_Masses[rigidID] = dmy_mass*dx3*RHO_particle[RigidID_Components[rigidID]];
		Rigid_IMasses[rigidID] = 1.0 / Rigid_Masses[rigidID];
		for (int d = 0; d < DIM; d++) xGs[rigidID][d] = dmy_com[d] / dmy_mass;
	}
}

void Make_phi_rigid_mass_OBL(const double *phi_sum, Particle* p) {
	const double dx = DX;
	const double dx3 = DX3;
	const int np_domain = NP_domain;
	int const* const* sekibun_cell = Sekibun_cell;
	int const* nlattice = Ns;

#pragma omp parallel for 
	for (int rigidID = 0; rigidID < Rigid_Number; rigidID++) {
		double dmy, dmy_phi, dmy_mass;
		int x_int[DIM], r_mesh[DIM];
		double dmy_com[DIM], residue[DIM], r[DIM], x[DIM];
		dmy_mass = dmy_com[0] = dmy_com[1] = dmy_com[2] = 0.0;

		for (int n = Rigid_Particle_Cumul[rigidID]; n < Rigid_Particle_Cumul[rigidID + 1]; n++) {
			double xp[DIM];
			for (int d = 0; d < DIM; d++) xp[d] = p[n].x[d];

			int sw_in_cell = Particle_cell(xp, dx, x_int, residue);
			sw_in_cell = 1;

			int sign, im;
			for (int mesh = 0; mesh < np_domain; mesh++) {
				sign = Relative_coord_check_stepover_Y(sekibun_cell[mesh], x_int, residue,
					sw_in_cell, nlattice, dx, r_mesh, r);

				for (int d = 0; d < DIM; d++) x[d] = r_mesh[d] * dx;

				im = (r_mesh[0] * NY*NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];

				dmy = Distance_OBL(x, xp);
				dmy_phi = Phi(dmy) / MAX(phi_sum[im], 1.0);

				dmy_mass += dmy_phi;
				for (int d = 0; d < DIM; d++) {
					dmy_com[d] += dmy_phi*(xp[d] + r[d]);
				}
			}
		}


		Rigid_Masses[rigidID] = dmy_mass*dx3*RHO_particle[RigidID_Components[rigidID]];
		Rigid_IMasses[rigidID] = 1.0 / Rigid_Masses[rigidID];
		for (int d = 0; d < DIM; d++) xGs[rigidID][d] = dmy_com[d] / dmy_mass;
	}
}

void Make_phi_rigid_inertia(const double *phi_sum, Particle* p) {
	const double dx = DX;
	const double dx3 = DX3;
	const int np_domain = NP_domain;
	int const* const* sekibun_cell = Sekibun_cell;
	int const* nlattice = Ns;

#pragma omp parallel for 
	for (int rigidID = 0; rigidID < Rigid_Number; rigidID++) {
		double dmy, dmy_phi;
		int x_int[DIM], r_mesh[DIM];
		double residue[DIM], r[DIM], x[DIM], dmy_inertia[DIM][DIM];
		double ri_x, ri_y, ri_z;

		dmy_inertia[0][0] = dmy_inertia[0][1] = dmy_inertia[0][2] = 0.0;
		dmy_inertia[1][0] = dmy_inertia[1][1] = dmy_inertia[1][2] = 0.0;
		dmy_inertia[2][0] = dmy_inertia[2][1] = dmy_inertia[2][2] = 0.0;

		for (int n = Rigid_Particle_Cumul[rigidID]; n < Rigid_Particle_Cumul[rigidID + 1]; n++) {
			double xp[DIM];
			for (int d = 0; d < DIM; d++) xp[d] = p[n].x[d];

			int sw_in_cell = Particle_cell(xp, dx, x_int, residue);
			sw_in_cell = 1;
			for (int mesh = 0; mesh < np_domain; mesh++) {
				Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, nlattice, dx, r_mesh, r);
				for (int d = 0; d < DIM; d++) x[d] = r_mesh[d] * dx;
				int im = (r_mesh[0] * NY*NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];

				ri_x = GRvecs[n][0] + r[0];
				ri_y = GRvecs[n][1] + r[1];
				ri_z = GRvecs[n][2] + r[2];

				dmy = Distance(x, xp);
				dmy_phi = Phi(dmy) / MAX(phi_sum[im], 1.0);
				dmy_inertia[0][0] += dmy_phi*(ri_y*ri_y + ri_z*ri_z);
				dmy_inertia[0][1] += dmy_phi*(-ri_x*ri_y);
				dmy_inertia[0][2] += dmy_phi*(-ri_x*ri_z);

				dmy_inertia[1][0] += dmy_phi*(-ri_y*ri_x);
				dmy_inertia[1][1] += dmy_phi*(ri_x*ri_x + ri_z*ri_z);
				dmy_inertia[1][2] += dmy_phi*(-ri_y*ri_z);

				dmy_inertia[2][0] += dmy_phi*(-ri_z*ri_x);
				dmy_inertia[2][1] += dmy_phi*(-ri_z*ri_y);
				dmy_inertia[2][2] += dmy_phi*(ri_x*ri_x + ri_y*ri_y);
			}
		}

		for (int d = 0; d < DIM; d++) {
			for (int e = 0; e < DIM; e++) {
				Rigid_Moments[rigidID][d][e] = dmy_inertia[d][e] * dx3*RHO_particle[RigidID_Components[rigidID]];
			}
		}
		Matrix_Inverse(Rigid_Moments[rigidID], Rigid_IMoments[rigidID], DIM);
		check_Inverse(Rigid_Moments[rigidID], Rigid_IMoments[rigidID], DIM);
	}
}

void Make_phi_rigid_inertia_OBL(const double *phi_sum, Particle* p) {
	const double dx = DX;
	const double dx3 = DX3;
	const int np_domain = NP_domain;
	int const* const* sekibun_cell = Sekibun_cell;
	int const* nlattice = Ns;

#pragma omp parallel for 
	for (int rigidID = 0; rigidID < Rigid_Number; rigidID++) {
		double dmy, dmy_phi;
		int x_int[DIM], r_mesh[DIM];
		double residue[DIM], r[DIM], x[DIM], dmy_inertia[DIM][DIM];
		double ri_x, ri_y, ri_z;

		dmy_inertia[0][0] = dmy_inertia[0][1] = dmy_inertia[0][2] = 0.0;
		dmy_inertia[1][0] = dmy_inertia[1][1] = dmy_inertia[1][2] = 0.0;
		dmy_inertia[2][0] = dmy_inertia[2][1] = dmy_inertia[2][2] = 0.0;

		for (int n = Rigid_Particle_Cumul[rigidID]; n < Rigid_Particle_Cumul[rigidID + 1]; n++) {
			double xp[DIM];
			for (int d = 0; d < DIM; d++) xp[d] = p[n].x[d];

			int sw_in_cell = Particle_cell(xp, dx, x_int, residue);
			sw_in_cell = 1;

			int im, sign;
			for (int mesh = 0; mesh < np_domain; mesh++) {
				sign = Relative_coord_check_stepover_Y(sekibun_cell[mesh], x_int, residue,
					sw_in_cell, nlattice, dx, r_mesh, r);

				for (int d = 0; d < DIM; d++) x[d] = r_mesh[d] * dx;

				im = (r_mesh[0] * NY*NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];

				ri_x = GRvecs[n][0] + r[0];
				ri_y = GRvecs[n][1] + r[1];
				ri_z = GRvecs[n][2] + r[2];

				dmy = Distance_OBL(x, xp);
				dmy_phi = Phi(dmy) / MAX(phi_sum[im], 1.0);

				dmy_inertia[0][0] += dmy_phi*(ri_y*ri_y + ri_z*ri_z);
				dmy_inertia[0][1] += dmy_phi*(-ri_x*ri_y);
				dmy_inertia[0][2] += dmy_phi*(-ri_x*ri_z);

				dmy_inertia[1][0] += dmy_phi*(-ri_y*ri_x);
				dmy_inertia[1][1] += dmy_phi*(ri_x*ri_x + ri_z*ri_z);
				dmy_inertia[1][2] += dmy_phi*(-ri_y*ri_z);

				dmy_inertia[2][0] += dmy_phi*(-ri_z*ri_x);
				dmy_inertia[2][1] += dmy_phi*(-ri_z*ri_y);
				dmy_inertia[2][2] += dmy_phi*(ri_x*ri_x + ri_y*ri_y);
			}
		}

		for (int d = 0; d < DIM; d++) {
			for (int e = 0; e < DIM; e++) {
				Rigid_Moments[rigidID][d][e] = dmy_inertia[d][e] * dx3*RHO_particle[RigidID_Components[rigidID]];
			}
		}

		Matrix_Inverse(Rigid_Moments[rigidID], Rigid_IMoments[rigidID], DIM);
		check_Inverse(Rigid_Moments[rigidID], Rigid_IMoments[rigidID], DIM);
	}
}

