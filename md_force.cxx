/*!
  \file md_force.cxx
  \brief Routines to compute MD forces on particles
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */

#include "md_force.h"


void Calc_f_Lennard_Jones_shear_cap_primitive_lnk(Particle *p, void(*distance0_func)(const double *x1, const double *x2, double &r12, double *x12), const double cap) {
	// Particle 変数の f に 
	// !! += 
	//で足す. f の初期値 が正しいと仮定している!!
	const double pair_cutoff = (!SW_PATCHY ? A_R_cutoff * LJ_dia : PATCHY_A_R_cutoff * SIGMA);
	double r_ij_vec[DIM] = { 0.0, 0.0, 0.0 };
	double r_ij = 0.0;
	double shear_stress[2] = { 0.0, 0.0 };
	double rigid_shear_stress[2] = { 0.0, 0.0 };

	// List Constructor
	int i, j;
	int *lscl;
	lscl = alloc_1d_int(Particle_Number);
	int lc[DIM];
	double lc_r[DIM];
	int mc[DIM];
	int lcyz, lcxyz;
	int ic[DIM], cn;
	lc[0] = int(NX / Cell_length);
	lc[1] = int(NY / Cell_length);
	lc[2] = int(NZ / Cell_length);
	for (int d = 0; d < DIM; d++) {
		lc_r[d] = Cell_length * DX;
	}
	lcyz = lc[1] * lc[2];
	lcxyz = lc[0] * lcyz;

	int *head;
	head = alloc_1d_int(lcxyz);

#ifdef _MPI
	double *dmy_mpi;
#endif
#ifdef _MPI
	Particle_Group_Communication(p, ONE_TO_ONE_LJ);
#endif

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(lcxyz, head)
#endif

	for (int d = 0; d < lcxyz; d++) head[d] = -1;

	for (int n = 0; n < Local_Particle_Number; n++) {
		Particle *p_n = &p[n];
		for (int d = 0; d < DIM; d++) {
			mc[d] = int((*p_n).x[d] / lc_r[d]);
		}
		cn = mc[0] * lcyz + mc[1] * lc[2] + mc[2];
		lscl[n] = head[cn];
		head[cn] = n;
	}

	for (ic[0] = 0; ic[0] < lc[0]; ic[0]++) {
		for (ic[1] = 0; ic[1] < lc[1]; ic[1]++) {
			for (ic[2] = 0; ic[2] < lc[2]; ic[2]++) {
				cn = ic[0] * lcyz + ic[1] * lc[2] + ic[2];
				// Scan the neighbor cells (including itself) of cell c
				int icl[DIM], icls[DIM];
				for (icl[0] = ic[0] - 1; icl[0] <= ic[0] + 1; icl[0]++) {
					for (icl[1] = ic[1] - 1; icl[1] <= ic[1] + 1; icl[1]++) {
						for (icl[2] = ic[2] - 1; icl[2] <= ic[2] + 1; icl[2]++) {
							// Consider periodic condition
							for (int d = 0; d < DIM; d++) {
								icls[d] = icl[d];
								if (icls[d] < 0) {
									icls[d] = lc[d] - 1;
								}
								if (icls[d] > (lc[d] - 1)) {
									icls[d] = 0;
								}
							}
							// Calculate the scalar cell index of the neighbor cell
							int cl = ((icls[0] + lc[0]) % lc[0])*lcyz + ((icls[1] + lc[1]) % lc[1])*lc[2] + ((icls[2] + lc[2]) % lc[2]);
							// Scan atom i in cell c
							i = head[cn];
							while (i != -1) {
								j = head[cl];
								double i_dir[DIM] = { 0.0, 0.0, 0.0 };
								int rigidID_i = -1;
								if (SW_PATCHY) rigid_body_rotation(i_dir, PATCHY_AXIS, p[i].q, BODY2SPACE);
								if (SW_PT == rigid) rigidID_i = Particle_RigidID[i];

								while (j != -1) {
									if (i > j && !rigid_chain(i, j) && !obstacle_chain(p[i].spec, p[j].spec)) {
										distance0_func(p[i].x, p[j].x, r_ij, r_ij_vec);

										if (r_ij < pair_cutoff) {
											double dmy_r, dmy_o;
											double j_dir[DIM] = { 0.0, 0.0, 0.0 };
											double n_ij_vec[DIM] = { 0.0, 0.0, 0.0 };
											dmy_r = dmy_o = 0.0;

											if (SW_PATCHY) {
												rigid_body_rotation(j_dir, PATCHY_AXIS, p[j].q, BODY2SPACE);
												n_ij_vec[0] = (j_dir[0] - i_dir[0]);
												n_ij_vec[1] = (j_dir[1] - i_dir[1]);
												n_ij_vec[2] = (j_dir[2] - i_dir[2]);

												patchy_janus_f(dmy_r, dmy_o, r_ij, n_ij_vec[0] * r_ij_vec[0] + n_ij_vec[1] * r_ij_vec[1] + n_ij_vec[2] * r_ij_vec[2], SIGMA);
												dmy_r = MIN(cap / r_ij, dmy_r);
											} else {
												dmy_r = MIN(cap / r_ij, Lennard_Jones_f(r_ij, LJ_dia));
											}

											//spherical particle forces
											double dmy_fi[DIM] = { 0.0, 0.0, 0.0 };
											for (int d = 0; d < DIM; d++) {
												dmy_fi[d] = (dmy_r)*(-r_ij_vec[d]) + (dmy_o)*(-n_ij_vec[d]);

												p[i].fr[d] += dmy_fi[d];
												p[j].fr[d] -= dmy_fi[d];
											}

											//spherical particle torques
											double dmy_ti[DIM] = { 0.0, 0.0, 0.0 };
											double dmy_tj[DIM] = { 0.0, 0.0, 0.0 };

											if (SW_PATCHY) {
												dmy_ti[0] = (dmy_o) * (r_ij_vec[1] * i_dir[2] - r_ij_vec[2] * i_dir[1]);
												dmy_ti[1] = (dmy_o) * (r_ij_vec[2] * i_dir[0] - r_ij_vec[0] * i_dir[2]);
												dmy_ti[2] = (dmy_o) * (r_ij_vec[0] * i_dir[1] - r_ij_vec[1] * i_dir[0]);

												dmy_tj[0] = (-dmy_o) * (r_ij_vec[1] * j_dir[2] - r_ij_vec[2] * j_dir[1]);
												dmy_tj[1] = (-dmy_o) * (r_ij_vec[2] * j_dir[0] - r_ij_vec[0] * j_dir[2]);
												dmy_tj[2] = (-dmy_o) * (r_ij_vec[0] * j_dir[1] - r_ij_vec[1] * j_dir[0]);
											}

											for (int d = 0; d < DIM; d++) {
												p[i].torque_r[d] += dmy_ti[d];
												p[j].torque_r[d] += dmy_tj[d];
											}

											//stress
											shear_stress[0] += (dmy_fi[0] * r_ij_vec[1]);
											shear_stress[1] += ((dmy_ti[2] + dmy_tj[2]) / 2.0);

											//rigid body forces & torques
											//particles treated as additive LJ centers: overlaps are not corrected
											if (SW_PT == rigid) {
												int rigidID_j = Particle_RigidID[j];

												for (int d = 0; d < DIM; d++) {
													forceGrs[rigidID_i][d] += dmy_fi[d];
													forceGrs[rigidID_j][d] -= dmy_fi[d];
												}

												torqueGrs[rigidID_i][0] += (dmy_ti[0] + (GRvecs[i][1] * dmy_fi[2] - GRvecs[i][2] * dmy_fi[1]));
												torqueGrs[rigidID_i][1] += (dmy_ti[1] + (GRvecs[i][2] * dmy_fi[0] - GRvecs[i][0] * dmy_fi[2]));
												torqueGrs[rigidID_i][2] += (dmy_ti[2] + (GRvecs[i][0] * dmy_fi[1] - GRvecs[i][1] * dmy_fi[0]));

												torqueGrs[rigidID_j][0] += (dmy_tj[0] - (GRvecs[j][1] * dmy_fi[2] - GRvecs[j][2] * dmy_fi[1]));
												torqueGrs[rigidID_j][1] += (dmy_tj[1] - (GRvecs[j][2] * dmy_fi[0] - GRvecs[j][0] * dmy_fi[2]));
												torqueGrs[rigidID_j][2] += (dmy_tj[2] - (GRvecs[j][0] * dmy_fi[1] - GRvecs[j][1] * dmy_fi[0]));

												double R_IJ_vec[DIM];
												double R_IJ;
												distance0_func(xGs[rigidID_i], xGs[rigidID_j], R_IJ, R_IJ_vec);
												rigid_shear_stress[0] += (dmy_fi[0] * R_IJ_vec[1]);
											}
										}
									}
									j = lscl[j];
								}
								i = lscl[i];
							}
						}
					}
				}
			}
		}
	}
	free_1d_int(lscl);
	free_1d_int(head);
#ifdef _MPI
	dmy_mpi = alloc_1d_double(procs);
	ierr = MPI_Allgather(&shear_stress[0], 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	shear_stress[0] = 0.0;
	for (int n = 0; n < procs; n++) {
		shear_stress[0] += dmy_mpi[n];
	}
	ierr = MPI_Allgather(&shear_stress[1], 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	shear_stress[1] = 0.0;
	for (int n = 0; n < procs; n++) {
		shear_stress[1] += dmy_mpi[n];
	}
	free_1d_double(dmy_mpi);
	Particle_Group_Communication(p, ONE_TO_ONE_SEKIBUN);
#endif

	dev_shear_stress_lj += shear_stress[0];
	dev_shear_stress_rot += shear_stress[1];

	if (SW_PT == rigid) {
		double dmy_shear = 0.0;
#pragma omp parallel for reduction(+:dmy_shear)
		for (int rigidID = 0; rigidID < Rigid_Number; rigidID++) {
			double IinvN[DIM] = { 0.0, 0.0, 0.0 };
			M_v_prod(IinvN, Rigid_IMoments[rigidID][0], torqueGrs[rigidID]);
			const double Jyy = (Rigid_Moments[rigidID][2][2] + Rigid_Moments[rigidID][0][0] - Rigid_Moments[rigidID][1][1]) / 2.0;
			const double Jzy = (-Rigid_Moments[rigidID][2][1]);
			dmy_shear += (Jyy*IinvN[2] - Jzy*IinvN[1]);
		}
		rigid_shear_stress[1] = dmy_shear;

		rigid_dev_shear_stress_lj += rigid_shear_stress[0];
		rigid_dev_shear_stress_rot += rigid_shear_stress[1];
	}
}

void Calc_f_Lennard_Jones_shear_cap_primitive(Particle *p, void(*distance0_func)(const double *x1, const double *x2, double &r12, double *x12), const double cap) {

#ifdef _MPI
	double *dmy_mpi;
#endif
#ifdef _MPI
	Particle_Group_Communication(p, ONE_TO_ONE_LJ);
#endif    

	const double pair_cutoff = (!SW_PATCHY ? A_R_cutoff * LJ_dia : PATCHY_A_R_cutoff * SIGMA);

	double shear_stress[2] = { 0.0, 0.0 };
	double rigid_shear_stress[2] = { 0.0, 0.0 };

	for (int n = 0; n < Update_Particle_Number; n++) {
		Particle *p_n = &p[n];
		double n_dir[DIM] = { 0.0, 0.0, 0.0 };
		int rigidID_n = -1;
		if (SW_PT == rigid) rigidID_n = Particle_RigidID[n];
		if (SW_PATCHY) rigid_body_rotation(n_dir, PATCHY_AXIS, p_n->q, BODY2SPACE);

		for (int m = n + 1; m < Update_Particle_Number; m++) {
			double m_dir[DIM] = { 0.0, 0.0, 0.0 };
			double n_ij_vec[DIM] = { 0.0, 0.0, 0.0 };
			double r_ij_vec[DIM] = { 0.0, 0.0, 0.0 };
			double r_ij = 0.0;

			distance0_func((*p_n).x, p[m].x, r_ij, r_ij_vec);

			if (r_ij < pair_cutoff && !rigid_chain(n, m) && !obstacle_chain(p[n].spec, p[m].spec)) {
				double dmy_r, dmy_o;
				dmy_r = dmy_o = 0.0;

				if (SW_PATCHY) {
					rigid_body_rotation(m_dir, PATCHY_AXIS, p[m].q, BODY2SPACE);
					n_ij_vec[0] = (m_dir[0] - n_dir[0]);
					n_ij_vec[1] = (m_dir[1] - n_dir[1]);
					n_ij_vec[2] = (m_dir[2] - n_dir[2]);

					patchy_janus_f(dmy_r, dmy_o, r_ij, n_ij_vec[0] * r_ij_vec[0] + n_ij_vec[1] * r_ij_vec[1] + n_ij_vec[2] * r_ij_vec[2], SIGMA);
				} else {
					dmy_r = MIN(cap / r_ij, Lennard_Jones_f(r_ij, LJ_dia));
					dmy_o = 0.0;
				}

				//forces
				double dmy_fn[DIM] = { 0.0, 0.0, 0.0 };
				for (int d = 0; d < DIM; d++) {
					dmy_fn[d] = (dmy_r) * (-r_ij_vec[d]) + (dmy_o) * (-n_ij_vec[d]);
					(*p_n).fr[d] += dmy_fn[d];
					p[m].fr[d] -= dmy_fn[d];
				}

				//torques	  
				double dmy_tn[DIM] = { 0.0, 0.0, 0.0 };
				double dmy_tm[DIM] = { 0.0, 0.0, 0.0 };
				if (SW_PATCHY) {
					dmy_tn[0] = (dmy_o)*(r_ij_vec[1] * n_dir[2] - r_ij_vec[2] * n_dir[1]);
					dmy_tn[1] = (dmy_o)*(r_ij_vec[2] * n_dir[0] - r_ij_vec[0] * n_dir[2]);
					dmy_tn[2] = (dmy_o)*(r_ij_vec[0] * n_dir[1] - r_ij_vec[1] * n_dir[0]);

					dmy_tm[0] = (-dmy_o)*(r_ij_vec[1] * m_dir[2] - r_ij_vec[2] * m_dir[1]);
					dmy_tm[1] = (-dmy_o)*(r_ij_vec[2] * m_dir[0] - r_ij_vec[0] * m_dir[2]);
					dmy_tm[2] = (-dmy_o)*(r_ij_vec[0] * m_dir[1] - r_ij_vec[1] * m_dir[0]);
				}
				for (int d = 0; d < DIM; d++) {
					(*p_n).torque_r[d] += dmy_tn[d];
					p[m].torque_r[d] += dmy_tm[d];
				}

				//stress
				shear_stress[0] += (dmy_fn[0] * r_ij_vec[1]);
				shear_stress[1] += ((dmy_tn[2] + dmy_tm[2]) / 2.0);

				// rigid body forces & torques
				if (SW_PT == rigid) {
					int rigidID_m = Particle_RigidID[m];

					for (int d = 0; d < DIM; d++) {
						forceGrs[rigidID_n][d] += dmy_fn[d];
						forceGrs[rigidID_m][d] -= dmy_fn[d];
					}

					torqueGrs[rigidID_n][0] += (dmy_tn[0] + (GRvecs[n][1] * dmy_fn[2] - GRvecs[n][2] * dmy_fn[1]));
					torqueGrs[rigidID_n][1] += (dmy_tn[1] + (GRvecs[n][2] * dmy_fn[0] - GRvecs[n][0] * dmy_fn[2]));
					torqueGrs[rigidID_n][2] += (dmy_tn[2] + (GRvecs[n][0] * dmy_fn[1] - GRvecs[n][1] * dmy_fn[0]));

					torqueGrs[rigidID_m][0] += (dmy_tm[0] - (GRvecs[m][1] * dmy_fn[2] - GRvecs[m][2] * dmy_fn[1]));
					torqueGrs[rigidID_m][1] += (dmy_tm[1] - (GRvecs[m][2] * dmy_fn[0] - GRvecs[m][0] * dmy_fn[2]));
					torqueGrs[rigidID_m][2] += (dmy_tm[2] - (GRvecs[m][0] * dmy_fn[1] - GRvecs[m][1] * dmy_fn[0]));

					double R_IJ_vec[DIM];
					double R_IJ;
					distance0_func(xGs[rigidID_n], xGs[rigidID_m], R_IJ, R_IJ_vec);
					rigid_shear_stress[0] += (dmy_fn[0] * R_IJ_vec[1]);
				}
			}
		}
	}

	for (int n = 0; n < Update_Particle_Number; n++) {
		Particle *p_n = &p[n];
		double n_dir[DIM] = { 0.0, 0.0, 0.0 };
		int rigidID_n = -1;
		if (SW_PT == rigid) rigidID_n = Particle_RigidID[n];
		if (SW_PATCHY) rigid_body_rotation(n_dir, PATCHY_AXIS, p_n->q, BODY2SPACE);

		for (int m = Update_Particle_Number; m < Particle_Number; m++) {
			double m_dir[DIM] = { 0.0, 0.0, 0.0 };
			double n_ij_vec[DIM] = { 0.0, 0.0, 0.0 };
			double r_ij_vec[DIM] = { 0.0, 0.0, 0.0 };
			double r_ij = 0.0;

			distance0_func((*p_n).x, p[m].x, r_ij, r_ij_vec);

			if (r_ij < pair_cutoff && !rigid_chain(n, m) && !obstacle_chain(p[n].spec, p[m].spec)) {
				double dmy_r, dmy_o;
				dmy_r = dmy_o = 0.0;

				if (SW_PATCHY) {
					rigid_body_rotation(m_dir, PATCHY_AXIS, p[m].q, BODY2SPACE);
					n_ij_vec[0] = (m_dir[0] - n_dir[0]);
					n_ij_vec[1] = (m_dir[1] - n_dir[1]);
					n_ij_vec[2] = (m_dir[2] - n_dir[2]);

					patchy_janus_f(dmy_r, dmy_o, r_ij, n_ij_vec[0] * r_ij_vec[0] + n_ij_vec[1] * r_ij_vec[1] + n_ij_vec[2] * r_ij_vec[2], SIGMA);
				} else {
					dmy_r = MIN(cap / r_ij, Lennard_Jones_f(r_ij, LJ_dia));
					dmy_o = 0.0;
				}

				//forces
				double dmy_fn[DIM] = { 0.0, 0.0, 0.0 };
				for (int d = 0; d < DIM; d++) {
					dmy_fn[d] = (dmy_r) * (-r_ij_vec[d]) + (dmy_o) * (-n_ij_vec[d]);
					(*p_n).fr[d] += dmy_fn[d];
					//p[m].fr[d]   -= dmy_fn[d];
				}

				//torques	  
				double dmy_tn[DIM] = { 0.0, 0.0, 0.0 };
				double dmy_tm[DIM] = { 0.0, 0.0, 0.0 };
				if (SW_PATCHY) {
					dmy_tn[0] = (dmy_o)*(r_ij_vec[1] * n_dir[2] - r_ij_vec[2] * n_dir[1]);
					dmy_tn[1] = (dmy_o)*(r_ij_vec[2] * n_dir[0] - r_ij_vec[0] * n_dir[2]);
					dmy_tn[2] = (dmy_o)*(r_ij_vec[0] * n_dir[1] - r_ij_vec[1] * n_dir[0]);

					dmy_tm[0] = (-dmy_o)*(r_ij_vec[1] * m_dir[2] - r_ij_vec[2] * m_dir[1]);
					dmy_tm[1] = (-dmy_o)*(r_ij_vec[2] * m_dir[0] - r_ij_vec[0] * m_dir[2]);
					dmy_tm[2] = (-dmy_o)*(r_ij_vec[0] * m_dir[1] - r_ij_vec[1] * m_dir[0]);
				}
				for (int d = 0; d < DIM; d++) {
					(*p_n).torque_r[d] += dmy_tn[d];
					//p[m].torque_r[d]   += dmy_tm[d];
				}

				//stress
				//avoid double count
				shear_stress[0] += (dmy_fn[0] * r_ij_vec[1]) * 0.5;
				shear_stress[1] += ((dmy_tn[2] + dmy_tm[2]) / 2.0) * 0.5;

				// rigid body forces & torques
				if (SW_PT == rigid) {
					int rigidID_m = Particle_RigidID[m];

					for (int d = 0; d < DIM; d++) {
						forceGrs[rigidID_n][d] += dmy_fn[d];
						forceGrs[rigidID_m][d] -= dmy_fn[d];
					}

					torqueGrs[rigidID_n][0] += (dmy_tn[0] + (GRvecs[n][1] * dmy_fn[2] - GRvecs[n][2] * dmy_fn[1]));
					torqueGrs[rigidID_n][1] += (dmy_tn[1] + (GRvecs[n][2] * dmy_fn[0] - GRvecs[n][0] * dmy_fn[2]));
					torqueGrs[rigidID_n][2] += (dmy_tn[2] + (GRvecs[n][0] * dmy_fn[1] - GRvecs[n][1] * dmy_fn[0]));

					torqueGrs[rigidID_m][0] += (dmy_tm[0] - (GRvecs[m][1] * dmy_fn[2] - GRvecs[m][2] * dmy_fn[1]));
					torqueGrs[rigidID_m][1] += (dmy_tm[1] - (GRvecs[m][2] * dmy_fn[0] - GRvecs[m][0] * dmy_fn[2]));
					torqueGrs[rigidID_m][2] += (dmy_tm[2] - (GRvecs[m][0] * dmy_fn[1] - GRvecs[m][1] * dmy_fn[0]));

					double R_IJ_vec[DIM];
					double R_IJ;
					distance0_func(xGs[rigidID_n], xGs[rigidID_m], R_IJ, R_IJ_vec);
					rigid_shear_stress[0] += (dmy_fn[0] * R_IJ_vec[1]);
				}
			}
		}
	}

#ifdef _MPI
	dmy_mpi = alloc_1d_double(procs);
	ierr = MPI_Allgather(&shear_stress[0], 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	shear_stress[0] = 0.0;
	for (int n = 0; n < procs; n++) {
		shear_stress[0] += dmy_mpi[n];
	}
	ierr = MPI_Allgather(&shear_stress[1], 1, MPI_DOUBLE, dmy_mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	shear_stress[1] = 0.0;
	for (int n = 0; n < procs; n++) {
		shear_stress[1] += dmy_mpi[n];
	}
	free_1d_double(dmy_mpi);
	//Particle_Group_Communication (p, ONE_TO_ONE_SEKIBUN);
#endif

	dev_shear_stress_lj += shear_stress[0];
	dev_shear_stress_rot += shear_stress[1];

	if (SW_PT == rigid) {
		double dmy_shear = 0.0;
#pragma omp parallel for reduction(+:dmy_shear)
		for (int rigidID = 0; rigidID < Rigid_Number; rigidID++) {
			double IinvN[DIM] = { 0.0, 0.0, 0.0 };
			M_v_prod(IinvN, Rigid_IMoments[rigidID][0], torqueGrs[rigidID]);
			const double Jyy = (Rigid_Moments[rigidID][2][2] + Rigid_Moments[rigidID][0][0] - Rigid_Moments[rigidID][1][1]) / 2.0;
			const double Jzy = (-Rigid_Moments[rigidID][2][1]);
			dmy_shear += (Jyy*IinvN[2] - Jzy*IinvN[1]);
		}
		rigid_shear_stress[1] = dmy_shear;

		rigid_dev_shear_stress_lj += rigid_shear_stress[0];
		rigid_dev_shear_stress_rot += rigid_shear_stress[1];
	}
}

void Add_f_gravity(Particle *p) {
	static const double Gravity_on_fluid = G*RHO * 4. / 3.*M_PI * SQ(RADIUS)* RADIUS;

	if (SW_PT != rigid) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(Update_Particle_Number, p, G_direction, Gravity_on_fluid, MASS_RATIOS)
#endif
		for (int n = 0; n < Update_Particle_Number; n++) {
			p[n].fr[G_direction] -= Gravity_on_fluid * (MASS_RATIOS[p[n].spec] - 1.0);
		}
	} else {
#pragma omp parallel for
		for (int rigidID = 0; rigidID < Rigid_Number; rigidID++) {
			const int rigid_spec = RigidID_Components[rigidID];
			const double rigid_volume = Rigid_Masses[rigidID] / RHO_particle[rigid_spec];
			forceGrs[rigidID][G_direction] -= G*RHO*rigid_volume*(MASS_RATIOS[rigid_spec] - 1.0);
		}
	}
}

void Calc_f_slip_correct_precision(Particle *p, double const* const* u, const CTime &jikan) {
	static const double dmy0 = DX3*RHO;
	double dmy = dmy0 / jikan.dt_fluid;
	int *nlattice;
	nlattice = Ns;
	double *buffer;
	double xp[DIM];
	int x_int[DIM];
	double residue[DIM];
	int sw_in_cell;
	double *force;
	double *torque;
	int r_mesh[DIM];
	int dmy_mesh[DIM];
	int position;
	int size;
	double r[DIM];
	double dmy_fp[DIM];
	double x[DIM];
	double dmyR;
	double dmy_phi;
	int pspec;
#if defined _MPI
	double **recvdata;
	double **senddata;
	int **uplist;
	int *idchk;
	int *recvbuf;
	int *recvranks;
	int *sendbuf;
	int *sendcounts;
	int *listcounts;
	int *sendranks;
	int x_id;
	int y_id;
	int id;
#endif

#if defined _MPI
	if (procs > 1) {
		size = DIM + 1;
		buffer = calloc_1d_double(Local_Particle_Number * size);
	}
#endif
	force = alloc_1d_double(Local_Particle_Number);
	if (ROTATION) {
		torque = alloc_1d_double(Local_Particle_Number);
	}
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (Local_Particle_Number, Update_Particle_Number, p, DX, NP_domain, Sekibun_cell, Ns, RADIUS, NPs, NZ_, Angular2v, ROTATION, u, procs, buffer, size, nlattice, dmy_mesh, position) private (xp,x_int,residue,sw_in_cell,force,torque,r_mesh,r,dmy_fp,x,dmyR,dmy_phi,pspec) 
#endif
	for (int n = 0; n < Local_Particle_Number; n++) {
		for (int d = 0; d < DIM; d++) {
			xp[d] = p[n].x[d];
			force[d] = torque[d] = 0.0;
		}

		sw_in_cell = Particle_cell(xp, DX, x_int, residue);
		sw_in_cell = 1;

		for (int mesh = 0; mesh < NP_domain; mesh++) {
			Relative_coord(Sekibun_cell[mesh], x_int, residue, sw_in_cell, nlattice, DX, r_mesh, r);
			if (Range_coord(r_mesh, dmy_mesh)) {
				double x[DIM];
				for (int d = 0; d < DIM; d++) {
					x[d] = r_mesh[d] * DX;
				}
				dmyR = Distance(x, xp);
				dmy_phi = Phi(dmyR, RADIUS);

				int im = REALMODE_ARRAYINDEX(dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
				for (int d = 0; d < DIM; d++) {
					dmy_fp[d] = u[d][im] * dmy_phi;
					force[d] += dmy_fp[d];
				}
				// torque
				torque[0] += (r[1] * dmy_fp[2] - r[2] * dmy_fp[1]);
				torque[1] += (r[2] * dmy_fp[0] - r[0] * dmy_fp[2]);
				torque[2] += (r[0] * dmy_fp[1] - r[1] * dmy_fp[0]);
			}
		}// mesh
		pspec = p[n].spec;

#ifdef _MPI
		if (procs > 1 && Update_Particle_Number <= n) {
			position = n * size;
			// バッファに粒子データを詰め込む
			buffer[position + 0] = force[0];
			buffer[position + 1] = force[1];
			buffer[position + 2] = force[2];
			if (ROTATION) {
				buffer[position + 3] = torque[0];
				buffer[position + 4] = torque[1];
				buffer[position + 5] = torque[2];
			}
		}
#endif
	}
#ifdef _MPI
	if (procs > 1) {
		idchk = alloc_1d_int(procs);
		recvbuf = alloc_1d_int(3 * procs);
		recvcounts = alloc_1d_int(procs);
		recvdata = alloc_2d_double(procs, Update_Particle_Number * size);
		recvranks = alloc_1d_int(procs);
		sendbuf = alloc_1d_int(3 * procs);
		sendcounts = alloc_1d_int(procs);
		senddata = alloc_2d_double(procs, Local_Particle_Number * size);
		uplist = alloc_2d_int(procs, Update_Particle_Number);
		listcounts = alloc_1d_int(procs);
		sendranks = alloc_1d_int(procs);
		// 初期化
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(Local_Particle_Number, Update_Particle_Number, procs, recvbuf, recvcounts, recvdata, recvranks, sendbuf, sendcounts, senddata, sendranks, size, idchk, uplist, listcounts)
#endif
		for (int prc = 0; prc < procs; prc++) {
			idchk[prc] = 0;
			recvbuf[2 * prc] = 0;
			recvbuf[2 * prc + 1] = 0;
			recvcounts[prc] = 0;
			recvranks[prc] = MPI_PROC_NULL;
			for (int n = 0; n < (Update_Particle_Number * size); n++) {
				recvdata[prc][n] = 0.0;
			}
			sendbuf[2 * prc] = 0;
			sendbuf[2 * prc + 1] = 0;
			sendcounts[prc] = 0;
			sendranks[prc] = MPI_PROC_NULL;
			for (int n = 0; n < (Local_Particle_Number * size); n++) {
				senddata[prc][n] = 0.0;
			}
			for (int n = 0; n < Update_Particle_Number; n++) {
				uplist[prc][n] = 0;
			}
			listcounts[prc] = 0;
		}
		// 足しこみ先リスト作成
		// プロセス担当粒子の粒子座標から、プロセス担当粒子あるいは参照用粒子と
		// して粒子が格納されているランク(送信元プロセスID)を求める
		// ※粒子格納順序が保障される(プロセス担当粒子、参照用粒子共に送信元
		//   プロセスの昇順に格納する。また各々の送信元プロセスが送信した粒子順序で格納)ため、
		for (int n = 0; n < Update_Particle_Number; n++) {
			for (int i = 0; i < sekibun_size; i++) {
				id = SEKIBUN_REFID(procid, i);
				uplist[id][listcounts[id]] = n;
				listcounts[id]++;
			}
			for (int i = 0; i < procs; i++) {
				idchk[i] = 0;
			}
		}
		// 粒子データ送信先の振り分け
		for (int n = Update_Particle_Number; n < Local_Particle_Number; n++) {
			position = n * size;
			x_id = (int)(p[n].x[0] * IDX);
			y_id = (int)(p[n].x[1] * IDX);
			// 粒子座標から、その粒子がプロセス担当粒子として存在するランクに送信する
			id = ID(xmesh[x_id], ymesh[y_id]);
			sendranks[id] = id;
			recvranks[id] = procid;
			for (int j = 0; j < size; j++) {
				senddata[id][sendcounts[id] + j] = buffer[position + j];
			}
			sendcounts[id] += size;
			recvcounts[id] += size;
		}
		// 通信先ランクおよび通信データ個数の通知
		for (int prc = 0; prc < procs; prc++) {
			sendbuf[3 * prc] = recvranks[prc];
			sendbuf[3 * prc + 1] = recvcounts[prc];
		}

		ierr = MPI_Alltoall(sendbuf, 3, MPI_INT, recvbuf, 3, MPI_INT,
			MPI_COMM_WORLD);

		MPI_ERRORCHECK(ierr);
		for (int prc = 0; prc < procs; prc++) {
			recvranks[prc] = recvbuf[3 * prc];
			recvcounts[prc] = recvbuf[3 * prc + 1];
		}
		// 粒子通信
		for (int prc = 0; prc < procs; prc++) {

			ierr = MPI_Irecv(recvdata[prc], recvcounts[prc], MPI_DOUBLE,
				recvranks[prc], TAG(procid, recvranks[prc]),
				MPI_COMM_WORLD, (ireq + prc));

			MPI_ERRORCHECK(ierr);
		}
		for (int prc = 0; prc < procs; prc++) {

			ierr = MPI_Isend(senddata[prc], sendcounts[prc], MPI_DOUBLE,
				sendranks[prc], TAG(sendranks[prc], procid),
				MPI_COMM_WORLD, (ireq + procs + prc));

			MPI_ERRORCHECK(ierr);
		}

		ierr = MPI_Waitall((2 * procs), ireq, ista);

		MPI_ERRORCHECK(ierr);
		for (int prc = 0; prc < procs; prc++) {
			recvcounts[prc] /= size;
		}
		assert(recvranks[procid] == MPI_PROC_NULL);
		// リストを元に足し合わせ
		for (int prc = 0; prc < procs; prc++) {
			if (recvranks[prc] == MPI_PROC_NULL) continue;
			assert(recvcounts[prc] == listcounts[prc]);
			for (int n = 0; n < recvcounts[prc]; n++) {
				position = n * size;
				id = uplist[prc][n];
				force[0] += recvdata[prc][position + 0];
				force[1] += recvdata[prc][position + 1];
				force[2] += recvdata[prc][position + 2];
				if (ROTATION) {
					torque[0] += recvdata[prc][position + 3];
					torque[1] += recvdata[prc][position + 4];
					torque[2] += recvdata[prc][position + 5];
				}
			}
		}
		free_1d_double(buffer);
		free_1d_int(idchk);
		free_1d_int(recvbuf);
		free_1d_int(recvcounts);
		free_2d_double(recvdata);
		free_1d_int(recvranks);
		free_1d_int(sendbuf);
		free_1d_int(sendcounts);
		free_2d_double(senddata);
		free_1d_int(sendranks);
		free_2d_int(uplist);
		free_1d_int(listcounts);
	}
#endif
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(Update_Particle_Number, force, torque, ROTATION, dmy, p)
#endif
	for (int n = 0; n < Update_Particle_Number; n++) {
		for (int d = 0; d < DIM; d++) {
			p[n].f_slip[d] = (dmy * force[d]);
			if (ROTATION) {
				p[n].torque_slip[d] = (dmy * torque[d]);
			}
		}
	}
	free_1d_double(force);
	free_1d_double(torque);
#ifdef _MPI
	Particle_Group_Communication(p, ONE_TO_ONE_SEKIBUN);
#endif
}

void Calc_f_hydro_correct_precision(Particle *p, double const *phi_sum, double const *const *u, const CTime &jikan) {
	static const double dmy0 = -DX3 * RHO;
	double dmy = dmy0 / jikan.dt_fluid;
	int *nlattice;
	nlattice = Ns;
#if defined _MPI
	double **force;
	double **torque;
	double *buffer;

	force = alloc_2d_double(Local_Particle_Number, DIM);
	if (ROTATION) {
		torque = alloc_2d_double(Local_Particle_Number, DIM);
	}
	int position;
	int size;
#else
	double force[DIM];
	double torque[DIM];
#endif
	double xp[DIM];
	double vp[DIM];
	double omega_p[DIM];
	double x[DIM];
	double r[DIM];
	double residue[DIM];
	double v_rot[DIM];
	double dmy_fp[DIM];
	double dmy_phi;
	double dmyR;
	int dmy_mesh[DIM];
	int r_mesh[DIM];
	int x_int[DIM];
	int sw_in_cell;

	int im;
	int pspec;
	double forceg[DIM];
	double torqueg[DIM];
#if defined _MPI
	double **recvdata;
	double **senddata;
	int **uplist;
	int *idchk;
	int *recvbuf;
	int *recvranks;
	int *sendbuf;
	int *sendcounts;
	int *listcounts;
	int *sendranks;
	int x_id;
	int y_id;
	int id;
#endif
	int rigidID;
	// initialize forceGs and torqueGs
	for (int rigidID = 0; rigidID < Rigid_Number; rigidID++) {
		for (int d = 0; d < DIM; d++) {
			forceGs[rigidID][d] = 0.0;
			torqueGs[rigidID][d] = 0.0;
		}
	}

#if defined _MPI
	if (procs > 1) {
		size = (ROTATION) ? (2 * DIM) + 1 : DIM + 1;
		buffer = calloc_1d_double(Local_Particle_Number * size);
	}
#endif


#ifndef _MPI
#ifdef _OPENMP
#pragma omp parallel for private(                                              \
    xp, vp, omega_p, x_int, residue, sw_in_cell, force, torque, r_mesh, r,     \
    dmy_fp, x, dmyR, dmy_phi, v_rot, pspec, rigidID, forceg, torqueg)
#endif
#endif
	for (int n = 0; n < Local_Particle_Number; n++) {
		if (SW_PT == rigid)
			rigidID = Particle_RigidID[n];
		for (int d = 0; d < DIM; d++) {
#ifdef _MPI
			force[n][d] = 0.0;
			if (ROTATION) {
				torque[n][d] = 0.0;
			}
#else
			force[d] = torque[d] = 0.0;
			forceg[d] = torqueg[d] = 0.0;
#endif
			xp[d] = p[n].x[d];
			vp[d] = p[n].v[d];
			omega_p[d] = p[n].omega[d];
		}
		sw_in_cell = Particle_cell(xp, DX, x_int, residue);
		sw_in_cell = 1;
		for (int mesh = 0; mesh < NP_domain; mesh++) {
			Relative_coord(Sekibun_cell[mesh], x_int, residue, sw_in_cell, Ns, DX,
				r_mesh, r);
#ifdef _MPI
			if (Range_coord(r_mesh, dmy_mesh)) {
#endif
				for (int d = 0; d < DIM; d++) {
					x[d] = r_mesh[d] * DX;
				}
				dmyR = Distance(x, xp); // vesion2.00 needs this value
				Angular2v(omega_p, r, v_rot);
#ifdef _MPI
				im = REALMODE_ARRAYINDEX(dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
#else
				im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
#endif
				//3.32 or 3.40?
				//dmy_phi = Phi(Distance(x, xp), RADIUS);
				dmy_phi = Phi(dmyR, RADIUS) / MAX(phi_sum[im], 1.0);
				for (int d = 0; d < DIM; d++) {
					dmy_fp[d] = ((vp[d] + v_rot[d]) - u[d][im]) * dmy_phi;
#ifdef _MPI
					force[n][d] += dmy_fp[d];
#else
					force[d] += dmy_fp[d];
#endif
				}
				//if (ROTATION) {
#ifdef _MPI
					torque[n][0] += (r[1] * dmy_fp[2] - r[2] * dmy_fp[1]);
					torque[n][1] += (r[2] * dmy_fp[0] - r[0] * dmy_fp[2]);
					torque[n][2] += (r[0] * dmy_fp[1] - r[1] * dmy_fp[0]);
#else
					torque[0] += (r[1] * dmy_fp[2] - r[2] * dmy_fp[1]);
					torque[1] += (r[2] * dmy_fp[0] - r[0] * dmy_fp[2]);
					torque[2] += (r[0] * dmy_fp[1] - r[1] * dmy_fp[0]);
#endif
				//}
				if (SW_PT == rigid) {
					for (int d = 0; d < DIM; d++)
						forceg[d] += dmy_fp[d];
					torqueg[0] += ((GRvecs[n][1] + r[1]) * dmy_fp[2] -
						(GRvecs[n][2] + r[2]) * dmy_fp[1]);
					torqueg[1] += ((GRvecs[n][2] + r[2]) * dmy_fp[0] -
						(GRvecs[n][0] + r[0]) * dmy_fp[2]);
					torqueg[2] += ((GRvecs[n][0] + r[0]) * dmy_fp[1] -
						(GRvecs[n][1] + r[1]) * dmy_fp[0]);
				}
#ifdef _MPI
			}
#endif
		} // mesh
		pspec = p[n].spec;


#ifndef _MPI
		for (int d = 0; d < DIM; d++) {
			p[n].f_hydro[d] = (dmy * force[d]);
			p[n].f_slip[d] = 0.0;
		}
		if (ROTATION) {
			for (int d = 0; d < DIM; d++) {
				p[n].torque_hydro[d] = (dmy * torque[d]);
				p[n].torque_slip[d] = 0.0;
			}
		}
#endif


		if (SW_PT == rigid) {
			for (int d = 0; d < DIM; d++) {
#ifdef _OPENMP
#pragma omp atomic
#endif
				forceGs[rigidID][d] += dmy * forceg[d];
#ifdef _OPENMP
#pragma omp atomic
#endif
				torqueGs[rigidID][d] += dmy * torqueg[d];
			}
		}

#ifdef _MPI
		if (procs > 1 && Update_Particle_Number <= n) {
			position = n * size;
			// バッファに粒子データを詰め込む
			buffer[position + 0] = force[n][0];
			buffer[position + 1] = force[n][1];
			buffer[position + 2] = force[n][2];
			if (ROTATION) {
				buffer[position + 3] = torque[n][0];
				buffer[position + 4] = torque[n][1];
				buffer[position + 5] = torque[n][2];
			}
		}
#endif
	}

#ifdef _MPI
	if (procs > 1) {
		idchk = alloc_1d_int(procs);
		recvbuf = alloc_1d_int(3 * procs);
		recvcounts = alloc_1d_int(procs);
		recvdata = alloc_2d_double(procs, Local_Particle_Number * size);
		recvranks = alloc_1d_int(procs);
		sendbuf = alloc_1d_int(3 * procs);
		sendcounts = alloc_1d_int(procs);
		senddata = alloc_2d_double(procs, Local_Particle_Number * size);
		uplist = alloc_2d_int(procs, Local_Particle_Number);
		listcounts = alloc_1d_int(procs);
		sendranks = alloc_1d_int(procs);
		// 初期化
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(Local_Particle_Number, Update_Particle_Number, procs, recvbuf,      \
           recvcounts, recvdata, recvranks, sendbuf, sendcounts, senddata,     \
           sendranks, size, idchk, uplist, listcounts)
#endif
		for (int prc = 0; prc < procs; prc++) {
			idchk[prc] = 0;
			recvbuf[2 * prc] = 0;
			recvbuf[2 * prc + 1] = 0;
			recvcounts[prc] = 0;
			recvranks[prc] = MPI_PROC_NULL;
			for (int n = 0; n < (Local_Particle_Number * size); n++) {
				recvdata[prc][n] = 0.0;
			}
			sendbuf[2 * prc] = 0;
			sendbuf[2 * prc + 1] = 0;
			sendcounts[prc] = 0;
			sendranks[prc] = MPI_PROC_NULL;
			for (int n = 0; n < (Local_Particle_Number * size); n++) {
				senddata[prc][n] = 0.0;
			}
			for (int n = 0; n < Local_Particle_Number; n++) {
				uplist[prc][n] = 0;
			}
			listcounts[prc] = 0;
		}
		// 足しこみ先リスト作成
		// プロセス担当粒子の粒子座標から、プロセス担当粒子あるいは参照用粒子と
		// して粒子が格納されているランク(送信元プロセスID)を求める
		// ※粒子格納順序が保障される(プロセス担当粒子、参照用粒子共に送信元
		//   プロセスの昇順に格納する。また各々の送信元プロセスが送信した粒子順序で格納)ため、
		for (int n = 0; n < Update_Particle_Number; n++) {
			for (int i = 0; i < sekibun_size; i++) {
				id = SEKIBUN_REFID(procid, i);
				uplist[id][listcounts[id]] = n;
				listcounts[id]++;
			}
			for (int i = 0; i < procs; i++) {
				idchk[i] = 0;
			}
		}

		// 粒子データ送信先の振り分け
		for (int n = Update_Particle_Number; n < Local_Particle_Number; n++) {
			position = n * size;
			x_id = (int)(p[n].x[0] * IDX);
			y_id = (int)(p[n].x[1] * IDX);
			// 粒子座標から、その粒子がプロセス担当粒子として存在するランクに送信する
			id = ID(xmesh[x_id], ymesh[y_id]);
			sendranks[id] = id;
			recvranks[id] = procid;
			for (int j = 0; j < size; j++) {
				senddata[id][sendcounts[id] + j] = buffer[position + j];
			}
			sendcounts[id] += size;
			recvcounts[id] += size;
		}
		// 通信先ランクおよび通信データ個数の通知
		for (int prc = 0; prc < procs; prc++) {
			sendbuf[3 * prc] = recvranks[prc];
			sendbuf[3 * prc + 1] = recvcounts[prc];
		}

		ierr =
			MPI_Alltoall(sendbuf, 3, MPI_INT, recvbuf, 3, MPI_INT, MPI_COMM_WORLD);

		MPI_ERRORCHECK(ierr);
		for (int prc = 0; prc < procs; prc++) {
			recvranks[prc] = recvbuf[3 * prc];
			recvcounts[prc] = recvbuf[3 * prc + 1];
		}
		// 粒子通信
		for (int prc = 0; prc < procs; prc++) {

			ierr =
				MPI_Irecv(recvdata[prc], recvcounts[prc], MPI_DOUBLE, recvranks[prc],
					TAG(procid, recvranks[prc]), MPI_COMM_WORLD, (ireq + prc));

			MPI_ERRORCHECK(ierr);
		}
		for (int prc = 0; prc < procs; prc++) {

			ierr = MPI_Isend(senddata[prc], sendcounts[prc], MPI_DOUBLE,
				sendranks[prc], TAG(sendranks[prc], procid),
				MPI_COMM_WORLD, (ireq + procs + prc));

			MPI_ERRORCHECK(ierr);
		}

		ierr = MPI_Waitall((2 * procs), ireq, ista);

		MPI_ERRORCHECK(ierr);
		for (int prc = 0; prc < procs; prc++) {
			recvcounts[prc] /= size;
		}
		assert(recvranks[procid] == MPI_PROC_NULL);
		// リストを元に足し合わせ
		for (int prc = 0; prc < procs; prc++) {
			if (recvranks[prc] == MPI_PROC_NULL)
				continue;
			assert(recvcounts[prc] == listcounts[prc]);
			for (int n = 0; n < recvcounts[prc]; n++) {
				position = n * size;
				id = uplist[prc][n];
				force[id][0] += recvdata[prc][position + 0];
				force[id][1] += recvdata[prc][position + 1];
				force[id][2] += recvdata[prc][position + 2];
				if (ROTATION) {
					torque[id][0] += recvdata[prc][position + 3];
					torque[id][1] += recvdata[prc][position + 4];
					torque[id][2] += recvdata[prc][position + 5];
				}
			}
		}
		free_1d_double(buffer);
		free_1d_int(idchk);
		free_1d_int(recvbuf);
		free_1d_int(recvcounts);
		free_2d_double(recvdata);
		free_1d_int(recvranks);
		free_1d_int(sendbuf);
		free_1d_int(sendcounts);
		free_2d_double(senddata);
		free_1d_int(sendranks);
		free_2d_int(uplist);
		free_1d_int(listcounts);
	}
#endif
#ifdef _MPI
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(Update_Particle_Number, force, torque, ROTATION, dmy, p)
#endif
	for (int n = 0; n < Update_Particle_Number; n++) {
		for (int d = 0; d < DIM; d++) {
			p[n].f_hydro[d] = (dmy * force[n][d]);
			p[n].f_slip[d] = 0.0;
			if (ROTATION) {
				p[n].torque_hydro[d] = (dmy * torque[n][d]);
				p[n].torque_slip[d] = 0.0;
			}
		}
	}
#endif
#ifdef _MPI
	free_2d_double(force);
	if (ROTATION) {
		free_2d_double(torque);
	}
#endif
#ifdef _MPI
	Particle_Group_Communication(p, ONE_TO_ONE_SEKIBUN);
#endif
}
void Calc_f_hydro_correct_precision_OBL(Particle *p, double const* phi_sum, double const* const* u, const CTime &jikan) {
#ifdef _MPI
	static const double dmy0 = -DX3*RHO;
	double dmy = dmy0 / jikan.dt_fluid;
	int *nlattice;
	nlattice = Ns;
	double xp[DIM], vp[DIM], omega_p[DIM];
	int x_int[DIM];
	double residue[DIM];
	int sw_in_cell;

	double r[DIM];
	double dmy_fp[DIM];
	double x[DIM];
	double dmyR;
	double dmy_phi;
	double dmy_rhop;
	double dmy_ry;
	double v_rot[DIM];
	int sign;
	int size;
	int im;
	int dmy_mesh[DIM];
	int r_mesh[DIM];
	int position;
	double sum_force = 0.0;
	double sum_volume = 0.0;

	double forceg[DIM];
	double torqueg[DIM];
	double dVg[Rigid_Number][DIM];
	double dWg[Rigid_Number][DIM];

	int rigidID;

	int skip_flag[Particle_Number];

	// initialize forceGs and torqueGs
	for (int rigidID = 0; rigidID < Rigid_Number; rigidID++) {
		for (int d = 0; d < DIM; d++) {
			forceGs[rigidID][d] = 0.0;
			torqueGs[rigidID][d] = 0.0;
		}
	}

	Reset_phi(Hydro_force);
	Reset_phi(Hydro_force_new);

	for (int n = 0; n < Particle_Number; n++) {
		skip_flag[n] = 1;
		for (int d = 0; d < DIM; d++) {
			xp[d] = p[n].x[d];
		}
		sw_in_cell = Particle_cell(xp, DX, x_int, residue);// {1,0}
		sw_in_cell = 1;
		for (int mesh = 0; mesh < NP_domain; mesh++) {

			// pre judge.
			r_mesh[1] = (x_int[1] + Sekibun_cell[mesh][1] + nlattice[1]) % nlattice[1];
			if (!(r_mesh[1] >= PREV_NPs[REAL][1] && NEXT_NPs[REAL][1] > r_mesh[1])) {
				p_info[n][mesh].im = -1;
				continue;
			}

			sign = Relative_coord_check_stepover_Y(Sekibun_cell[mesh], x_int, residue, sw_in_cell, nlattice, DX, r_mesh, r);
			p_info[n][mesh].sign = sign;
			p_info[n][mesh].dmy_ry = (r_mesh[1] + sign*L_particle[1]);
			memcpy(p_info[n][mesh].r, r, sizeof(double) * 3);
			if (Range_coord(r_mesh, dmy_mesh)) {
				im = REALMODE_ARRAYINDEX(dmy_mesh[0], dmy_mesh[1], dmy_mesh[2]);
				dmyR = 0.;
				for (int d = 0; d < DIM; d++) {
					dmyR += SQ(r[d]);
				}
				dmyR = sqrt(dmyR);
				dmy_phi = Phi(dmyR, RADIUS) / MAX(phi_sum[im], 1.0);
				p_info[n][mesh].dmyR = sqrt(dmyR);
				p_info[n][mesh].dmy_phi = Phi(dmyR, RADIUS) / MAX(phi_sum[im], 1.0);
				skip_flag[n] = 0;
			} else {
				im = -1;
			}
			p_info[n][mesh].dmyR = dmyR;
			p_info[n][mesh].dmy_phi = dmy_phi;
			p_info[n][mesh].im = im;
		}
	}

	for (int n = 0; n < Particle_Number; n++) {
		dmy_rhop = RHO_particle[p[n].spec];
		if (SW_PT == rigid) rigidID = Particle_RigidID[n];

		for (int d = 0; d < DIM; d++) {
			xp[d] = p[n].x[d];
			vp[d] = p[n].v[d];
			omega_p[d] = p[n].omega[d];

			h_obl_force_local[(3 * p[n].id) + d] = 0.0;
			if (ROTATION) {
				h_obl_torque_local[(3 * p[n].id) + d] = 0.0;
			}
			forceg[d] = torqueg[d] = 0.0;
		}

		h_obl_volume_local[p[n].id] = 0.0;
		h_obl_Itrace_local[p[n].id] = 0.0;
		sw_in_cell = 1;

		if (skip_flag[n] == 1) continue;

		for (int mesh = 0; mesh < NP_domain; mesh++) {
			sign = p_info[n][mesh].sign;
			memcpy(r, p_info[n][mesh].r, sizeof(double) * 3);
			dmyR = p_info[n][mesh].dmyR;
			dmy_phi = p_info[n][mesh].dmy_phi;
			im = p_info[n][mesh].im;
			dmy_ry = p_info[n][mesh].dmy_ry;

			if (im >= 0) {
				Angular2v(omega_p, r, v_rot);

				for (int d = 0; d < DIM; d++) {
					if (!(d == 0)) {
						dmy_fp[d] = ((vp[d] + v_rot[d]) - u[d][im])*dmy_phi;
					} else {
						dmy_fp[d] = ((vp[d] - sign*Shear_rate_eff*L_particle[1] + v_rot[d]) - u[d][im])*dmy_phi;
					}
					h_obl_force_local[(3 * p[n].id) + d] += dmy_fp[d];
				}

				if (ROTATION) {
					h_obl_torque_local[(3 * p[n].id) + 0] += (r[1] * dmy_fp[2] - r[2] * dmy_fp[1]);
					h_obl_torque_local[(3 * p[n].id) + 1] += (r[2] * dmy_fp[0] - r[0] * dmy_fp[2]);
					h_obl_torque_local[(3 * p[n].id) + 2] += (r[0] * dmy_fp[1] - r[1] * dmy_fp[0]);
				}

				if (SW_PT == rigid) {
					for (int d = 0; d < DIM; d++) forceg[d] += dmy_fp[d];
					torqueg[0] += ((GRvecs[n][1] + r[1]) * dmy_fp[2] - (GRvecs[n][2] + r[2]) * dmy_fp[1]);
					torqueg[1] += ((GRvecs[n][2] + r[2]) * dmy_fp[0] - (GRvecs[n][0] + r[0]) * dmy_fp[2]);
					torqueg[2] += ((GRvecs[n][0] + r[0]) * dmy_fp[1] - (GRvecs[n][1] + r[1]) * dmy_fp[0]);
				}

				Hydro_force[im] += dmy_fp[0] * dmy_ry*dmy_rhop;//viscocity
				Hydro_force_new[im] += dmy_fp[0] * dmy_ry*dmy_rhop;

				h_obl_volume_local[p[n].id] += dmy_phi;
				h_obl_Itrace_local[p[n].id] += dmy_phi*SQ(dmyR);
			}
		}//mesh
	}// Particle_Number

	ierr = MPI_Allreduce(h_obl_force_local, h_obl_force, (DIM*Particle_Number), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_ERRORCHECK(ierr);
	ierr = MPI_Allreduce(h_obl_torque_local, h_obl_torque, (DIM*Particle_Number), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_ERRORCHECK(ierr);
	ierr = MPI_Allreduce(h_obl_Itrace_local, h_obl_Itrace, Particle_Number, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_ERRORCHECK(ierr);
	ierr = MPI_Allreduce(h_obl_volume_local, h_obl_volume, Particle_Number, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_ERRORCHECK(ierr);


	for (int n = 0; n < Update_Particle_Number; n++) {
		for (int d = 0; d < DIM; d++) {
			p[n].f_hydro[d] = (dmy * h_obl_force[(3 * p[n].id) + d]);
			if (ROTATION) {
				p[n].torque_hydro[d] = (dmy * h_obl_torque[(3 * p[n].id) + d]);
			}
		}
	}

	for (int n = 0; n < Particle_Number; n++) {
		h_obl_Itrace[n] *= 2.0 / 3.0;
		sum_force += h_obl_force[3 * n];
		sum_volume += h_obl_volume[n];
	}

	sum_volume /= RHO;

	for (int n = 0; n < Particle_Number; n++) {
		if (skip_flag[n] == 1) continue;
		int rigidID;
		dmy_rhop = RHO_particle[p[n].spec];
		if (SW_PT == rigid) rigidID = Particle_RigidID[n];

		for (int d = 0; d < DIM; d++) {
			xp[d] = p[n].x[d];
		}

		sw_in_cell = 1;

		for (int mesh = 0; mesh < NP_domain; mesh++) {
			sign = p_info[n][mesh].sign;
			memcpy(r, p_info[n][mesh].r, sizeof(double) * 3);
			dmyR = p_info[n][mesh].dmyR;
			dmy_phi = p_info[n][mesh].dmy_phi;
			im = p_info[n][mesh].im;
			dmy_ry = p_info[n][mesh].dmy_ry;

			if (im >= 0) {
				if (fabs(h_obl_volume[p[n].id]) > 1.0e-14 && fabs(h_obl_Itrace[p[n].id]) > 1.0e-14) {
					double dmy_stress_p = p[n].momentum_depend_fr[0] / h_obl_volume[p[n].id];
					Hydro_force[im] -= dmy_stress_p*dmy_ry*dmy_phi;//viscocity
					if (SW_PT != rigid) {
						double dmy_stress_v = -RHO*h_obl_force[(3 * p[n].id)] / h_obl_volume[p[n].id];
						double dmy_stress_w = -RHO*(h_obl_torque[(3 * p[n].id) + 1] * r[2] - h_obl_torque[(3 * p[n].id) + 2] * r[1]) / h_obl_Itrace[p[n].id];
						Hydro_force_new[im] += (dmy_stress_v + dmy_stress_w)*dmy_ry*dmy_phi;
					}
				}
				Hydro_force[im] -= sum_force*dmy_ry*dmy_phi / sum_volume;//viscocity

				if (SW_PT == rigid) {
					rigidID = Particle_RigidID[n];
					double dmy_stress_v = dVg[rigidID][0];
					double dmy_stress_w = dWg[rigidID][1] * (GRvecs[n][2] + r[2]) - dWg[rigidID][2] * (GRvecs[n][1] + r[1]);
					Hydro_force_new[im] += (dmy_stress_v + dmy_stress_w)*dmy_ry*dmy_phi*dmy_rhop;
				}
			}
		}
	}

	if (SW_PT == rigid) {
		for (int rigidID = 0; rigidID < Rigid_Number; rigidID++) {
			for (int d = 0; d < DIM; d++) {
				dVg[rigidID][d] = jikan.dt_md * Rigid_IMasses[rigidID] * forceGs[rigidID][d];
			}

			int rigid_first = Rigid_Particle_Cumul[rigidID];
			double Ibody[DIM] = { Rigid_Moments_body[rigidID][0][0],
								 Rigid_Moments_body[rigidID][1][1],
								 Rigid_Moments_body[rigidID][2][2] };
			MD_solver_omega_Euler_update(dWg[rigidID], omegaGs[rigidID], torqueGs[rigidID], Ibody, p[rigid_first].q, jikan.dt_md);
		}
	}
#else
static const double dmy0 = -DX3*RHO;
double dmy = dmy0 / jikan.dt_fluid;

int *nlattice;
nlattice = Ns;

double xp[DIM], vp[DIM], omega_p[DIM];
int x_int[DIM];
double residue[DIM];
int sw_in_cell;
double force[DIM];
double torque[DIM];
int r_mesh[DIM];
double r[DIM];
double dmy_fp[DIM];
double x[DIM];
double dmyR;
double dmy_phi;
double dmy_rhop;
double dmy_ry;
double v_rot[DIM];
double volume[Particle_Number];
double Itrace[Particle_Number];
int sign;
int im;
double sum_force = 0.0;
double sum_volume = 0.0;

double forceg[DIM];
double torqueg[DIM];
double dVg[Rigid_Number][DIM];
double dWg[Rigid_Number][DIM];
// initialize forceGs and torqueGs
for (int rigidID = 0; rigidID<Rigid_Number; rigidID++) {
	for (int d = 0; d<DIM; d++) {
		forceGs[rigidID][d] = 0.0;
		torqueGs[rigidID][d] = 0.0;
	}
}

Reset_phi(Hydro_force);
Reset_phi(Hydro_force_new);
#pragma omp parallel for \
   private(xp,vp,omega_p,x_int,residue,sw_in_cell,force,torque,r_mesh,r,dmy_fp,x,\
           dmyR,dmy_phi,dmy_rhop,dmy_ry,v_rot,sign,im,forceg,torqueg) 
for (int n = 0; n < Particle_Number; n++) {
	int rigidID;
	dmy_rhop = RHO_particle[p[n].spec];
	if (SW_PT == rigid) rigidID = Particle_RigidID[n];

	for (int d = 0; d < DIM; d++) {
		xp[d] = p[n].x[d];
		vp[d] = p[n].v[d];
		omega_p[d] = p[n].omega[d];

		force[d] = torque[d] = 0.0;
		forceg[d] = torqueg[d] = 0.0;
	}

	volume[n] = 0.0;
	Itrace[n] = 0.0;
	sw_in_cell
		= Particle_cell(xp, DX, x_int, residue);
	sw_in_cell = 1;

	for (int mesh = 0; mesh < NP_domain; mesh++) {
		sign = Relative_coord_check_stepover_Y(Sekibun_cell[mesh], x_int, residue,
			sw_in_cell, nlattice, DX, r_mesh, r);
		dmyR = 0.;
		for (int d = 0; d<DIM; d++) {
			x[d] = r_mesh[d] * DX;
			dmyR += SQ(r[d]);
		}
		im = (r_mesh[0] * NY*NZ_ + r_mesh[1] * NZ_ + r_mesh[2]);

		dmyR = sqrt(dmyR);
		dmy_phi = Phi(dmyR, RADIUS) / MAX(phi_sum[im], 1.0);

		Angular2v(omega_p, r, v_rot);

		for (int d = 0; d < DIM; d++) {
			if (!(d == 0)) {
				dmy_fp[d] = ((vp[d] + v_rot[d]) - u[d][im])*dmy_phi;
			} else {
				dmy_fp[d] = ((vp[d] - sign*Shear_rate_eff*L_particle[1] + v_rot[d])
					- u[d][im])*dmy_phi;
			}
			force[d] += dmy_fp[d];
		}
		{// torque
			torque[0] += (r[1] * dmy_fp[2] - r[2] * dmy_fp[1]);
			torque[1] += (r[2] * dmy_fp[0] - r[0] * dmy_fp[2]);
			torque[2] += (r[0] * dmy_fp[1] - r[1] * dmy_fp[0]);
		}
		if (SW_PT == rigid) {
			for (int d = 0; d<DIM; d++) forceg[d] += dmy_fp[d];
			torqueg[0] += ((GRvecs[n][1] + r[1]) * dmy_fp[2] - (GRvecs[n][2] + r[2]) * dmy_fp[1]);
			torqueg[1] += ((GRvecs[n][2] + r[2]) * dmy_fp[0] - (GRvecs[n][0] + r[0]) * dmy_fp[2]);
			torqueg[2] += ((GRvecs[n][0] + r[0]) * dmy_fp[1] - (GRvecs[n][1] + r[1]) * dmy_fp[0]);
		}

		dmy_ry = (r_mesh[1] + sign*L_particle[1]);
#pragma omp atomic
		Hydro_force[im] += dmy_fp[0] * dmy_ry*dmy_rhop;//viscocity
#pragma omp atomic
		Hydro_force_new[im] += dmy_fp[0] * dmy_ry*dmy_rhop;

		volume[n] += dmy_phi;
		Itrace[n] += dmy_phi*SQ(dmyR);
	}//mesh
	Itrace[n] *= 2.0 / 3.0;

	for (int d = 0; d < DIM; d++) {
		p[n].f_hydro[d] = (dmy * force[d]);
	}
	if (ROTATION) {
		for (int d = 0; d < DIM; d++) {
			p[n].torque_hydro[d] = (dmy * torque[d]);
		}
	}
	if (SW_PT == rigid) {
		for (int d = 0; d<DIM; d++) {
#pragma omp atomic
			forceGs[rigidID][d] += dmy * forceg[d];

#pragma omp atomic
			torqueGs[rigidID][d] += dmy * torqueg[d];
		}
	}

	for (int mesh = 0; mesh < NP_domain; mesh++) {
		sign = Relative_coord_check_stepover_Y(Sekibun_cell[mesh], x_int, residue,
			sw_in_cell, nlattice, DX, r_mesh, r);
		im = (r_mesh[0] * NY*NZ_ + r_mesh[1] * NZ_ + r_mesh[2]);

		dmyR = 0.;
		for (int d = 0; d<DIM; d++) {
			dmyR += SQ(r[d]);
		}
		dmyR = sqrt(dmyR);
		dmy_phi = Phi(dmyR, RADIUS) / MAX(phi_sum[im], 1.0);

		dmy_ry = (r_mesh[1] + sign*L_particle[1]);
		double dmy_stress_p = p[n].momentum_depend_fr[0] / volume[n];

#pragma omp atomic
		Hydro_force[im] -= dmy_stress_p*dmy_ry*dmy_phi;//viscocity

		if (SW_PT != rigid) {
			double dmy_stress_v = -RHO*force[0] / volume[n];
			double dmy_stress_w = -RHO*(torque[1] * r[2] - torque[2] * r[1]) / Itrace[n];
#pragma omp atomic
			Hydro_force_new[im] += (dmy_stress_v + dmy_stress_w)*dmy_ry*dmy_phi;

		}
	}

# pragma omp critical
	{
		sum_force += force[0];
		sum_volume += volume[n];
	}
}// Particle_Number
sum_volume /= RHO;

if (SW_PT == rigid) {
#pragma omp parallel for 
	for (int rigidID = 0; rigidID < Rigid_Number; rigidID++) {
		for (int d = 0; d < DIM; d++) {
			dVg[rigidID][d] = jikan.dt_md * Rigid_IMasses[rigidID] * forceGs[rigidID][d];
		}

		int rigid_first = Rigid_Particle_Cumul[rigidID];
		double Ibody[DIM] = { Rigid_Moments_body[rigidID][0][0],
			Rigid_Moments_body[rigidID][1][1],
			Rigid_Moments_body[rigidID][2][2] };
		MD_solver_omega_Euler_update(dWg[rigidID], omegaGs[rigidID], torqueGs[rigidID],
			Ibody, p[rigid_first].q, jikan.dt_md);
	}
}

#pragma omp parallel for \
   private(xp,x_int,residue,sw_in_cell,r_mesh,r,x,dmyR,dmy_phi,dmy_ry,sign,im,dmy_rhop) 
for (int n = 0; n < Particle_Number; n++) {
	int rigidID;
	dmy_rhop = RHO_particle[p[n].spec];
	for (int d = 0; d < DIM; d++) {
		xp[d] = p[n].x[d];
	}
	sw_in_cell = Particle_cell(xp, DX, x_int, residue);// {1,0} が返ってくる
	sw_in_cell = 1;

	for (int mesh = 0; mesh < NP_domain; mesh++) {
		sign = Relative_coord_check_stepover_Y(Sekibun_cell[mesh], x_int, residue, \
			sw_in_cell, nlattice, DX, r_mesh, r);
		im = (r_mesh[0] * NY*NZ_ + r_mesh[1] * NZ_ + r_mesh[2]);

		dmyR = 0;
		for (int d = 0; d<DIM; d++) {
			x[d] = r_mesh[d] * DX;
			dmyR += SQ(r[d]);
		}
		dmyR = sqrt(dmyR);
		dmy_phi = Phi(dmyR, RADIUS) / MAX(phi_sum[im], 1.0);

		dmy_ry = (r_mesh[1] + sign*L_particle[1]);
#pragma omp atomic      
		Hydro_force[im] -= sum_force*dmy_ry*dmy_phi / sum_volume;//viscocity

		if (SW_PT == rigid) {
			rigidID = Particle_RigidID[n];
			double dmy_stress_v = dVg[rigidID][0];
			double dmy_stress_w = dWg[rigidID][1] * (GRvecs[n][2] + r[2]) - dWg[rigidID][2] * (GRvecs[n][1] + r[1]);
#pragma omp atomic
			Hydro_force_new[im] += (dmy_stress_v + dmy_stress_w)*dmy_ry*dmy_phi*dmy_rhop;
		}
	}
}//Particle_Number
#endif
}

