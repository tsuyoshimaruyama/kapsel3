/*!
  \file rigid.h
  \brief Rigid particles as stiff chains
  \authors T. Kobiki
  \date 2012/12/29
  \version 1.0
*/
#ifndef RIGID_H
#define RIGID_H

#include "input.h"
#include "Matrix_Inverse.h"
#include "matrix_diagonal.h"

/*!
  \brief Initialize the geometry and center of mass position for each
  of the rigid particles
 */
inline void init_set_xGs(Particle *p){
	// initialize xGs[][]
#pragma omp parallel for schedule(dynamic, 1)
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		for(int d=0; d<DIM; d++) xGs[rigidID][d] = 0.0;
	}
	
	int rigidID;
	
        //center of mass
#pragma omp parallel for schedule(dynamic, 1) private(rigidID)
	for(int n=0; n<Particle_Number; n++){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++) {
#pragma omp atomic
                  xGs[rigidID][d] += p[n].x[d];
                }
	}

#pragma omp parallel for schedule(dynamic, 1)
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		for(int d=0; d<DIM; d++) xGs[rigidID][d] /= Rigid_Particle_Numbers[rigidID];
	}

        //position vectors from center of mass to individual beads
#pragma omp parallel for schedule(dynamic, 1) private(rigidID)
	for(int n=0; n<Particle_Number; n++){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++) {
                  GRvecs[n][d] = p[n].x[d] - xGs[rigidID][d];
                }
	}
}

/*!
  \brief Update center of mass position after all beads have been
  rigidly displaced
 */
inline void set_xGs(Particle *p){
	int rigid_first_n = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:rigid_first_n)
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		for(int d=0; d<DIM; d++){
			xGs[rigidID][d] = p[rigid_first_n].x[d] - GRvecs[rigid_first_n][d];
			xGs[rigidID][d] = fmod(xGs[rigidID][d] + L_particle[d], L_particle[d]);
		}
		rigid_first_n += Rigid_Particle_Numbers[rigidID];
	}
	
	//check(for debug)
	if(rigid_first_n != Particle_Number) fprintf(stderr, "debug: set_xGs() error");
}

//diagonalize inertia tensors
inline void set_Rigid_Coordinates(Particle *p){
#pragma omp parallel for schedule(dynamic, 1)
  for(int rigidID = 0; rigidID < Rigid_Number; rigidID++){
    double **eigen_vector;
    double eigen_value[DIM];
    int iter;
    double dmy_R[DIM][DIM];
    quaternion dmy_q;
    eigen_vector = alloc_2d_double(DIM, DIM);
    jacobi(Rigid_Moments[rigidID], eigen_vector, eigen_value, iter, DIM);
    M_coordinate_frame(eigen_vector[2], eigen_vector[0], eigen_vector[1]);
    for(int i = 0; i < DIM; i++){
      for(int j = 0; j < DIM; j++){
        dmy_R[j][i] = eigen_vector[i][j];
      }
    }
    rm_rqtn(dmy_q, dmy_R);

    for(int n = Rigid_Particle_Cumul[rigidID]; n < Rigid_Particle_Cumul[rigidID+1]; n++){
      qtn_init(p[n].q, dmy_q);
    }
  }

#pragma omp parallel for schedule(dynamic, 1)
  for(int n = 0; n < Particle_Number; n++){
    rigid_body_rotation(GRvecs_body[n], GRvecs[n], p[n].q, SPACE2BODY);
  }
}

/*!
  \brief Compute the change in relative bead positions (from COM) due to rotation of
  the rigid body
 */
inline void solver_GRvecs(const CTime &jikan, string CASE){
	int rigidID;
	double bufvec[DIM];
        double dwt[DIM];
	if(CASE == "Euler"){
#pragma omp parallel for schedule(dynamic, 1) private(rigidID, bufvec, dwt)
		for(int n=0; n<Particle_Number; n++){
			rigidID = Particle_RigidID[n];
                        dwt[0] = omegaGs[rigidID][0];
                        dwt[1] = omegaGs[rigidID][1];
                        dwt[2] = omegaGs[rigidID][2];

			bufvec[0] = (dwt[1]*GRvecs[n][2] - dwt[2]*GRvecs[n][1]) * jikan.dt_md;
			bufvec[1] = (dwt[2]*GRvecs[n][0] - dwt[0]*GRvecs[n][2]) * jikan.dt_md;
			bufvec[2] = (dwt[0]*GRvecs[n][1] - dwt[1]*GRvecs[n][0]) * jikan.dt_md;
			for(int d=0; d<DIM; d++) GRvecs[n][d] += bufvec[d];
		}
	}
	else if(CASE == "AB2"){
#pragma omp parallel for schedule(dynamic, 1) private(rigidID, bufvec, dwt)
		for(int n=0; n<Particle_Number; n++){
			rigidID = Particle_RigidID[n];
                        dwt[0] = 3.0*omegaGs[rigidID][0] - omegaGs_old[rigidID][0];
                        dwt[1] = 3.0*omegaGs[rigidID][1] - omegaGs_old[rigidID][1];
                        dwt[2] = 3.0*omegaGs[rigidID][2] - omegaGs_old[rigidID][2];

			bufvec[0] = (dwt[1]*GRvecs[n][2] - dwt[2]*GRvecs[n][1]) * jikan.hdt_md;
			bufvec[1] = (dwt[2]*GRvecs[n][0] - dwt[0]*GRvecs[n][2]) * jikan.hdt_md;
			bufvec[2] = (dwt[0]*GRvecs[n][1] - dwt[1]*GRvecs[n][0]) * jikan.hdt_md;
			for(int d=0; d<DIM; d++) GRvecs[n][d] += bufvec[d];
		}
	}
	else{
		fprintf(stderr, "error, string CASE in solver_GRvecs()");
		exit_job(EXIT_FAILURE);
	}
}

/*! 
  \brief set Rigid_Masses and Rigid_IMasses and Rigid_Moments and Rigid_IMoments
  \note Dont use it before (init_set_xGs) or (set_xGs and solver_GRvecs)
*/
inline void set_Rigid_MMs(Particle *p){
	// initialize Rigid_Masses and Rigid_Moments and Rigid_IMoments
#pragma omp parallel for schedule(dynamic, 1)
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		Rigid_Masses[rigidID] = 0.0;
		for(int row=0; row<DIM; row++){
			for(int column=0; column<DIM; column++){
				Rigid_Moments[rigidID][row][column] = 0.0;
				Rigid_IMoments[rigidID][row][column] = 0.0;
			}
		}
	}
	
	int rigidID;
	double dmy_mass, dmy_moi;
#pragma omp parallel for schedule(dynamic, 1) private(rigidID, dmy_mass, dmy_moi)
	for(int n=0; n<Particle_Number; n++){
		rigidID = Particle_RigidID[n];
                dmy_mass = MASS[ RigidID_Components[rigidID] ];
                dmy_moi = MOI[ RigidID_Components[rigidID] ];
                double &ri_x = GRvecs[n][0];
                double &ri_y = GRvecs[n][1];
                double &ri_z = GRvecs[n][2];
                
		// mass of p[n] adds to Rigid_Masses[rigidID];
#pragma omp atomic
		Rigid_Masses[rigidID] += dmy_mass;

		// moment of p[n] adds to Rigid_Moments[rigidID]
#pragma omp atomic
		Rigid_Moments[rigidID][0][0] += dmy_moi + dmy_mass * ( ri_y*ri_y + ri_z*ri_z );
#pragma omp atomic
		Rigid_Moments[rigidID][0][1] += dmy_mass * ( -ri_x*ri_y );
#pragma omp atomic
		Rigid_Moments[rigidID][0][2] += dmy_mass * ( -ri_x*ri_z );

#pragma omp atomic
		Rigid_Moments[rigidID][1][0] += dmy_mass * ( -ri_y*ri_x );
#pragma omp atomic
		Rigid_Moments[rigidID][1][1] += dmy_moi + dmy_mass * ( ri_z*ri_z + ri_x*ri_x );
#pragma omp atomic
		Rigid_Moments[rigidID][1][2] += dmy_mass * ( -ri_y*ri_z );

#pragma omp atomic
		Rigid_Moments[rigidID][2][0] += dmy_mass * ( -ri_z*ri_x );
#pragma omp atomic
		Rigid_Moments[rigidID][2][1] += dmy_mass * ( -ri_z*ri_y );
#pragma omp atomic
		Rigid_Moments[rigidID][2][2] += dmy_moi + dmy_mass * ( ri_x*ri_x + ri_y*ri_y );
	}
	
	char str[256];
	// inverse
#pragma omp parallel for schedule(dynamic, 1)
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		Rigid_IMasses[rigidID] = 1. / Rigid_Masses[rigidID];
		Matrix_Inverse(Rigid_Moments[rigidID], Rigid_IMoments[rigidID], DIM);
		check_Inverse(Rigid_Moments[rigidID], Rigid_IMoments[rigidID], DIM);
	}
}

/*!
  \brief Set rigid body velocities for individual beads
 */
inline void set_particle_vomegas(Particle *p){

	int rigidID;
#pragma omp parallel for schedule(dynamic, 1) private(rigidID)
	for(int n=0; n<Particle_Number; n++){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++){
			p[n].omega[d] = omegaGs[rigidID][d];
		}
		p[n].v[0] = velocityGs[rigidID][0] + omegaGs[rigidID][1]*GRvecs[n][2] - omegaGs[rigidID][2]*GRvecs[n][1];
		p[n].v[1] = velocityGs[rigidID][1] + omegaGs[rigidID][2]*GRvecs[n][0] - omegaGs[rigidID][0]*GRvecs[n][2];
		p[n].v[2] = velocityGs[rigidID][2] + omegaGs[rigidID][0]*GRvecs[n][1] - omegaGs[rigidID][1]*GRvecs[n][0];
	}
}

/*!
  \brief Update rigid velocities (VelocityGs) and angular velocities (OmegaGs)
  \note set_Rigid_VOGs() after calculating xGs, Rigid_IMoments, forceGs and torqueGs!!
*/
inline void calc_Rigid_VOGs(Particle *p, const CTime &jikan, string CASE){
	set_Rigid_MMs(p);
	
	int rigidID;
	//calc forceGrs, forceGvs
#pragma omp parallel for schedule(dynamic, 1) private(rigidID)
	for(int n=0; n<Particle_Number; n++){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++){
#pragma omp atomic
			forceGrs[rigidID][d] += p[n].fr[d];
			//forceGvs[rigidID][d] += p[n].fv[d];
			//torqueGrs[rigidID][d] += p[n].torquer[d];
			//torqueGvs[rigidID][d] += p[n].torquev[d];	//fv, torquer and torquev are constant zero.
		}
	}
	//set olds
#pragma omp parallel for schedule(dynamic, 1)
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		for(int d=0; d<DIM; d++){
			velocityGs_old[rigidID][d] = velocityGs[rigidID][d];
			omegaGs_old[rigidID][d] = omegaGs[rigidID][d];
		}
	}
	
	//calc velocityGs and omegaGs
	if(CASE == "Euler"){
#pragma omp parallel for schedule(dynamic, 1)
		for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
			if(Rigid_Motions[ RigidID_Components[rigidID] ] == 0) continue;	// if "fix"
			for(int d1=0; d1<DIM; d1++){
				velocityGs[rigidID][d1] += jikan.dt_md * Rigid_IMasses[rigidID]
                                  * ( forceGs[rigidID][d1] + forceGrs[rigidID][d1] );
				for(int d2=0; d2<DIM; d2++) omegaGs[rigidID][d1] += jikan.dt_md * Rigid_IMoments[rigidID][d1][d2]
                                                              * ( torqueGs[rigidID][d2] );
			}
		}
	
	}else if(CASE == "AB2"){
#pragma omp parallel for schedule(dynamic, 1)
		for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
			if(Rigid_Motions[ RigidID_Components[rigidID] ] == 0) continue;	// if "fix"
			for(int d1=0; d1<DIM; d1++){
				velocityGs[rigidID][d1] += jikan.hdt_md * Rigid_IMasses[rigidID]
                                  * ( 2. * forceGs[rigidID][d1] + forceGrs[rigidID][d1] + forceGrs_previous[rigidID][d1]
                                      );
				for(int d2=0; d2<DIM; d2++) omegaGs[rigidID][d1] += jikan.hdt_md * Rigid_IMoments[rigidID][d1][d2]
                                                              * ( 2. * torqueGs[rigidID][d2] );
			}
		}
	
	}else{
		fprintf(stderr, "error, string CASE in calc_Rigid_VOGs()");
		exit_job(EXIT_FAILURE);
	}


		
	// renew old previous and initialize
#pragma omp parallel for schedule(dynamic, 1)
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		for(int d=0; d<DIM; d++){
			 forceGrs_previous[rigidID][d] = forceGrs[rigidID][d];
			 forceGrs[rigidID][d] = 0.0;
			 forceGs[rigidID][d] = 0.0;
			 torqueGs[rigidID][d] = 0.0;
		}
	}
}
#endif
