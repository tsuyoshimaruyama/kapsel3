//
//	rigid.h 12/12/29 T.K
//
#ifndef RIGID_H
#define RIGID_H

#include "input.h"
#include "Matrix_Inverse.h"


inline void set_xGs(Particle *p){
	// initialize xGs[][]
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		for(int d=0; d<DIM; d++) xGs[rigidID][d] = 0.0;
	}
	
	int rigidID;
	
	for(int n=0; n<Particle_Number; n++){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++) xGs[rigidID][d] += p[n].x[d];
	}
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		for(int d=0; d<DIM; d++) xGs[rigidID][d] /= Rigid_Particle_Numbers[rigidID];
	}
}


// set Rigid_Masses and Rigid_IMasses and Rigid_Moments and Rigid_IMoments
inline void set_Rigid_MMs(Particle *p){
	set_xGs(p);
	// initialize Rigid_Masses and Rigid_Moments and Rigid_IMoments
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		Rigid_Masses[rigidID] = 0.0;
		for(int row=0; row<DIM; row++){
			for(int column=0; column<DIM; column++){
				Rigid_Moments[rigidID][row][column] = 0.0;
				Rigid_IMoments[rigidID][row][column] = 0.0;
			}
		}
	}
	
	double GRvec[DIM];
	int rigidID;
	
	for(int n=0; n<Particle_Number; n++){
		rigidID = Particle_RigidID[n];
		// set GRvec
		for(int d=0; d<DIM; d++){
			GRvec[d] = p[n].x[d] - xGs[rigidID][d];
		}
		// mass of p[n] adds to Rigid_Masses[rigidID];
		Rigid_Masses[rigidID] += MASS[ RigidID_Components[rigidID] ];
		// moment of p[n] adds to Rigid_Moments[rigidID]
		Rigid_Moments[rigidID][0][0] += MOI[ RigidID_Components[rigidID] ]
													+ MASS[ RigidID_Components[rigidID] ]
														* ( GRvec[1]*GRvec[1] + GRvec[2]*GRvec[2] );
		Rigid_Moments[rigidID][0][1] += MASS[ RigidID_Components[rigidID] ]
														* ( - GRvec[0]*GRvec[1] );
		Rigid_Moments[rigidID][0][2] += MASS[ RigidID_Components[rigidID] ]
														* ( - GRvec[0]*GRvec[2] );
		Rigid_Moments[rigidID][1][0] += MASS[ RigidID_Components[rigidID] ]
														* ( - GRvec[1]*GRvec[0] );
		Rigid_Moments[rigidID][1][1] += MOI[ RigidID_Components[rigidID] ]
													+ MASS[ RigidID_Components[rigidID] ]
														* ( GRvec[2]*GRvec[2] + GRvec[0]*GRvec[0] );
		Rigid_Moments[rigidID][1][2] += MASS[ RigidID_Components[rigidID] ]
														* ( - GRvec[1]*GRvec[2] );
		Rigid_Moments[rigidID][2][0] += MASS[ RigidID_Components[rigidID] ]
														* ( - GRvec[2]*GRvec[0] );
		Rigid_Moments[rigidID][2][1] += MASS[ RigidID_Components[rigidID] ]
														* ( - GRvec[2]*GRvec[1] );
		Rigid_Moments[rigidID][2][2] += MOI[ RigidID_Components[rigidID] ]
													+ MASS[ RigidID_Components[rigidID] ]
														* ( GRvec[0]*GRvec[0] + GRvec[1]*GRvec[1] );
	}
	
	char str[256];
	// inverse
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		Rigid_IMasses[rigidID] = 1. / Rigid_Masses[rigidID];
		Matrix_Inverse(Rigid_Moments[rigidID], Rigid_IMoments[rigidID], DIM);
		
//		sprintf(str, "Rigid_Moments[%d].csv", rigidID);
//		Matrix_output(Rigid_Moments[rigidID], DIM, str); 
//		sprintf(str, "Rigid_IMoments[%d].csv", rigidID);
//		Matrix_output(Rigid_IMoments[rigidID], DIM, str); 
		
		check_Inverse(Rigid_Moments[rigidID], Rigid_IMoments[rigidID], DIM);
	}

}

inline void set_particle_vomegas(Particle *p){
	set_xGs(p);
	double GRvec[DIM];
	int rigidID;
	for(int n=0; n<Particle_Number; n++){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++){
			p[n].omega[d] = omegaGs[rigidID][d];
			GRvec[d] = p[n].x[d] - xGs[rigidID][d];
		}
		p[n].v[0] = velocityGs[rigidID][0] + omegaGs[rigidID][1]*GRvec[2] - omegaGs[rigidID][2]*GRvec[1];
		p[n].v[1] = velocityGs[rigidID][1] + omegaGs[rigidID][2]*GRvec[0] - omegaGs[rigidID][0]*GRvec[2];
		p[n].v[2] = velocityGs[rigidID][2] + omegaGs[rigidID][0]*GRvec[1] - omegaGs[rigidID][1]*GRvec[0];
	}
}

// set VelocityGs and OmegaGs
// ### set_Rigid_VOGs() after calculating xGs, Rigid_IMoments, forceGs and torqueGs!! ###
inline void calc_Rigid_VOGs(Particle *p, const CTime &jikan, string CASE){
	set_Rigid_MMs(p);
	
	int rigidID;
	//calc forceGrs, forceGvs
	for(int n=0; n<Particle_Number; n++){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++){
			forceGrs[rigidID][d] += p[n].fr[d];
			forceGvs[rigidID][d] += p[n].fv[d];
			//torqueGrs[rigidID][d] += p[n].torquer[d]; 	//検討が必要。
			//torqueGvs[rigidID][d] += p[n].torqueGvs[d];	//おそらくtorquerとtorquevを計算しているところで計算して足し合わせるのが確実だが、
															//その辺りは理解できていないので未実装
		}
	}
	
	//calc velocityGs and omegaGs
	if(CASE == "Euler"){
		for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
			if(Rigid_Motions[rigidID] == 0) continue;	// if "fix"
			for(int d1=0; d1<DIM; d1++){
				velocityGs[rigidID][d1] += jikan.dt_md * Rigid_IMasses[rigidID]
											* ( forceGs[rigidID][d1] + forceGrs[rigidID][d1] + forceGvs[rigidID][d1] );
				for(int d2=0; d2<DIM; d2++) omegaGs[rigidID][d1] += jikan.dt_md * Rigid_IMoments[rigidID][d1][d2]
																	* ( torqueGs[rigidID][d2] + torqueGrs[rigidID][d2] + torqueGvs[rigidID][d2] );
			}
		}
	
	}else if(CASE == "Euler_hydro"){
		for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
			if(Rigid_Motions[rigidID] == 0) continue;	// if "fix"
			for(int d1=0; d1<DIM; d1++){
				velocityGs[rigidID][d1] += jikan.dt_md * Rigid_IMasses[rigidID]
											* ( forceGs[rigidID][d1] + forceGvs[rigidID][d1] 
												+ 0.5*(forceGrs[rigidID][d1] + forceGrs_previous[rigidID][d1])
											);
				for(int d2=0; d2<DIM; d2++) omegaGs[rigidID][d1] += jikan.dt_md * Rigid_IMoments[rigidID][d1][d2]
																	* ( torqueGs[rigidID][d2] + torqueGvs[rigidID][d2]
																		+ 0.5*(torqueGrs[rigidID][d2] + torqueGrs_previous[rigidID][d2])
																	);
			}
		}
	
	}else if(CASE == "AB2_hydro"){
		for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
			if(Rigid_Motions[rigidID] == 0) continue;	// if "fix"
			for(int d1=0; d1<DIM; d1++){
				velocityGs[rigidID][d1] += jikan.hdt_md * Rigid_IMasses[rigidID]
											* ( 2. * forceGs[rigidID][d1]
												+ 3. * forceGvs[rigidID][d1] - forceGvs_previous[rigidID][d1]
												+ forceGrs[rigidID][d1] + forceGrs_previous[rigidID][d1]
											);
				for(int d2=0; d2<DIM; d2++) omegaGs[rigidID][d1] += jikan.hdt_md * Rigid_IMoments[rigidID][d1][d2]
																	* ( 2. * torqueGs[rigidID][d2]
																		+ 3. * torqueGvs[rigidID][d2] - torqueGvs_previous[rigidID][d2]
																		+ torqueGrs[rigidID][d2] + torqueGrs_previous[rigidID][d2]
																	);
			}
		}
	
	}else{
		fprintf(stderr, "error, string CASE in calc_Rigid_VOGs()");
		exit_job(EXIT_FAILURE);
	}

    
	// renew old previous and initialize
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		for(int d=0; d<DIM; d++){
			 forceGrs_previous[rigidID][d] = forceGrs[rigidID][d];
			 forceGrs[rigidID][d] = 0.0;
			 forceGvs_previous[rigidID][d] = forceGvs[rigidID][d];
			 forceGvs[rigidID][d] = 0.0;
			 forceGs[rigidID][d] = 0.0;
			 torqueGrs_previous[rigidID][d] = torqueGrs[rigidID][d];
			 torqueGrs[rigidID][d] = 0.0;
			 torqueGvs_previous[rigidID][d] = torqueGvs[rigidID][d];
			 torqueGvs[rigidID][d] = 0.0;
			 torqueGs[rigidID][d] = 0.0;
		}
	}
}
			


#endif