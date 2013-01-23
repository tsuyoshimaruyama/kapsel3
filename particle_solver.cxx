///
// $Id: particle_solver.cxx,v 1.1 2006/06/27 18:41:28 nakayama Exp $
//
#include "particle_solver.h"

void MD_solver_position_Euler(Particle *p, const CTime &jikan){
#pragma omp parallel for schedule(dynamic, 1)
    for(int n=0; n<Particle_Number; n++){
	for(int d=0; d<DIM; d++){
	    p[n].x[d] += jikan.dt_md * p[n].v[d];
	    p[n].x[d] = fmod(p[n].x[d]+L_particle[d] , L_particle[d]);
	    
	    assert(p[n].x[d] >= 0);
	    assert(p[n].x[d] < L[d]);
	}
    }
    set_xGs(p);		// T.K 13/01/20
}
void MD_solver_position_AB2(Particle *p, const CTime &jikan){
#pragma omp parallel for schedule(dynamic, 1)
    for(int n=0; n<Particle_Number; n++){
	for(int d=0; d<DIM; d++){
	    p[n].x[d] += jikan.hdt_md * (3.*p[n].v[d] - p[n].v_old[d]);
	    p[n].x[d] = fmod(p[n].x[d]+L_particle[d] , L_particle[d]);
	    
	    assert(p[n].x[d] >= 0);
	    assert(p[n].x[d] < L[d]);
	}
    }
    set_xGs(p);		// T.K 13/01/20
}
void MD_solver_velocity_Euler(Particle *p, const CTime &jikan){
  Force(p);
	calc_Rigid_VOGs(p, jikan, "Euler");	// T.K 13/01/04
    double dmy;
    double dmy_rot;
    // T.K 13/01/04
    int rigidID;
    double GRvec[DIM];
    #pragma omp parallel for schedule(dynamic,1) private(dmy, dmy_rot)
  for(int n=0; n< Particle_Number; n++){
    //double dmy = jikan.dt_md * IMASS[p[n].spec];
    //double dmy_rot = jikan.dt_md * IMOI[p[n].spec];
    dmy = jikan.dt_md * IMASS[p[n].spec];
    dmy_rot = jikan.dt_md * IMOI[p[n].spec];
    
    // T.K 13/01/04
	if(SW_PT == rigid){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++) GRvec[d] = p[n].x[d] - xGs[rigidID][d];
	}
	
    for(int d=0; d<DIM; d++){  
      {
	p[n].v_old[d] = p[n].v[d];
	
	// T.K 13/01/04
	if(SW_PT == rigid){
		p[n].v[0] = velocityGs[rigidID][0] + omegaGs[rigidID][1]*GRvec[2] - omegaGs[rigidID][2]*GRvec[1];
		p[n].v[1] = velocityGs[rigidID][1] + omegaGs[rigidID][2]*GRvec[0] - omegaGs[rigidID][0]*GRvec[2];
		p[n].v[2] = velocityGs[rigidID][2] + omegaGs[rigidID][0]*GRvec[1] - omegaGs[rigidID][1]*GRvec[0];
	}
	else{
		p[n].v[d] += ( dmy * 
			       ( p[n].f_hydro[d] + p[n].fr[d] + p[n].fv[d] )
			       );
	}
	
	p[n].fr_previous[d] = p[n].fr[d];
	p[n].fr[d] = 0.0;
	p[n].fv_previous[d] = p[n].fv[d];
	p[n].fv[d] = 0.0;
	//p[n].f_hydro_previous[d] = p[n].f_hydro[d];
	p[n].f_hydro[d] = 0.0;
      }
      {
	p[n].omega_old[d] = p[n].omega[d];
	
	// T.K 13/01/04
	if(SW_PT == rigid){
		p[n].omega[d] = omegaGs[rigidID][d];
	}
	else{
		p[n].omega[d] += ( dmy_rot 
				   * ( p[n].torque_hydro[d] 
				       + p[n].torquer[d] 
				       + p[n].torquev[d] )
				   );
	}
	
	p[n].torquer_previous[d] = p[n].torquer[d];
	p[n].torquer[d] = 0.0;
	p[n].torquev_previous[d] = p[n].torquev[d];
	p[n].torquev[d] = 0.0;
	//p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
	p[n].torque_hydro[d] = 0.0;
      }
    }
  }
}
void MD_solver_velocity_Euler_hydro(Particle *p, const CTime &jikan){
  Force(p);
  calc_Rigid_VOGs(p, jikan, "Euler_hydro");	// T.K 13/01/04
    double dmy;
    double dmy_rot;
    // T.K 13/01/04
    int rigidID;
    double GRvec[DIM];
#pragma omp parallel for schedule(dynamic,1) private(dmy, dmy_rot)
  for(int n=0; n< Particle_Number; n++){
    //static double dmy = jikan.dt_md * IMASS[p[n].spec];
    //static double dmy_rot = jikan.dt_md * IMOI[p[n].spec];
    dmy = jikan.dt_md * IMASS[p[n].spec];
    dmy_rot = jikan.dt_md * IMOI[p[n].spec];
    
    // T.K 13/01/04
	if(SW_PT == rigid){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++) GRvec[d] = p[n].x[d] - xGs[rigidID][d];
	}
	
    for(int d=0; d<DIM; d++){  
      {
	p[n].v_old[d] = p[n].v[d];
	
	// T.K 13/01/04
	if(SW_PT == rigid){
		p[n].v[0] = velocityGs[rigidID][0] + omegaGs[rigidID][1]*GRvec[2] - omegaGs[rigidID][2]*GRvec[1];
		p[n].v[1] = velocityGs[rigidID][1] + omegaGs[rigidID][2]*GRvec[0] - omegaGs[rigidID][0]*GRvec[2];
		p[n].v[2] = velocityGs[rigidID][2] + omegaGs[rigidID][0]*GRvec[1] - omegaGs[rigidID][1]*GRvec[0];
	}
	else{
		p[n].v[d] += ( dmy * 
			       ( p[n].f_hydro[d] + p[n].fv[d] 
				 + 0.5*(p[n].fr[d] + p[n].fr_previous[d] )
				 )
			       );
	}
	
	p[n].fr_previous[d] = p[n].fr[d];
	p[n].fr[d] = 0.0;
	p[n].fv_previous[d] = p[n].fv[d];
	p[n].fv[d] = 0.0;
	//p[n].f_hydro_previous[d] = p[n].f_hydro[d];
	p[n].f_hydro[d] = 0.0;
      }
      {
	p[n].omega_old[d] = p[n].omega[d];
	
	// T.K 13/01/04
	if(SW_PT == rigid){
		p[n].omega[d] = omegaGs[rigidID][d];
	}
	else{
		p[n].omega[d] += ( dmy_rot 
				   * ( p[n].torque_hydro[d] 
				       + 0.5 * (p[n].torquer[d] + p[n].torquer[d] )
				       + p[n].torquev[d] )
				   );
	}
	
	p[n].torquer_previous[d] = p[n].torquer[d];
	p[n].torquer[d] = 0.0;
	p[n].torquev_previous[d] = p[n].torquev[d];
	p[n].torquev[d] = 0.0;
	//p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
	p[n].torque_hydro[d] = 0.0;
      }
    }
  }
}
void MD_solver_velocity_AB2_hydro(Particle *p, const CTime &jikan){
    Force(p);
    calc_Rigid_VOGs(p, jikan, "AB2_hydro");	// T.K 13/01/04
    double dmy;
    double dmy_rot;
    // T.K 13/01/04
    int rigidID;
    double GRvec[DIM];
#pragma omp parallel for schedule(dynamic,1) private(dmy, dmy_rot)
    for(int n=0; n< Particle_Number; n++){
	//double dmy = jikan.hdt_md * IMASS[p[n].spec];
	//double dmy_rot = jikan.hdt_md * IMOI[p[n].spec];
	dmy = jikan.hdt_md * IMASS[p[n].spec];
	dmy_rot = jikan.hdt_md * IMOI[p[n].spec];
	
    // T.K 13/01/04
	if(SW_PT == rigid){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++) GRvec[d] = p[n].x[d] - xGs[rigidID][d];
	}
	
	for(int d=0; d<DIM; d++){  
	    {
		p[n].v_old[d] = p[n].v[d];
		// T.K 13/01/04
		if(SW_PT == rigid){
			p[n].v[0] = velocityGs[rigidID][0] + omegaGs[rigidID][1]*GRvec[2] - omegaGs[rigidID][2]*GRvec[1];
			p[n].v[1] = velocityGs[rigidID][1] + omegaGs[rigidID][2]*GRvec[0] - omegaGs[rigidID][0]*GRvec[2];
			p[n].v[2] = velocityGs[rigidID][2] + omegaGs[rigidID][0]*GRvec[1] - omegaGs[rigidID][1]*GRvec[0];
		}
		else{
			p[n].v[d] += 
			    dmy * (2.*p[n].f_hydro[d]
				   + 3.* p[n].fv[d] - p[n].fv_previous[d] // AB2
				   + p[n].fr[d] + p[n].fr_previous[d] // CN
				);
		}
		
		p[n].fr_previous[d] = p[n].fr[d];
		p[n].fr[d] = 0.0;
		p[n].fv_previous[d] = p[n].fv[d];
		p[n].fv[d] = 0.0;
		//p[n].f_hydro_previous[d] = p[n].f_hydro[d];
		p[n].f_hydro[d] = 0.0;
	    }
	    {
		p[n].omega_old[d] = p[n].omega[d];
		
		// T.K 13/01/04
		if(SW_PT == rigid){
			p[n].omega[d] = omegaGs[rigidID][d];
		}
		else{
			p[n].omega[d] += dmy_rot 
			    *( 2.* p[n].torque_hydro[d]
			       + 3.* p[n].torquev[d] -  p[n].torquev_previous[d] // AB2
			       + p[n].torquer[d] + p[n].torquer_previous[d] // CN
				);
		}
		
		p[n].torquer_previous[d] = p[n].torquer[d];
		p[n].torquer[d] = 0.0;
		p[n].torquev_previous[d] = p[n].torquev[d];
		p[n].torquev[d] = 0.0;
		//p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
		p[n].torque_hydro[d] = 0.0;
	    }
	}
    }
}

//OBL function
void MD_solver_position_Euler_OBL(Particle *p, const CTime &jikan){
#pragma omp parallel for schedule(dynamic, 1)
    for(int n=0; n<Particle_Number; n++){
	for(int d=0; d<DIM; d++){
	    p[n].x_previous[d] = p[n].x[d];
	    p[n].x[d] += jikan.dt_md * p[n].v[d];
	    
	}
	
	double signY = p[n].x[1];
	p[n].x[1] = fmod(p[n].x[1] + L_particle[1] , L_particle[1]);
	signY -= p[n].x[1];
	int sign = (int) signY;
	if (!(sign == 0)) {
	    sign = sign/abs(sign);
	}
	p[n].x[0] += -sign*degree_oblique*L_particle[1];
	p[n].v[0] -= sign*Shear_rate_eff*L_particle[1];
	p[n].v_old[0] -= sign*Shear_rate_eff*L_particle[1];
	
	p[n].x[0] = fmod(p[n].x[0] + L_particle[0] , L_particle[0]);
	
	p[n].x[2] = fmod(p[n].x[2] + L_particle[2] , L_particle[2]);

	for(int d = 0; d < DIM; d++){
	    assert(p[n].x[d] >= 0);
	    assert(p[n].x[d] < L[d]);
	}
    }
    set_xGs(p);		// T.K 13/01/20
}
void MD_solver_position_AB2_OBL(Particle *p, const CTime &jikan){
#pragma omp parallel for schedule(dynamic, 1)
    for(int n=0; n<Particle_Number; n++){
	for(int d=0; d<DIM; d++){
	    p[n].x_previous[d] = p[n].x[d];
	    p[n].x[d] += jikan.hdt_md * (3.*p[n].v[d] - p[n].v_old[d]);
	    
	}
	double signY = p[n].x[1];
	p[n].x[1] = fmod(p[n].x[1] + L_particle[1] , L_particle[1]);
	signY -= p[n].x[1];
	int sign = (int) signY;
	if (!(sign == 0)) {
	    sign  = sign/abs(sign);
	}
	p[n].x[0] -= (double)sign*degree_oblique*L_particle[1];
	p[n].v[0] -= (double)sign*Shear_rate_eff*L_particle[1];
	p[n].v_old[0] -= (double)sign*Shear_rate_eff*L_particle[1];
	
	p[n].x[0] = fmod(p[n].x[0] + L_particle[0] , L_particle[0]);
	
	p[n].x[2] = fmod(p[n].x[2] + L_particle[2] , L_particle[2]);
	
	for(int d = 0; d < DIM; d++){
	    assert(p[n].x[d] >= 0);
	    assert(p[n].x[d] < L[d]);
	}
    }
    set_xGs(p);		// T.K 13/01/20
}
void MD_solver_velocity_Euler_OBL(Particle *p, const CTime &jikan){
    Force_OBL(p);
    calc_Rigid_VOGs(p, jikan, "Euler");	// T.K 13/01/04
    double dmy;
    double dmy_rot;
    // T.K 13/01/04
    int rigidID;
    double GRvec[DIM];
#pragma omp parallel for schedule(dynamic,1) private(dmy, dmy_rot)
    for(int n=0; n< Particle_Number; n++){
	//static double dmy = jikan.dt_md * IMASS[p[n].spec];
	//static double dmy_rot = jikan.dt_md * IMOI[p[n].spec];
	dmy = jikan.dt_md * IMASS[p[n].spec];
	dmy_rot = jikan.dt_md * IMOI[p[n].spec];
	
    // T.K 13/01/04
	if(SW_PT == rigid){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++) GRvec[d] = p[n].x[d] - xGs[rigidID][d];
	}
	
	for(int d=0; d<DIM; d++){  
	    {
		p[n].v_old[d] = p[n].v[d];
		
		// T.K 13/01/04
		if(SW_PT == rigid){
			p[n].v[0] = velocityGs[rigidID][0] + omegaGs[rigidID][1]*GRvec[2] - omegaGs[rigidID][2]*GRvec[1];
			p[n].v[1] = velocityGs[rigidID][1] + omegaGs[rigidID][2]*GRvec[0] - omegaGs[rigidID][0]*GRvec[2];
			p[n].v[2] = velocityGs[rigidID][2] + omegaGs[rigidID][0]*GRvec[1] - omegaGs[rigidID][1]*GRvec[0];
		}
		else{
			p[n].v[d] += ( dmy * 
				       ( p[n].f_hydro[d] + p[n].fr[d] + p[n].fv[d]) 
			    );
		}
		
		p[n].momentum_depend_fr[d] = jikan.dt_md*p[n].fr[d];
		
		p[n].fr_previous[d] = p[n].fr[d];
		p[n].fr[d] = 0.0;
		p[n].fv_previous[d] = p[n].fv[d];
		p[n].fv[d] = 0.0;
		//p[n].f_hydro_previous[d] = p[n].f_hydro[d];
		p[n].f_hydro[d] = 0.0;
	    }
	    {
		p[n].omega_old[d] = p[n].omega[d];
		
		// T.K 13/01/04
		if(SW_PT == rigid){
			p[n].omega[d] = omegaGs[rigidID][d];
		}
		else{
			p[n].omega[d] += ( dmy_rot 
					   * ( p[n].torque_hydro[d] 
					       + p[n].torquer[d]
					       + p[n].torquev[d] )
			    );
		}
		
		p[n].torquer_previous[d] = p[n].torquer[d];
		p[n].torquer[d] = 0.0;
		p[n].torquev_previous[d] = p[n].torquev[d];
		p[n].torquev[d] = 0.0;
		//p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
		p[n].torque_hydro[d] = 0.0;
	    }
	}
    }
}
void MD_solver_velocity_AB2_hydro_OBL(Particle *p, const CTime &jikan){
    Force_OBL(p);
    calc_Rigid_VOGs(p, jikan, "AB2_hydro");	// T.K 13/01/04
    double dmy;
    double dmy_rot;
    // T.K 13/01/04
    int rigidID;
    double GRvec[DIM];
#pragma omp parallel for schedule(dynamic,1) private(dmy, dmy_rot)
    for(int n=0; n< Particle_Number; n++){
	//double dmy = jikan.hdt_md * IMASS[p[n].spec];
	//double dmy_rot = jikan.hdt_md * IMOI[p[n].spec];
	dmy = jikan.hdt_md * IMASS[p[n].spec];
	dmy_rot = jikan.hdt_md * IMOI[p[n].spec];
	
    // T.K 13/01/04
	if(SW_PT == rigid){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++) GRvec[d] = p[n].x[d] - xGs[rigidID][d];
	}
	
	for(int d=0; d<DIM; d++){  
	    {
		p[n].v_old[d] = p[n].v[d];
		
	// T.K 13/01/04
	if(SW_PT == rigid){
		p[n].v[0] = velocityGs[rigidID][0] + omegaGs[rigidID][1]*GRvec[2] - omegaGs[rigidID][2]*GRvec[1];
		p[n].v[1] = velocityGs[rigidID][1] + omegaGs[rigidID][2]*GRvec[0] - omegaGs[rigidID][0]*GRvec[2];
		p[n].v[2] = velocityGs[rigidID][2] + omegaGs[rigidID][0]*GRvec[1] - omegaGs[rigidID][1]*GRvec[0];
	}
	else{
			p[n].v[d] += 
			    dmy * (2.*p[n].f_hydro[d]
				   + 3.* p[n].fv[d] - p[n].fv_previous[d] // AB2
				   + p[n].fr[d] + p[n].fr_previous[d] // CN
				);
	}

		p[n].momentum_depend_fr[d] = jikan.hdt_md*(p[n].fr[d] + p[n].fr_previous[d]);
		
		p[n].fr_previous[d] = p[n].fr[d];
		p[n].fr[d] = 0.0;
		p[n].fv_previous[d] = p[n].fv[d];
		p[n].fv[d] = 0.0;
		//p[n].f_hydro_previous[d] = p[n].f_hydro[d];
		p[n].f_hydro[d] = 0.0;
	    }
	    {
		p[n].omega_old[d] = p[n].omega[d];
		
		// T.K 13/01/04
		if(SW_PT == rigid){
			p[n].omega[d] = omegaGs[rigidID][d];
		}
		else{
			p[n].omega[d] += dmy_rot 
			    *( 2.* p[n].torque_hydro[d]
			       + 3.* p[n].torquev[d] -  p[n].torquev_previous[d] // AB2
			       + p[n].torquer[d] + p[n].torquer_previous[d] // CN
				);
		}
		
		p[n].torquer_previous[d] = p[n].torquer[d];
		p[n].torquer[d] = 0.0;
		p[n].torquev_previous[d] = p[n].torquev[d];
		p[n].torquev[d] = 0.0;
		//p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
		p[n].torque_hydro[d] = 0.0;
	    }
	}
    }
}
