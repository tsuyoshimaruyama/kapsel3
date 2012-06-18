///
// $Id: particle_solver.cxx,v 1.1 2006/06/27 18:41:28 nakayama Exp $
//
#include "particle_solver.h"

inline void MD_solver_orientation_Euler(Particle &p, const double &dt){
  if(ROTATION){
    quaternion dqdt;
    qtn_init(p.q_old, p.q);
    qdot(dqdt, p.q, p.omega, SPACE_FRAME);
    qtn_add(p.q, dqdt, dt);
    qtn_normalize(p.q);
  }
}

inline void MD_solver_orientation_AB2(Particle &p, const double &dt){
  if(ROTATION){
    double wb[DIM];
    double wb_old[DIM];
    rigid_body_rotation(wb, p.omega, p.q, SPACE2BODY);
    rigid_body_rotation(wb_old, p.omega_old, p.q_old, SPACE2BODY);
    for(int d = 0; d < DIM; d++){
      wb[d] = 3.0*wb[d] - wb_old[d];
    }
    
    quaternion dqdt;
    qtn_init(p.q_old, p.q);
    qdot(dqdt, p.q, wb, BODY_FRAME);
    qtn_add(p.q, dqdt, dt);
    qtn_normalize(p.q);
  }
}

inline void MD_solver_orientation_SB2(Particle &p, const double &dt){
  if(ROTATION){
    double wb[DIM];
    for(int d = 0; d < DIM; d++){
      wb[d] = p.omega[d] + dt/2.0 * IMOI[p.spec]*torque_hydro[d];
    }
    
    quaternion dqdt;
    qtn_init(p.q_old, p.q);
    qdot(dqdt, p.q, wb, SPACE_FRAME);
    qtn_add(p.q, dqdt, dt);
    qtn_normalize(p.q);
  }
}

void MD_solver_position_Euler(Particle *p, const CTime &jikan)
{
#pragma omp parallel for schedule(dynamic, 1)
  for(int n = 0; n < Particle_Number; n++) {
    for(int d = 0; d < DIM; d++) {
      p[n].x[d] += jikan.dt_md * p[n].v[d];
      p[n].x[d] = fmod(p[n].x[d] + L_particle[d] , L_particle[d]);

      assert(p[n].x[d] >= 0);
      assert(p[n].x[d] < L[d]);
    }
    
    MD_solver_orientation_Euler(p[n], jikan.dt_md);
  }

}

void MD_solver_position_AB2(Particle *p, const CTime &jikan)
{

#pragma omp parallel for schedule(dynamic, 1)
  for(int n = 0; n < Particle_Number; n++) {
    for(int d = 0; d < DIM; d++) {
      p[n].x[d] += jikan.hdt_md * (3.*p[n].v[d] - p[n].v_old[d]);
      p[n].x[d] = fmod(p[n].x[d] + L_particle[d] , L_particle[d]);

      assert(p[n].x[d] >= 0);
      assert(p[n].x[d] < L[d]);
    }

    MD_solver_orientation_AB2(p[n], jikan.hdt_md);
  }
}

inline void self_propulsion(Particle &p, double *force, double *torque){
  for(int d = 0; d < DIM; d++){
    force[d] = 0.0;
    torque[d] = 0.0;
  }
  if(SW_JANUS_MOTOR){
    int &spec = p.spec;
    int self_flag = janus_propulsion[spec];
    if(self_flag == swimmer || self_flag == tumbler){
      rigid_body_rotation(force, janus_force[spec], p.q, BODY2SPACE);
    }
    if(ROTATION && (self_flag == rotator || self_flag == tumbler)){
      rigid_body_rotation(torque, janus_torque[spec], p.q, BODY2SPACE);
    }
  }
}

void MD_solver_velocity_Euler(Particle *p, const CTime &jikan)
{
  Force(p);

  double dmy;
  double dmy_rot;
  double self_force[DIM];
  double self_torque[DIM];
#pragma omp parallel for schedule(dynamic,1) private(dmy, dmy_rot, self_force, self_torque)
  for(int n = 0; n < Particle_Number; n++) {
    //double dmy = jikan.dt_md * IMASS[p[n].spec];
    //double dmy_rot = jikan.dt_md * IMOI[p[n].spec];
    dmy = jikan.dt_md * IMASS[p[n].spec];
    dmy_rot = jikan.dt_md * IMOI[p[n].spec];
    self_propulsion(p[n], self_force, self_torque);
    for(int d = 0; d < DIM; d++) {
      {
	p[n].v_old[d] = p[n].v[d];
	p[n].v[d] += ( dmy *
		       (p[n].f_hydro[d] + p[n].f_slip[d] + 
			p[n].fr[d] + self_force[d]) );

	p[n].fr_previous[d] = p[n].fr[d];
	p[n].fr[d] = 0.0;

	p[n].f_hydro_previous[d] = p[n].f_hydro[d];
	p[n].f_hydro[d] = 0.0;

	p[n].f_slip_previous[d] = p[n].f_slip[d];
	p[n].f_slip[d] = 0.0;

      }
      {
	p[n].omega_old[d] = p[n].omega[d];
	p[n].omega[d] += ( dmy_rot * (p[n].torque_hydro[d] + p[n].torque_slip[d] + self_torque[d]));

	p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
	p[n].torque_hydro[d] = 0.0;

	p[n].torque_slip_previous[d] = p[n].torque_slip[d];
	p[n].torque_slip[d] = 0.0;
      }
    }
  }
}

void MD_solver_velocity_iter(Particle *p, const CTime &jikan, 
				   const ITER &iter_flag){
  if(iter_flag != start_iter && iter_flag != new_iter && 
     iter_flag != end_iter && iter_flag != reset_iter){
    fprintf(stderr, "Error: wrong Euler iter_flag\n");
    exit_job(EXIT_FAILURE);
  }
  if(iter_flag == start_iter){//only compute forces once
    Force(p);
  }

  if(iter_flag == start_iter || iter_flag == new_iter){
    double dmy;
    double dmy_rot;
    double self_force[DIM];
    double self_torque[DIM];

    if(jikan.ts == 0){//Euler update
      #pragma omp parallel for schedule(dynamic, 1) private(dmy, dmy_rot, self_force, self_torque)
      for(int n = 0; n < Particle_Number; n++){
	dmy = jikan.dt_md * IMASS[p[n].spec];
	dmy_rot = jikan.dt_md * IMOI[p[n].spec];
	self_propulsion(p[n], self_force, self_torque);
	
	for(int d = 0; d < DIM; d++){
	  p[n].v_old[d] = p[n].v[d];
	  p[n].v[d] += (dmy * (p[n].f_hydro[d] + p[n].f_slip[d] + p[n].fr[d] + self_force[d]));
	  
	  p[n].omega_old[d] = p[n].omega[d];
	  p[n].omega[d] += (dmy_rot * (p[n].torque_hydro[d] + p[n].torque_slip[d] + self_torque[d]));
	}
      }//Euler
    }else{//AB2_hydro update
      #pragma omp parallel for schedule(dynamic, 1) private(dmy, dmy_rot, self_force, self_torque)
      for(int n = 0; n < Particle_Number; n++){
	dmy = jikan.hdt_md * IMASS[p[n].spec];
	dmy_rot = jikan.hdt_md * IMOI[p[n].spec];
	self_propulsion(p[n], self_force, self_torque);

	for(int d = 0; d < DIM; d++){
	  p[n].v_old[d] = p[n].v[d];
	  p[n].v[d] += 
	    dmy * (2.0 * p[n].f_hydro[d] + 2.0 * p[n].f_slip[d] + 2.0 * self_force[d]
		   + p[n].fr[d] + p[n].fr_previous[d]);
	  
	  p[n].omega_old[d] = p[n].omega[d];
	  p[n].omega[d] +=
	    dmy_rot * (2.0 * p[n].torque_hydro[d] + 2.0 * p[n].torque_slip[d] + 2.0 * self_torque[d]);
	}
      }//AB2
    }

  }//start_iter || new_iter
  else if(iter_flag == reset_iter){
    #pragma omp parallel for schedule(dynamic, 1)
    for(int n = 0; n < Particle_Number; n++){
      for(int d = 0; d < DIM; d++){
	p[n].v[d] = p[n].v_old[d];
	p[n].omega[d] = p[n].omega_old[d];
      }
    }
    
  }// reset_iter
  else if(iter_flag == end_iter){
    #pragma omp parallel for schedule(dynamic, 1)
    for(int n = 0; n < Particle_Number; n++){
      for(int d = 0; d < DIM; d++){
	{
	  p[n].fr_previous[d] = p[n].fr[d];
	  p[n].fr[d] = 0.0;
	  
	  p[n].f_hydro_previous[d] = p[n].f_hydro[d];
	  p[n].f_hydro[d] = 0.0;
	  
	  p[n].f_slip_previous[d] = p[n].f_slip[d];
	  p[n].f_slip[d] = 0.0;
	}
	{
	  p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
	  p[n].torque_hydro[d] = 0.0;
	  
	  p[n].torque_slip_previous[d] = p[n].torque_slip[d];
	  p[n].torque_slip[d] = 0.0;
	}
      }
    }

  }//end_iter

}

void MD_solver_velocity_Euler_hydro(Particle *p, const CTime &jikan){
  Force(p);

  double dmy;
  double dmy_rot;
  double self_force[DIM];
  double self_torque[DIM];
#pragma omp parallel for schedule(dynamic,1) private(dmy, dmy_rot, self_force, self_torque)
  for(int n = 0; n < Particle_Number; n++) {
    //static double dmy = jikan.dt_md * IMASS[p[n].spec];
    //static double dmy_rot = jikan.dt_md * IMOI[p[n].spec];
    dmy = jikan.dt_md * IMASS[p[n].spec];
    dmy_rot = jikan.dt_md * IMOI[p[n].spec];
    self_propulsion(p[n], self_force, self_torque);

    for(int d = 0; d < DIM; d++) {
      {
	p[n].v_old[d] = p[n].v[d];
	p[n].v[d] += ( dmy *
		       ( p[n].f_hydro[d] +
			 0.5*(p[n].fr[d] + p[n].fr_previous[d] ) +
			 self_force[d]
			 )
		       );

	p[n].fr_previous[d] = p[n].fr[d];
	p[n].fr[d] = 0.0;
	p[n].f_hydro_previous[d] = p[n].f_hydro[d];
	p[n].f_hydro[d] = 0.0;
      }
      {
	p[n].omega_old[d] = p[n].omega[d];
	p[n].omega[d] += ( dmy_rot * (p[n].torque_hydro[d] + self_torque[d]));

	p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
	p[n].torque_hydro[d] = 0.0;
      }
    }
  }
}

void MD_solver_velocity_AB2_hydro(Particle *p, const CTime &jikan){
  Force(p);

  double dmy;
  double dmy_rot;
  double self_force[DIM];
  double self_torque[DIM];

#pragma omp parallel for schedule(dynamic,1) private(dmy, dmy_rot, self_force, self_torque)
  for(int n = 0; n < Particle_Number; n++) {
	//double dmy = jikan.hdt_md * IMASS[p[n].spec];
	//double dmy_rot = jikan.hdt_md * IMOI[p[n].spec];
    dmy = jikan.hdt_md * IMASS[p[n].spec];
    dmy_rot = jikan.hdt_md * IMOI[p[n].spec];
    self_propulsion(p[n], self_force, self_torque);
    for(int d = 0; d < DIM; d++) {
      {
	p[n].v_old[d] = p[n].v[d];
	p[n].v[d] +=
	  dmy * (2.*p[n].f_hydro[d]
		 + 2.*p[n].f_slip[d]
		 + p[n].fr[d] + p[n].fr_previous[d] // CN
		 + 2.0 * self_force[d]
		 );

	p[n].fr_previous[d] = p[n].fr[d];
	p[n].fr[d] = 0.0;

	p[n].f_hydro_previous[d] = p[n].f_hydro[d];
	p[n].f_hydro[d] = 0.0;

	p[n].f_slip_previous[d] = p[n].f_slip[d];
	p[n].f_slip[d] = 0.0;
      }
      {
	p[n].omega_old[d] = p[n].omega[d];
	p[n].omega[d] += dmy_rot * ( 2.* p[n].torque_hydro[d] + 
				     2.* p[n].torque_slip[d] +
				     2.* self_torque[d]);

	p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
	p[n].torque_hydro[d] = 0.0;

	p[n].torque_slip_previous[d] = p[n].torque_slip[d];
	p[n].torque_slip[d] = 0.0;
      }
    }

  }
}

//OBL function
void MD_solver_position_Euler_OBL(Particle *p, const CTime &jikan){
#pragma omp parallel for schedule(dynamic, 1)
  for(int n = 0; n < Particle_Number; n++) {
    for(int d = 0; d < DIM; d++) {
      p[n].x_previous[d] = p[n].x[d];
      p[n].x[d] += jikan.dt_md * p[n].v[d];

    }

    double signY = p[n].x[1];
    p[n].x[1] = fmod(p[n].x[1] + L_particle[1] , L_particle[1]);
    signY -= p[n].x[1];
    int sign = (int) signY;
    if (!(sign == 0)) {
      sign = sign / abs(sign);
    }
    p[n].x[0] += -sign * degree_oblique * L_particle[1];
    p[n].v[0] -= sign * Shear_rate_eff * L_particle[1];
    p[n].v_old[0] -= sign * Shear_rate_eff * L_particle[1];

    p[n].x[0] = fmod(p[n].x[0] + L_particle[0] , L_particle[0]);

    p[n].x[2] = fmod(p[n].x[2] + L_particle[2] , L_particle[2]);

    for(int d = 0; d < DIM; d++) {
      assert(p[n].x[d] >= 0);
      assert(p[n].x[d] < L[d]);
    }

    MD_solver_orientation_Euler(p[n], jikan.dt_md);
  }
}

void MD_solver_position_AB2_OBL(Particle *p, const CTime &jikan){
#pragma omp parallel for schedule(dynamic, 1)
  for(int n = 0; n < Particle_Number; n++) {
    for(int d = 0; d < DIM; d++) {
      p[n].x_previous[d] = p[n].x[d];
      p[n].x[d] += jikan.hdt_md * (3.*p[n].v[d] - p[n].v_old[d]);

    }
    double signY = p[n].x[1];
    p[n].x[1] = fmod(p[n].x[1] + L_particle[1] , L_particle[1]);
    signY -= p[n].x[1];
    int sign = (int) signY;
    if (!(sign == 0)) {
      sign  = sign / abs(sign);
    }
    p[n].x[0] -= (double)sign * degree_oblique * L_particle[1];
    p[n].v[0] -= (double)sign * Shear_rate_eff * L_particle[1];
    p[n].v_old[0] -= (double)sign * Shear_rate_eff * L_particle[1];

    p[n].x[0] = fmod(p[n].x[0] + L_particle[0] , L_particle[0]);

    p[n].x[2] = fmod(p[n].x[2] + L_particle[2] , L_particle[2]);

    for(int d = 0; d < DIM; d++) {
      assert(p[n].x[d] >= 0);
      assert(p[n].x[d] < L[d]);
    }
    
    MD_solver_orientation_AB2(p[n], jikan.hdt_md);
  }
}
void MD_solver_velocity_Euler_OBL(Particle *p, const CTime &jikan){
  Force_OBL(p);

  double dmy;
  double dmy_rot;
  double self_force[DIM];
  double self_torque[DIM];
#pragma omp parallel for schedule(dynamic,1) private(dmy, dmy_rot, self_force, self_torque)
  for(int n = 0; n < Particle_Number; n++) {
	//static double dmy = jikan.dt_md * IMASS[p[n].spec];
	//static double dmy_rot = jikan.dt_md * IMOI[p[n].spec];
    dmy = jikan.dt_md * IMASS[p[n].spec];
    dmy_rot = jikan.dt_md * IMOI[p[n].spec];
    self_propulsion(p[n], self_force, self_torque);
    for(int d = 0; d < DIM; d++) {
      {
	p[n].v_old[d] = p[n].v[d];
	p[n].v[d] += ( dmy *
		       ( p[n].f_hydro[d] + p[n].fr[d] + self_force[d])
		       );

	p[n].momentum_depend_fr[d] = jikan.dt_md * p[n].fr[d];

	p[n].fr_previous[d] = p[n].fr[d];
	p[n].fr[d] = 0.0;
	p[n].f_hydro_previous[d] = p[n].f_hydro[d];
	p[n].f_hydro[d] = 0.0;
      }
      {
	p[n].omega_old[d] = p[n].omega[d];
	p[n].omega[d] += ( dmy_rot * (p[n].torque_hydro[d] + self_torque[d]));

	p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
	p[n].torque_hydro[d] = 0.0;
      }
    }
  }
}
void MD_solver_velocity_AB2_hydro_OBL(Particle *p, const CTime &jikan){
  Force_OBL(p);

  double dmy;
  double dmy_rot;
  double self_force[DIM];
  double self_torque[DIM];
#pragma omp parallel for schedule(dynamic,1) private(dmy, dmy_rot, self_force, self_torque)
  for(int n = 0; n < Particle_Number; n++) {
	//double dmy = jikan.hdt_md * IMASS[p[n].spec];
	//double dmy_rot = jikan.hdt_md * IMOI[p[n].spec];
    dmy = jikan.hdt_md * IMASS[p[n].spec];
    dmy_rot = jikan.hdt_md * IMOI[p[n].spec];
    self_propulsion(p[n], self_force, self_torque);
    for(int d = 0; d < DIM; d++) {
      {
	p[n].v_old[d] = p[n].v[d];
	p[n].v[d] +=
	  dmy * (2.*p[n].f_hydro[d]
		 + p[n].fr[d] + p[n].fr_previous[d] // CN
		 + 2.0 * self_force[d]
		 );

	p[n].momentum_depend_fr[d] = jikan.hdt_md * (p[n].fr[d] + p[n].fr_previous[d]);

	p[n].fr_previous[d] = p[n].fr[d];
	p[n].fr[d] = 0.0;
	p[n].f_hydro_previous[d] = p[n].f_hydro[d];
	p[n].f_hydro[d] = 0.0;
      }
      {
	p[n].omega_old[d] = p[n].omega[d];
	p[n].omega[d] += dmy_rot * ( 2.* p[n].torque_hydro[d] + 2.* self_torque[d]);

	p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
	p[n].torque_hydro[d] = 0.0;
      }
    }
  }
}
