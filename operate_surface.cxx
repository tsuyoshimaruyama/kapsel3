#include "operate_surface.h"

void Make_particle_momentum_factor(double const* const* u, Particle *p){
  //////////////////////////////
  const double dx = DX;
  const double dx3 = DX3;
  const int np_domain = NP_domain;
  int const* const* sekibun_cell = Sekibun_cell;
  const int* Nlattice = Ns;
  const double radius = RADIUS;
  //////////////////////////////
  int sw_in_cell, pspec;
  int x_int[DIM], r_mesh[DIM];
  double dmy_r, dmy_phi, dmy_region;
  double xp[DIM], r[DIM], x[DIM], residue[DIM], u_fluid[DIM];

  double dmy_xi, dmy_theta, dmy_tau;
  double dmy_sin, dmy_sin2;
  double n_r[DIM], n_theta[DIM], n_tau[DIM];
  double M0, SM0;
  double M1[DIM], SM1[DIM], dv_s[DIM], domega_s[DIM];
  double M2[DIM][DIM], SM2[DIM][DIM], ST[DIM][DIM], SU[DIM][DIM], SV[DIM][DIM];
  double slip_mode, slip_vel, slip_magnitude;
  double r_x_theta[DIM], r_x_us[DIM], us[DIM];
#pragma omp parallel for schedule(dynamic, 1) \
  private(sw_in_cell, pspec, x_int, r_mesh, dmy_r, dmy_phi, dmy_region, xp, \
	  r, x, residue, u_fluid, dmy_xi, dmy_theta, dmy_tau, dmy_sin, dmy_sin2,	\
	  n_r, n_theta, n_tau, M0, SM0, M1, SM1, dv_s, domega_s, M2, SM2, ST, SU, SV, \
	  slip_mode, slip_vel, slip_magnitude, r_x_theta, r_x_us, us)
  for(int n = 0; n < Particle_Number; n++){
    pspec = p[n].spec;

    for(int d = 0; d < DIM; d++){
      xp[d] = p[n].x[d];
      
      M1[d] = SM1[d] = dv_s[d] = domega_s[d] = 0.0;
      for(int l = 0; l < DIM; l++){
	M2[d][l] = SM2[d][l] = ST[d][l] = SU[d][l] = SV[d][l] = 0.0;
      }
    }
    M0 = SM0 = 0.0;
    
    sw_in_cell = Particle_cell(xp, dx, x_int, residue);
    sw_in_cell = 1;
    for(int mesh = 0; mesh < np_domain; mesh++){
      Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
      int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
      for(int d = 0; d < DIM; d++){
	x[d] = r_mesh[d] * dx;
	u_fluid[d] = u[d][im];
      }
      dmy_r = Distance(x, xp);
      dmy_phi = Phi(dmy_r, radius);
      dmy_xi = ABS(dmy_r - radius);
      
      //particle properties
      M0 += dmy_phi;
      for(int d = 0; d < DIM; d++){
	M1[d] += dmy_phi * r[d];
	M2[d][d] += dmy_phi * dmy_r * dmy_r;
	for(int l = 0; l < DIM; l++){
	  M2[d][l] -= dmy_phi * r[d] * r[l];
	}
      }

      //surface properties
      if(janus_propulsion[pspec] == slip &&  less_than_mp(dmy_xi, HXI)){

	slip_vel = janus_slip_vel[pspec] * janus_slip_scale;
	slip_mode = janus_slip_mode[pspec];
	if(janus_slip_region == surface_slip){
	  dmy_region = (1.0 - dmy_phi) * DPhi_compact_sin_norm(dmy_r, radius);
	}else{
	  dmy_region = (1.0 - dmy_phi);
	}
	
	SM0 += dmy_region;
	Squirmer_coord(r, n_r, n_theta, n_tau, dmy_r, dmy_theta, dmy_tau, p[n]);
	dmy_sin = sin(dmy_theta);
	dmy_sin2 = sin(2.0 * dmy_theta);
	slip_magnitude = slip_vel * (dmy_sin + slip_mode * dmy_sin2);
	if(janus_slip_boundary == full_boundary){
	  for(int d = 0; d < DIM; d++){
	    us[d] = n_theta[d] * slip_magnitude - u_fluid[d];
	  }
	}else{ // partial tangential boundary
	  double dmy_udot = u_fluid[0] * n_theta[0] + u_fluid[1] * n_theta[1] + u_fluid[2] * n_theta[2];
	  for(int d = 0; d < DIM; d++){
	    us[d] = n_theta[d] * (slip_magnitude - dmy_udot);
	  }
	}
	
	v_cross(r_x_theta, r, n_theta);
	v_cross(r_x_us, r, us);
	for(int d = 0; d < DIM; d++){
	  dv_s[d] += dmy_region * us[d];
	  domega_s[d] += dmy_region * r_x_us[d];
	  
	  SM1[d] += dmy_region * r[d];
	  SM2[d][d] += dmy_region * dmy_r * dmy_r;
	  for(int l = 0; l < DIM; l++){
	    SM2[d][l] -= dmy_region * r[d] * r[l];
	    ST[d][l] += dmy_region * n_theta[d] * n_theta[l];
	    SU[d][l] += dmy_region * n_theta[d] * r_x_theta[l];
	    SV[d][l] += dmy_region * r_x_theta[d] * r_x_theta[l];
	  }
	}
      }//slip surface

    }//mesh

    p[n].mass = M0;
    p[n].surface_mass = SM0;
    for(int d = 0; d < DIM; d++){
      
      p[n].mass_center[d] = M1[d];
      p[n].surface_mass_center[d] = SM1[d];
      p[n].surface_dv[d] = dv_s[d];
      p[n].surface_domega[d] = domega_s[d];
      
      for(int l = 0; l < DIM; l++){
	p[n].inertia[d][l] = M2[d][l];
	p[n].surface_inertia[d][l] = SM2[d][l];
	p[n].surfaceT[d][l] = ST[d][l];
	p[n].surfaceU[d][l] = SU[d][l];
	p[n].surfaceV[d][l] = SV[d][l];
      }
    }

  }//Particle_Number
}

inline void slip_droplet(double *vp, double *wp, double *delta_v, double *delta_omega, Particle p){
  double LL[DIM][DIM];
  double dP[DIM], dL[DIM], L0[DIM], mc[DIM];
  double imass, mc2;

  for(int d = 0; d < DIM; d++){
    vp[d] = p.v_slip[d];
    wp[d] = p.omega_slip[d];
    mc[d] = p.mass_center[d];
    LL[d][0] = LL[d][1] = LL[d][2] = 0.0;
  }
  imass = 1.0 / p.mass;

  //momentum change due to slip at boundary
  if(janus_slip_boundary == full_boundary){
    double w_x_r[DIM];
    double r_x_v[DIM];
    v_cross(w_x_r, wp, p.surface_mass_center);
    v_cross(r_x_v, p.surface_mass_center, vp);
    for(int d = 0; d < DIM; d++){
      dP[d] = -(p.surface_mass * vp[d] + w_x_r[d] + p.surface_dv[d]);
      dL[d] = -(r_x_v[d] + 
		p.surface_inertia[d][0] * wp[0] + p.surface_inertia[d][1] * wp[1] + p.surface_inertia[d][2] * wp[2] +
		p.surface_domega[d]);
    }
  }else{//partial_boundary
    for(int d = 0; d < DIM; d++){
      dP[d] = -(p.surfaceT[d][0] * vp[0] + p.surfaceT[d][1] * vp[1] + p.surfaceT[d][2] * vp[2] +
		p.surfaceU[d][0] * wp[0] + p.surfaceU[d][1] * wp[1] + p.surfaceU[d][2] * wp[2] +
		p.surface_dv[d]);
      
      dL[d] = -(p.surfaceU[0][d] * vp[0] + p.surfaceU[1][d] * vp[1] + p.surfaceU[2][d] * vp[2] +
		p.surfaceV[d][0] * wp[0] + p.surfaceV[d][1] * wp[1] + p.surfaceV[d][2] * wp[2] +
		p.surface_domega[d]);
    }
  }

  // Compute angular velocity of droplet
  v_cross(L0, mc, dP, imass);
  mc2 = (mc[0] * mc[0] + mc[1] * mc[1] + mc[2] * mc[2]) * imass;
  for(int d = 0; d < DIM; d++){
    L0[d] = dL[d] - L0[d];
    LL[d][d] = -mc2;

    for(int l = 0; l < DIM; l++){
      LL[d][l] += (p.inertia[d][l] + imass * mc[d] * mc[l]);
    }
  }
  M_inv(LL);
  M_v_prod(delta_omega, LL, L0);

  // Compute linear velocity of droplet
  v_cross(L0, delta_omega, mc);
  for(int d = 0; d < DIM; d++){
    delta_v[d] = imass * (dP[d] - L0[d]);
  }
}

void Make_force_u_slip_particle(double **up, double const* const* u, Particle *p, const CTime &jikan){
  //////////////////////// System parameters
  const double dx = DX;
  const double dx3 = DX3;
  const int np_domain = NP_domain;
  int const* const* sekibun_cell = Sekibun_cell;
  const int* Nlattice = Ns;
  const double radius = RADIUS;

  ////////////////////////  Function variables
  int sw_in_cell, pspec;
  int x_int[DIM], r_mesh[DIM];
  double dmy_r, dmy_phi, dmy_region;
  double xp[DIM], vp[DIM], omega_p[DIM], v_rot[DIM];
  double delta_v[DIM], delta_omega[DIM], delta_v_rot[DIM];
  double r[DIM], x[DIM], residue[DIM], u_fluid[DIM];
  
  double dmy_xi, dmy_theta, dmy_tau;
  double dmy_sin, dmy_sin2, dmy_vslip, slip_mode, slip_vel, dmy_us;
  double n_r[DIM], n_theta[DIM], n_tau[DIM];
  double dmy_fv[DIM], force_s[DIM], torque_s[DIM], force_p[DIM], torque_p[DIM];

#pragma omp parallel for schedule(dynamic, 1) \
  private(sw_in_cell, pspec, x_int, r_mesh, dmy_r, dmy_phi, dmy_region, xp, vp, omega_p, v_rot, \
	  delta_v, delta_omega, delta_v_rot, r, x, residue, u_fluid, \
	  dmy_xi, dmy_theta, dmy_tau, dmy_sin, dmy_sin2, dmy_vslip, slip_mode, slip_vel, dmy_us,  \
	  n_r, n_theta, n_tau, dmy_fv, force_s, torque_s, force_p, torque_p)
  for(int n = 0; n < Particle_Number; n++){
    pspec = p[n].spec;

    for(int d = 0; d < DIM; d++){
      p[n].f_slip[d] = 0.0;
      p[n].torque_slip[d] = 0.0;
    }
    
    if(janus_propulsion[pspec] == slip){
      
      slip_vel = janus_slip_vel[pspec] * janus_slip_scale;
      slip_mode = janus_slip_mode[pspec];
      slip_droplet(vp, omega_p, delta_v, delta_omega, p[n]);
      for(int d = 0; d < DIM; d++){
	xp[d] = p[n].x[d];
	force_s[d] = torque_s[d] = 0.0;
	force_p[d] = torque_p[d] = 0.0;
      }
      sw_in_cell = Particle_cell(xp, dx, x_int, residue);
      sw_in_cell = 1;

      for(int mesh = 0; mesh < np_domain; mesh++){
	Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
	int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
	for(int d = 0; d < DIM; d++){
	  u_fluid[d] = u[d][im];
	  x[d] = r_mesh[d] * dx;
	}
	dmy_r = Distance(x, xp);
	dmy_phi = Phi(dmy_r, radius);	 
	dmy_xi = ABS(dmy_r - radius);

	{//fluid particle domain (droplet with counter flow)
          Angular2v(delta_omega, r, delta_v_rot);
          for(int d = 0; d < DIM; d++){
            dmy_fv[d] = dmy_phi * (delta_v[d] + delta_v_rot[d]);
            force_p[d] += dmy_fv[d];
            
#pragma omp atomic
            up[d][im] += dmy_fv[d];
          }
          {
            torque_p[0] += (r[1] * dmy_fv[2] - r[2] * dmy_fv[1]);
            torque_p[1] += (r[2] * dmy_fv[0] - r[0] * dmy_fv[2]);
            torque_p[2] += (r[0] * dmy_fv[1] - r[1] * dmy_fv[0]);
          }
        }
	
	if(janus_slip_region == surface_slip){
	  dmy_region = (1.0 - dmy_phi) * DPhi_compact_sin_norm(dmy_r, radius);
	}else{
	  dmy_region = (1.0 - dmy_phi);
	}
	if(less_than_mp(dmy_xi, HXI)){//interface domain
	  Angular2v(omega_p, r, v_rot);
	  Squirmer_coord(r, n_r, n_theta, n_tau, dmy_r, dmy_theta, dmy_tau, p[n]);
	  dmy_sin = sin(dmy_theta);
	  dmy_sin2 = sin(2.0 * dmy_theta);
	  
	  for(int d = 0; d < DIM; d++){
	    dmy_fv[d] = (vp[d] + v_rot[d]);
	  }
	  dmy_vslip = slip_vel * (dmy_sin + slip_mode * dmy_sin2);
	  if(janus_slip_boundary == full_boundary){
	    for(int d = 0; d < DIM; d++){
	      dmy_fv[d] = dmy_region * (dmy_fv[d] + n_theta[d]*dmy_vslip - u_fluid[d]);
	    }
	  }else { //partial tangential boundary
	    dmy_us = dmy_vslip + 
	      (dmy_fv[0] - u_fluid[0])*n_theta[0] + (dmy_fv[1] - u_fluid[1])*n_theta[1] + (dmy_fv[2] - u_fluid[2])*n_theta[2]; 
	    dmy_us *= dmy_region;
	    for(int d = 0; d < DIM; d++){
	      dmy_fv[d] = n_theta[d] * dmy_us;
	    }
	  }
	  for(int d = 0; d < DIM; d++){
	    force_s[d] += dmy_fv[d];
#pragma omp atomic
	    up[d][im] += dmy_fv[d];
	  }
	  {
	    torque_s[0] += (r[1] * dmy_fv[2] - r[2] * dmy_fv[1]);
	    torque_s[1] += (r[2] * dmy_fv[0] - r[0] * dmy_fv[2]);
	    torque_s[2] += (r[0] * dmy_fv[1] - r[1] * dmy_fv[0]);
	  }
	}//interface_domain
      }//mesh
      for(int d = 0; d < DIM; d++){
	force_p[d] = -force_p[d];
	torque_p[d] = -torque_p[d];
      }

      if((v_rms(force_p, force_s) > LARGE_TOL_MP || 
          v_rms(torque_p, torque_s) > LARGE_TOL_MP)){
        fprintf(stderr, "###############################");
	fprintf(stderr, "# Momentum Conservation Warning : %10.8E %10.8E\n",
		v_rms(force_p, force_s), v_rms(torque_p, torque_s));
        fprintf(stderr, "# Force  : %10.8E %10.8E %10.8E %10.8E %10.8E %10.8E\n", force_p[0], force_p[1], force_p[2], force_s[0], force_s[1], force_s[2]);
        fprintf(stderr, "# Torque : %10.8E %10.8E %10.8E %10.8E %10.8E %10.8E\n", torque_p[0], torque_p[1], torque_p[2], torque_s[0], torque_s[1], torque_s[2]);
        fprintf(stderr, "###############################");
	
      }
    }// slip_particle ?
  }// Particle_Number
}

