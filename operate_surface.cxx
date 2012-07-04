#include "operate_surface.h"
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
  double dmy_r, dmy_phi;
  double xp[DIM], vp[DIM], omega_p[DIM], v_rot[DIM];
  double r[DIM], x[DIM], residue[DIM];
  
  double dmy_xi, dmy_theta, dmy_tau;
  double dmy_vdot, dmy_udot;
  double dmy_mass, dmy_sin, dmy_sin2;
  double n_r[DIM], n_theta[DIM], n_tau[DIM];
  double Vv[DIM], Uv[DIM], Sv[DIM], dmy_fv[DIM], polar_axis[DIM];
  double force[DIM], torque[DIM];
  double slip_mode, slip_vel;

  for(int n = 0; n < Particle_Number; n++){
    pspec = p[n].spec;
    if(janus_propulsion[pspec] == slip){

      slip_vel = janus_slip_vel[pspec] * janus_slip_scale;
      slip_mode = janus_slip_mode[pspec];
      Janus_direction(polar_axis, p[n]);
      for(int d = 0; d < DIM; d++){
	xp[d] = p[n].x[d];
	vp[d] = p[n].v[d];
	omega_p[d] = p[n].omega[d];

	Vv[d] = Sv[d] = Uv[d] = 0.0;
	force[d] = torque[d] = 0.0;
      }
      dmy_mass = 0.0;
      sw_in_cell = Particle_cell(xp, dx, x_int, residue);
      sw_in_cell = 1;

      { // compute hydrodynamic force 
	for(int mesh = 0; mesh < np_domain; mesh++){
	  Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
	  int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
	  for(int d = 0; d < DIM; d++){
	    x[d] = r_mesh[d] * dx;
	  }
	  dmy_r = Distance(x, xp);
	  dmy_phi = 1.0 - Phi(dmy_r, radius);
	  dmy_mass += (1.0 - dmy_phi);
	  dmy_xi = ABS(dmy_r - radius);

	  if(dmy_xi <= HXI && dmy_phi > 0.0){//interface domain
	    Angular2v(omega_p, r, v_rot);
	    Spherical_coord(r, n_r, n_theta, n_tau, dmy_r, dmy_theta, dmy_tau, p[n]);
	    dmy_sin = sin(dmy_theta);
	    dmy_sin2 = sin(2.0 * dmy_theta);
	    
	    for(int d = 0; d < DIM; d++){
	      dmy_fv[d] = (vp[d] + v_rot[d]);
	    }
	    dmy_udot = u[0][im]*n_theta[0] + u[1][im]*n_theta[1] + u[2][im]*n_theta[2];
	    dmy_vdot = dmy_fv[0]*n_theta[0] + dmy_fv[1]*n_theta[1] + dmy_fv[2]*n_theta[2];
	    if(janus_slip_debug == full_tangent){
	      for(int d = 0; d < DIM; d++){
		Uv[d] = n_theta[d] * dmy_udot;
		Vv[d] = n_theta[d] * dmy_vdot;
		Sv[d] = n_theta[d] * (slip_vel * (dmy_sin + slip_mode * dmy_sin2));
	      }
	    }else if(janus_slip_debug == particle_tangent){
	      for(int d = 0; d < DIM; d++){
		Uv[d] = u[d][im];
		Vv[d] = n_theta[d] * dmy_vdot;
		Sv[d] = n_theta[d] * (slip_vel * (dmy_sin + slip_mode * dmy_sin2));
	      }
	    }else if(janus_slip_debug == no_tangent){
	      for(int d = 0; d < DIM; d++){
		Uv[d] = u[d][im];
		Vv[d] = dmy_fv[d];
		Sv[d] = n_theta[d] * (slip_vel * (dmy_sin + slip_mode * dmy_sin2));
	      }
	    }else{
	      fprintf(stderr, "OPERATE_SURFACE tangent error\n");
	      exit_job(EXIT_FAILURE);
	    }
	    /*	    fprintf(stderr, "# %8.6E %8.6E %8.6E\n",
		    dmy_fv[0]*n_r[0] + dmy_fv[1]*n_r[1] + dmy_fv[2]*n_r[2],
		    dmy_fv[0]*n_theta[0] + dmy_fv[1]*n_theta[1] + dmy_fv[2]*n_theta[2],
		    dmy_fv[0]*n_tau[0] + dmy_fv[1]*n_tau[1] + dmy_fv[2]*n_tau[2]);*/

	    for(int d = 0; d < DIM; d++){
	      dmy_fv[d] = dmy_phi * (Sv[d] + (Vv[d] - Uv[d]));
	      up[d][im] += dmy_fv[d];
	      force[d] += dmy_fv[d];
	    }
	    {
	      torque[0] += (r[1] * dmy_fv[2] - r[2] * dmy_fv[1]);
	      torque[1] += (r[2] * dmy_fv[0] - r[0] * dmy_fv[2]);
	      torque[2] += (r[0] * dmy_fv[1] - r[1] * dmy_fv[0]);
	    }
	  }
	}//mesh
      } // Hydrodynamic force

      { // particle force	
	dmy_mass = -IMASS_RATIOS[pspec] / dmy_mass; //eff_mass_ratio ?
	for(int d = 0; d < DIM; d++){
	  force[d] *= dmy_mass;
	  torque[d] *= dmy_mass;
	}

	for(int mesh = 0; mesh < np_domain; mesh++){
	  Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
	  int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
	  for(int d = 0; d < DIM; d++){
	    x[d] = r_mesh[d] * dx;
	  }
	  dmy_r = Distance(x, xp);
	  dmy_phi = Phi(dmy_r, radius);

	  for(int d = 0; d < DIM; d++){
	    dmy_fv[d] = dmy_phi * force[d];
	    up[d][im] += dmy_fv[d];
	  }
	}
	
	for(int d = 0; d < DIM; d++){
	  p[n].f_slip[d] = 0.0;
	  p[n].torque_slip[d] = 0.0;
	}
      }// particle force

    }// slip_particle ?
  }// Particle_Number
}

void Make_force_u_slip_particle_norm(double **up, double const* const* u, Particle *p, const CTime &jikan){
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
  double dmy_r, dmy_phi;
  double xp[DIM], vp[DIM], omega_p[DIM], v_rot[DIM];
  double r[DIM], x[DIM], residue[DIM];
  
  double dmy_xi, dmy_theta, dmy_tau, dmy_slip;
  double dmy_mass, dmy_sin, dmy_sin2, dmy_dir;
  double n_r[DIM], n_theta[DIM], n_tau[DIM];
  double Vv[DIM], Uv[DIM], Sv[DIM], dmy_fv[DIM], polar_axis[DIM];
  double force[DIM], torque[DIM];

  //////////////////////// Normalization variables
  int dmy_id;
  double dmy_cs[2], dmy_ds[2], dmy_ca[2];
  double slip_mode, mode_factor, mode_scale[2];
  double stick_factor, stick_scale[2];
  double slip_vel, slip_factor, slip_scale[2];
  const double momentum_flux = 8.0/3.0 * M_PI * radius * radius * dx;
  const double surface_area = 4.0 * M_PI * radius * radius * dx;

  for(int n = 0; n < Particle_Number; n++){
    pspec = p[n].spec;
    if(janus_propulsion[pspec] == slip){

      slip_vel = janus_slip_vel[pspec];
      slip_mode = janus_slip_mode[pspec];
      Janus_direction(polar_axis, p[n]);
      for(int d = 0; d < DIM; d++){
	xp[d] = p[n].x[d];
	vp[d] = p[n].v[d];
	omega_p[d] = p[n].omega[d];

	Vv[d] = Sv[d] = Uv[d] = 0.0;
	force[d] = torque[d] = 0.0;
      }
      sw_in_cell = Particle_cell(xp, dx, x_int, residue);
      sw_in_cell = 1;

      /*
      { // compute surface normalization factors
	dmy_mass = 0.0;
	dmy_cs[0] = dmy_cs[1] = 0.0;
	dmy_ca[0] = dmy_ca[1] = 0.0;
	dmy_ds[0] = dmy_ds[1] = 0.0;

	for(int mesh = 0; mesh < np_domain; mesh++){
	  Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
	  int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
	  for(int d = 0; d < DIM; d++){
	    x[d] = r_mesh[d] * dx;
	  }
	  dmy_r = Distance(x, xp);
	  dmy_phi = 1.0 - Phi(dmy_r, radius);
	  dmy_mass += (1.0 - dmy_phi);
	  dmy_xi = ABS(dmy_r - radius);
	  
	  if(dmy_xi <= HXI && dmy_phi > 0.0){//interface domain
	    Angular2v(omega_p, r, v_rot);
	    Spherical_coord(r, n_r, n_theta, n_tau, dmy_r, dmy_theta, dmy_tau, p[n]);
	    dmy_sin = sin(dmy_theta);
	    dmy_sin2 = dmy_phi * dmy_sin * sin(2.0 * dmy_theta);
	    dmy_sin *= dmy_phi * dmy_sin;
	    dmy_id = ((dmy_sin2 >= 0.0) ? 1 : 0);

	    dmy_cs[dmy_id] += dmy_sin;
	    dmy_ca[dmy_id] += dmy_sin2;
	    dmy_ds[dmy_id] += dmy_phi;
	  }
	}//mesh
      } // compute normalization
      */

      { // compute hydrodynamic force 
	dmy_mass = 0.0;
	for(int d = 0; d < DIM; d++){
	  force[d] = torque[d] = 0.0;
	}
	for(int d = 0; d < 2; d++){
	  slip_scale[d] = momentum_flux / (2.0 * dx3 * ABS(dmy_cs[d]));
	  mode_scale[d] = surface_area / (4.0 * dx3 * ABS(dmy_ca[d]));
	  stick_scale[d] = surface_area / (2.0 * dx3 * ABS(dmy_ds[d]));
	}
	dmy_cs[0] = dmy_cs[1] = 0.0;
	dmy_ca[0] = dmy_ca[1] = 0.0;

	for(int mesh = 0; mesh < np_domain; mesh++){
	  Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
	  int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
	  for(int d = 0; d < DIM; d++){
	    x[d] = r_mesh[d] * dx;
	  }
	  dmy_r = Distance(x, xp);
	  dmy_phi = 1.0 - Phi(dmy_r, radius);
	  dmy_mass += (1.0 - dmy_phi);
	  dmy_xi = ABS(dmy_r - radius);

	  if(dmy_xi <= HXI && dmy_phi > 0.0){//interface domain
	    Angular2v(omega_p, r, v_rot);
	    Spherical_coord(r, n_r, n_theta, n_tau, dmy_r, dmy_theta, dmy_tau, p[n]);
	    dmy_sin = sin(dmy_theta);
	    dmy_sin2 = sin(2.0 * dmy_theta);
	    //	    dmy_id = (dmy_sin2 >= 0.0 ? 1 : 0);

	    //	    slip_factor = slip_scale[dmy_id];
	    //	    mode_factor = 0.5 * (mode_scale[0] + mode_scale[1]);
	    //	    stick_factor = stick_scale[dmy_id];
	    slip_factor = 0.847;
	    mode_factor = 0.847;
	    stick_factor = 1.0;

	    //	    dmy_cs[dmy_id] += slip_factor * dmy_phi * dmy_sin * dmy_sin;
	    //	    dmy_ca[dmy_id] += mode_factor * dmy_phi * dmy_sin * dmy_sin2;

	    for(int d = 0; d < DIM; d++){
	      dmy_fv[d] = (vp[d] + v_rot[d]);
	    }

	    slip_factor *= slip_vel;
	    mode_factor *= slip_vel * slip_mode;
	    for(int d = 0; d < DIM; d++){
	      //	      Uv[d] = (u[0][im] * n_theta[0] + u[1][im] * n_theta[1] + u[2][im] * n_theta[2]);
	      Uv[d] = u[d][im];
	      Vv[d] = stick_factor * (dmy_fv[0] * n_theta[0] + dmy_fv[1] * n_theta[1] + dmy_fv[2] * n_theta[2]);
	      Sv[d] = (slip_factor * dmy_sin) + (mode_factor * sin(2.0 * dmy_theta));
	    }

	    for(int d = 0; d < DIM; d++){
	      dmy_dir = dmy_phi * n_theta[d];
	      dmy_fv[d] = dmy_dir * (Vv[d] + Sv[d]) - dmy_phi * Uv[d];

	      up[d][im] += dmy_fv[d];
	      force[d] += dmy_fv[d];
	    }
	    {
	      torque[0] += (r[1] * dmy_fv[2] - r[2] * dmy_fv[1]);
	      torque[1] += (r[2] * dmy_fv[0] - r[0] * dmy_fv[2]);
	      torque[2] += (r[0] * dmy_fv[1] - r[1] * dmy_fv[0]);
	    }
	  }
	}//mesh
      }
      { // particle force	
	dmy_mass = -IMASS_RATIOS[pspec] / dmy_mass; //eff_mass_ratio ?
	for(int d = 0; d < DIM; d++){
	  force[d] *= dmy_mass;
	  torque[d] *= dmy_mass;
	}

	for(int mesh = 0; mesh < np_domain; mesh++){
	  Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
	  int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
	  for(int d = 0; d < DIM; d++){
	    x[d] = r_mesh[d] * dx;
	  }
	  dmy_r = Distance(x, xp);
	  dmy_phi = Phi(dmy_r, radius);

	  for(int d = 0; d < DIM; d++){
	    dmy_fv[d] = dmy_phi * force[d];
	    up[d][im] += dmy_fv[d];
	  }
	}

	/*
	dmy_mass = (MASS[pspec] / jikan.dt_md); 
	for(int d = 0; d < DIM; d++){
	  p[n].f_slip[d] = (dmy_mass * force[d]);
	}
	if(ROTATION){
	  for(int d = 0; d < DIM; d++){
	    p[n].torque_slip[d] = (dmy_mass * force[d]);
	  }
	}
	fprintf(stderr, "%10.8E %10.8E %10.8E\n", 
		torque[0], torque[1], torque[2]);
	*/
      }

    }// slip_particle ?
  }// Particle_Number
}
