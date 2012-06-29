#include "operate_surface.h"
void Make_force_u_slip_particle(double **up, double const* const* u, Particle *p, const CTime &jikan){
  ////////////////////////
  const double dx = DX;
  const double dx3 = DX3;
  const int np_domain = NP_domain;
  int const* const* sekibun_cell = Sekibun_cell;
  const int* Nlattice = Ns;
  const double radius = RADIUS;
  ////////////////////////
  double xp[DIM], vp[DIM], omega_p[DIM], v_rot[DIM];
  double r[DIM], x[DIM], residue[DIM];
  int x_int[DIM], r_mesh[DIM];
  double dmy_r, dmy_phi;
  int sw_in_cell, pspec;

  double n_r[DIM], n_theta[DIM], n_tau[DIM];
  double slip_vel, slip_mode, dmy_xi, dmy_theta, dmy_tau, dmy_slip;
  double Vv[DIM], Uv[DIM], Sv[DIM], force[DIM], torque[DIM], dmy_fv[DIM], polar_axis[DIM];
  double dmy_cs, dmy_ca, dmy_cs0, dmy_cs1, dmy_ca0, dmy_ca1, dmy_ds;
  double dmy_mass, dmy_sin, dmy_sin2, dmy_dir;
  double slip_scale, stick_scale;
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

      { // compute normalization factors
	dmy_cs = dmy_ca  = dmy_ds = dmy_mass = 0.0;
	dmy_cs0 = dmy_cs1 = dmy_ca0 = dmy_ca1 = 0.0;

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
	    
	    dmy_ds += dmy_phi;
	    dmy_cs += dmy_sin;
	    dmy_ca += dmy_sin2;
	    if(dmy_sin2 >= 0.0){
	      dmy_cs1 += dmy_sin;
	      dmy_ca1 += dmy_sin2;
	    }else{
	      dmy_cs0 += dmy_sin;
	      dmy_ca0 += dmy_sin2;
	    }
	  }
	}//mesh
	slip_scale = momentum_flux / (dx3 * dmy_cs);
	stick_scale = surface_area / (dx3 * dmy_ds);
	
	fprintf(stderr, "%5.3f %5.3f %5.3f\n",
		slip_scale, stick_scale, dmy_ca, 
		momentum_flux / (2.0 * dmy_cs1), momentum_flux / (2.0 * dmy_cs0),
		surface_area / (4.0 * dmy_ca1), surface_area / (4.0 * dmy_ca0)
		);
      } // compute normalization
      
      { // compute hydrodynamic force 
	for(int d = 0; d < DIM; d++){
	  force[d] = torque[d] = 0.0;
	}
	stick_scale = 1.0;
	slip_scale *= slip_vel;
	for(int mesh = 0; mesh < np_domain; mesh++){
	  Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
	  int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
	  for(int d = 0; d < DIM; d++){
	    x[d] = r_mesh[d] * dx;
	  }
	  dmy_r = Distance(x, xp);
	  dmy_phi = 1.0 - Phi(dmy_r, radius);
	  dmy_xi = ABS(dmy_r - radius);

	  if(dmy_xi <= HXI && dmy_phi > 0.0){//interface domain
	    Angular2v(omega_p, r, v_rot);
	    Spherical_coord(r, n_r, n_theta, n_tau, dmy_r, dmy_theta, dmy_tau, p[n]);
	    dmy_sin = sin(dmy_theta);

	    for(int d = 0; d < DIM; d++){
	      dmy_fv[d] = (vp[d] + v_rot[d]);
	    }
	    for(int d = 0; d < DIM; d++){
	      Uv[d] = (u[0][im] * n_theta[0] + u[1][im] * n_theta[1] + u[2][im] * n_theta[2]);
	      Vv[d] = (dmy_fv[0] * n_theta[0] + dmy_fv[1] * n_theta[1] + dmy_fv[2] * n_theta[2]);
	      Sv[d] = (dmy_sin + slip_mode * sin(2.0 * dmy_theta));
	    }
	    for(int d = 0; d < DIM; d++){
	      dmy_dir = dmy_phi * n_theta[d];
	      dmy_fv[d] = dmy_dir * (stick_scale * Vv[d] + slip_scale * Sv[d] - Uv[d]);

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

	dmy_mass = (MASS[pspec] / jikan.dt_md); 
	for(int d = 0; d < DIM; d++){
	  p[n].f_slip[d] = (dmy_mass * force[d]);
	}
	if(ROTATION){
	  for(int d = 0; d < DIM; d++){
	    p[n].torque_slip[d] = (dmy_mass * force[d]);
	  }
	}
      }

    }// slip_particle ?
  }// Particle_Number
}

