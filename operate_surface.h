/*!
  \file operate_surface.h
  \brief Routines to control slip velocity at particle fluid boundaries
  \author J. Molina
  \date 22/06/2012
  \version 1.0
 */
#ifndef OPERATE_SURFACE_H
#define OPERATE_SURFACE_H

#include "variable.h"
#include "make_phi.h"
#include "particle_solver.h"

/*void Make_f_slip_particle_primitive(double const* const * u,
				    double **up,
				    Particle *p,
				    const int &SW_MC,
				    const double &dx,
				    const int &np_domain,
				    int **sekibun_cell,
				    const int Nlattice[DIM],
				    const double radius = RADIUS);
*/

void Make_force_u_slip_particle(double **up, double const* const* u,
				Particle *p, const CTime &jikan);

inline void update_angular(double* l, const double *x, const double *v, const double alpha=1.0){
  l[0] += alpha * (x[1] * v[2] + x[2] * v[1]);
  l[1] += alpha * (x[2] * v[0] + x[0] * v[2]);
  l[2] += alpha * (x[0] * v[1] + x[1] * v[0]);
}
inline void momentum_check(double const* const* up, Particle *p, const CTime &jikan){
  ////////////////////////
  const double dx = DX;
  const int np_domain = NP_domain;
  int const* const* sekibun_cell = Sekibun_cell;
  const int* Nlattice = Ns;
  const double radius = RADIUS;
  static const double dmy0 = DX3 * RHO;
  /////////////////////////
  double xp[DIM], vp[DIM], omega_p[DIM], v_rot[DIM], r[DIM], x[DIM], fv[DIM],residue[DIM], fu[DIM];
  int x_int[DIM], r_mesh[DIM];
  int sw_in_cell, pspec;
  double dmy_r, dmy_phi, dmy_phic, dmy_mass, dmy_xi;
  /////////////////////////
  double tot_p[DIM], part_p[DIM], fluid_p[DIM], md_p[DIM], full_p[DIM];
  double tot_int_p[DIM], part_int_p[DIM], fluid_int_p[DIM], md_int_p[DIM];

  double tot_w[DIM], part_w[DIM], fluid_w[DIM], md_w[DIM], full_w[DIM];
  double tot_int_w[DIM], part_int_w[DIM], fluid_int_w[DIM], md_int_w[DIM];
  for(int d = 0; d < DIM; d++){
    tot_p[d] = part_p[d] = fluid_p[d] = md_p[d] = 0.0;
    tot_int_p[d] = part_int_p[d] = fluid_int_p[d] = md_int_p[d] = 0.0;

    tot_w[d] = part_w[d] = fluid_w[d] = md_w[d] = 0.0;
    tot_int_w[d] = part_int_w[d] = fluid_int_w[d] = md_int_w[d] = 0.0;
  }
  dmy_mass = 0.0;
  for(int n = 0; n < Particle_Number; n++){
    pspec = p[n].spec;
    for(int d = 0; d < DIM; d++){
      xp[d] = p[n].x[d];
      vp[d] = p[n].v[d];
      omega_p[d] = p[n].omega[d];
    }

    sw_in_cell = Particle_cell(xp, dx, x_int, residue);
    sw_in_cell = 1;
    for(int mesh = 0; mesh < np_domain; mesh++){
      Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
      Angular2v(omega_p, r, v_rot);
      int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
      for(int d = 0; d < DIM; d++){
	x[d] = r_mesh[d] * dx;
	fu[d] = up[d][im];
      }
      dmy_r = Distance(x, xp);
      dmy_phi = Phi(dmy_r, radius);
      dmy_phic = 1.0 - dmy_phi;
      dmy_mass += (1.0 - dmy_phi);
      dmy_xi = ABS(dmy_r - radius);
      for(int d = 0; d < DIM; d++){
	fv[d] = (vp[d] + v_rot[d]);
      }

      for(int d = 0; d < DIM; d++){
	md_p[d] += dmy_phi * fv[d];
	part_p[d] += dmy_phi * fu[d];
	fluid_p[d] += dmy_phic * fu[d];
      }
      {
	update_angular(md_w, x, fv, dmy_phi);
	update_angular(part_w, x, fu, dmy_phi);
	update_angular(fluid_w, x, fu, dmy_phic);
      }

      if(dmy_r <= radius + dx){
	for(int d = 0; d < DIM; d++){
	  tot_p[d] += fu[d];
	}
	{
	  update_angular(tot_w, x, fu);
	}
	
	if(dmy_xi <= HXI){
	  for(int d = 0; d < DIM; d++){
	    tot_int_p[d] += fu[d];
	    md_int_p[d] += dmy_phi * fv[d];
	    part_int_p[d] += dmy_phi * fu[d];
	    fluid_int_p[d] += dmy_phic * fu[d];
	  }
	  {
	    update_angular(tot_int_w, x, fu);
	    update_angular(md_int_w, x, fv, dmy_phi);
	    update_angular(part_int_w, x, fu, dmy_phi);
	    update_angular(fluid_int_w, x, fu, dmy_phic);
	  }
	}
      }
    }
  }

  for(int d = 0; d < DIM; d++){
    full_p[d] = fluid_p[d] = 0.0;
    full_w[d] = fluid_w[d] = 0.0;
  }
  int outside = 1;
  for(int i = 0; i < NX; i++){
    x[0] = (double)i * dx;
    for(int j = 0; j < NY; j++){
      x[1] = (double)j * dx;
      for(int k = 0; k < NZ; k++){
	x[2] = (double)k * dx;

	int im = (i * NY * NZ_) + (j * NZ_) + k;
	for(int d = 0; d < DIM; d++){
	  fu[d] = up[d][im];
	  full_p[d] += fu[d];
	}
	update_angular(full_w, x, fu);

	for(int n = 0; n < Particle_Number; n++){
	  pspec = p[n].spec;
	  for(int d = 0; d < DIM; d++){
	    xp[d] = p[n].x[d];
	    vp[d] = p[n].v[d];
	    omega_p[d] = p[n].omega[d];
	  }
	  dmy_r = Distance(x, xp);
	  dmy_phi = Phi(dmy_r, radius);
	  dmy_phic = 1.0 - dmy_phi;
	  if(dmy_phi > 0.){
	    outside = 0.;
	    break;
	  }
	}

	if(outside){
	  for(int d = 0; d < DIM; d++){
	    fluid_p[d] += up[d][im];
	  }
	  update_angular(fluid_w, x, fu);
	}
      }
    }
  }

  fprintf(stderr, "Momentum : %10.8E %10.8E %10.8E\n", full_p[0], full_p[1], full_p[2]);
}


inline void shit_slip_particle(double **up,
			       double const* const* u,
			       Particle *p){
  //////////////////////////////////
  const double dx = DX;
  const int np_domain = NP_domain;
  int const* const* sekibun_cell = Sekibun_cell;
  const int * Nlattice = Ns;
  const double radius = RADIUS;
  static const double dmy0 = DX3*RHO;
  //////////////////////////////////

  double xp[DIM], vp[DIM], vp_old[DIM], omega_p[DIM], v_rot[DIM], fv[DIM];
  double r[DIM], x[DIM], residue[DIM];
  int x_int[DIM], r_mesh[DIM];
  double dmy_r, dmy_phi;
  int sw_in_cell, pspec;

  double SM[DIM][DIM], SMI[DIM][DIM], Vv[DIM], Sv[DIM], Uv[DIM];
  double n_r[DIM], n_theta[DIM], n_tau[DIM];
  double polar_axis[DIM], force[DIM], torque[DIM];
  double slip_vel, slip_mode, dmy_xi, dmy_theta, dmy_tau, dmy_slip, dmy_dir;
  double dmy_cs, dmy_ds, slip_scale, stick_scale, f_scale, cmass;

#pragma omp parallel for schedule(dynamic, 1)\
  private(xp, vp, omega_p, v_rot, fv, r, x, residue, x_int, r_mesh, dmy_r, dmy_phi, \
	  sw_in_cell, pspec, SM, Vv, Sv, n_r, n_theta, n_tau, polar_axis, force, torque, \
	  slip_vel, slip_mode, dmy_xi, dmy_theta, dmy_tau, dmy_slip, dmy_cs, dmy_ds, slip_scale, f_scale, cmass)
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
      }
      sw_in_cell = Particle_cell(xp, dx, x_int, residue);
      sw_in_cell = 1;

      {// Compute particle slip velocity correction
	dmy_cs = dmy_ds = 0.0;
	for(int d = 0; d < DIM; d++){
	  cmass = 0.0;
	  Vv[d] = Sv[d] = Uv[d] = 0.0;
	  SM[d][0] = SM[d][1] = SM[d][2] = 0.0;
	}

	for(int mesh = 0; mesh < np_domain; mesh++){
	  Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
	  int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
	  for(int d = 0; d < DIM; d++){
	    x[d] = r_mesh[d] * dx;
	  }
	  dmy_r = Distance(x, xp);
	  dmy_phi = 1.0 - Phi(dmy_r, radius);
	  cmass += (1.0 - dmy_phi);
	  dmy_xi = ABS(dmy_r - radius);

	  if(dmy_xi <= HXI && dmy_phi > 0.0){ // interface domain only
	    Spherical_coord(r, n_r, n_theta, n_tau, dmy_r, dmy_theta, dmy_tau, p[n]);	  
	    Angular2v(omega_p, r, v_rot);
	    
	    dmy_cs += dmy_phi * SQ(sin(dmy_theta));
	    dmy_ds += dmy_phi;
	    for(int d = 0; d < DIM; d++){
	      fv[d] = (vp[d] + v_rot[d]);
	    }
	    for(int i = 0; i < DIM; i++){
	      dmy_dir = dmy_phi * n_theta[i];
	      Uv[i] += dmy_dir * (u[0][im] * n_theta[0] + u[1][im] * n_theta[1] + u[2][im] * n_theta[2]);
	      Vv[i] += dmy_dir * (fv[0] * n_theta[0] + fv[1] * n_theta[1] + fv[2] * n_theta[2]);
	      Sv[i] += dmy_dir * (sin(dmy_theta) + slip_mode * sin(2.0 * dmy_theta));
	      for(int j = 0; j < DIM; j++){
		SM[i][j] += dmy_dir * n_theta[j];  
	      }
	    }
	  } // interface domain
	}// mesh

	{
	  //eff_mass_ratio not globally defined, use MASS_RATIO ?
	  cmass *= dmy0;
	  slip_scale = ((8.0*M_PI/3.0*radius*radius) / (dx*dx*dmy_cs));
	  stick_scale = (4.0*M_PI*radius*radius) / (dx*dx*dmy_ds);
	  M_scale(SM, dmy0);
	  M_copy(SMI, SM);

	  // Vv : momentum exchange needed to enforce slip with current velocity
	  // force: particle force needed to ensure momentum conservation
	  for(int d = 0; d < DIM; d++){
	    Vv[d] = dmy0 * (stick_scale * Vv[d] + slip_scale * slip_vel * Sv[d] - Uv[d]);  
	    SMI[d][d] += cmass;
	  }
	  M_inv(SMI);
	  M_v_prod(force, SMI, Vv, -1.0);
	}
      }// slip velocity correction

      {// Compute fluid slip velocity 
	for(int d = 0; d < DIM; d++){ // read updated (slip) particle velocities
	  vp_old[d] = vp[d];
	  vp[d] += force[d];
	  Vv[d] = Sv[d] = 0.0;
	}
	for(int mesh = 0; mesh < np_domain; mesh++){
	  Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
	  int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
	  for(int d = 0; d < DIM; d++){
	    x[d] = r_mesh[d] * dx;
	  }
	  dmy_r = Distance(x, xp);
	  dmy_phi = 1.0 - Phi(dmy_r, radius);
	  dmy_xi = ABS(dmy_r - radius);

	  if(dmy_xi <= HXI && dmy_phi > 0.0){ //interface domain
	    Angular2v(omega_p, r, v_rot);
	    Spherical_coord(r, n_r, n_theta, n_tau, dmy_r, dmy_theta, dmy_tau, p[n]);
	    for(int d = 0; d < DIM; d++){
	      fv[d] = stick_scale*(vp[d] + v_rot[d]) - u[d][im];
	    }
	    dmy_slip = slip_scale * slip_vel * (sin(dmy_theta) + slip_mode * sin(2.0 * dmy_theta))
	      + (fv[0] * n_theta[0] + fv[1] * n_theta[1] + fv[2] * n_theta[2]);

	    for(int d = 0; d < DIM; d++){ // careful with parallel update
	      up[d][im] += dmy_phi * (slip_scale * slip_vel * n_theta[d]);
	    }
	  }// interface domain

	  //Check momentum conservation
	  for(int d = 0; d < DIM; d++){
	    Vv[d] += (1.0 - dmy_phi) * force[d];
	    Sv[d] += up[d][im];
	  }
	} // mesh
      }
      {
	f_scale = -(Sv[0]*polar_axis[0] + Sv[1]*polar_axis[1] + Sv[2]*polar_axis[2]);
	f_scale /= (Vv[0]*polar_axis[0] + Vv[1]*polar_axis[1] + Vv[2]*polar_axis[2]);
	assert(f_scale > 0);
	for(int d = 0; d < DIM; d++){
	  force[d] *= f_scale;
	  p[n].v[d] += force[d];
	  p[n].f_slip_previous[d] = force[d];
	  p[n].torque_slip_previous[d] = 0.0;
	}

	for(int d = 0; d < DIM; d++){
	  fv[d] = f_scale * Vv[d] + Sv[d];
	}
	fprintf(stderr, "slip: %8.6g %8.6g %8.6g %8.6g %8.6g %8.6g %8.6g\n",
		stick_scale, slip_scale, f_scale,
		sqrt(fv[0]*fv[0] + fv[1]*fv[1] + fv[2]*fv[2]),
		fv[0], fv[1], fv[2]);
      }
    }// slip particle ?
  }// Particle Number
}
#endif
