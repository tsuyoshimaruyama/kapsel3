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
  double xp[DIM], vp[DIM], omega_p[DIM], v_rot[DIM];
  double r[DIM], x[DIM], residue[DIM];

  double dmy_xi, dmy_theta, dmy_tau;
  double dmy_sin, dmy_sin2;
  double n_r[DIM], n_theta[DIM], n_tau[DIM], polar_axis[DIM];
  double M0, M1[DIM], M2[DIM][DIM], ST[DIM][DIM], SU[DIM][DIM], SV[DIM][DIM], dv_s[DIM], domega_s[DIM];
  double slip_mode, slip_vel, slip_magnitude;
  double r_x_theta[DIM], r_x_us[DIM], us[DIM];
#pragma omp parallel for schedule(dynamic, 1) \
  private(sw_in_cell, pspec, x_int, r_mesh, dmy_r, dmy_phi, dmy_region, xp, vp, omega_p, \
	  v_rot, r, x, residue, dmy_xi, dmy_theta, dmy_tau, dmy_sin, dmy_sin2, dmy_vdot, \
	  n_r, n_theta, n_tau, dmy_fv, polar_axis, M0, M1, M2, ST, SU, SV, dv_s, domega_s, \
	  slip_mode, slip_vel, slip_magnitude, r_x_theta, r_x_us, us)
  for(int n = 0; n < Particle_Number; n++){
    pspec = p[n].spec;

    for(int i = 0; i < DIM; i++){
      xp[i] = p[n].x[i];
      vp[i] = p[n].v[i];
      omega_p[i] = p[n].omega[i];
      
      M1[i] = dv_s[i] = domega_s[i] = 0.0;
      for(int j = 0; j < DIM; j++){
	M2[i][j] = ST[i][j] = SU[i][j] = SV[i][j] = 0.0;
      }
    }
    M0 = 0.0;
    
    sw_in_cell = Particle_cell(xp, dx, x_int, residue);
    sw_in_cell = 1;
    for(int mesh = 0; mesh < np_domain; mesh++){
      Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
      int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
      for(int d = 0; d < DIM; d++){
	x[d] = r_mesh[d] * dx;
      }
      dmy_r = Distance(x, xp);
      dmy_phi = Phi(dmy_r, radius);
      dmy_xi = ABS(dmy_r - radius);
      
      //particle properties
      M0 += dmy_phi;
      for(int i = 0; i < DIM; i++){
	M1[i] += dmy_phi * r[i];
	M2[i][i] += dmy_phi * dmy_r * dmy_r;
	for(int j = 0; j < DIM; j++){
	  M2[i][j] -= dmy_phi * r[i] * r[j];
	}
      }

      if(janus_propulsion[pspec] == slip){	//surface properties
	slip_vel = janus_slip_vel[pspec] * janus_slip_scale;
	slip_mode = janus_slip_mode[pspec];
	
	if(janus_slip_region == surface_slip){
	  dmy_region = (1.0 - dmy_phi) * DPhi_compact_sin_norm(dmy_r, radius);
	}else{
	  dmy_region = (1.0 - dmy_phi);
	}
	if(dmy_xi <= HXI && dmy_region > 0.0){
	  Angular2v(omega_p, r, v_rot);
	  Spherical_coord(r, n_r, n_theta, n_tau, dmy_r, dmy_theta, dmy_tau, p[n]);
	  dmy_sin = sin(dmy_theta);
	  dmy_sin2 = sin(2.0 * dmy_theta);
	  
	  slip_magnitude = slip_vel * (dmy_sin + slip_mode * dmy_sin2);
	  for(int i = 0; i < DIM; i++){
	    us[i] = n_theta[i] * slip_magnitude - u[i][im];
	  }
	  
	  v_cross(r_x_theta, r, n_theta);
	  v_cross(r_x_us, r, us);
	  for(int i = 0; i < DIM; i++){
	    dv_s[i] += dmy_region * us[i];
	    domega_s[i] += dmy_region * r_x_us[i];
	    
	    for(int j = 0; j < DIM; j++){
	      ST[i][j] += dmy_region * n_theta[i] * n_theta[j];
	      SU[i][j] += dmy_region * n_theta[i] * r_x_theta[j];
	      SV[i][j] += dmy_region * r_x_theta[i] * r_x_theta[j];
	    }
	  }
	}//interface_region
      }//slip_particle
    }//mesh

    p[n].mass = M0;
    for(int i = 0; i < DIM; i++){
      
      p[n].mass_center[i] = M1[i];
      p[n].surface_dv[i] = dv_s[i];
      p[n].surface_domega[i] = domega_s[i];
      
      for(int j = 0; j < DIM; j++){
	p[n].inertia[i][j] = M2[i][j];
	p[n].surfaceT[i][j] = ST[i][j];
	p[n].surfaceU[i][j] = SU[i][j];
	p[n].surfaceV[i][j] = SV[i][j];
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


  for(int d = 0; d < DIM; d++){
    dP[d] = -(p.surfaceT[d][0] * vp[0] + p.surfaceT[d][1] * vp[1] + p.surfaceT[d][2] * vp[2] +
	      p.surfaceU[d][0] * wp[0] + p.surfaceU[d][1] * wp[1] + p.surfaceU[d][2] * wp[2] +
	      p.surface_dv[d]);

    dL[d] = -(p.surfaceU[0][d] * vp[0] + p.surfaceU[1][d] * vp[1] + p.surfaceU[2][d] * vp[2] +
	      p.surfaceV[d][0] * wp[0] + p.surfaceV[d][1] * wp[1] + p.surfaceV[d][2] * wp[2] +
	      p.surface_domega[d]);
  }

  // Compute angular velocity of droplet
  v_cross(L0, mc, dP, imass);
  mc2 = (mc[0] * mc[0] + mc[1] * mc[1] + mc[2] * mc[2]) * imass;
  for(int i = 0; i < DIM; i++){
    L0[i] = dL[i] - L0[i];
    LL[i][i] = -mc2;

    for(int j = 0; j < DIM; j++){
      LL[i][j] += (p.inertia[i][j] + imass * mc[i] * mc[j]);
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

void Make_force_u_slip_particle_noscale(double **up, double const* const* u, Particle *p, const CTime &jikan){
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
  double r[DIM], x[DIM], residue[DIM];
  
  double dmy_xi, dmy_theta, dmy_tau;
  double dmy_sin, dmy_sin2, dmy_vdot, dmy_vslip, slip_mode, slip_vel;
  double n_r[DIM], n_theta[DIM], n_tau[DIM];
  double dmy_fv[DIM], polar_axis[DIM];
  double force_s[DIM], torque_s[DIM], force_p[DIM], torque_p[DIM];

#pragma omp parallel for schedule(dynamic, 1) \
  private(sw_in_cell, pspec, x_int, r_mesh, dmy_r, dmy_phi, dmy_region, xp, vp, omega_p, v_rot, r, x, residue, \
	  dmy_xi, dmy_theta, dmy_tau, dmy_vdot, dmy_sin, dmy_sin2, n_r, n_theta, n_tau, \
	  dmy_fv, polar_axis, force, torque, slip_mode, slip_vel)
  for(int n = 0; n < Particle_Number; n++){
    pspec = p[n].spec;
    
    if(janus_propulsion[pspec] == slip){
      
      slip_vel = janus_slip_vel[pspec] * janus_slip_scale;
      slip_mode = janus_slip_mode[pspec];
      Janus_direction(polar_axis, p[n]);
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
	  x[d] = r_mesh[d] * dx;
	}
	dmy_r = Distance(x, xp);
	dmy_phi = Phi(dmy_r, radius);	 
	dmy_xi = ABS(dmy_r - radius);

	{//fluid particle domain
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
	if(dmy_xi <= HXI && dmy_region > 0.0 ){//interface domain
	  Angular2v(omega_p, r, v_rot);
	  Spherical_coord(r, n_r, n_theta, n_tau, dmy_r, dmy_theta, dmy_tau, p[n]);
	  dmy_sin = sin(dmy_theta);
	  dmy_sin2 = sin(2.0 * dmy_theta);
	  
	  for(int d = 0; d < DIM; d++){
	    dmy_fv[d] = (vp[d] + v_rot[d]);
	  }
	  dmy_vdot = dmy_fv[0]*n_theta[0] + dmy_fv[1]*n_theta[1] + dmy_fv[2]*n_theta[2];
	  dmy_vslip = slip_vel * (dmy_sin + slip_mode * dmy_sin2);
	  for(int d = 0; d < DIM; d++){
	    dmy_fv[d] = dmy_region * (n_theta[d] * (dmy_vdot + dmy_vslip) - u[d][im]);
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
	
	p[n].f_slip[d] = 0.0;
	p[n].torque_slip[d] = 0.0;
      }

      if(v_rms(force_p, force_s) > 1.0e-12 ||
	 v_rms(torque_p, torque_s) > 1.0e-12){
	fprintf(stderr, "# Momentum Conservation Warning : %10.8E %10.8E\n",
		v_rms(force_p, force_s), v_rms(torque_p, torque_s));
	
      }
    }// slip_particle ?
  }// Particle_Number
}

/*
 */
void Make_force_u_slip_particle_scale(double **up, double const* const* u, Particle *p, const CTime &jikan){
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
  double r[DIM], x[DIM], residue[DIM];
  
  double dmy_xi, dmy_theta, dmy_tau;
  double dmy_mass, dmy_sin, dmy_sin2, dmy_vdot;
  double n_r[DIM], n_theta[DIM], n_tau[DIM];
  double Vv[DIM], Uv[DIM], Sv[DIM], dmy_fv[DIM], polar_axis[DIM];
  double force[DIM], torque[DIM];
  double slip_mode, slip_vel;

  //////////////////////// Normalization variables
  int dmy_id;
  double dmy_cs[2], dmy_ds[2], dmy_ca[2];
  double mode_factor, mode_scale[2];
  double area_factor, area_scale[2];
  double slip_factor, slip_scale[2];
  const double momentum_flux = 8.0/3.0 * M_PI * radius * radius * dx;
  const double surface_area = 4.0 * M_PI * radius * radius * dx;

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
      sw_in_cell = Particle_cell(xp, dx, x_int, residue);
      sw_in_cell = 1;

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

	  if(janus_slip_region == surface_slip){
	    dmy_region = DPhi_compact_sin_norm(dmy_r, radius);
	  }else{
	    dmy_region = 1.0;
	  }
	  
	  if(dmy_xi <= HXI && dmy_phi > 0.0){//interface domain
	    Angular2v(omega_p, r, v_rot);
	    Spherical_coord(r, n_r, n_theta, n_tau, dmy_r, dmy_theta, dmy_tau, p[n]);
	    dmy_sin = sin(dmy_theta);
	    dmy_sin2 = sin(2.0 * dmy_theta);
	    dmy_id = ((dmy_sin2 >= 0.0) ? 1 : 0);

	    dmy_cs[dmy_id] += dmy_region * dmy_phi * dmy_sin * dmy_sin;
	    dmy_ca[dmy_id] += dmy_region * dmy_phi * dmy_sin * dmy_sin2;
	    dmy_ds[dmy_id] += dmy_region * dmy_phi;
	  }
	}//mesh
      } // compute normalization

      { // compute hydrodynamic force 
	dmy_mass = 0.0;
	for(int d = 0; d < DIM; d++){
	  force[d] = torque[d] = 0.0;
	}
	for(int d = 0; d < 2; d++){
	  slip_scale[d] = momentum_flux / (2.0 * dx3 * ABS(dmy_cs[d]));
	  mode_scale[d] = surface_area / (4.0 * dx3 * ABS(dmy_ca[d]));
	  area_scale[d] = surface_area / (2.0 * dx3 * ABS(dmy_ds[d]));
	}
	fprintf(stderr, "%8.6E %8.6E %8.6E %8.6E %8.6E %8.6E\n",
		slip_scale[0], slip_scale[1], 
		mode_scale[0], mode_scale[1],
		area_scale[0], area_scale[1]);
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

	  if(janus_slip_region == surface_slip){
	    dmy_region = DPhi_compact_sin_norm(dmy_r, radius);
	  }else{
	    dmy_region = 1.0;
	  }

	  if(dmy_xi <= HXI && dmy_phi > 0.0){//interface domain
	    Angular2v(omega_p, r, v_rot);
	    Spherical_coord(r, n_r, n_theta, n_tau, dmy_r, dmy_theta, dmy_tau, p[n]);
	    dmy_sin = sin(dmy_theta);
	    dmy_sin2 = sin(2.0 * dmy_theta);
	    dmy_id = (dmy_sin2 >= 0.0 ? 1 : 0);

	    slip_factor = slip_scale[dmy_id];
	    mode_factor = mode_scale[dmy_id];
	    area_factor = area_scale[dmy_id];

	    dmy_cs[dmy_id] += slip_factor * dmy_region * dmy_phi * dmy_sin * dmy_sin;
	    dmy_ca[dmy_id] += mode_factor * dmy_region * dmy_phi * dmy_sin * dmy_sin2;
	    dmy_ds[dmy_id] += area_factor * dmy_region * dmy_phi;

	    for(int d = 0; d < DIM; d++){
	      dmy_fv[d] = (vp[d] + v_rot[d]);
	    }
	    dmy_vdot = dmy_fv[0]*n_theta[0] + dmy_fv[1]*n_theta[1] + dmy_fv[2]*n_theta[2];

	    slip_factor *= slip_vel;
	    mode_factor *= slip_vel * slip_mode;
	    for(int d = 0; d < DIM; d++){
	      Uv[d] = u[d][im];
	      Vv[d] = n_theta[d] * dmy_vdot;
	      Sv[d] = n_theta[d] * (slip_factor*dmy_sin + mode_factor * dmy_sin2);
	    }

	    for(int d = 0; d < DIM; d++){
	      dmy_fv[d] = dmy_phi * dmy_region * (Sv[d] + (Vv[d] - Uv[d]));
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

      for(int d = 0; d < 2; d++){
	  slip_scale[d] = momentum_flux / (2.0 * dx3 * ABS(dmy_cs[d]));
	  mode_scale[d] = surface_area / (4.0 * dx3 * ABS(dmy_ca[d]));
	  area_scale[d] = surface_area / (2.0 * dx3 * ABS(dmy_ds[d]));
	}
	dmy_cs[0] = dmy_cs[1] = 0.0;
	dmy_ca[0] = dmy_ca[1] = 0.0;
	dmy_ds[0] = dmy_ds[1] = 0.0;

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
