//
// $Id: make_phi.cxx,v 1.3 2006/10/20 00:24:36 iwashita Exp $
//

#include "make_phi.h"

void (*Angular2v)(const double *omega, const double *r, double *v);

int NP_domain;
int NP_domain_interface;
int **Sekibun_cell;
int **Sekibun_cell_interface;

/////////////
void Make_surface_normal(double **surface_normal
			   ,const Particle *p){
  
  for(int d=0;d<DIM;d++){
    Reset_phi(surface_normal[d]);
  }

    double xp[DIM];
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell;
    int r_mesh[DIM];
    double r[DIM];
    double dmy_r;
    double dmy;
	double ir;
#pragma omp parallel for schedule(dynamic, 1) private(xp,x_int,residue,sw_in_cell,r_mesh,r,dmy_r,dmy,ir)
  for(int n=0; n < Particle_Number; n++){
    for(int d=0;d<DIM;d++){
      xp[d] = p[n].x[d];
      
      assert(xp[d] >= 0);
      assert(xp[d] < L[d]);
    }
    sw_in_cell 
      = Particle_cell(xp, DX, x_int, residue);// {1,0} $B$,JV$C$F$/$k(B
    sw_in_cell = 1;
    for(int mesh=0; mesh < NP_domain; mesh++){
      Relative_coord(Sekibun_cell[mesh]
		     ,x_int, residue, sw_in_cell, Ns, DX, r_mesh, r);
      dmy_r= sqrt(SQ(r[0])+SQ(r[1])+SQ(r[2]));
      dmy = ABS(dmy_r - RADIUS);
      if(dmy < HXI){
	  ir = 1./dmy_r;
	for(int d=0;d<DIM;d++){
	  surface_normal[d][(r_mesh[0]*NY*NZ_)+(r_mesh[1]*NZ_)+r_mesh[2]] 
	    += r[d]*ir;
	}
      }
    }
  }
}
/////////////
inline void Make_rho_field_primitive(double *phi
				     ,Particle *p
				     ,const double &dx
				     ,const int &np_domain
				     ,int **sekibun_cell
				     ,int Nlattice[DIM]
    ){
    
    double drho;
    double xp[DIM];
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell;
    int r_mesh[DIM];
    double r[DIM];
    double x[DIM];
    double dmy;
    double dmy_phi;
    int im;
#pragma omp parallel for schedule(dynamic, 1) private(drho,xp,x_int,residue,sw_in_cell,r_mesh,r,x,dmy,dmy_phi)
    for(int n=0; n < Particle_Number; n++){
	drho = RHO_particle[p[n].spec] - RHO;
	for(int d=0;d<DIM;d++){
	    xp[d] = p[n].x[d];
	}
	
	sw_in_cell 
	    = Particle_cell(xp, dx, x_int, residue);// {1,0} $B$,JV$C$F$/$k(B
	sw_in_cell = 1;
	for(int mesh=0; mesh < np_domain; mesh++){
	    Relative_coord(sekibun_cell[mesh]
			   , x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
	    for(int d=0;d<DIM;d++){
		x[d] = r_mesh[d] * dx;
	    }
	    dmy=Distance(x,xp);
	    dmy_phi=Phi(dmy)*drho;
	    phi[(r_mesh[0]*NY*NZ_)+(r_mesh[1]*NZ_)+r_mesh[2]] += dmy_phi;
	}
    }
#pragma omp parallel for schedule(dynamic, 1) private(im)
    for(int i=0;i<Nlattice[0];i++){
	for(int j=0;j<Nlattice[1];j++){
	    for(int k=0;k<Nlattice[2];k++){
		im=(i*NY*NZ_)+(j*NZ_)+k;
		phi[im] += RHO; 
	    }
	}
    }
}
void Make_rho_field(double *phi
		    ,Particle *p
		    ){
  int *nlattice;
  nlattice = Ns;
  Make_rho_field_primitive(phi, p, DX, NP_domain,Sekibun_cell,nlattice);
}

//
inline double janus_geometry(const Particle &p, const double normal[DIM]){
  double body_normal[DIM];
  double cos_theta;

  rigid_body_rotation(body_normal, normal, p.q, SPACE2BODY);
  if(janus_axis[p.spec] == x_axis){
    cos_theta = body_normal[0];
  }else if(janus_axis[p.spec] == y_axis){
    cos_theta = body_normal[1];
  }else if(janus_axis[p.spec] == z_axis){
    cos_theta = body_normal[2];
  }
  else if(janus_axis[p.spec] == no_axis){
    cos_theta = 1.0;
  }else{
    fprintf(stderr, "Error: %d not a janus particle\n", p.spec);
    exit_job(EXIT_FAILURE);
  }

  return ((cos_theta >= 0.0) ? 1.0 : -1.0);
}

void Make_phi_janus_particle(double *phi, double *id_phi, Particle *p){
  double xp[DIM], residue[DIM], r[DIM];
  double dmy_r;
  int x_int[DIM], r_mesh[DIM];
  int *nlattice;
  int im, sw_in_cell;
  nlattice = Ns;

  //Reset janus labels for grid points
  Reset_phi(id_phi, -1.0);

  for(int n = 0; n < Particle_Number; n++){
    
    for(int d = 0; d < DIM; d++){
      xp[d] = p[n].x[d];
    }
    
    sw_in_cell = Particle_cell(xp, DX, x_int, residue);
    sw_in_cell = 1;
    
    for(int mesh = 0; mesh < NP_domain; mesh++){
      Relative_coord(Sekibun_cell[mesh], x_int, residue, sw_in_cell,
		     nlattice, DX, r_mesh, r);
      dmy_r = sqrt(SQ(r[0]) + SQ(r[1]) + SQ(r[2]));
      for(int d = 0; d < DIM; d++){
	r[d] /= dmy_r;
      }
      im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
      
      //New mesh point or closer particle
      if(id_phi[im] < 0 || dmy_r < id_phi[im]){
	phi[im] = janus_geometry(p[n], r) * ABS(phi[im]);
	id_phi[im] = dmy_r;
      }
      
    }//mesh
  }//Particle Number
}

void Make_phi_janus_particle_OBL(double *phi, double *id_phi, Particle *p){
  double xp[DIM], residue[DIM], r[DIM], x[DIM];
  double dmy_r;
  int x_int[DIM], r_mesh[DIM];
  int *nlattice;
  int im, sign, sw_in_cell;
  nlattice = Ns;

  // Reset janus labels for grid points
  Reset_phi(id_phi, -1.0);

  for(int n = 0; n < Particle_Number; n++){
    for(int d = 0; d < DIM; d++){
      xp[d] = p[n].x[d];
    }
    
    sw_in_cell = Particle_cell(xp, DX, x_int, residue);
    sw_in_cell = 1;
    
    for(int mesh = 0; mesh < NP_domain; mesh++){
      sign = Relative_coord_check_stepover_Y(Sekibun_cell[mesh], x_int, 
					     residue, sw_in_cell, nlattice, 
					     DX, r_mesh, r);
      dmy_r = sqrt(SQ(r[0]) + SQ(r[1]) + SQ(r[2]));
      for(int d = 0; d < DIM; d++){
	r[d] /= dmy_r;
      }
      im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
      
      //New mesh point or closer particle
      if(id_phi[im] < 0 || dmy_r < id_phi[im]){
	phi[im] = janus_geometry(p[n], r) * ABS(phi[im]);
	id_phi[im] = dmy_r;
      }
      
    }//mesh
  }//Particle Number
}


inline void Make_phi_u_primitive(double *phi
				 ,double **up
				 ,Particle *p
				 ,const int &SW_UP
				 ,const double &dx
				 ,const int &np_domain
				 ,int **sekibun_cell
				 ,const int Nlattice[DIM]
				 ,const double radius = RADIUS
    ){
    
    double xp[DIM],vp[DIM],omega_p[DIM];
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell; 
    int r_mesh[DIM];
    double r[DIM];
    double x[DIM];
    double dmy;
    double dmy_phi;
    double v_rot[DIM];
    int im;
#pragma omp parallel for schedule(dynamic, 1) private(xp,vp,omega_p,x_int,residue,sw_in_cell,r_mesh,r,x,dmy,dmy_phi,v_rot)
    for(int n=0; n < Particle_Number; n++){
	for(int d=0;d<DIM;d++){
	    xp[d] = p[n].x[d];
	    vp[d] = p[n].v[d];
	    omega_p[d] = p[n].omega[d];
	}
	
	sw_in_cell 
	    = Particle_cell(xp, dx, x_int, residue);// {1,0} $B$,JV$C$F$/$k(B
	sw_in_cell = 1;
	for(int mesh=0; mesh < np_domain; mesh++){
	    Relative_coord(sekibun_cell[mesh]
			   , x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
	    //dmy = 0.;
	    for(int d=0;d<DIM;d++){
		x[d] = r_mesh[d] * dx;
		//dmy += SQ(r[d]);
	    }
	    dmy = Distance(x, xp);
	    //dmy = sqrt(dmy);
	    dmy_phi= Phi(dmy,radius);
	    im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2]; 
	    phi[im] += dmy_phi;
	    
	    if(SW_UP){
		Angular2v(omega_p, r, v_rot);
		for(int d=0;d<DIM;d++){
		    up[d][im] += ( (vp[d]+v_rot[d]) * dmy_phi);
		}
	    }
	}
	p[n].eff_mass_ratio 
	    = 1.0;
    }
    
    // koba code //
    if(SW_UP){
	double idmy_phi;
#pragma omp parallel for schedule(dynamic, 1) private(im,idmy_phi)
	for (int i = 0; i < NX; i++) {
	    for (int j = 0; j < NY; j++) {
		for (int k = 0; k < NZ; k++) {
		    im=(i*NY*NZ_)+(j*NZ_)+k;
		    if (phi[im] > 1.){
			idmy_phi=1./phi[im];
			up[0][im] *= idmy_phi;
			up[1][im] *= idmy_phi;
			up[2][im] *= idmy_phi;
			phi[im] = 1.;
		    }
		}
	    }
	}
    }
}
inline void Make_phi_u_primitive_OBL(double *phi
				     ,double **up
				     ,Particle *p
				     ,const int &SW_UP
				     ,const double &dx
				     ,const int &np_domain
				     ,int **sekibun_cell
				     ,const int Nlattice[DIM]
				     ,const double radius = RADIUS
    ){
    
    double xp[DIM],vp[DIM],omega_p[DIM];
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell; 
    int r_mesh[DIM];
    double r[DIM];
    double x[DIM];
    double dmy;
    double dmy_phi;
    double v_rot[DIM];
    int sign;
    int im;
#pragma omp parallel for schedule(dynamic, 1) private(xp,vp,omega_p,x_int,residue,sw_in_cell,r_mesh,r,x,dmy,dmy_phi,v_rot,sign)
    for(int n=0; n < Particle_Number; n++){
	for(int d=0;d<DIM;d++){
	    xp[d] = p[n].x[d];
	    vp[d] = p[n].v[d];
	    omega_p[d] = p[n].omega[d];
	}
	
	sw_in_cell 
	    = Particle_cell(xp, dx, x_int, residue);// {1,0} $B$,JV$C$F$/$k(B
	sw_in_cell = 1;
	for(int mesh=0; mesh < np_domain; mesh++){
	    sign = Relative_coord_check_stepover_Y(sekibun_cell[mesh]
						   , x_int
						   , residue
						   , sw_in_cell
						   , Nlattice
						   , dx
						   , r_mesh
						   , r);

	    dmy = 0.;
	    for(int d=0;d<DIM;d++){
		x[d] = r_mesh[d] * dx;
		dmy+= SQ(r[d]);
	    }
	    dmy = sqrt(dmy);

	    dmy_phi= Phi(dmy,radius);
	    im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
	    phi[im] += dmy_phi;
	    
	    if(SW_UP){
		Angular2v(omega_p, r, v_rot);
		for(int d=0;d<DIM;d++){
		    if (!(d == 0)) {
			up[d][im] += (vp[d] + v_rot[d])*dmy_phi;
		    } else {
			up[0][im] += ( (vp[d] - sign*Shear_rate_eff*L_particle[1] + v_rot[d])*dmy_phi);
		    }
		}
	    }
	}
	p[n].eff_mass_ratio 
	    = 1.0;
    }
    
    // koba code //
    if(SW_UP){
	double idmy_phi;
#pragma omp parallel for schedule(dynamic, 1) private(im,idmy_phi)
	for (int i = 0; i < NX; i++) {
	    for (int j = 0; j < NY; j++) {
		for (int k = 0; k < NZ; k++) {
		    im=(i*NY*NZ_)+(j*NZ_)+k;
		    if (phi[im] > 1.){
			idmy_phi=1./phi[im];
			up[0][im] *= idmy_phi;
			up[1][im] *= idmy_phi;
			up[2][im] *= idmy_phi;
			phi[im] = 1.;
		    }
		}
	    }
	}
    }
}

void Make_u_slip_particle(double const* const* u,
			  double **up,
			  Particle *p,
			  const CTime &jikan){
  const double dx = DX;
  const int np_domain = NP_domain;
  int const * const * sekibun_cell = Sekibun_cell;
  const int * Nlattice = Ns;
  const double radius = RADIUS;
  static const double dmy0 = -DX3*RHO;
  double dmy = dmy0 / jikan.dt_fluid;
  //////

  double xp[DIM], vp[DIM], omega_p[DIM];
  int x_int[DIM];
  double residue[DIM];
  int sw_in_cell;
  int r_mesh[DIM];
  double r[DIM];
  double x[DIM];
  double dmy_r;
  double dmy_phi;
  double v_rot[DIM];

  double n_r[DIM], n_theta[DIM], n_tau[DIM];
  double delta_v[DIM];
  double slip_vel, slip_mode, dmy_xi, dmy_theta, dmy_tau, dmy_slip;
  double polar_axis[DIM], force[DIM], torque[DIM], dmy_fp[DIM];
  double dmy_cs, slip_scale;

  double MS[DIM][DIM], QS[DIM];

  int slip_flag = 0;
#pragma omp parallel for schedule(dynamic, 1) \
  private(xp, vp, omega_p, x_int, residue, sw_in_cell, r_mesh, r, x, dmy_r, dmy_phi, \
	  v_rot, n_r, n_theta, n_tau, delta_v, slip_vel, slip_mode, dmy_xi, dmy_theta, dmy_tau, dmy_slip, \
	  polar_axis, force, torque, dmy_fp, dmy_cs, slip_scale)
  for(int n = 0; n < Particle_Number; n++){

    for(int d = 0; d < DIM; d++){
      force[d] = 0.0;
      torque[d] = 0.0;
      QS[d] = 0.0;
      MS[d][0] = MS[d][1] = MS[d][2] = 0.0;
    }

    if(janus_propulsion[p[n].spec] == slip){
      slip_flag = 1;
      slip_vel = janus_slip_vel[p[n].spec];
      slip_mode = janus_slip_mode[p[n].spec];
      
      for(int d = 0; d < DIM; d++){
	xp[d] = p[n].x[d];
	vp[d] = p[n].v[d];
	omega_p[d] = p[n].omega[d];
      }
      
      sw_in_cell = Particle_cell(xp, dx, x_int, residue);
      sw_in_cell = 1;

      dmy_cs = 0.0;
      Janus_direction(polar_axis, p[n]);
      double dmy_v[DIM];
      dmy_v[0] = dmy_v[1] = dmy_v[2] = 0.0;
      // precompute momentum flux along janus axis
      for(int mesh = 0; mesh < np_domain; mesh++){
	Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
	int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
	
	for(int d = 0; d < DIM; d++){
	  x[d] = r_mesh[d] * dx;
	}
	dmy_r = Distance(x, xp);
	dmy_phi = 1.0 - Phi(dmy_r, radius);
	dmy_xi = ABS(dmy_r - radius);

	if(dmy_xi < HXI && dmy_phi > 0.0){
	  Spherical_coord(r, n_r, n_theta, n_tau, dmy_r, dmy_theta, dmy_tau, p[n]);
	  dmy_cs -= slip_vel * sin(dmy_theta) * sin(dmy_theta) * dmy_phi;
	  for(int i = 0; i < DIM; i++){
	    for(int j = 0; j < DIM; j++){
	      MS[i][j] = dmy_phi * n_theta[i] * n_theta[j];
	    }
	  }
	}	
      }//raw mesh
      fprintf(stderr, "%10.8g   %10.8g   %10.8g   %10.8g   %10.8g   %10.8g\n",
	      MS[0][0] * p[n].v[0] + MS[0][1] * p[n].v[1] + MS[0][2] * p[n].v[2],
	      MS[1][0] * p[n].v[0] + MS[1][1] * p[n].v[1] + MS[1][2] * p[n].v[2],
	      MS[2][0] * p[n].v[0] + MS[2][1] * p[n].v[1] + MS[2][2] * p[n].v[2],
	      dmy_v[0], dmy_v[1], dmy_v[2]);
      MS[0][0] += 1.0;
      MS[1][1] += 1.0;
      MS[2][2] += 1.0;

      // normalized slip velocity
      slip_scale = (-8.0*M_PI/3.0*slip_vel*RADIUS*RADIUS) / (DX*DX*dmy_cs);
      for(int mesh = 0; mesh < np_domain; mesh++){
	Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
	int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];

	for(int d = 0; d < DIM; d++){
	  x[d] = r_mesh[d] * dx;
	}
	dmy_r = Distance(x,xp);
	dmy_phi = 1.0 - Phi(dmy_r, radius);
	dmy_xi = ABS(dmy_r - radius);

	if(dmy_xi < HXI && dmy_phi > 0.0){
	  Angular2v(omega_p, r, v_rot);
	  Spherical_coord(r, n_r, n_theta, n_tau, dmy_r, dmy_theta, dmy_tau, p[n]);

	  for(int d = 0; d < DIM; d++){
	    delta_v[d] = (vp[d] + v_rot[d]) - u[d][im];
	  }
	  dmy_slip = slip_scale * slip_vel * (sin(dmy_theta) + slip_mode * sin(2.0 * dmy_theta))
	    + (delta_v[0] * n_theta[0] + delta_v[1] * n_theta[1] + delta_v[2] * n_theta[2]);
	  for(int d = 0; d < DIM; d++){
	    dmy_fp[d] = (dmy_slip * n_theta[d] * dmy_phi);
	    force[d] += dmy_fp[d];
	  }

	  {//torque
	    torque[0] += (r[1] * dmy_fp[2] - r[2] * dmy_fp[1]);
	    torque[1] += (r[2] * dmy_fp[0] - r[0] * dmy_fp[2]);
	    torque[2] += (r[0] * dmy_fp[1] - r[1] * dmy_fp[0]);
	  }
	}//interface domain
      }//normalized mesh
    }//slip particle?

    M_inv(MS);
    M_v_prod(force, MS);
    for(int d = 0; d < DIM; d++){
      vp[d] += jikan.hdt_md * IMASS[p[n].spec] * 2.0 * force[d];
      p[n].v[d] += jikan.hdt_md * IMASS[p[n].spec] * 2.0 * force[d];
      p[n].f_slip_previous = p[n].f_slip[d];
      p[n].f_slip[d] = (dmy * force[d] * p[n].eff_mass_ratio);
    }
    if(ROTATION){
      for(int d = 0; d < DIM; d++){
	p[n].torque_slip_previous[d] = p[n].torque_slip[d];
	p[n].torque_slip[d] = (dmy * torque[d] * p[n].eff_mass_ratio);
      }
    }

    for(int mesh = 0; mesh < np_domain; mesh++){
      Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
      int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];

      for(int d = 0; d < DIM; d++){
	x[d] = r_mesh[d] * dx;
      }
      dmy_r = Distance(x,xp);
      dmy_phi = 1.0 - Phi(dmy_r, radius);
      dmy_xi = ABS(dmy_r - radius);

      if(dmy_xi < HXI && dmy_phi > 0.0){
	for(int d = 0; d < DIM; d++){
	  delta_v[d] = (vp[d] + v_rot[d]) - u[d][im];
	  dmy_slip = slip_scale * slip_vel * (sin(dmy_theta) + slip_mode * sin(2.0 * dmy_theta))
	    + (delta_v[0] * n_theta[0] + delta_v[1] * n_theta[1] + delta_v[2] * n_theta[2]);
	}
	for(int d = 0; d < DIM; d++){
	  dmy_fp[d] = (dmy_slip * n_theta[d] * dmy_phi);
	  up[d][im] += dmy_fp[d];
	}
      }
    //    fprintf(stderr, "slip: %10.8g %10.8g %10.8g %10.8g\n",
    //	    slip_scale, force[0], force[1], force[2]);
  }//Particle_number

  //  if(slip_flag){
  //    Solenoidal_u(up);
  //  }
  }
}

void Make_phi_particle(double *phi
		       ,Particle *p
		       ,const double radius
		       ){
  const int SW_UP = 0;
  double **dmy_up;
  int *nlattice;
  nlattice = Ns;
  Make_phi_u_primitive(phi, dmy_up, p, SW_UP,DX,NP_domain
		       ,Sekibun_cell
		       ,nlattice
		       ,radius);
}
void Make_phi_u_particle(double *phi
			 ,double **up
			 ,Particle *p
			 ){
  const int SW_UP = 1;
  int *nlattice;
  nlattice = Ns;
  Make_phi_u_primitive(phi, up, p, SW_UP,DX,NP_domain,Sekibun_cell,nlattice);
  }

void Make_phi_particle_OBL(double *phi
			   ,Particle *p
			   ,const double radius
    ){
    const int SW_UP = 0;
    double **dmy_up;
    int *nlattice;
    nlattice = Ns;
    Make_phi_u_primitive_OBL(phi, dmy_up, p, SW_UP,DX,NP_domain
			     ,Sekibun_cell
			     ,nlattice
			     ,radius);
}
void Make_phi_u_particle_OBL(double *phi
			     ,double **up
			     ,Particle *p
    ){
    const int SW_UP = 1;
    int *nlattice;
    nlattice = Ns;
    Make_phi_u_primitive_OBL(phi, up, p, SW_UP,DX,NP_domain,Sekibun_cell,nlattice);
}

void Make_phi_u_advection(double *phi, double **up, Particle *p){
  // map only V_p, excepting \Omega_p to the field up
  int *nlattice;
  if(SW_EQ == Shear_Navier_Stokes || SW_EQ == Shear_Navier_Stokes_Lees_Edwards){
    nlattice = Ns_shear;
  }else {
    nlattice = Ns;
  } 

    double xp[DIM],vp[DIM];
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell;
    int r_mesh[DIM];
    double r[DIM];
    double x[DIM];
    double dmy;
    double dmy_phi;
#pragma omp parallel for schedule(dynamic, 1) private(xp,vp,x_int,residue,sw_in_cell,r_mesh,r,x,dmy,dmy_phi)
  for(int n=0; n < Particle_Number; n++){
    for(int d=0;d<DIM;d++){
      xp[d] = p[n].x[d];
      vp[d] = p[n].v[d];
      {
	assert(p[n].x[d] >= 0);
	assert(p[n].x[d] < L[d]);
      }
    }
    
    sw_in_cell 
      = Particle_cell(xp, DX, x_int, residue);// {1,0} $B$,JV$C$F$/$k(B
    sw_in_cell = 1;
    for(int mesh=0; mesh < NP_domain; mesh++){
      Relative_coord(Sekibun_cell[mesh], x_int, residue, sw_in_cell, nlattice, DX, r_mesh, r);
      for(int d=0;d<DIM;d++){
	x[d] = r_mesh[d] * DX;
      }
      dmy = Distance(x, xp);
      dmy_phi= Phi(dmy);
      phi[(r_mesh[0]*NY*NZ_)+(r_mesh[1]*NZ_)+r_mesh[2]] += dmy_phi;
      
      for(int d=0;d<DIM;d++){
	up[d][(r_mesh[0]*NY*NZ_)+(r_mesh[1]*NZ_)+r_mesh[2]] += ( vp[d] * dmy_phi);
      }
    }
    p[n].eff_mass_ratio 
      = 1.0;
  }
}
