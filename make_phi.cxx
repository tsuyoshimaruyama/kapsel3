/*!
  \file make_phi.cxx
  \brief Routines to compute smooth profile and grid particle properties
  \author T. Iwashita
  \date 2006/10/20
  \version 1.3
 */

#include "make_phi.h"

void (*Angular2v)(const double *omega, const double *r, double *v);

int NP_domain;
int **Sekibun_cell;

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
      = Particle_cell(xp, DX, x_int, residue);// {1,0} が返ってくる
    sw_in_cell = 1;
    for(int mesh=0; mesh < NP_domain; mesh++){
      Relative_coord(Sekibun_cell[mesh]
		     ,x_int, residue, sw_in_cell, Ns, DX, r_mesh, r);
      dmy_r= sqrt(SQ(r[0])+SQ(r[1])+SQ(r[2]));
      dmy = ABS(dmy_r - RADIUS);
      if(dmy < HXI){
	  ir = 1./dmy_r;
          int im = (r_mesh[0]*NY*NZ_)+(r_mesh[1]*NZ_)+r_mesh[2];
#pragma omp atomic
	  surface_normal[0][im] += r[0]*ir;
#pragma omp atomic
	  surface_normal[1][im] += r[1]*ir;
#pragma omp atomic
	  surface_normal[2][im] += r[2]*ir;
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
#pragma omp parallel for schedule(dynamic, 1) private(drho,xp,x_int,residue,sw_in_cell,r_mesh,r,x,dmy,dmy_phi)
    for(int n=0; n < Particle_Number; n++){
	drho = RHO_particle[p[n].spec] - RHO;
	for(int d=0;d<DIM;d++){
	    xp[d] = p[n].x[d];
	}
	
	sw_in_cell 
	    = Particle_cell(xp, dx, x_int, residue);// {1,0} が返ってくる
	sw_in_cell = 1;
	for(int mesh=0; mesh < np_domain; mesh++){
	    Relative_coord(sekibun_cell[mesh]
			   , x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
	    for(int d=0;d<DIM;d++){
		x[d] = r_mesh[d] * dx;
	    }
	    dmy=Distance(x,xp);
	    dmy_phi=Phi(dmy)*drho;
#pragma omp atomic
	    phi[(r_mesh[0]*NY*NZ_)+(r_mesh[1]*NZ_)+r_mesh[2]] += dmy_phi;
	}
    }

    int im;
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

#pragma omp parallel for schedule(dynamic, 1) \
  private(xp,vp,omega_p,x_int,residue,sw_in_cell,r_mesh,r,x,dmy,dmy_phi,v_rot,im)
    for(int n=0; n < Particle_Number; n++){
	for(int d=0;d<DIM;d++){
	    xp[d] = p[n].x[d];
	    vp[d] = p[n].v[d];
	    omega_p[d] = p[n].omega[d];
	}
	
	sw_in_cell 
	    = Particle_cell(xp, dx, x_int, residue);// {1,0} が返ってくる
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
#pragma omp atomic
	    phi[im] += dmy_phi;
	    
	    if(SW_UP){
		Angular2v(omega_p, r, v_rot);
#pragma omp atomic
                up[0][im] += ( (vp[0]+v_rot[0]) * dmy_phi);
#pragma omp atomic
                up[1][im] += ( (vp[1]+v_rot[1]) * dmy_phi);
#pragma omp atomic
                up[2][im] += ( (vp[2]+v_rot[2]) * dmy_phi);
	    }
	}
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

inline void Make_phi_particle_sum_primitive(double *phi,
                                            double *phi_sum,
                                            Particle *p,
                                            const double &dx,
                                            const int &np_domain,
                                            int **sekibun_cell,
                                            const int Nlattice[DIM],
                                            const double radius){
#pragma omp parallel for schedule(dynamic, 1)
  for(int n = 0; n < Particle_Number; n++){
    double xp[DIM];
    for(int d = 0; d < DIM; d++) xp[d] = p[n].x[d];

    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell = Particle_cell(xp, dx, x_int, residue);
    sw_in_cell = 1;
    
    int im;
    int r_mesh[DIM];
    double dmy, dmy_phi;
    double r[DIM], x[DIM];
    for(int mesh = 0; mesh < np_domain; mesh++){
      Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);

      for(int d = 0; d < DIM; d++) x[d] = r_mesh[d]*DX;

      dmy     = Distance(x, xp);
      dmy_phi = Phi(dmy, radius);

#pragma omp atomic
      phi_sum[(r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2]] += dmy_phi;
    }
  }

  {
    int im;
#pragma omp parallel for schedule(dynamic, 1)
    for(int i = 0; i < NX; i++){
      for(int j = 0; j < NY; j++){
        for(int k = 0; k < NZ; k++){
          im = (i * NY * NZ_) + (j * NZ_) + k;
          phi[im] = MIN(phi_sum[im], 1.0);
        }
      }
    }
  }
}

inline void Make_u_particle_sum_primitive(double **up,
                                          double const* phi_sum,
                                          Particle *p,
                                          const double &dx,
                                          const int &np_domain,
                                          int const* const* sekibun_cell,
                                          const int Nlattice[DIM],
                                          const double radius){
#pragma omp parallel for schedule(dynamic, 1)
  for(int n = 0; n < Particle_Number; n++){
    double xp[DIM], vp[DIM], omega_p[DIM];
    for(int d = 0; d < DIM; d++){
      xp[d] = p[n].x[d];
      vp[d] = p[n].v[d];
      omega_p[d] = p[n].omega[d];
    }

    int im, sw_in_cell;
    int x_int[DIM], r_mesh[DIM];
    double residue[DIM], r[DIM], x[DIM], v_rot[DIM];
    double dmy, dmy_phi;
    sw_in_cell = Particle_cell(xp, dx, x_int, residue);
    sw_in_cell = 1;
    for(int mesh = 0; mesh < np_domain; mesh++){
      Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
      for(int d = 0; d < DIM; d++) x[d] = r_mesh[d]*DX;

      im = (r_mesh[0]*NY*NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];      

      dmy = Distance(x, xp);
      dmy_phi = (phi_sum[im] <= 1.0 ? Phi(dmy, radius) : Phi(dmy,radius)/phi_sum[im]);

      Angular2v(omega_p, r, v_rot);
#pragma omp atomic
      up[0][im] += ((vp[0] + v_rot[0]) * dmy_phi);
#pragma omp atomic
      up[1][im] += ((vp[1] + v_rot[1]) * dmy_phi);
#pragma omp atomic
      up[2][im] += ((vp[2] + v_rot[2]) * dmy_phi);
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
#pragma omp parallel for schedule(dynamic, 1) private(xp,vp,omega_p,x_int,residue,sw_in_cell,r_mesh,r,x,dmy,dmy_phi,v_rot,sign,im)
    for(int n=0; n < Particle_Number; n++){
	for(int d=0;d<DIM;d++){
	    xp[d] = p[n].x[d];
	    vp[d] = p[n].v[d];
	    omega_p[d] = p[n].omega[d];
	}
	
	sw_in_cell 
	    = Particle_cell(xp, dx, x_int, residue);// {1,0} が返ってくる
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
#pragma omp atomic
	    phi[im] += dmy_phi;
	    
	    if(SW_UP){
		Angular2v(omega_p, r, v_rot);
#pragma omp atomic
                up[0][im] += ( (vp[0] - sign*Shear_rate_eff*L_particle[1] + v_rot[0])*dmy_phi);
#pragma omp atomic
                up[1][im] += (vp[1] + v_rot[1])*dmy_phi;
#pragma omp atomic
                up[2][im] += (vp[2] + v_rot[2])*dmy_phi;
	    }
	}
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
void Make_phi_particle_sum(double *phi, double* phi_sum, Particle *p, const double radius){
  int *nlattice;
  nlattice = Ns;
  Make_phi_particle_sum_primitive(phi, phi_sum, p, DX, NP_domain, Sekibun_cell, nlattice, radius);
}
void Make_u_particle_sum(double **up,  double const* phi_sum,  Particle *p, const double radius){
  int *nlattice;
  nlattice = Ns;
  Make_u_particle_sum_primitive(up, phi_sum, p, DX, NP_domain, Sekibun_cell, nlattice, radius);
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
      = Particle_cell(xp, DX, x_int, residue);// {1,0} が返ってくる
    sw_in_cell = 1;
    for(int mesh=0; mesh < NP_domain; mesh++){
      Relative_coord(Sekibun_cell[mesh], x_int, residue, sw_in_cell, nlattice, DX, r_mesh, r);
      for(int d=0;d<DIM;d++){
	x[d] = r_mesh[d] * DX;
      }
      dmy = Distance(x, xp);
      dmy_phi= Phi(dmy);
      int im = (r_mesh[0]*NY*NZ_)+(r_mesh[1]*NZ_)+r_mesh[2];
#pragma omp atomic
      phi[im] += dmy_phi;
#pragma omp atomic
      up[0][im] += ( vp[0] * dmy_phi);
#pragma omp atomic
      up[1][im] += ( vp[1] * dmy_phi);
#pragma omp atomic
      up[2][im] += ( vp[2] * dmy_phi);
    }
  }
}
