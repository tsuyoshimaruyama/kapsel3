//
// $Id: make_phi.h,v 1.3 2006/07/28 17:14:55 nakayama Exp $
//
#ifndef MAKE_PHI_H
#define MAKE_PHI_H

#include <assert.h> 
#include "input.h"
#include "variable.h"
#include "profile.h"
#include "rigid_body.h"
#include "operate_omega.h"

extern void (*Angular2v)(const double *omega, const double *r, double *v);
extern int NP_domain;
extern int NP_domain_interface;
extern int **Sekibun_cell;
extern int **Sekibun_cell_interface;

void Make_phi_u_advection(double *phi, double **up, Particle *p);

void Make_phi_janus_particle(double *phi,
			     double *id_phi,
			     Particle *p);
void Make_phi_janus_particle_OBL(double *phi,
				 double *id_phi,
				 Particle *p);

void Make_phi_particle(double *phi
		       ,Particle *p
		       ,const double radius = RADIUS
		       );
void Make_phi_u_particle(double *phi, double **up, Particle *p);
void Make_phi_particle_OBL(double *phi
			   ,Particle *p
			   ,const double radius = RADIUS
    );
void Make_phi_u_particle_OBL(double *phi, double **up, Particle *p);
void Make_surface_normal(double **surface_normal
			   ,const Particle *p);
void Make_rho_field(double *phi,Particle *p);

inline int Particle_cell(const double *xp
			 ,const double &dx
			 ,int *x_int
			 ,double *residue){
  const double idx = 1./dx;
  {
    assert(xp[0] >= 0);
    assert(xp[0] < L_particle[0]);

    assert(xp[1] >= 0);
    assert(xp[1] < L_particle[1]);

    assert(xp[2] >= 0);
    assert(xp[2] < L_particle[2]);
  }
  
  // r_i の左下のメッシュ座標を計算
  int sw_in_cell=0;
  for(int d=0; d<DIM; d++){
    double dmy = (xp[d] *idx);
    x_int[d] = (int)dmy;

    residue[d] = (dmy - (double)x_int[d]) * dx;
    if( dmy - (double)x_int[d] > 0.0 ){
      sw_in_cell = 1;
    }
  }
  return sw_in_cell;
}
inline void Relative_coord(const int *cell
			   ,const int *x_int
			   ,const double *residue
			   ,const int &sw_in_cell
			   ,const int Nlattice[DIM]
			   ,const double dx
			   ,int *r_mesh
			   ,double *r){

    for(int d=0;d<DIM;d++){
	r_mesh[d] = (x_int[d] + cell[d] + Nlattice[d] ) % Nlattice[d]; 
	{
	    assert( r_mesh[d] < Nlattice[d] );
	    assert( r_mesh[d] >= 0 );
	}
	r[d] = (double)cell[d] * dx - residue[d]; 
    }
}
inline void Relative_coord_OBL(const int *cell
			       ,const int *x_int
			       ,const double *residue
			       ,const int &sw_in_cell
			       ,const int Nlattice[DIM]
			       ,const double dx
			       ,int *r_mesh
			       ,double *r){
 
    double signY = x_int[1] + cell[1];
    r_mesh[1] = (x_int[1] + cell[1] + Nlattice[1]) % Nlattice[1];
    signY -= r_mesh[1];
    int sign = (int)signY;
    if (!(sign == 0)) {
	sign  = sign/abs(sign);
    }

    r_mesh[0] = (int)fmod(x_int[0] 
			  - sign*degree_oblique*Nlattice[1] + residue[0]/dx +
			  cell[0] + 4*Nlattice[0], Nlattice[0]);
    double x_oblique_residue =
	fmod(x_int[0] - 
	     sign*degree_oblique*Nlattice[1] + residue[0]/dx +
	     cell[0] + 4*Nlattice[0], Nlattice[0]) - r_mesh[0];
    
    r_mesh[2] = (x_int[2] + cell[2] + Nlattice[2]) % Nlattice[2];   

    for(int d=0;d<DIM;d++){
	{
	    assert( r_mesh[d] < Nlattice[d] );
	    assert( r_mesh[d] >= 0 );
	}
    }
    r[0] = (double)cell[0] * dx - x_oblique_residue*dx;
    r[1] = (double)cell[1] * dx - residue[1];
    r[2] = (double)cell[2] * dx - residue[2];   
}
inline int Relative_coord_check_stepover_Y(const int *cell
					   ,const int *x_int
					   ,const double *residue
					   ,const int &sw_in_cell
					   ,const int Nlattice[DIM]
					       ,const double dx
					   ,int *r_mesh
					   ,double *r){
    
    double signY = x_int[1] + cell[1];
    r_mesh[1] = (x_int[1] + cell[1] + Nlattice[1]) % Nlattice[1];
    signY -= r_mesh[1];
    int sign = (int)signY;
    if (!(sign == 0)) {
	sign  = sign/abs(sign);
    }

    r_mesh[0] = (int)fmod(x_int[0] 
			  - sign*degree_oblique*Nlattice[1] + residue[0]/dx +
			  cell[0] + 4*Nlattice[0], Nlattice[0]);
    double x_oblique_residue =
	fmod(x_int[0] - 
	     sign*degree_oblique*Nlattice[1] + residue[0]/dx +
	     cell[0] + 4*Nlattice[0], Nlattice[0]) - r_mesh[0];
    
    r_mesh[2] = (x_int[2] + cell[2] + Nlattice[2]) % Nlattice[2];   

    for(int d=0;d<DIM;d++){
	{
	    assert( r_mesh[d] < Nlattice[d] );
	    assert( r_mesh[d] >= 0 );
	}
    }
    r[0] = (double)cell[0] * dx - x_oblique_residue*dx;
    r[1] = (double)cell[1] * dx - residue[1];
    r[2] = (double)cell[2] * dx - residue[2];

    return sign;
}


inline void Janus_direction(double *polar_axis, const Particle &p){
  double ex[DIM] = {1.0, 0.0, 0.0};
  double ey[DIM] = {0.0, 1.0, 0.0};
  double ez[DIM] = {0.0, 0.0, 1.0};
  double *e3;
  int dmy_axis = janus_axis[p.spec];

  if(dmy_axis == z_axis){
    e3 = ez;
  }else if(dmy_axis == y_axis){
    e3 = ey;
  }else if(dmy_axis == x_axis){
    e3 = ex;
  }else{
    fprintf(stderr, "Error: Invalid Janus axis for Janus_direction");
    exit_job(EXIT_FAILURE);
  }
  rigid_body_rotation(polar_axis, e3, p.q, BODY2SPACE);
}
inline void Spherical_coord(const double *x, double *r, double *theta, double *phi, 
			    double &norm_r, double &theta_angle, double &phi_angle, const Particle &p){
  double ex[DIM] = {1.0, 0.0, 0.0};
  double ey[DIM] = {0.0, 1.0, 0.0};
  double ez[DIM] = {0.0, 0.0, 1.0};
  double *e1, *e2, *e3;
  double r12[DIM];
  double dot_e3_r, dot_e1_r12, dmy_norm;

  // determine janus (e3) axis
  int dmy_axis = janus_axis[p.spec];
  if(dmy_axis == z_axis){
    e1 = ex;
    e2 = ey;
    e3 = ez;
  }else if(dmy_axis == y_axis){
    e1 = ez;
    e2 = ex;
    e3 = ey;
  }else if(dmy_axis == x_axis){
    e1 = ey;
    e2 = ez;
    e3 = ex;
  }else{
    fprintf(stderr, "Error: Invalid Janus axis for Spherical_coord\n");
    exit_job(EXIT_FAILURE);
  }

  //space -> janus body coordinates
  rigid_body_rotation(r, x, p.q, SPACE2BODY);

  //r normal vector
  norm_r = sqrt( SQ(r[0]) + SQ(r[1]) + SQ(r[2]) );
  assert(norm_r > 0.);
  dmy_norm = 1.0/norm_r;
  for(int d = 0; d < DIM; d++){
    r[d] *= dmy_norm;
  }

  // phi normal vector
  v_cross(phi, e3, r);  
  dmy_norm = sqrt(SQ(phi[0]) + SQ(phi[1]) + SQ(phi[2]));
  if(dmy_norm > 0.){// r NOT parallel to janus axis
    dmy_norm = 1.0/dmy_norm;
    for(int d = 0; d < DIM; d++){
      phi[d] *= dmy_norm;
    }
    
    //theta normal vector
    v_cross(theta, phi, r);
    dmy_norm = 1.0/sqrt(SQ(theta[0]) + SQ(theta[1]) + SQ(theta[2]));
    for(int d = 0; d < DIM; d++){
      theta[d] *= dmy_norm;
    }
    
    //theta angle
    dot_e3_r = v_inner_prod(e3, r);
    theta_angle = acos(dot_e3_r);
    
    //phi angle
    dmy_norm = 0.0;
    for(int d = 0; d < DIM; d++){
      r12[d] = r[d] - e3[d] * dot_e3_r;
      dmy_norm += r12[d]*r12[d];
    }
    dot_e1_r12 = v_inner_prod(e1, r12) / sqrt(dmy_norm);
    phi_angle = acos(dot_e1_r12);
  }else{ // r parallel to janus axis (phi / theta not uniquely defined)
    // define such that surface integral of theta is parallel to z (of phi is null)
    if(e3[0]*r[0] + e3[1]*r[1] + e3[2]*r[2] > 0){ // parallel
      for(int d = 0; d < DIM; d++){
	phi[d] = -e1[d];
	theta[d] = e2[d];
      }
      theta_angle = 0.0;
      phi_angle = PI_half;
    }else{ // anti-parallel (mirror image - left-handed frame !)
      for(int d = 0; d < DIM; d++){
	phi[d] = e1[d];
	theta[d] = -e2[d];
      }
      theta_angle = M_PI;
      phi_angle = 3.0*PI_half;
    }
  }
    
  //body -> space coordinates
  rigid_body_rotation(r, p.q, BODY2SPACE);
  rigid_body_rotation(theta, p.q, BODY2SPACE);
  rigid_body_rotation(phi, p.q, BODY2SPACE);
}


inline void Angular2v_rot_off(const double *omega, const double *r, double *v){
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
}
inline void Angular2v_rot_on(const double *omega, const double *r, double *v){
    v[0] = omega[1] * r[2] - omega[2] * r[1];
    v[1] = omega[2] * r[0] - omega[0] * r[2];
    v[2] = omega[0] * r[1] - omega[1] * r[0];
}
inline void Reset_phi_primitive(double *phi
				,const int nx
				,const int ny
				,const int nz
				,const double &value
    ){
    int im;
#pragma omp parallel for schedule(dynamic, 1) private(im)
    for(int i=0;i<nx;i++){
	for(int j=0;j<ny;j++){
	    for(int k=0;k<nz;k++){
		im =(i*NY*NZ_)+(j*NZ_)+k; 
		phi[im] = value;
	    }
	}
    }
}
inline void Reset_phi(double *phi, const double value = 0.0){
  Reset_phi_primitive(phi, NX, NY, NZ_, value);
}
inline void Reset_phi_u(double *phi, double **up){
    int im;
//#pragma omp parallel for schedule(dynamic, 8) private(im)
#pragma omp parallel private(im)
    {
#pragma omp for nowait schedule(dynamic, 1)
	for(int i=0;i<NX;i++){
	    for(int j=0;j<NY;j++){
		for(int k=0;k<NZ_;k++){
		    im =(i*NY*NZ_)+(j*NZ_)+k; 
		    phi[im] = 0.; 
		}
	    }
	}
#pragma omp for nowait schedule(dynamic, 1)
	for(int i=0;i<NX;i++){
	    for(int j=0;j<NY;j++){
		for(int k=0;k<NZ_;k++){
		    im =(i*NY*NZ_)+(j*NZ_)+k; 
		    up[0][im] = 0.;
		}
	    }
	}
#pragma omp for nowait schedule(dynamic, 1)
	for(int i=0;i<NX;i++){
	    for(int j=0;j<NY;j++){
		for(int k=0;k<NZ_;k++){
		    im =(i*NY*NZ_)+(j*NZ_)+k; 
		    up[1][im] = 0.;
		}
	    }
	}
#pragma omp for schedule(dynamic, 1)
	for(int i=0;i<NX;i++){
	    for(int j=0;j<NY;j++){
		for(int k=0;k<NZ_;k++){
		    im =(i*NY*NZ_)+(j*NZ_)+k; 
		    up[2][im] = 0.;
		}
	    }
	}
    }
}
inline void Reset_u(double **up){
  int im;
#pragma omp parallel private(im)
  {
#pragma omp for nowait schedule(dynamic, 1)
    for(int i = 0; i < NX; i++){
      for(int j = 0; j < NY; j++){
	for(int k = 0; k < NZ_; k++){
	  im = (i * NY * NZ_) + (j * NZ_) + k;
	  up[0][im] = 0.0;
	}
      }
    }/* end omp for up[0] */

#pragma omp for nowait schedule(dynamic, 1)
    for(int i = 0; i < NX; i++){
      for(int j = 0; j < NY; j++){
	for(int k = 0; k < NZ_; k++){
	  im = (i * NY * NZ_) + (j * NZ_) + k;
	  up[1][im] = 0.0;
	}
      }
    }
  }/* end omp for up[1] */

#pragma omp for nowait schedule(dynamic, 1)
  for(int i = 0; i < NX; i++){
    for(int j = 0; j < NY; j++){
      for(int k = 0; k < NZ_; k++){
	im = (i * NY * NZ_) + (j * NZ_) + k;
	up[2][im] = 0.0;
      } 
    }
  }/* end omp for up[2] */
}

#endif
