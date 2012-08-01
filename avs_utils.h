/*
  Read all avs data
 */
#ifndef AVS_UTILS_H
#define AVS_UTILS_H

#include <string>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <assert.h>
#include "alloc.h"
#include "macro.h"
#include "parameter_define.h"
#include "udfmanager.h"
#include "rigid_body.h"
#include "quaternion.h"

enum JAX{x_axis, y_axis, z_axis};
extern double ex[DIM]; //
extern double ey[DIM]; //
extern double ez[DIM]; //
extern double *e1; //
extern double *e2; //
extern double *e3; //

//GLOBAL DATA
const int cbuffer = 128;
extern char AVS_dir[cbuffer];   // 
extern char Out_dir[cbuffer];   // 
extern char Out_name[cbuffer];  // 
extern char Outp_dir[cbuffer];  // 
extern char Outp_name[cbuffer]; // 

extern FILE *fluid_field;
extern FILE *particle_field;
extern FILE *fluid_cod;
extern FILE *particle_cod;

extern int fluid_veclen;
extern int particle_veclen;


extern int GTS;  //
extern int Num_snap; //

extern double A; //
extern double A_XI; //
extern double DX; //

extern int Ns[DIM]; //
extern int &NX; //
extern int &NY; //
extern int &NZ; //
extern int lbox[DIM]; //

extern int Nspec; //
extern int Nparticles; //
extern int *p_spec; //
extern double *p_slip; //
extern double *p_slipmode; //
extern JAX *p_axis; //

// Per Frame Data
//fluid
extern FILE *fluid_data;
extern double **u; //
//particle
extern FILE *particle_data;
extern double UP; // 
extern double B1; //
extern double B2; //
extern double B1_app; //
extern double r0[DIM];  //
extern int r0_int[DIM]; //
extern double res0[DIM];//
extern double v0[DIM]; //
extern double w0[DIM]; //
extern double f0[DIM]; //
extern double t0[DIM]; //
extern double Q0[DIM][DIM]; //
extern double n0[DIM];
extern quaternion q0; //


// Get grid point closest to particle
inline void Particle_cell(const double *xp, const double &dx,
			 int *x_int,  double *residue){
  const double idx = 1.0/dx;
  {
    for(int d = 0; d < DIM; d++){
      assert(xp[d] >= 0 && xp[d] < lbox[0]);
    }
    for(int d = 0; d < DIM; d++){
      double dmy = (xp[d] * idx);
      x_int[d] = (int) dmy;

      residue[d] = (dmy - (double)x_int[d]) * dx;
    }
  }
}

// Profile functions
inline double H(const double x){
  return x > 0 ? exp(-SQ(DX/x)) : 0;
}
inline double Phi(const double &x, const double radius){
  double HXI = A_XI / 2.0;
  double dmy = H(radius + HXI - x);
  return dmy / (dmy + H(x - radius + HXI));
}

// Aplly periodic boundary conditions to grid coordinates
inline void PBC_ip(int &ip, const int &Ns){
  ip = (ip + Ns) % Ns;
  assert(ip >= 0 && ip < Ns);
}
inline void PBC_ip(int *ip){
  PBC_ip(ip[0], NX);
  PBC_ip(ip[1], NY);
  PBC_ip(ip[2], NZ);
}
inline void PBC_ip(int *ip, const int &NX, const int &NY, const int &NZ){
  PBC_ip(ip[0], NX);
  PBC_ip(ip[1], NY);
  PBC_ip(ip[2], NZ);
}

// Get (linear) array position for multi-dimensional arrays
inline int linear_id(const int &i, const int &j, const int &k){
  return (i * NY * NZ) + (j * NZ) + k;
}
inline int linear_id(const int &i, const int &j, const int &k, 
		     const int &NX, const int &NY, const int &NZ){
  return (i * NY * NZ) + (j * NZ) + k;
}

// PBC Distance function
inline void Distance(const double *x1, const double *x2,
		     double &r12, double *x12){
  double dmy = 0.0;
  for(int d = 0; d < DIM; d++){
    x12[d] = x2[d] - x1[d];
    x12[d] -= (double)Nint(x12[d]/lbox[d]) * lbox[d];
    dmy += SQ(x12[d]);    
  }
  r12 = sqrt(dmy);
  if(dmy > 0){
    for(int d = 0; d < DIM; d++){
      x12[d] /= r12;
    }
  }
}

// read system data from udf file used to launch KAPSEL
void get_system_data(UDFManager *ufin);

// read headers from avs *.fld files
void initialize_avs();

// close *.fld files
void close_avs();

// restart avs
void reset_avs();

// initialize new avs frame from *.fld files
void setup_avs_frame();

// clar frame data
void clear_avs_frame();

// initialize particle data
void initialize_avs_p(const int &pid);

// intialize avs mem
void initialize_avs_mem();

// read fluid velocity data
void read_avs_u();

// read data for pid particle
void read_avs_p(const int &pid);



#endif
