#ifndef AVS_PROCESS_H
#define AVS_PROCESS_H
#include <string.h>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <assert.h>
#include "alloc.h"
#include "udfmanager.h"
#include "macro.h"
#include "parameter_define.h"
#include "rigid_body.h"
#include "quaternion.h"

extern int Cs[DIM];
extern int &CX; //
extern int &CY; //
extern int &CZ; //
extern int HCs[DIM];
extern int &HCX; //
extern int &HCY; //
extern int &HCZ; //




void init_io(int argc, char *argv[], UDFManager* &udf_in, int &pid){
  
  char* fname;
  int pset, lset, fset, vset;
  pset = lset = fset = vset = 0;
  for(int i = 1; i < argc; i++){
    if(strcmp(argv[i], "-p") == 0){
      if(i+1 == argc || argv[i+1][0] == '-'){
	fprintf(stderr, "Error: No particle selected\n");
	wrong_invocation();
      }
      pid = atoi(argv[++i]);
      if(--pid < 0){
	fprintf(stderr, "Incorrect id\n");
	wrong_invocation();
      }
      pset = 1;
    }else if(strcmp(argv[i],"-l") == 0){
      if(i+1 == argc || argv[i+1][0] == '-'){
	fprintf(stderr, "Error: No box size selected\n");
	wrong_invocation();
      }
      Cs[0] = atoi(argv[++i]);
      if(Cs[0] <= 0){
	fprintf(stderr, "Incorrect size\n");
	wrong_invocation();
      }      
      if(Cs[0]%2 == 0){
	Cs[0]++;
      }
      Cs[1] = Cs[2] = Cs[0];
      HCs[0] = HCs[1] = HCs[2] = Cs[0] / 2;
      lset = 1;
    }else if(strcmp(argv[i], "-i") == 0){
      if(i+1 == argc || argv[i+1][0] == '-'){
	fprintf(stderr, "Error: No udf file selected\n");
	wrong_invocation();
      }
      i++;
      fname = argv[i];
      if(file_check(argv[i])){
	udf_in = new UDFManager(argv[i]);
      }else{
	fprintf(stderr, "Error: Cannot open file\n");
      }
      fset = 1;
    }
  }

  if(!pset){
    fprintf(stderr, "Error: No particle selected\n");
    wrong_invocation();
  }
  if(!lset){
    fprintf(stderr, "Error: No box size selected\n");
    wrong_invocation();
  }
  if(!fset){
    fprintf(stderr, "Error: No udf file selected\n");
    wrong_invocation();
  }

  fprintf(stderr, " # Centered Particle: %8d\n", pid + 1);
  fprintf(stderr, " # Box Size         : (%d,%d,%d)\n", 
	  Cs[0], Cs[1], Cs[2]);
  fprintf(stderr, " # HBox Size        : (%d,%d,%d)\n", 
	  HCs[0], HCs[1], HCs[2]);
  fprintf(stderr, " # UDF File         : %s\n", fname);

}

inline void read_ux(double *ux, FILE *fstream){
  float dmy;
  for(int k = 0; k < NZ; k++){
    for(int j = 0; j < NY; j++){
      for(int i = 0; i < NX; i++){
	int im = linear_id(i, j, k);
	fread(&dmy, sizeof(float), 1, fstream);
	ux[im] = (double)dmy;
      }
    }
  }
}

inline void read_pid(const int &pid, double x[DIM], FILE *fstream){
  float dmy;
  for(int d = 0; d < DIM; d++){
    for(int i = 0; i < Nparticles; i++){
      fread(&dmy, sizeof(float), 1, fstream);
      if(i == pid){
	x[DIM - 1 - d] = (double)dmy;
      }
    }
  }
}
inline void read_pid(const int &pid, double &a, FILE *fstream){
  float dmy;
  for(int i = 0; i < Nparticles; i++){
    fread(&dmy, sizeof(float), 1, fstream);
    if(i == pid){
      a = (double)dmy;
    }
  }
}
inline void read_pid(const int &pid, double Q[DIM][DIM], FILE *fstream){
  float dmy;
  for(int m = 0; m < DIM; m++){
    for(int n = 0; n < DIM; n++){
      for(int i = 0; i < Nparticles; i++){
	fread(&dmy, sizeof(float), 1, fstream);
	if(i == pid){
	  Q[n][m] = (double) dmy;
	}
      }
    }
  }
}
inline void read_u(){
  read_ux(u[0], fluid_data);
  read_ux(u[1], fluid_data);
  read_ux(u[2], fluid_data);
}
inline void read_p(const int &pid){
  double A0;
  float dmy;

  read_pid(pid, r0, particle_cod);
  read_pid(pid, A0, particle_data);
  read_pid(pid, v0, particle_data);
  read_pid(pid, w0, particle_data);
  read_pid(pid, f0, particle_data);
  read_pid(pid, t0, particle_data);
  read_pid(pid, Q0, particle_data);
  rm_rqtn(q0, Q0);

  double jax[DIM];
  rigid_body_rotation(jax, e3, q0, BODY2SPACE);
  B1_app = 3.0/2.0 * (jax[0]*v0[0] + jax[1]*v0[1] + jax[2]*v0[2]);

  Particle_cell(r0, DX, r0_int, res0);
}


void wrong_invocation(){
  fprintf(stderr, "usage: pavs -p PID -l Ns -i UDF\n");
  fprintf(stderr, "       -p PID\t ID (1...N) of centered particle\n");
  fprintf(stderr, "       -l Ns \t Box size length\n");
  fprintf(stderr, "       -i UDF\t Input file\n");
  exit_job(EXIT_FAILURE);
}

#endif
