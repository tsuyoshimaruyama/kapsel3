#include "avs_process.h"
void process_avs_frame(const V_OPTION &vflag);
int main(int argc, char* argv[]){
  int pid;                // id of centered particle
  int *p_spec;            // species id for all particles
  JAX *sp_axis;           // janus axis for each species

  V_OPTION vflag;         // output flag for velocites
  UDFManager *udf_in;     // UDF input file
  
  initialize(argc, argv, udf_in, pid, vflag);
  get_system_data(udf_in, p_spec, sp_axis);
  initialize_avs();

  u = (double **)malloc(sizeof(double*)* DIM);
  post_u = (double **)malloc(sizeof(double*)*DIM);
  for(int d = 0; d < DIM; d++){
    u[d] = alloc_1d_double(NX*NY*NZ);
    post_u[d] = alloc_1d_double(CX*CY*CZ);
  }
  assert(pid < Nparticles);

  for(int i = 0; i <= Num_snap; i++){
    setup_avs_frame();
    read_u();
    read_p(pid);
    process_avs_frame(vflag);
    clear_avs_frame();
  }
  fclose(fluid_field);
  fclose(particle_field);
  fclose(post_field);
}


void init_interpol(int pid[8][DIM], double pval[8], double pcoeff[64], const double *field, const int*p0){
  const int ei[8][DIM] = 
    {{0, 0, 0}, 
     {1, 0, 0},
     {0, 1, 0},
     {1, 1, 0},
     {0, 0, 1},
     {1, 0, 1},
     {0, 1, 1},
     {1, 1, 1}};
  double pdx[8], pdy[8], pdz[8];
  double pdxdy[8], pdxdz[8], pdydz[8], pdxdydz[8];
  int im;

  //get function values at cube edges
  for(int l = 0; l < 8; l++){
    for(int d = 0; d < DIM; d++){
      pid[l][d] = p0[d] + ei[l][d];
    }
    PBC_ip(pid[l]);
    pval[l] = field[linear_id(pid[l][0], pid[l][1], pid[l][2])];
  }
  for(int l = 0; l < 64; l++){
    pcoeff[l] = 0.0;
  }

  //Use finite difference approximations for partial derivatives
  //Assuming unit cube (using scaled coordinates)
  if(SW_INTERPOL == cubic_interpol){
    for(int l = 0; l < 8; l++){
      pdx[l] = pdy[l] = pdz[l] = pdxdy[l] = pdxdz[l] = pdydz[l] = pdxdydz[l] = 0.0;
      int &x = pid[l][0];
      int &y = pid[l][1];
      int &z = pid[l][2];

      for(int i = -1; i <= 1; i = i + 2){
	int dmy_xi = x + i; PBC_ip(dmy_xi, NX);
	int dmy_yi = y + i; PBC_ip(dmy_yi, NY);
	int dmy_zi = z + i; PBC_ip(dmy_zi, NZ);

	pdx[l] += i * field[linear_id(dmy_xi, y, z)];
	pdy[l] += i * field[linear_id(x, dmy_yi, z)];
	pdz[l] += i * field[linear_id(x, y, dmy_zi)];

	for(int j = -1; j <= 1; j = j + 2){
	  int dmy_yj = y + j; PBC_ip(dmy_yj, NY);
	  int dmy_zj = z + j; PBC_ip(dmy_zj, NZ);

	  pdxdy[l] += i * j * field[linear_id(dmy_xi, dmy_yj, z)];
	  pdxdz[l] += i * j * field[linear_id(dmy_xi, y, dmy_zj)];
	  pdydz[l] += i * j * field[linear_id(x, dmy_yi, dmy_zj)];

	  for(int k = -1; k <= 1; k = k + 2){
	    int dmy_zk = z + k; PBC_ip(dmy_zk, NZ);

	    pdxdydz[l] += i * j * k * field[linear_id(dmy_xi, dmy_yj, dmy_zk)];
	  }//k
	}//j
      }//i
      pdx[l] /= 2.0;
      pdy[l] /= 2.0;
      pdz[l] /= 2.0;
      pdxdy[l] /= 4.0;
      pdxdz[l] /= 4.0;
      pdydz[l] /= 4.0;
      pdxdydz[l] /= 8.0;
    }//l
    tricubic_get_coeff(pcoeff, pval, pdx, pdy, pdz, pdxdy, pdxdz, pdydz, pdxdydz);
  }
}

inline double interpol(const int pid[8][DIM], const double pval[8], double pcoeff[64],
		       const double *r){
  //scaled distance from edge to r
  double x = r[0] / DX;
  double y = r[1] / DX;
  double z = r[2] / DX;
  assert(x >= 0 && y >= 0 && z >= 0);
  assert(x <= 1 && y <= 1 && z <= 1);
  return tricubic_eval(pcoeff, x, y, z);
}
inline void binary_write(double *a, FILE *fstream){
  for(int k = 0; k < CZ; k++){
    for(int j = 0; j < CY; j++){
      for(int i = 0; i < CX; i++){
	int im = linear_id(i, j, k, CX, CY, CZ);
	float dmy = (float)a[im];
	fwrite(&dmy, sizeof(float), 1, fstream);
      }
    }
  }
}
void process_avs_frame(const V_OPTION &vflag){
  double rp[DIM];
  double sp[DIM];
  double cp[DIM];
  int ei[DIM];
  int edge[DIM];
  int im;
  
  int interpol_pid[8][DIM];
  double interpol_pval[8];
  double interpol_coeff[64];
  double interu;

  for(int i = -HCX; i <= HCX; i++){
    ei[0] = i;
    for(int j = -HCY; j <= HCY; j++){
      ei[1] = j;
      for(int k = -HCZ; k <= HCZ; k++){
	ei[2] = k;

	im = ((i + HCX) * CY * CZ) + ((j + HCY) * CZ) + (k + HCZ);

	for(int d = 0; d < DIM; d++){
	  rp[d] = (double)ei[d] * DX;     //particle frame
	  sp[d] = r0[d] + rp[d];          //global frame
	  sp[d] = fmod(sp[d] + lbox[d], lbox[d]);
	  edge[d] = (int)(sp[d]/DX);
	  cp[d] = sp[d] - edge[d]*DX;

	  if(!(edge[d] >= 0 && edge[d] < Ns[d] && cp[d] >= 0)){
	    fprintf(stderr, "rp  : %f\n", rp[d]);
	    fprintf(stderr, "sp  : %f\n", sp[d]);
	    fprintf(stderr, "edge: %d\n", edge[d]);
	    fprintf(stderr, "cp  : %g\n", cp[d]);
	    assert(false);
	  }
	}
	for(int d = 0; d < DIM; d++){
	  init_interpol(interpol_pid, interpol_pval, interpol_coeff, u[d], edge);
	  post_u[d][im] = interpol(interpol_pid, interpol_pval, interpol_coeff, cp);
	  if(vflag == relative_velocity){
	    post_u[d][im] -= v0[d];
	  }
	}
      }
    }
  }
  
  binary_write(u[0], post_data);
  binary_write(u[1], post_data);
  binary_write(u[2], post_data);
}


