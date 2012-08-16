#ifndef AVS_PROCESS_H
#define AVS_PROCESS_H

#include "avs_utils.h"
// 
extern double ** u_avg;
extern double ** u_gold;
extern double ** u_phi_avg;
extern double ** u_phi_gold;
extern int *u_cnt;
extern int *u_phi_cnt;
extern int nskip;

void wrong_invocation(){
  fprintf(stderr, "usage: pavs -p PID -l Ns -i UDF [-err fr]\n");
  fprintf(stderr, "       -p   PID\t ID (1...N) of centered particle\n");
  fprintf(stderr, "       -i   UDF\t Input file\n");
  exit_job(EXIT_FAILURE);
}

inline void Cylindrical_coord(const double *x, 
                              double *n_h, double *n_s, double *n_phi,
                              double &h, double &s, double &phi, const JAX &ej){
  double r[DIM], dmy_r;
  dmy_r = sqrt(SQ(x[0]) + SQ(x[1]) + SQ(x[2]));
  for(int d = 0; d < DIM; d++){
    r[d] = x[d] / dmy_r;
  }

  double dot_e3_r = v_inner_prod(e3, r);
  h = dot_e3_r * dmy_r;
  for(int d = 0; d < DIM; d++){
    n_h[d] = e3[d];
  }

  if(less_than_mp(ABS(dot_e3_r), 1.0)){
    for(int d = 0; d < DIM; d++){
      n_s[d] = x[d] - h*n_h[d];
    }
    s = sqrt(SQ(n_s[0]) + SQ(n_s[1]) + SQ(n_s[2]));
    for(int d = 0; d < DIM; d++){
      n_s[d] /= s;
    }
    
    double dot_e1_s = v_inner_prod(e1, n_s);
    double dot_e2_s = v_inner_prod(e2, n_s);
    phi = acos(dot_e1_s);
    if(negative_mp(dot_e2_s) || (zero_mp(dot_e2_s) && negative_mp(dot_e1_s))){
      for(int d = 0; d < DIM; d++){
        n_s[d] *= -1.0;
      }
      phi += M_PI;
    }
    v_cross(n_phi, n_h, n_s);
  }else{
    s = 0.0;
    phi = 0.0;
    for(int d = 0; d < DIM; d++){
      n_s[d] = e2[d];
      n_phi[d] = -e1[d];
    }
  }
  for(int d = 0; d < DIM; d++){
    assert(n_h[d] == n_h[d]);
    assert(n_s[d] == n_s[d]);
    assert(n_phi[d] == n_phi[d]);
  }
  assert(h == h);
  assert(s == s);
  assert(phi == phi);
}

inline void smooth_v_2d(const int &j, const int &k, double const* const* u, double *avg){
  const int ndir = 5;
  int md[ndir][DIM-1] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};
  int e0[DIM-1] = {j,k};
  int ei[DIM-1];
  for(int d = 0; d < DIM; d++){
    avg[d] = 0.0;
  }
  for(int i = 0; i < ndir; i++){
    for(int d = 0; d < DIM -1; d++){
      ei[d] = e0[d] + md[i][d];
    }
    PBC_ip(ei[0], JNY);
    PBC_ip(ei[1], JNZ);

    int im = (ei[0] * JNZ) + ei[1];
    for(int d = 0; d < DIM; d++){
      avg[d] += u[d][im];
    }
  }

  for(int d = 0; d < DIM; d++){
    avg[d] /= ((double)ndir);
  }
}

inline void smooth_v(const int &i, const int &j, const int &k, double const* const* u,
             double *avg){
  const int ndir = 27;
  int md[ndir][DIM] ;
  int e0[DIM] = {i, j, k};
  int ei[DIM];
  int dir[DIM] = {-1, 0, 1};
  int dmy = 0;

  for(int d = 0; d < DIM; d++){
    avg[d] = 0.0;
    for(int l = 0; l < DIM; l++){
      for(int m = 0; m < DIM; m++){
        md[dmy][0] = dir[d];
        md[dmy][1] = dir[l];
        md[dmy][2] = dir[m];
        dmy++;
      }
    }
  }
  assert(dmy == ndir);

  for(int i = 0; i < ndir; i++){
    for(int d = 0; d < DIM; d++){
      ei[d] = e0[d] + md[i][d];
    }
    PBC_ip(ei);
    int im = (ei[0] * NY * NZ) + (ei[1] * NZ) + ei[2];
    for(int d = 0; d < DIM; d++){
      avg[d] += u[d][im];
    }
  }

  for(int d = 0; d < DIM; d++){
    avg[d] /= ((double) ndir);
  }
}

void v_gold(double ww[DIM], const double nr[DIM], const double &rdist, 
	    const double &B1, const double &B2, const double &RADIUS, const COORD_SYSTEM &coord){
  double dmy, cos_theta;
  double jax[DIM], nsin_theta[DIM];
  double nd1[DIM], nd2[DIM], nd3[DIM];
  double ar, ar2, ar3;

  if(coord == SPACE_FRAME){
    rigid_body_rotation(jax, e3, q0, BODY2SPACE);
  }else{
    for(int d = 0; d < DIM; d++){
      jax[d] = e3[d];
    }
  }
  
  if(rdist >= RADIUS){
    cos_theta = jax[0]*nr[0] + jax[1]*nr[1] + jax[2]*nr[2];
    ar = RADIUS/rdist;
    ar2 = ar*ar;
    ar3 = ar*ar2;

    for(int d = 0; d < DIM; d++){
      nsin_theta[d] = cos_theta * nr[d] - jax[d];
    }

    for(int d = 0; d < DIM; d++){
      nd1[d] = cos_theta * nr[d] - 1.0/3.0 * jax[d];
      nd2[d] = 2.0 * cos_theta * nsin_theta[d];
      nd3[d] = (3.0 * cos_theta * cos_theta - 1.0) * nr[d];
    }
    for(int d = 0; d < DIM; d++){
      ww[d] = B1 * ar3 * nd1[d] + 
	B2/2.0 * ar2 * ( ar2 * nd2[d] + (ar2 - 1.0) * nd3[d] );
    }
  }else{
    for(int d = 0; d < DIM; d++){
      ww[d] = jax[d] * 2./3. * B1;
    }
  }
}

void init_io(int argc, char *argv[], UDFManager* &udf_in, int &pid){
  
  char* fname;
  int pset, fset;
  pset = fset = 0;
  nskip = 0;
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
    }else if(strcmp(argv[i], "-s") == 0){
      if(i+1 == argc || argv[i+1][0] == '-'){
        fprintf(stderr, "Error: No skip steps specified\n");
        wrong_invocation();
      }
      nskip = atoi(argv[++i]);
      fprintf(stderr, "Skipping first %d frames...\n", nskip);
    }
  }
  if(!pset){
    fprintf(stderr, "Error: No particle selected\n");
    wrong_invocation();
  }
  if(!fset){
    fprintf(stderr, "Error: No udf file selected\n");
    wrong_invocation();
  }
}

void init_avg_mem(const double &RADIUS){
  
  u_avg = (double**) malloc(sizeof(double*) * DIM);
  u_phi_avg = (double**) malloc(sizeof(double*) * DIM);
  u_phi_gold = (double **) malloc(sizeof(double*) * DIM);
  u_gold = (double**) malloc(sizeof(double*) * DIM);
  for(int d = 0; d < DIM; d++){
    u_avg[d] = alloc_1d_double(NX*NY*NZ);
    u_phi_avg[d] = alloc_1d_double(JNY*JNZ);
    u_phi_gold[d] = alloc_1d_double(JNY*JNZ);
    u_gold[d] = alloc_1d_double(NX*NY*NZ);
  }
  u_cnt = alloc_1d_int(NX*NY*NZ);
  u_phi_cnt = alloc_1d_int(JNY*JNZ);
  fprintf(stderr, "Allocated: %d %d\n", DIM, NX*NY*NZ);

  for(int i = 0; i < NX; i++){
    for(int j = 0; j < NY; j++){
      for(int k = 0; k < NZ; k++){
        int im = (i * NY * NZ) + (j * NZ) + k;
        u_cnt[im] = 0;
        u_avg[0][im] = 0.0;
        u_avg[1][im] = 0.0;
        u_avg[2][im] = 0.0;
      }
    }
  }
  for(int j = 0; j < JNY; j++){
    for(int k = 0; k < JNZ; k++){
      int im = (j * JNZ) + k;
      u_phi_cnt[im] = 0;
      u_phi_avg[0][im] = 0.0;
      u_phi_avg[1][im] = 0.0;
      u_phi_avg[2][im] = 0.0;
      u_phi_gold[0][im] = 0.0;
      u_phi_gold[1][im] = 0.0;
      u_phi_gold[2][im] = 0.0;
    }
  }

  int ei[DIM], gi[DIM];
  double rp[DIM], ww[DIM], dmy_r;
  for(int ii = -HNX; ii <= HNX; ii++){
    ei[0] = ii;
    gi[0] = ei[0] + HNs[0];
    for(int jj = -HNY; jj <= HNY; jj++){
      ei[1] = jj;
      gi[1] = ei[1] + HNs[1];
      for(int kk = -HNZ; kk <= HNZ; kk++){
        ei[2] = kk;
        gi[2] = ei[2] + HNs[2];
        PBC_ip(gi);
        int im = (gi[0] * NY * NZ) + (gi[1] * NZ) + gi[2];

        dmy_r = 0;
        for(int d = 0; d < DIM; d++){
          rp[d] = (double)ei[d] * DX;
          dmy_r += rp[d]*rp[d];
        }
        dmy_r = sqrt(dmy_r);
        if(ii != 0 || jj != 0 || kk != 0){
          for(int d = 0; d < DIM; d++){
            rp[d] /= dmy_r;
          }
        }
        v_gold(ww, rp, dmy_r, B1, B2, RADIUS, BODY_FRAME);
        for(int d = 0; d < DIM; d++){
          u_gold[d][im] = ww[d];
        }
      }
    }
  }
}

void clean_exit(){
  for(int d = 0; d < DIM; d++){
    free_1d_double(u_avg[d]);
    free_1d_double(u_gold[d]);
    free_1d_double(u_phi_avg[d]);
    free_1d_double(u_phi_gold[d]);
  }
  free(u_avg);
  free(u_phi_avg);
  free(u_phi_gold);
  free(u_gold);
  free_1d_int(u_cnt);
}


#endif
