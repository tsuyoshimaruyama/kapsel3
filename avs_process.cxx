#include "avs_process.h"
void process_avs_frame(const V_OPTION &vflag, const int &frame);
void v_gold(double ww[DIM], const double nr[DIM], const double &rdist, 
	    const double &B1, const double &B2);
void get_image(const int &frame, const int &iso);
void get_image_rotate(const int &frame, const int &npi);
int main(int argc, char* argv[]){
  int pid;                // id of centered particle
  int *p_spec;            // species id for all particles
  JAX *sp_axis;           // janus axis for each species
  double *sp_slip;        // slip velocity for each species
  double *sp_slipmode;    // B2 coefficient

  V_OPTION vflag;         // output flag for velocites
  UDFManager *udf_in;     // UDF input file
  
  initialize(argc, argv, udf_in, pid, vflag);
  get_system_data(udf_in, p_spec, sp_axis, sp_slip, sp_slipmode);
  initialize_avs();

  u = (double **)malloc(sizeof(double*)* DIM);
  post_u = (double **)malloc(sizeof(double*)*DIM);
  for(int d = 0; d < DIM; d++){
    u[d] = alloc_1d_double(NX*NY*NZ);
    post_u[d] = alloc_1d_double(CX*CY*CZ);
  }
  assert(pid < Nparticles);
  //janus axis of centered particle
  if(sp_axis[p_spec[pid]] == x_axis){
    e3 = ex;
    e1 = ey;
    e2 = ez;
  }else if(sp_axis[p_spec[pid]] == y_axis){
    e3 = ey;
    e1 = ez;
    e2 = ex;
  }else if(sp_axis[p_spec[pid]] == z_axis){
    e3 = ez;
    e1 = ex;
    e2 = ey;
  }
  //slip velocity of centered particle
  B1_real = sp_slip[p_spec[pid]];
  B2 = sp_slipmode[p_spec[pid]];

  for(int i = 0; i <= Num_snap; i++){
    setup_avs_frame();
    read_u();
    read_p(pid);
    process_avs_frame(vflag, i);
    get_image(i, 0);
    get_image(i, 1);
    get_image(i, 2);
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
void process_avs_frame(const V_OPTION &vflag, const int &frame){
  double rp[DIM], sp[DIM], cp[DIM];
  int ei[DIM], edge[DIM];
  int im;
  
  int interpol_pid[8][DIM];
  double interpol_pval[8];
  double interpol_coeff[64];

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
  double jax[DIM];
  rigid_body_rotation(jax, e3, q0, BODY2SPACE);
  fprintf(stderr, "# %d: B= %7.5g %7.5g %7.5g\n", 
	  frame, 
	  B1_app, B1_real, B2);
  
  binary_write(u[0], post_data);
  binary_write(u[1], post_data);
  binary_write(u[2], post_data);
}
void v_gold(double ww[DIM], const double nr[DIM], const double &rdist, 
	    const double &B1, const double &B2){
  double dmy, cos_theta;
  double jax[DIM];
  rigid_body_rotation(jax, e3, q0, BODY2SPACE);
  double nd1[DIM], nd2[DIM], nd3[DIM];
  double ar, ar2, ar3;

  if(rdist > A){
    cos_theta = jax[0]*nr[0] + jax[1]*nr[1] + jax[2]*nr[2];
    ar = (A)/rdist;
    ar2 = ar*ar;
    ar3 = ar*ar2;

    for(int d = 0; d < DIM; d++){
      nd1[d] = cos_theta * nr[d] - 1.0/3.0 * jax[d];
      nd2[d] = cos_theta * (cos_theta * nr[d] - jax[d]);
      nd3[d] = (3.0 * cos_theta * cos_theta - 1.0) * nr[d];
    }
    for(int d = 0; d < DIM; d++){
      ww[d] = B1 * ar3 * nd1[d] + 
	B2/2.0 * ar2 * ( 2.0 * ar2 * nd2[d] + (ar2 - 1.0) * nd3[d] );
    }
  }else{
    for(int d = 0; d < DIM; d++){
      ww[d] = 0.0;
    }
  }
}

void point_id(const int &iso, const int &ii, const int &jj, const int &kk, int *ei){
  switch(iso){
  case 0:
    ei[0] = ii;
    ei[1] = jj;
    ei[2] = kk;
    break;
  case 1:
    ei[0] = jj;
    ei[1] = ii;
    ei[2] = kk;
    break;
  case 2:
    ei[0] = jj;
    ei[1] = kk;
    ei[2] = ii;
    break;
  }
}
void get_image(const int &frame, const int &iso){
  double rp[DIM], sp[DIM], cp[DIM], vv[DIM], ww[DIM], dmy;
  int ei[DIM], edge[DIM];
  int Nj, Nk;
  int interpol_pid[8][DIM];
  double interpol_pval[8];
  double interpol_coeff[64];

  char dmy_path[256];
  FILE *fstream;

  double *eu, *ev, *ew;

  switch(iso){
  case 0:
    sprintf(dmy_path, "%s/plots/yz_%d.dat", AVS_dir, frame);
    eu = ey;
    ev = ex;
    ew = ez;

    Nj = HCY;
    Nk = HCZ;
    break;
  case 1:
    sprintf(dmy_path, "%s/plots/xz_%d.dat", AVS_dir, frame);
    eu = ex;
    ev = ey;
    ew = ez;

    Nj = HCX;
    Nk = HCZ;
    break;
  case 2:
    sprintf(dmy_path, "%s/plots/xy_%d.dat", AVS_dir, frame);
    eu = ex;
    ev = ez;
    ew = ey;

    Nj = HCX;
    Nk = HCY;
    break;
  }
  fstream = filecheckopen(dmy_path, "w");
  double v01 = v_inner_prod(v0, ev);
  double v02 = v_inner_prod(v0, eu);
  double v03 = v_inner_prod(v0, ew);
  double fdomain;

  int ii = 0;
  for(int jj = -Nj; jj <= Nj; jj++){
    for(int kk = -Nk; kk <= Nk; kk++){
      point_id(iso, ii, jj, kk, ei);

      dmy = 0.0;
      for(int d = 0; d < DIM; d++){
	rp[d] = (double)ei[d] * DX;
	sp[d] = r0[d] + rp[d];
	sp[d] = fmod(sp[d] + lbox[d], lbox[d]);
	edge[d] = (int)(sp[d]/DX);
	cp[d] = sp[d] - edge[d]*DX;

	dmy += rp[d]*rp[d];
      }
      dmy = sqrt(dmy);
      for(int d = 0; d < DIM; d++){
	rp[d] = (dmy > 0 ? rp[d] / dmy : 0.0);
      }
      v_gold(ww, rp, dmy, B1_real, B2);

      if(dmy > A){
	for(int d = 0; d < DIM; d++){
	  init_interpol(interpol_pid, interpol_pval, interpol_coeff, u[d], edge);
	  vv[d] = interpol(interpol_pid, interpol_pval, interpol_coeff, cp);
	}
	fdomain = 1.0;
      }else{
	for(int d = 0; d < DIM; d++){
	  vv[d] = 0.0;
	}
	fdomain = 0.0;
      }

      double vv1 = v_inner_prod(vv, ev);
      double vv2 = v_inner_prod(vv, eu);
      double vv3 = v_inner_prod(vv, ew);

      double ww1 = v_inner_prod(ww, ev);
      double ww2 = v_inner_prod(ww, eu);
      double ww3 = v_inner_prod(ww, ew);

      fprintf(fstream, 
	      "%3d %3d %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	      jj, kk, dmy, vv1, vv2, vv3, ww1, ww2, ww3,
	      fdomain*v01, fdomain*v02, fdomain*v03);
    }
  }
  fclose(fstream);
}

void get_image_rotate(const int &frame, const int &npi){
  double rp[DIM], sp[DIM], cp[DIM], vv[DIM], ww[DIM], dmy;
  int edge[DIM];
  int interpol_pid[8][DIM];
  double ei[DIM];
  double interpol_pval[8];
  double interpol_coeff[64];
  double theta, cs, sn;

  char dmy_path[256];
  FILE *fstream;
  sprintf(dmy_path, "%s/plots/wz-PI%d_%d.dat", AVS_dir, npi, frame);
  fstream = filecheckopen(dmy_path, "w");
  theta = M_PI / (double)npi;
  cs = cos(theta);
  sn = sin(theta);

  double eu[DIM] = {cs, sn, 0.0};
  double ev[DIM] = {-sn, cs, 0.0};
  double ew[DIM] = {0.0, 0.0, 1.0};

  int HN_UV = min(HCX, min(HCY, HCZ));
  for(int jj = -HN_UV; jj <= HN_UV; jj++){
    ei[0] = (double)jj * cs;
    ei[1] = (double)jj * sn;
    for(int kk = -HN_UV; kk <= HN_UV; kk++){
      ei[2] = (double) kk;

      dmy = 0.0;
      for(int d = 0; d < DIM; d++){
	rp[d] = ei[d] * DX;
	sp[d] = r0[d] + rp[d];
	sp[d] = fmod(sp[d] + lbox[d], lbox[d]);
	edge[d] = (int)(sp[d]/DX);
	cp[d] = sp[d] - edge[d]*DX;

	dmy += rp[d]*rp[d];
      }
      dmy = sqrt(dmy);
      for(int d = 0; d < DIM; d++){
	rp[d] = (dmy > 0 ? rp[d] / dmy : 0.0);
      }
      v_gold(ww, rp, dmy, B1_real, B2);

      if(dmy > A){
	for(int d = 0; d < DIM; d++){
	  init_interpol(interpol_pid, interpol_pval, interpol_coeff, u[d], edge);
	  vv[d] = interpol(interpol_pid, interpol_pval, interpol_coeff, cp);
	}
      }else{
	for(int d = 0; d < DIM; d++){
	  vv[d] = 0.0;
	}
      }

      fprintf(fstream, "%3d %3d %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	      jj, kk, dmy,
	      v_inner_prod(vv, ev), v_inner_prod(vv, eu), v_inner_prod(vv, ew),
	      v_inner_prod(ww, ev), v_inner_prod(ww, eu), v_inner_prod(ww, ew),
	      v_inner_prod(v0, ev), v_inner_prod(v0, eu), v_inner_prod(v0, ew));
    }
  }
  fclose(fstream);
}



