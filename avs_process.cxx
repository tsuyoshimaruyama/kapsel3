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
    get_image(i, 0);
    get_image(i, 1);
    get_image(i, 2);
    clear_avs_frame();
  }
  fclose(fluid_field);
  fclose(particle_field);
  fclose(post_field);
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
	      "%3d %3d %8.5g  %8.5g  %8.5g  %8.5g  %8.5g  %8.5g  %8.5g  %8.5g  %8.5g  %8.5g \n",
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



