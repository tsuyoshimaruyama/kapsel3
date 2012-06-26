#include "avs_process.h"
void v_gold(double ww[DIM], const double nr[DIM], const double &rdist, 
	    const double &B1, const double &B2);
void get_image(const int &frame, const int &iso);
void get_rms(const int &frame);

int main(int argc, char* argv[]){
  int pid;                // id of centered particle
  int *p_spec;            // species id for all particles
  JAX *sp_axis;           // janus axis for each species
  double *sp_slip;        // slip velocity for each species
  double *sp_slipmode;    // B2 coefficient

  UDFManager *udf_in;     // UDF input file
  
  rms[0] = rms[1] = rms[2] = rms[3] =  0.0;
  initialize(argc, argv, udf_in, pid);
  get_system_data(udf_in, p_spec, sp_axis, sp_slip, sp_slipmode);
  initialize_avs();

  u = (double **)malloc(sizeof(double*)* DIM);
  for(int d = 0; d < DIM; d++){
    u[d] = alloc_1d_double(NX*NY*NZ);
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
  B2 = B2*B1_real;
  UP = 2./3. * B1_real;
  fprintf(stderr, "# UP = %g (%g)\n",UP,B1_real);

  for(int i = 0; i <= Num_snap; i++){
    setup_avs_frame();
    read_u();
    read_p(pid);
    get_rms(i);
    get_image(i, 0);
    get_image(i, 2);
    clear_avs_frame();
  }
  fprintf(stderr, "# Error: %8.5e %8.5e %8.5e %8.5e\n", 
	  rms[0]/(double)Num_snap, 
	  rms[1]/(double)Num_snap, 
	  rms[2]/(double)Num_snap,
	  rms[3]/(double)Num_snap);
  fclose(fluid_field);
  fclose(particle_field);
}

/////
void v_gold(double ww[DIM], const double nr[DIM], const double &rdist, 
	    const double &B1, const double &B2){
  double dmy, cos_theta;
  double jax[DIM], nsin_theta[DIM];
  rigid_body_rotation(jax, e3, q0, BODY2SPACE);
  double nd1[DIM], nd2[DIM], nd3[DIM];
  double ar, ar2, ar3;

  if(rdist >= A){
    cos_theta = jax[0]*nr[0] + jax[1]*nr[1] + jax[2]*nr[2];
    ar = (A)/rdist;
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
      ww[d] = 0.0;
    }
  }
  for(int d = 0; d < DIM; d++){
    assert(ww[d] == ww[d]);
  }
}

void get_rms(const int &frame){
  int ei[DIM], count;
  double rgrid[DIM], nr[DIM], ww[DIM], df[DIM];
  double e0, e1, e2, e3, dmy_r, scale, dmy_e;
  e0 = e1 = e2 = e3 = 0.0;
  count = 0;
  for(int i = 0; i < NX; i++){
    ei[0] = i;
    for(int j = 0; j < NY; j++){
      ei[1] = j;
      for(int k = 0; k < NZ; k++){
	ei[2] = k;

	for(int d = 0; d < DIM; d++){
	  rgrid[d] = (double)ei[d] * DX;
	}
	Distance(r0, rgrid, dmy_r, nr);
	if(dmy_r >= A){
	  int im = (i * NY * NZ) + (j * NZ) + k;
	  scale = dmy_r/(A);
	  scale *= scale;

	  v_gold(ww, nr, dmy_r, B1_real, B2);
	  for(int d = 0; d < DIM; d++){
	    df[d] = (u[d][im] - ww[d])/UP;
	  }
	  dmy_e = v_inner_prod(df, df);

	  e0 += dmy_e;
	  e1 += dmy_e * scale;
	  e2 += dmy_e * scale*scale;
	  e3 += dmy_e * scale*scale*scale;
	  count ++;
	}
      }
    }
  }
  if(count > 0){
    e0 = sqrt(e0 / (double)count);
    e1 = sqrt(e1 / (double)count);
    e2 = sqrt(e2 / (double)count);
    e3 = sqrt(e3 / (double)count);
  }
  fprintf(stderr, "%d %8.5g %8.5g %8.5g %8.5g\n", frame, e0, e1, e2, e3);
  if(frame > 0){
    rms[0] += e0;
    rms[1] += e1;
    rms[2] += e2;
    rms[3] += e3;
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
  double rp[DIM], vv[DIM], ww[DIM], dmy;
  int ei[DIM], gi[DIM];
  int Nj, Nk;
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
  double rr1 = v_inner_prod(res0, ev);
  double rr2 = v_inner_prod(res0, eu);
  double rr3 = v_inner_prod(res0, ew);
  double pdomain, fdomain;

  int ii = 0;
  for(int jj = -Nj; jj <= Nj; jj++){
    for(int kk = -Nk; kk <= Nk; kk++){
      point_id(iso, ii, jj, kk, ei);
      for(int d = 0; d < DIM; d++){
	gi[d] = ei[d] + r0_int[d];
      }
      PBC_ip(gi);
      int im = (gi[0] * NY * NZ) + (gi[1] * NZ) + gi[2];

      dmy = 0.0;
      for(int d = 0; d < DIM; d++){
	rp[d] = (double)ei[d] * DX - res0[d];
	dmy += rp[d]*rp[d];
      }
      dmy = sqrt(dmy);
      for(int d = 0; d < DIM; d++){
	rp[d] = (dmy > 0 ? rp[d] / dmy : 0.0);
      }
      v_gold(ww, rp, dmy, B1_real, B2);
      pdomain = Phi(dmy, A);
      if(dmy >= A){
	fdomain = 1.0;
      }else{
	fdomain = 0.0;
      }
      

      for(int d = 0; d < DIM; d++){
	vv[d] = fdomain*u[d][im];
      }


      double vv1 = v_inner_prod(vv, ev);
      double vv2 = v_inner_prod(vv, eu);
      double vv3 = v_inner_prod(vv, ew);

      double ww1 = v_inner_prod(ww, ev);
      double ww2 = v_inner_prod(ww, eu);
      double ww3 = v_inner_prod(ww, ew);

      fprintf(fstream, 
	      "%3d %3d %8.5g  %8.5g  %8.5g  %8.5g  %8.5g  %8.5g  %8.5g  %8.5g  %8.5g  %8.5g %8.5g %8.5g %8.5g %8.5f\n",
	      jj, kk, dmy, 
	      vv1, vv2, vv3, 
	      fdomain*ww1, fdomain*ww2, fdomain*ww3,
	      fdomain*v01, fdomain*v02, fdomain*v03,
	      rr1, rr2, rr3, 
	      pdomain);
    }
  }
  fclose(fstream);
}

