#include "avs_process.h"

int Cs[DIM];
int &CX = Cs[0];
int &CY = Cs[1];
int &CZ = Cs[2];
int HCs[DIM];
int &HCX = HCs[0];
int &HCY = HCs[1];
int &HCZ = HCs[2];

void v_gold(double ww[DIM], const double nr[DIM], const double &rdist, 
	    const double &B1, const double &B2, const double &RADIUS);
void get_image(const int &frame, const int &iso, const double &RADIUS);
void get_rms(const int &frame, const double &RADIUS);

int main(int argc, char* argv[]){
  int pid;                // id of centered particle
  UDFManager *udf_in;     // UDF input file
  
  init_io(argc, argv, udf_in, pid);
  get_system_data(udf_in);
  initialize(pid);

  for(int i = 0; i <= Num_snap; i++){
    setup_avs_frame();
    read_u();
    read_p(pid);
    get_image(i, 0, A);
    get_image(i, 1, A);
    get_image(i, 2, A);
    clear_avs_frame();
  }
  fclose(fluid_field);
  fclose(particle_field);
}

/////
void v_gold(double ww[DIM], const double nr[DIM], const double &rdist, 
	    const double &B1, const double &B2, const double &RADIUS){
  double dmy, cos_theta;
  double jax[DIM], nsin_theta[DIM];
  rigid_body_rotation(jax, e3, q0, BODY2SPACE);
  double nd1[DIM], nd2[DIM], nd3[DIM];
  double ar, ar2, ar3;

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
      ww[d] = 0.0;
    }
  }
}

void get_rms(const int &frame, const double &RADIUS){
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
	if(dmy_r >= RADIUS){
	  int im = (i * NY * NZ) + (j * NZ) + k;
	  scale = dmy_r/RADIUS;
	  scale *= scale;

	  v_gold(ww, nr, dmy_r, B1, B2, RADIUS);
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
}

// plot (jj,kk) 
void point_id(const int &iso, const int &ii, const int &jj, const int &kk, int *ei){
  switch(iso){
  case 0:
    ei[0] = ii; //skip
    ei[1] = jj;
    ei[2] = kk;
    break;
  case 1:
    ei[0] = jj;
    ei[1] = ii; //skip
    ei[2] = kk;
    break;
  case 2:
    ei[0] = jj;
    ei[1] = kk;
    ei[2] = ii; //skip
    break;
  }
}
void get_image(const int &frame, const int &iso, const double &RADIUS){
  double rp[DIM], vv[DIM], ww[DIM], dmy;
  int ei[DIM], gi[DIM];
  int Nj, Nk;
  char dmy_path[256];
  char dmy_idpath[256];
  FILE *fstream;
  FILE *idstream;

  double *eu, *ev, *ew;

  switch(iso){
  case 0:
    sprintf(dmy_path, "%s/plots/yz_%05.0f.dat", 
            AVS_dir, (double)frame);
    sprintf(dmy_idpath, "%s/plots/p_yz_%05.0f.dat", 
            AVS_dir, (double) frame);
    eu = ey;
    ew = ez;
    ev = ex; //skip

    Nj = HCY;
    Nk = HCZ;
    break;
  case 1:
    sprintf(dmy_path, "%s/plots/xz_%05.0f.dat", 
            AVS_dir,(double) frame);
    sprintf(dmy_idpath, "%s/plots/p_xz_%05.0f.dat",
            AVS_dir, (double) frame);
    eu = ex;
    ew = ez;
    ev = ey; //skip


    Nj = HCX;
    Nk = HCZ;
    break;
  case 2:
    sprintf(dmy_path, "%s/plots/xy_%05.0f.dat", 
            AVS_dir, (double) frame);
    sprintf(dmy_idpath, "%s/plots/p_xy_%05.0f.dat", 
            AVS_dir, (double) frame);
    eu = ex;
    ew = ey;
    ev = ez; //skip

    Nj = HCX;
    Nk = HCY;
    break;
  }
  fstream = filecheckopen(dmy_path, "w");


  double r01, r02, r03;
  double v01, v02, v03;
  double pdomain, fdomain;
  double vscale = 1.0/UP;
  int ii = 0;
  {   // particle data
    r01 = v_inner_prod(res0, ev); // skip
    r02 = v_inner_prod(res0, eu);
    r03 = v_inner_prod(res0, ew);

    v01 = v_inner_prod(v0, ev);  //skip
    v02 = v_inner_prod(v0, eu);
    v03 = v_inner_prod(v0, ew);

    idstream = filecheckopen(dmy_idpath, "w");
    fprintf(idstream, "%6.4g %6.4g %6.4g %6.4g\n", 
            r02, r03, v02*vscale, v03*vscale);
    fclose(idstream);
  }
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
      v_gold(ww, rp, dmy, B1, B2, RADIUS);
      pdomain = Phi(dmy, RADIUS);
      if(dmy >= RADIUS){
	fdomain = 1.0;
      }else{
	fdomain = 0.0;
      }

      for(int d = 0; d < DIM; d++){
	vv[d] = u[d][im];
      }

      // fluid data
      // simulation
      double xx1 = ii * DX; // ignore
      double xx2 = jj * DX;
      double xx3 = kk * DX;
      double vv1 = v_inner_prod(vv, ev); //ignore
      double vv2 = v_inner_prod(vv, eu);
      double vv3 = v_inner_prod(vv, ew);

      // analytical 
      double ww1 = v_inner_prod(ww, ev); //ignore
      double ww2 = v_inner_prod(ww, eu);
      double ww3 = v_inner_prod(ww, ew);

      double rscale = dmy/RADIUS;

      // particle data
      if(jj == 0 && kk == 0){
        xx2 = r02;
        xx3 = r03;
        //simulation
        vv2 = v02;
        vv3 = v03;
        //analytical
        ww2 = UP * v_inner_prod(n0, eu);
        ww3 = UP * v_inner_prod(n0, ew);
        rscale = 1.0;
      }

      fprintf(fstream, 
	      "%6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g\n",
	      xx2, xx3, 
	      vv2*vscale, vv3*vscale, 
	      ww2*vscale, ww3*vscale,
              rscale, fdomain);
    }
  }
  fclose(fstream);
}

