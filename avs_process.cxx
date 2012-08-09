#include "avs_process.h"

int Cs[DIM];
int &CX = Cs[0];
int &CY = Cs[1];
int &CZ = Cs[2];

int HCs[DIM];
int &HCX = HCs[0];
int &HCY = HCs[1];
int &HCZ = HCs[2];

bool error_calc;
double RMIN;             // minimum distance for error analysis (= A + DX)
double RMAX;             // maximum distance for error analysis

int nshells;             // number of shells [cte]
int npoints;             // total number of points (all shells)
int* shell_points;       // number of points in shell
double* shell_chi2;      // relative chi2
double* shell_rms;       // rms error
double* shell_avg;       // average velocity
double* shell_avg2;      // average squared velocity
double* shell_delta;     // max relative error
double* shell_abs_delta; // max error

void v_gold(double ww[DIM], const double nr[DIM], const double &rdist, 
	    const double &B1, const double &B2, const double &RADIUS);
void get_image(const int &frame, const int &iso, const double &RADIUS);
void get_error_shells(const int &frame, const double &RADIUS);

inline void point_id(const int &iso, const int &ii, const int &jj, const int &kk, int *ei){
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

int main(int argc, char* argv[]){
  int pid;                // id of centered particle
  UDFManager *udf_in;     // UDF input file
  double RADIUS;
  
  init_io(argc, argv, udf_in, pid);

  get_system_data(udf_in);
  init_avs();
  init_avs_p(pid);
  init_avs_mem();

  init_error();  
  for(int i = 0; i <= Num_snap; i++){
    RADIUS = A;
    setup_avs_frame();
    read_avs_u();
    read_avs_p(pid);
    get_error_shells(i, RADIUS);
    get_image(i, 0, RADIUS);
    //    get_image(i, 1, RADIUS);
    //    get_image(i, 2, RADIUS);
    clear_avs_frame();
  }
  close_avs();
}

void get_error_shells(const int &frame, const double &RADIUS){
  int im, dmy_shell;
  int ei[DIM], gi[DIM];
  double avw, wnorm, wnorm2, dmy_rms, dmy_r, dmy_r2;
  double vv[DIM], ww[DIM], nr[DIM], rgrid[DIM], rp[DIM];

  zero_shells();
  for(int ii = -HCX; ii <= HCX; ii++){
    ei[0] = ii;
    for(int jj = -HCY; jj <= HCY; jj++){
      ei[1] = jj;
      for(int kk = -HCZ; kk <= HCZ; kk++){
        ei[2] = kk;

        dmy_r2 = 0.0;
        for(int d = 0; d < DIM; d++){
          rgrid[d] = (double)ei[d] * DX;
          rp[d] = rgrid[d] - res0[d];
          dmy_r2 += rp[d]*rp[d];
        }
        dmy_r2 = sqrt(dmy_r2);
        Distance(r0, rgrid, dmy_r, nr);
        equal_mp(dmy_r, dmy_r2);
        
        if(dmy_r >= RMIN && dmy_r <= RMAX){
          dmy_shell = (int) dmy_r / DX;
          for(int d = 0; d < DIM; d++){
            gi[d] = r0_int[d] + ei[d];
          }
          PBC_ip(gi);
          v_gold(ww, nr, dmy_r, B1, B2, RADIUS);
          wnorm2 = SQ(ww[0]) + SQ(ww[1]) + SQ(ww[2]);
          wnorm = sqrt(wnorm2);

          im = (gi[0] * NY * NZ) + (gi[1] * NZ) + gi[2];
          for(int d = 0; d < DIM; d++){
            vv[d] = u[d][im];
          }

          dmy_rms = 0.0;
          for(int d = 0; d < DIM; d++){
            dmy_rms += SQ(vv[d] - ww[d]);
          }
          avw = sqrt(dmy_rms);

          shell_delta[dmy_shell] = MAX(shell_delta[dmy_shell], avw / wnorm);
          shell_abs_delta[dmy_shell] = MAX(shell_abs_delta[dmy_shell], avw);

          shell_chi2[dmy_shell] += dmy_rms / wnorm2;
          shell_rms[dmy_shell] += dmy_rms;
          shell_avg[dmy_shell] += wnorm;
          shell_avg2[dmy_shell] += wnorm2;

          shell_points[dmy_shell]++;
          npoints++;
        }
        
      }//kk
    }//jj
  }//ii
  double dmy_norm;
  for(int i = 0; i < nshells; i++){
    if(shell_points[i] > 0){
      dmy_norm = 1.0 / (double) shell_points[i];
      shell_chi2[i] *= dmy_norm;
      shell_rms[i] *= dmy_norm;
      shell_avg[i] *= dmy_norm;
      shell_avg2[i] *= dmy_norm;
    }
  }

  char dmy_path[256];
  FILE *fstream;
  sprintf(dmy_path, "%s/plots/err_%05.0f.dat", AVS_dir, (double)frame);
  fstream = filecheckopen(dmy_path, "w");
  for(int i = 0; i < nshells; i++){
    if(shell_points[i] > 0){
      fprintf(fstream, "%6.4g %d %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g\n", 
              ((double)i + 0.5) * DX,
              shell_points[i],
              shell_chi2[i],
              shell_rms[i] / shell_avg2[i],
              shell_rms[i] / SQ(shell_avg[i]),
              shell_abs_delta[i] / sqrt(shell_avg2[i]),
              shell_abs_delta[i] / shell_avg[i],
              shell_delta[i],
              sqrt(shell_avg2[i]),
              shell_avg[i]
              );
    }else{
      fprintf(fstream, "%6.4g %d %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g\n",
              ((double)i + 0.5) * DX,
              0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
  }
  fclose(fstream);
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
      fdomain = 1.0 - pdomain;

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

      // error data
      double rscale = dmy/RADIUS;
      double rms, wmag;
      if(error_calc && dmy >= RMIN){
        wmag = SQ(ww1) + SQ(ww2) + SQ(ww3);
        rms = (SQ(vv1 - ww1) + SQ(vv2 - ww2) + SQ(vv3 - ww3)) / wmag;
      }else{
        rms = 0.0;
        wmag = 0.0;
      }
      fprintf(fstream, 
	      "%6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g\n",
	      xx2, xx3, 
	      vv2*vscale, vv3*vscale, 
	      ww2*vscale, ww3*vscale,
              rscale, fdomain,
              rms, sqrt(wmag)*vscale);
    }
    fprintf(fstream, "\n");
  }
  fclose(fstream);
}

