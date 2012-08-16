#include "avs_process.h"

double **u_avg;
double **u_phi_avg;
double **u_gold;
double **u_phi_gold;
int *u_cnt;
int *u_phi_cnt;
int nskip;

void update_average();
void print_average(const double &RADIUS);
void compute_phi_average();
int main(int argc, char* argv[]){
  int pid;                // id of centered particle
  UDFManager *udf_in;     // UDF input file

  
  init_io(argc, argv, udf_in, pid);
  get_system_data(udf_in);
  init_avs();
  init_avs_p(pid);
  init_avs_mem();


  double RADIUS = A;
  init_avg_mem(RADIUS);
  for(int i = 0; i <= Num_snap; i++){
    setup_avs_frame();
    if(i >= nskip){
      read_avs_u();
      read_avs_p(pid);
      update_average();
    }

    clear_avs_frame();
  }
  compute_phi_average();
  print_average(RADIUS);

  close_avs();
  clean_exit();
  exit(EXIT_SUCCESS);
  
}

void update_average(){
  double old_cnt, new_cnt, dmy_r, dmy;
  double rg[DIM], rp[DIM], vv[DIM], sv[DIM];
  int ei[DIM], di[DIM], im, im_avg, i0;

  i0 = HNs[0];
  dmy = 3.0/(2.0 * B1);
  for(int i = 0; i < NX; i++){
    ei[0] = i;
    for(int j = 0; j < NY; j++){
      ei[1] = j;
      for(int k = 0; k < NZ; k++){
        ei[2] = k;
        int im = (i * NY * NZ) + (j * NZ) + k;
        for(int d = 0; d < DIM; d++){
          rg[d] = ei[d] * DX;
          vv[d] = u[d][im];
        }
        rigid_body_rotation(sv, vv, q0, SPACE2BODY);
        Distance(r0, rg, dmy_r, rp);
        for(int d = 0; d < DIM; d++){
          di[d] = Nint(rp[d]/DX) + HNs[d];
        }
        PBC_ip(di);
        im_avg = (di[0] * NY * NZ) + (di[1] * NZ) + di[2];
        old_cnt = (double) u_cnt[im_avg]++;
        new_cnt = (double) u_cnt[im_avg];
        for(int d = 0; d < DIM; d++){
          u_avg[d][im_avg] = (u_avg[d][im_avg] * old_cnt + sv[d]) / new_cnt;
        }

      }
    }
  }
}

void compute_phi_average(){
  int im, im_avg;
  int gi[DIM], di[DIM-1];
  double rg[DIM], vv[DIM], ww[DIM];
  double n_h[DIM], n_s[DIM], n_phi[DIM];
  double h_val, s_val, phi_val;
  double v_phi[DIM], w_phi[DIM];
  double dmy_r, r_y, r_z, old_cnt, new_cnt;

  for(int i = 0; i < NX; i++){
    gi[0] = i - HNX;

    for(int j = 0; j < NY; j++){
      gi[1] = j - HNY;
      
      for(int k = 0; k < NZ; k++){
        gi[2] = k - HNY;
        im = (i * NY * NZ) + (j * NZ) + k;

        dmy_r = 0.0;
        for(int d = 0; d < DIM; d++){
          rg[d] = ((double)gi[d]) * DX;
          vv[d] = u_avg[d][im];
          ww[d] = u_gold[d][im];
          dmy_r += rg[0]*rg[0];
        }
        dmy_r = sqrt(dmy_r);
        if(non_zero_mp(dmy_r)){
          Cylindrical_coord(rg, n_h, n_s, n_phi, h_val, s_val, phi_val, ax0);
        }else{
          for(int d = 0; d < DIM; d++){
            n_h[d] = e3[d];
            n_s[d] = e2[d];
            n_phi[d] = -e1[d];
          }
        }
        r_y = v_inner_prod(rg, n_s);
        r_z = v_inner_prod(rg, n_h);
        for(int d = 0; d < DIM; d++){
          equal_mp(r_y * n_s[d] + r_z * n_h[d], rg[d]);
        }
        di[0] = Nint(r_y / DX) + HJNY;
        di[1] = Nint(r_z / DX) + HJNZ;
        PBC_ip(di[0], JNY);
        PBC_ip(di[1], JNZ);
        im_avg = (di[0] * JNZ) + di[1];
        old_cnt = (double) u_phi_cnt[im_avg]++;
        new_cnt = (double) u_phi_cnt[im_avg];

        v_phi[1] = v_inner_prod(n_s, vv);
        v_phi[2] = v_inner_prod(n_h, vv);
        w_phi[1] = v_inner_prod(n_s, ww);
        w_phi[2] = v_inner_prod(n_h, ww);
        for(int d = 0; d < DIM; d++){
          equal_mp(w_phi[1] * n_s[d] + w_phi[2] * n_h[d], ww[d]);
          w_phi[0] = 0.0;
        }
        v_phi[0] = v_inner_prod(n_phi, vv);

        for(int d = 0; d < DIM; d++){
          u_phi_avg[d][im_avg] = (u_phi_avg[d][im_avg] * old_cnt + v_phi[d]) / new_cnt;
          u_phi_gold[d][im_avg] = (u_phi_gold[d][im_avg] * old_cnt + w_phi[d]) / new_cnt;
        }

      }//k
    }//j
  }//i
}

void print_average(const double &RADIUS){
  const char* dmy_3d = "u_avg_3d.dat";
  const char* dmy_yz = "u_avg_yz.dat";
  const char* dmy_phi= "u_avg_phi.dat";
  FILE *f3d, *fyz, *fphi;
  f3d = filecheckopen(dmy_3d, "w");
  fyz = filecheckopen(dmy_yz, "w");
  fphi= filecheckopen(dmy_phi, "w");

  double dmy = 3.0 / (2.0 * B1);
  double dist; 
  double gi[DIM];
  double uu[DIM], ww[DIM];
  int i0 = HNX;
  double fdomain;

  for(int i = 0; i < NX; i++){
    gi[0] = (double)(i - HNX) * DX;

    for(int j = 0; j < NY; j++){
      gi[1] = (double)(j - HNY) * DX;
      
      for(int k = 0; k < NZ; k++){
        gi[2] = (double)(k - HNZ) * DX;

        smooth_v(i, j, k, u_avg, uu);
        dist = sqrt(SQ(gi[0]) + SQ(gi[1]) + SQ(gi[2]));
        fdomain = Phi(dist, RADIUS);
        fdomain = 1.0 - fdomain;

        int im = (i * NY * NZ) + (j * NZ) + k;
        fprintf(f3d, "%8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E\n",
                gi[0], gi[1], gi[2],
                dist/RADIUS, fdomain,
                u_gold[0][im]*dmy, u_gold[1][im]*dmy, u_gold[2][im]*dmy,
                u_avg[0][im]*dmy, u_avg[1][im]*dmy, u_avg[2][im]*dmy,
                uu[0]*dmy, uu[1]*dmy, uu[2]*dmy);
        if(i == i0){
          fprintf(fyz, "%8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E\n",
                  gi[1], gi[2],
                  dist/RADIUS, fdomain,
                  u_gold[1][im]*dmy, u_gold[2][im]*dmy,
                  u_avg[1][im]*dmy, u_avg[2][im]*dmy,
                  uu[1]*dmy, uu[2]*dmy);
        }
      } // k
      if(i == i0){
        fprintf(fyz, "\n");
      }
    }// j
    fprintf(f3d, "\n");
  }// i
  fclose(f3d);
  fclose(fyz);

  int i = HJNX;
  gi[0] = 0.0;
  for(int j = 0; j < JNY; j++){
    gi[1] = (double)(j - HJNY)*DX;

    for(int k = 0; k < JNZ; k++){
      gi[2] = (double)(k - HJNZ)*DX;

      smooth_v_2d(j, k, u_phi_avg, uu);
      smooth_v_2d(j, k, u_phi_gold, ww);
      dist = sqrt(SQ(gi[1]) + SQ(gi[2]));
      fdomain = Phi(dist, RADIUS);
      fdomain = 1.0 - fdomain;


      int im_phi = (j * JNZ) + k;
      fprintf(fphi, "%8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E %8.6E\n",
              gi[0], gi[1], gi[2],
              dist/RADIUS, fdomain,
              u_phi_gold[0][im_phi]*dmy, u_phi_gold[1][im_phi]*dmy, u_phi_gold[2][im_phi]*dmy,
              u_phi_avg[0][im_phi]*dmy, u_phi_avg[1][im_phi]*dmy, u_phi_avg[2][im_phi]*dmy,
              ww[0]*dmy, ww[1]*dmy, ww[2]*dmy,
              uu[0]*dmy, uu[1]*dmy, uu[2]*dmy);
      
    }//k
    fprintf(fphi, "\n");
  }//j
  fclose(fphi);
}
