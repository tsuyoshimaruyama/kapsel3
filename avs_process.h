#ifndef AVS_PROCESS_H
#define AVS_PROCESS_H

#include "avs_utils.h"
extern int Cs[DIM];
extern int &CX; //
extern int &CY; //
extern int &CZ; //
extern int HCs[DIM];
extern int &HCX; //
extern int &HCY; //
extern int &HCZ; //

// 
extern bool error_calc;
extern double RMIN;
extern double RMAX;

extern int nshells;
extern int npoints;
extern int* shell_points;
extern double* shell_chi2;
extern double* shell_rms;
extern double* shell_avg;
extern double* shell_avg2;
extern double* shell_delta;
extern double* shell_abs_delta;

void wrong_invocation(){
  fprintf(stderr, "usage: pavs -p PID -l Ns -i UDF [-err fr]\n");
  fprintf(stderr, "       -p   PID\t ID (1...N) of centered particle\n");
  fprintf(stderr, "       -l   Ns \t Box size length\n");
  fprintf(stderr, "       -i   UDF\t Input file\n");
  fprintf(stderr, "       -err    \t Error calculation [optional]\n");
  exit_job(EXIT_FAILURE);
}



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


void init_io(int argc, char *argv[], UDFManager* &udf_in, int &pid){
  
  char* fname;
  int pset, lset, fset, vset, eset;
  pset = lset = fset = vset = 0;
  error_calc = false;
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
      HCs[0] = Cs[0] / 2;
      Cs[1] = Cs[2] = Cs[0];
      HCs[1] = HCs[2] = HCs[0];
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
    }else if(strcmp(argv[i], "-err") == 0){
      error_calc = true;
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
  fprintf(stderr, "# Box Size: %d %d %d\n", Cs[0], HCs[0], 2*HCs[0] + 1);
}

void zero_shells(){
  for(int i = 0; i < nshells; i++){
    shell_points[i] = 0;
    shell_chi2[i] = 0.0;
    shell_rms[i] = 0.0;
    shell_avg[i] = 0.0;
    shell_avg2[i] = 0.0;
    shell_delta[i] = 0.0;
    shell_abs_delta[i] = 0.0;
  }
  npoints = 0;
}
void set_particle_box(){
  
}
void init_error(){
  if(error_calc && Nparticles > 1){
    fprintf(stderr, "Error calculation only valid for single squirmer\n");
    exit_job(EXIT_FAILURE);
  }
  // reset particle centered box size
  HCs[0] = Cs[0] / 2;
  while((double) HCs[0] * DX > maxlbox ||
        2*HCs[0] + 1 >= MIN(NX, MIN(NY, NZ))){
    Cs[0] = Cs[0] - 2;
    HCs[0] = Cs[0] / 2;
  }
  Cs[1] = Cs[2] = Cs[0];
  HCs[1] = HCs[2] = HCs[0];
  fprintf(stderr, "#EBOX Size: %d %d %d\n", Cs[0], HCs[0], 2*HCs[0] + 1);

  RMIN = A + HXI;
  nshells = HCs[0] - 1;        //x_int != x
  RMAX = nshells * DX;
  assert(RMAX > RMIN + DX);

  shell_points = alloc_1d_int(nshells);
  shell_chi2 = alloc_1d_double(nshells);
  shell_rms = alloc_1d_double(nshells);
  shell_avg = alloc_1d_double(nshells);
  shell_avg2 = alloc_1d_double(nshells);
  shell_delta = alloc_1d_double(nshells);
  shell_abs_delta = alloc_1d_double(nshells);
  zero_shells();
  fprintf(stderr, "#Error Calculation: %d %g %g\n\n", nshells, RMIN, RMAX);
}
double mean_shell(const double *&xshell, const int &imax){
  assert(imax >= 0 && imax <= nshells);
  double avg = 0.0;
  int count = 0;

  for(int i = 0; i < imax; i++){
    if(shell_points[i] > 0){
      avg += xshell[i];
      count++;
    }
  }
  return avg/((double)count);
}
double total_shell(const double *&xshell, const int &imax){
  assert(imax >= 0 && imax <= nshells);
  double sum = 0.0;

  for(int i = 0; i < imax; i++){
    sum += (double)shell_points[i] * xshell[i];
  }
  return sum / npoints;
}
double max_shell(const double *&xshell, const int &imax){
  assert(imax >= 0 && imax <= nshells);
  int i0 = 0;
  while(shell_points[i0] == 0) i0++;

  double xmax = xshell[i0];
  for(int i = i0 + 1; i < imax; i++){
    if(shell_points[i] > 0){
      xmax = MAX(xmax, xshell[i]);
    }
  }
  return xmax;
}
double min_shell(const double *&xshell, const int &imax){
  assert(imax >= 0 && imax <= nshells);
  int i0 = 0;
  while(shell_points[i0] == 0) i0++;

  double xmin = xshell[i0];
  for(int i = i0 + 1; i < imax; i++){
    if(shell_points[i] > 0){
      xmin = MIN(xmin, xshell[i]);
    }
  }
  return xmin;
}


#endif
