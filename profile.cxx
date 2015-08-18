/*!
  \file profile.cxx
  \brief Smooth particle profile routines
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */

#include "profile.h"

inline void Copy_sekibun_cell(const int &np_domain
			 ,int** &src
			 ,int** &dest
			 ){
  {
    dest = alloc_2d_int(np_domain,DIM);
    for(int n = 0;n<np_domain;n++){
      for(int d=0;d<DIM;d++){
	dest[n][d] = src[n][d];
      }
    }
    free_2d_int(src);
  }
}
bool sekibun_comparison(int* a, int* b){
  return SQ(a[0]) + SQ(a[1]) + SQ(a[2]) < SQ(b[0]) + SQ(b[1]) + SQ(b[2]);
}
void Particle_domain(double (*profile_func)(const double &x, const double radius)
		     ,int &np_domain
		     ,int** &sekibun_cell
		     ){
  // sekibun_cell の初期化関数
  
  double dmy_xp[DIM];

  for(int d= 0; d< DIM; d++){
    dmy_xp[d] = 0.;
  }


  int n_mesh = 8*((int)A + (int)A_XI + 3)*((int)A+(int)A_XI + 3)*((int)A+(int)A_XI+3);
  int nh = (int)A + (int)A_XI + 2;
  int **dmy_sekibun_cell = alloc_2d_int(n_mesh,DIM);
  double MAX_RADIUS = RADIUS + HXI + 1.7321*DX;
  assert(MAX_RADIUS < nh * DX);

  int mesh = 0;
  for(int i= -nh; i<=nh; i++){
    for(int j= -nh; j<=nh; j++){
      for(int k= -nh; k<=nh; k++){
	double dmy_x[DIM];
	dmy_x[0] = (double)i*DX; 
	dmy_x[1] = (double)j*DX; 
	dmy_x[2] = (double)k*DX; 
	
	double dmy = sqrt(
			  SQ(dmy_x[0]-dmy_xp[0])
			  +SQ(dmy_x[1]-dmy_xp[1])
			  +SQ(dmy_x[2]-dmy_xp[2])
			  );
	//double dmy_phi=profile_func(dmy,RADIUS+DX);

	//if(dmy_phi > 0.0){
	//if(dmy_phi > 1.e-10){
	if(dmy < MAX_RADIUS){
	  dmy_sekibun_cell[mesh][0] = i;
	  dmy_sekibun_cell[mesh][1] = j;
	  dmy_sekibun_cell[mesh][2] = k;
	  mesh++;
	}
      }//k
    }//j
  }//i
  np_domain = mesh;

  Copy_sekibun_cell(np_domain, dmy_sekibun_cell, sekibun_cell);
  std::sort(sekibun_cell, sekibun_cell+np_domain, sekibun_comparison);
  for(int i = 0; i < np_domain; i++){
    fprintf(stderr, "%5d %5d %5d %5d %12.6f\n",
	    sekibun_cell[i][0], sekibun_cell[i][1], sekibun_cell[i][2],
	    SQ(sekibun_cell[i][0]) + SQ(sekibun_cell[i][1]) + SQ(sekibun_cell[i][2]),
	    sqrt(SQ(sekibun_cell[i][0]) + SQ(sekibun_cell[i][1]) + SQ(sekibun_cell[i][2]))	    
	    );
  }
}
