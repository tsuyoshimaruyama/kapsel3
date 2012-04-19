
#include "quaternion.h"

inline void qtn_print(quaternion &q){
  printf("     %.6f %.6f %.6f %.6f\n",
	 q.s, q.v[0], q.v[1], q.v[2]);
}

/*!
  \brief Matrix multiplication C = A.B
  \param[in] ha heigth of A (number of rows)
  \param[in] wa width of A (number of columns)
  \param[in] wb width of B 
 */
void mat_mul(double *c, double *a, double *b, 
		    const int &ha, const int &wa, const int &wb)
{
  int ima, imb, imc;
  double ab;
  for(int i = 0; i < ha; i++){
    for(int j = 0; j < wb; j++){
      imc = i*wb + j;

      ab = 0.0;
      for(int k = 0; k < wa; k++){
	ima = i*wa + k;
	imb = k*wb + j;
	ab += a[ima] * b[imb];
      }
      
      c[imc] = ab;
    }
  }
}
/*!
  \brief Matrix transpose At = Transpose(A)
  \param[in] h height of A (number of rows)
  \param[in] w width of A (number of columns)
 */
void mat_trans(double *at, double *a, const int &h, const int &w){
  int im, imt;
  for(int i = 0; i < h; i++){
    for(int j = 0; j < w; j++){
      im = i*w + j;
      imt = j*h + i;
      at[imt] = a[im];
    }
  }
}

/*
  \brief Compare to test if both matrices are equal
 */
int mat_cmp(double *a, double *b, const int &h, const int &w){
  double tol = 5.0e-5;
  double dmax = 0.0;
  int im;
  for(int i = 0; i < h; i++){
    for(int j = 0; j < w; j++){
      im = i*w + j;
      dmax = (b[im]!=0.0 ? 
	      MAX(dmax, ABS((a[im] - b[im])/b[im])) :
	      MAX(dmax, ABS(a[im] - b[im])));
    }
  }
  return (dmax <= tol ? 1 : 0);
}
inline void mat_print(double *a, const int &h, const int &w){
  int im;
  for(int i = 0; i < h; i++){
    for(int j = 0; j < w; j++){
      im = i*w + j;
      printf("%.5f ", a[im]);
    }
    printf("\n");
  }
  printf("\n");
}


//Test quaternion algebra implementation
int main(int argc, char**argv){
  quaternion q1,q2,q3,q4,q5,q_one;
  quaternion qc1, qc2;
  quaternion qg1, qg2;
  double tol = 1.0e-5;
  
  qtn_init(q_one, 1.0, 0.0, 0.0, 0.0);
  qtn_init(q1, 0.800829, 0.906111, 0.628903, -0.0367839);
  qtn_init(q2, -0.894365, -0.658605, -0.908624, 0.646626);
  qtn_init(q3, -0.22935, 0.888158, 0.395134, -0.748591);
  qtn_init(q4, -0.708922, 0.71681, -0.0131132, 0.630038);
  qtn_init(q5, -0.212223, 0.00565613, 0.375573, -0.232433);
  
  printf("Testing Quaternion Albegra\n");
  // q1 + q2
  qtn_add(qc1, q1, q2);
  qtn_init(qg1, -0.093536, 0.247506, -0.279721, 0.609842);
  qtn_init(qc2, q1);
  qtn_add(qc2, q2);
  if(qtn_cmp(qc1, qg1, tol) && qtn_cmp(qc2, qg1, tol)){
    printf("test 1: ok\n");
  }  else{
    printf("test 1: failed!\n");
    qtn_print(qc1);
    qtn_print(qg1);
  }
  
  
  // q1.(-4.7*q3)
  qtn_prod(qc1, q1, q3, -4.7);
  qtn_init(qg1, 5.94304, -0.221791, -3.84381, 3.72046);
  qtn_init(qc2, q1);
  qtn_prod(qc2, q3, -4.7);
  if(qtn_cmp(qc1, qg1, tol) && qtn_cmp(qc2, qg1, tol)){
    printf("test 2: ok\n");
  }  else{
    printf("test 2: failed!\n");
    qtn_print(qc2);
    qtn_print(qg2);
  }

  // norm q4
  double qval = qtn_sqnorm(q4);
  double qref = 1.41351;
  if(ABS((qval - qref)/qref) < tol){
    printf("test 3: ok\n");
  }  else{
    printf("test 3: failed!\n");
    qtn_print(q4);
    printf(" %.5f %.5f %.10f\n", qval, qref,ABS((qval - qref)/qref));
  }

  // quaternion inverse
  qtn_inv(qc1, q5);
  qtn_init(qg1, -0.883707, -0.0235524, -1.56391, 0.967863);
  qtn_prod(qc2,qc1,q5);
  if(qtn_cmp(qc2, q_one) && qtn_cmp(qc1, qg1)){
    printf("test 4: ok\n");
  } else{
    printf("test 4: failed!\n");
    qtn_print(qc1);
    qtn_print(qg1);
    qtn_print(qc2);
    qtn_print(q_one);
  }

  // rotation quaternion
  printf("Testing Quaternions Rotations\n");
  double phi = 1.10326;
  double rn[DIM] = {-0.578009, -0.166735, 0.798815};
  double rv[DIM] = {-0.637696, -0.183952, 0.881302};
  qtn_init(qg1, 0.851671, -0.302921, -0.0873818, 0.41864);
  
  
  rv_qtn(qc1, rv);
  rv_qtn(qc2, phi, rn);
  qtn_scalar(phi, qc1);
  qtn_vector(rv, qc1);
  printf("Rotation : %.5f %.5f %.5f %.5f\n",
	 phi, rv[0], rv[1], rv[2]);
  if(qtn_cmp(qc1, qg1) && qtn_cmp(qc2, qg1)){
    printf("test 5: ok\n");
  }else{
    printf("test 5: failed!\n");
    qtn_print(qc1);
    qtn_print(qc2);
    qtn_print(qg1);
  }


  //rotation matrices 
  double QR_C[DIM][DIM];
  double QE_C[DIM][DIM+1];
  double QG_C[DIM][DIM+1];
  double QQ[DIM][DIM];

  double QRT_C[DIM][DIM];
  double QET_C[DIM+1][DIM];
  double QGT_C[DIM+1][DIM];


  double QR_G[DIM][DIM] = {{0.63421,-0.660148,-0.402471},
			   {0.766027,0.46596,0.442815},
			   {-0.104788,-0.589141,0.801207}};
  double QE_G[DIM][DIM+1] = {{0.302921,0.851671,-0.41864,-0.0873818},
			     {0.0873818,0.41864,0.851671,0.302921},
			     {-0.41864,0.0873818,-0.302921,0.851671}};
  double QG_G[DIM][DIM+1] = {{0.302921,0.851671,0.41864,0.0873818},
			     {0.0873818,-0.41864,0.851671,-0.302921},
			     {-0.41864,-0.0873818,0.302921,0.851671}};
  double ID[DIM][DIM] = {{1.0, 0.0, 0.0},
			 {0.0, 1.0, 0.0},
			 {0.0, 0.0, 1.0}};

  double QRT_G[DIM][DIM] = {{0.63421,0.766027,-0.104788},
			    {-0.660148,0.46596,-0.589141},
			    {-0.402471,0.442815,0.801207}};
  double QET_G[DIM+1][DIM] = {{0.302921,0.0873818,-0.41864},
			      {0.851671,0.41864,0.0873818},
			      {-0.41864,0.851671,-0.302921},
			      {-0.0873818,0.302921,0.851671}};
  double QGT_G[DIM+1][DIM] = {{0.302921,0.0873818,-0.41864},
			      {0.851671,-0.41864,-0.0873818},
			      {0.41864,0.851671,0.302921},
			      {0.0873818,-0.302921,0.851671}};


  qtn_rm(QR_C, qc1);
  qtn_rmt(QRT_C, qc1);
  qtn_rm_e(QE_C, qc1);
  qtn_rm_et(QET_C, qc1);
  qtn_rm_g(QG_C, qc1);
  qtn_rm_gt(QGT_C, qc1);

  // Make sure matrices are propertly computed
  if(mat_cmp(&QR_C[0][0], &QR_G[0][0], DIM, DIM)){//QR
    printf("test 6a: ok\n");
  }else{
    printf("test 6a: failed\n");
    mat_print(&QR_C[0][0], DIM, DIM);
    mat_print(&QR_G[0][0], DIM, DIM);
  }
  if(mat_cmp(&QRT_C[0][0], &QRT_G[0][0], DIM, DIM)){//QRT
    printf("test 6b: ok\n");
  }else{
    printf("test 6b: failed\n");
    mat_print(&QRT_C[0][0], DIM, DIM);
    mat_print(&QRT_G[0][0], DIM, DIM);
  }
  if(mat_cmp(&QE_C[0][0], &QE_G[0][0], DIM, DIM+1)){//QE
    printf("test 6c: ok\n");
  }else{
    printf("test 6c: failed\n");
    mat_print(&QE_C[0][0], DIM, DIM+1);
    mat_print(&QE_G[0][0], DIM, DIM+1);
  }
  if(mat_cmp(&QET_C[0][0], &QET_G[0][0], DIM+1, DIM)){//QET
    printf("test 6d: ok\n");
  }else{
    printf("test 6d: failed\n");
    mat_print(&QET_C[0][0], DIM+1, DIM);
    mat_print(&QET_G[0][0], DIM+1, DIM);
  }
  if(mat_cmp(&QG_C[0][0], &QG_G[0][0], DIM, DIM+1)){//QG
    printf("test 6e: ok\n");
  }else{
    printf("test 6e: failed\n");
    mat_print(&QG_C[0][0], DIM, DIM+1);
    mat_print(&QG_G[0][0], DIM, DIM+1);
  }
  if(mat_cmp(&QGT_C[0][0], &QGT_G[0][0], DIM+1, DIM)){//QGT
    printf("test 6f: ok\n");
  }else{
    printf("test 6f: failed\n");
    mat_print(&QGT_C[0][0], DIM+1, DIM);
    mat_print(&QGT_G[0][0], DIM+1, DIM);
  }
  
  // R*RT = ID
  mat_mul(&QQ[0][0], &QR_C[0][0], &QRT_C[0][0], DIM, DIM, DIM);
  if(mat_cmp(&QQ[0][0], &ID[0][0], DIM, DIM)){
    printf("test 7: ok\n");
  }else{
    printf("test 7: failed\n");
    mat_print(&QQ[0][0], DIM, DIM);
  }
  // E*ET = ID
  mat_mul(&QQ[0][0], &QE_C[0][0], &QET_C[0][0], DIM, DIM+1,DIM);
  if(mat_cmp(&QQ[0][0], &ID[0][0], DIM, DIM)){
    printf("test 8: ok\n");
  }else{
    printf("test 8: failed\n");
    mat_print(&QQ[0][0], DIM, DIM);
  }
  // G*GT = ID
  mat_mul(&QQ[0][0], &QG_C[0][0], &QGT_C[0][0], DIM, DIM+1,DIM);
  if(mat_cmp(&QQ[0][0], &ID[0][0], DIM, DIM)){
    printf("test 9: ok\n");
  }else{
    printf("test 9: failed\n");
    mat_print(&QQ[0][0], DIM, DIM);
  }
  // R = E GT
  mat_mul(&QQ[0][0], &QE_C[0][0], &QGT_C[0][0], DIM, DIM+1, DIM);
  if(mat_cmp(&QQ[0][0], &QR_C[0][0], DIM, DIM)){
    printf("test 10: ok\n");
  }else{
    printf("test 10: failed\n");
    mat_print(&QQ[0][0], DIM, DIM);
    mat_print(&QR_C[0][0], DIM, DIM);
  }
  // RT = G ET
  mat_mul(&QQ[0][0], &QG_C[0][0], &QET_C[0][0], DIM, DIM+1, DIM);
  if(mat_cmp(&QQ[0][0], &QRT_C[0][0], DIM, DIM)){
    printf("test 11: ok\n");
  }else{
    printf("test 11: failed\n");
    mat_print(&QQ[0][0], DIM, DIM);
    mat_print(&QRT_C[0][0], DIM, DIM);
  }

  return 0;
}
