#ifndef LAD3_H
#define LAD3_H

//////////////////////////////////////////////////
//////////////////////////////////////////////////
// Basic Matrix / Vector routines 
//////////////////////////////////////////////////
//////////////////////////////////////////////////
inline double v_rms(const double a[DIM], const double b[DIM]){
  double c[DIM];
  for(int d = 0; d < DIM; d++){
    c[d] = a[d] - b[d];
  }
  return sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
}

//
// Copy routines
inline void v_copy(double copy[DIM], const double original[DIM]){
  for(int i = 0; i < DIM; i++)
    copy[i] = original[i];
}
inline void M_copy(double copy[DIM][DIM], const double original[DIM][DIM]){
  for(int i = 0; i < DIM; i++)
    for(int j = 0; j < DIM; j++)
      copy[i][j] = original[i][j];
}

// Compare two matrices (user specified relative error)
inline int M_cmp(const double A[DIM][DIM], const double B[DIM][DIM],
		 const double tol=LARGE_TOL_MP){
  bool mequal = true;
  for(int i = 0; i < DIM; i++){
    for(int j = 0; j < DIM; j++){
      mequal = mequal && equal_tol(A[i][j], B[i][j], tol);
    }
  }
  return (mequal ? 1 : 0);
}

//
// Scaling routines
inline void v_scale(double v[DIM], const double scale){
  for(int i = 0; i < DIM; i++)
    v[i] *= scale;
}
inline void M_scale(double A[DIM][DIM], const double scale){
  for(int i = 0; i < DIM; i++)
    for(int j = 0; j < DIM; j++)
      A[i][j] *= scale;
}

// 
// Vector add : c = a + alpha * b
inline void v_add(double c[DIM], const double a[DIM], const double b[DIM],
		  const double alpha=1.0){
  for(int i = 0; i < DIM; i++)
    c[i] = a[i] + alpha*b[i]; 
}
inline void v_add(double c[DIM], const double b[DIM], 
		  const double alpha=1.0){
  double a[DIM];
  v_copy(a, c);
  v_add(c, a, b, alpha);
}

//
// Matrix add : C = A + alpha*B
inline void M_add(double C[DIM][DIM], const double A[DIM][DIM], 
		  const double B[DIM][DIM], const double alpha=1.0){
  for(int i = 0; i < DIM; i++)
    for(int j = 0; j < DIM; j++)
      C[i][j] = A[i][j] + alpha * B[i][j];
}
inline void M_add(double C[DIM][DIM], const double B[DIM][DIM],
		  const double alpha=1.0){
  double A[DIM][DIM];
  M_copy(A, C);
  M_add(C, A, B, alpha);
}

//
// Matrix Multiply : C = alpha*(A.B)
inline void M_prod(double AB[DIM][DIM], 
		   const double A[DIM][DIM],
		   const double B[DIM][DIM],
		   const double alpha=1.0){
  double dmy;
  for(int i = 0; i < DIM; i++){
    for(int j = 0; j < DIM; j++){
      dmy = 0.0;
      for(int k = 0; k < DIM; k++){
	dmy += A[i][k]*B[k][j];
      }
      AB[i][j] = alpha*dmy;
    }
  }
}
inline void M_prod(double AB[DIM][DIM], const double B[DIM][DIM],
		   const double alpha=1.0){
  double A[DIM][DIM];
  M_copy(A, AB);
  M_prod(AB, A, B, alpha);
}

//
// Matrix determinant
inline double M_det(const double A[DIM][DIM]){
  assert(DIM == 3);
  return -A[0][2]*A[1][1]*A[2][0] + A[0][1]*A[1][2]*A[2][0] + 
    A[0][2]*A[1][0]*A[2][1] - A[0][0]*A[1][2]*A[2][1] - 
    A[0][1]*A[1][0]*A[2][2] + A[0][0]*A[1][1]*A[2][2];
}

//
// Matrix Inverse : B = alpha*(A^-1)
inline void M_inv(double B[DIM][DIM], const double A[DIM][DIM],
		  const double alpha=1.0){
  assert(DIM == 3);
  double idetA = M_det(A);
  assert(non_zero_mp(idetA));
  idetA = alpha/idetA;

  B[0][0] = (A[1][1]*A[2][2] - A[1][2]*A[2][1])*idetA;
  B[0][1] = (A[0][2]*A[2][1] - A[0][1]*A[2][2])*idetA;
  B[0][2] = (A[0][1]*A[1][2] - A[0][2]*A[1][1])*idetA;

  B[1][0] = (A[1][2]*A[2][0] - A[1][0]*A[2][2])*idetA;
  B[1][1] = (A[0][0]*A[2][2] - A[0][2]*A[2][0])*idetA;
  B[1][2] = (A[0][2]*A[1][0] - A[0][0]*A[1][2])*idetA;

  B[2][0] = (A[1][0]*A[2][1] - A[1][1]*A[2][0])*idetA;
  B[2][1] = (A[0][1]*A[2][0] - A[0][0]*A[2][1])*idetA;
  B[2][2] = (A[0][0]*A[1][1] - A[0][1]*A[1][0])*idetA;
}
inline void M_inv(double B[DIM][DIM], const double alpha=1.0){
  double A[DIM][DIM];
  M_copy(A, B);
  M_inv(B, A, alpha);
}

//
// Matrix Transpose
inline void M_trans(double B[DIM][DIM], const double A[DIM][DIM],
		    const double alpha=1.0){
  for(int i = 0; i < DIM; i++)
    for(int j = 0; j < DIM; j++)
      B[i][j] = alpha*A[j][i];
}
inline void M_trans(double B[DIM][DIM], const double alpha=1.0){
  double A[DIM][DIM];
  M_copy(A, B);
  M_trans(B, A, alpha);
}

//
// Verify rotation matrix
inline void M_isValidRotation(double QR[DIM][DIM]){
  double ID[DIM][DIM] = {{1.0, 0.0, 0.0},
		      {0.0, 1.0, 0.0},
		      {0.0, 0.0, 1.0}};
  double iQR[DIM][DIM];
  double tQR[DIM][DIM];
  int det, test_det, test_inv;

  test_det = (equal_tol(M_det(QR), 1.0, LARGE_TOL_MP) ? 1 : 0);
  M_trans(tQR, QR);
  M_prod(iQR, tQR, QR);
  test_inv = M_cmp(iQR, ID, LARGE_TOL_MP);
  assert(test_det && test_inv);
}

//
// Vector cross product: c = a x (alpha * b)
inline void v_cross(double c[DIM], 
		    const double a[DIM],
		    const double b[DIM], 
		    const double alpha=1.0){
  assert(DIM == 3);
  c[0] = alpha*(a[1]*b[2] - a[2]*b[1]);
  c[1] = alpha*(a[2]*b[0] - a[0]*b[2]);
  c[2] = alpha*(a[0]*b[1] - a[1]*b[0]);
}
inline void v_cross(double c[DIM],
		    const double b[DIM],
		    const double alpha=1.0){
  double a[DIM];
  v_copy(a, c);
  v_cross(c, a, b, alpha);
}

//
// Vector inner product: ab = alpha(a.b)
inline double v_inner_prod(const double a[DIM],
		    const double b[DIM]){
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
inline double v_norm(const double a[DIM]){
  return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

//
// Vector outer product: (AB)_ij = alpha*(a_i*b_j)
inline void v_outer_prod(double AB[DIM][DIM],
		  const double a[DIM],
		  const double b[DIM],
		  const double alpha=1.0){
  for(int i = 0; i < DIM; i++)
    for(int j = 0; j < DIM; j++)
      AB[i][j] = alpha*a[i]*b[j];
}

//
// Matrix - Vector product : x' = alpha*(A.x)
inline void M_v_prod(double Ax[DIM],
		     const double A[DIM][DIM],
		     const double x[DIM],
		     const double alpha=1.0){
  double dmy;
  for(int i = 0; i < DIM; i++){
    dmy = 0.0;
    for(int j = 0; j < DIM; j++){
      dmy += A[i][j] * x[j];
    }
    Ax[i] = alpha*dmy;
  }
}
inline void M_v_prod_add(double Ax[DIM],
			 const double A[DIM][DIM],
			 const double x[DIM],
			 const double alpha=1.0){
  double dmy;
  for(int i = 0; i < DIM; i++){
    dmy = 0.0;
    for(int j = 0; j < DIM; j++){
      dmy += A[i][j] * x[j];
    }
    Ax[i] += alpha*dmy;
  }
}
inline void M_v_prod(double Ax[DIM],
		     const double A[DIM][DIM],
		     const double alpha=1.0){
  double x[DIM];
  v_copy(x, Ax);
  M_v_prod(Ax, A, x, alpha);
}

//
// Vector - Matrix product : x' = alpha*(x.A) = alpha*(Transpose(A).x)
inline void v_M_prod(double xA[DIM],
		     const double x[DIM],
		     const double A[DIM][DIM],
		     const double alpha=1.0){
  double dmy;  
  for(int i = 0; i < DIM; i++){
    dmy = 0.0;
    for(int j = 0; j < DIM; j++){
      dmy += x[j]*A[j][i];
    }
    xA[i] = alpha*dmy;
  }
}
inline void v_M_prod(double xA[DIM],
		     const double A[DIM][DIM],
		     const double alpha=1.0){
  double x[DIM];
  v_copy(x, xA);
  v_M_prod(xA, x, A, alpha);
}

#endif
