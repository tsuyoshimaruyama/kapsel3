/*!
  \file ribid_body.h
  \brief Routines to solve equations of motion for orientation of a 
  symmetrical sphere
  \details Details on the quaternion calculations can be found in
  Graf, B. (2008). Quaternions and dynamics. 
  http://arxiv.org/abs/0811.2889
  Details on the Lie-Group representation of rotations and the use of simplectic integrators for rigid body motion are given in
  Marsden, J. E., & Ratiu, T. S. (2010). Introduction to Mechanics and Symmetry: A Basic Exposition of Classical Mechanical Systems (Texts in Applied Mathematics) (p. 285). Springer.
  and
  Dullweber, A., Leimkuhler, B., & McLachlan, R. (1997). Symplectic splitting methods for rigid body molecular dynamics. The Journal of Chemical Physics, 107(15), 5840–5851.
  and
  Buss, Samuel (2001). Accurate and efficient simulation of rigid body rotations. 
  \author J. Molina
  \date 2012/03/29
  \version 1.0
 */
#ifndef RIBID_BODY_H
#define RIBID_BODY_H

#include "quaternion.h"

enum COORD_SYSTEM {BODY_FRAME, SPACE_FRAME};
enum COORD_TRANS {BODY2SPACE, SPACE2BODY};

static double QRTOL=1.0e-12;
static double QRTOL_LARGE=1.0e-6;
//////////////////////////////////////////////////
//////////////////////////////////////////////////
// Basic Matrix / Vector routines 
//////////////////////////////////////////////////
//////////////////////////////////////////////////

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

// Compare two matrices
inline int M_cmp(const double A[DIM][DIM], const double B[DIM][DIM],
		 const double tol = QRTOL_LARGE){
  double dmax=0.0;
  for(int i = 0; i < DIM; i++)
    for(int j = 0; j < DIM; j++)
      dmax = ((B[i][j] != 0.0) ? 
	      MAX(dmax, ABS((A[i][j] - B[i][j]) / B[i][j])) :
	      MAX(dmax, ABS(A[i][j]))
	      );
  return ((dmax <= tol) ? 1 : 0);
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
  assert(idetA != 0.0);
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
inline void M_isValidRotation(double QR[DIM][DIM], 
			    const double tol=QRTOL){
  double ID[DIM][DIM] = {{1.0, 0.0, 0.0},
		      {0.0, 1.0, 0.0},
		      {0.0, 0.0, 1.0}};
  double iQR[DIM][DIM];
  double tQR[DIM][DIM];
  int det, test_det, test_inv;

  test_det = (ABS(M_det(QR) - 1.0) < tol) ? 1 : 0;
  M_trans(tQR, QR);
  M_prod(iQR, tQR, QR);
  test_inv = M_cmp(iQR, ID, tol);
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

/*!
  \brief Compute random rotation matrix
 */
inline void random_rotation(double QR[DIM][DIM]){
  quaternion dmy_q;
  random_rqtn(dmy_q);
  rqtn_rm(QR, dmy_q);
  M_isValidRotation(QR);
}


//////////////////////////////////////////////////
//////////////////////////////////////////////////
// Vector Transformation: Space <--> Body Coordinate Frames
//////////////////////////////////////////////////
//////////////////////////////////////////////////


/*!
  \brief Transform between body/space and space/body frames given the
  current orientation QUATERNION
 */
void rigid_body_rotation(double rotated[DIM], 
				const double original[DIM], 
				const quaternion &q, 
				const COORD_TRANS &transform);

/*!
  \brief Transform between body/space and space/body frames given the
  current orientation QUATERNION (in place)
 */
inline void rigid_body_rotation(double rotated[DIM], 
				const quaternion &q, 
				const COORD_TRANS &transform){
  double original[DIM];
  v_copy(original, rotated);
  rigid_body_rotation(rotated, original, q, transform);
}


/*!
  \brief Transform between body/space and space/body frames given the 
  current orientation MATRIX
 */
void rigid_body_rotation(double rotated[DIM], 
			 const double original[DIM],
			 const double QR[DIM][DIM], 
			 const COORD_TRANS &transform);
/*!
  \brief Transform between body/space and space/body frames given the 
  current orientation MATRIX (in place)
 */
inline void rigid_body_rotation(double rotated[DIM], 
				const double QR[DIM][DIM], 
				const COORD_TRANS &transform){
  double original[DIM];
  v_copy(original, rotated);
  rigid_body_rotation(rotated, original, QR, transform);
}

//////////////////////////////////////////////////
//
// Kinematics Calculations
//

// Quaternion Representation
/*!
  \brief Compute time derivative of orientation quaternion
  \param[out] dqdt time derivative of q
  \param[in] q current orientation quaternion
  \param[in] w current angular velocity
  \param[in] coord coordinate system of given angular velocity
 */
void qdot(quaternion &dqdt, 
	  const quaternion &q, 
	  const double omega[DIM], 
	  const COORD_SYSTEM &coord);
/*!
  \brief Compute time derivative of orientation quaternion (in place)
 */
inline void qdot(quaternion &dqdt,
		 const double omega[DIM],
		 const COORD_SYSTEM &coord);

/*!
  \brief Compute time derivative of orientation matrix
  \param[out] dQRdt time derivative of QR
  \param[in] QR current orientation matrix
  \param[in] omega current angular velocity
  \param[in] coord coordinate system for given angular velocity
 */
void Qdot(double dQRdt[DIM][DIM],
	  const double QR[DIM][DIM],
	  const double omega[DIM],
	  const COORD_SYSTEM &coord);
/*!
  \brief Compute time derivative of orientation matrix (in place)
 */
inline void Qdot(double dQRdt[DIM][DIM],
		 const double omega[DIM],
		 const COORD_SYSTEM &coord);

/*!
  \brief Compute time derivative of angular velocity given
  the current velocity and torque on the body, by solving the Euler
  Equations for rigid body motion
  \param[out] dwdt time derivative of angular velocity (body frame)
  \param[in] w angular velocity (body frame)
  \param[in] tau torque (body frame)
  \param[in] QR orientation matrix
  \param[in] I inertia tensor (body frame)
 */
void wdot(double dwdt[DIM], //rotation matrix represenation
	  const double w[DIM],
	  const double tau[DIM],
	  const double I[DIM][DIM]);
/*!
  \brief Compute time derivative of angular velocity (in place)
 */
inline void wdot(double dwdt[DIM],
		 const double tau[DIM],
		 const double I[DIM][DIM]);

/*!
  \brief Update orientation/ angular velocity using fourth order RK scheme
  of a free rigid body (torque=0)
  \param[in,out] q orientation quaternion
  \param[in,out] omega angular velocity (BODY FRAME)
  \param[in] I intertia tensor
  \param[in] dt time step
 */
void propagate_wq_RK4(quaternion &q, 
		      double omega[DIM],
		      const double I[DIM][DIM],
		      const double &dt);


//Rotation Matrix - Lie Group representation
/*!
  \brief Update orientation/ angular velocity using fourth order RK scheme 
  for omega and first order Magnus series expansion for orientation matrix
  for free rigid body (torque = 0)
  \param[in,out] QR orientation matrix
  \param[in,out] omega angular velocity
  \param[in] I inertia tensor
  \param[in] dt time step
 */
void propagate_w_RK4_Q_Euler(double QR[DIM][DIM], 
			     double omega[DIM],
			     const double I[DIM][DIM],
			     const double &dt);
				    
#endif
