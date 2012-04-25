/*!
  \file quaternion.h
  \brief Implements simple quaternion algebra
  \details Details on the notation can be found in 
  Graf, B. (2008). Quaternions and dynamics. 
  http://arxiv.org/abs/0811.2889
  \author J. Molina
  \date 2012/03/28
  \version 1.0
*/
#ifndef QUATERNION_H
#define QUATERNION_H

#include <assert.h>
#include <math.h>
#include "macro.h"
#include "parameter_define.h"

static double QTOL=1.0e-15;
static double QTOL_LARGE=1.0e-5;

typedef struct quaternion {
  double s; //scalar part
  double v[DIM]; //vector part
} quaternion;


inline double qtn_q0(const quaternion &q){
  return q.s;
}
inline double qtn_q1(const quaternion &q){
  return q.v[0];
}
inline double qtn_q2(const quaternion &q){
  return q.v[1];
}
inline double qtn_q3(const quaternion &q){
  return q.v[2];
}
inline void qtn_scalar(double &a, const quaternion &q){
  a = q.s;
}
inline void qtn_vector(double v[DIM], const quaternion &q){
  for(int i = 0; i < DIM; i++){
    v[i] = q.v[i];
  }
}

/*!
  \brief Initialize by specifying all four components
 */
inline void qtn_init(quaternion &q, const double &a0, const double &a1,
		     const double &a2, const double &a3){
  q.s = a0;
  q.v[0] = a1;
  q.v[1] = a2;
  q.v[2] = a3;
}

/*!
 \brief Initialize by specyfying scalar, vector parts
 */
inline void qtn_init(quaternion &q, const double &si, const double vi[DIM]){
  q.s = si;
  for(int i = 0; i < DIM; i++){
    q.v[i] = vi[i];
  }
}

/*!
  \brief Initialize by copy
 */
inline void qtn_init(quaternion &q, const quaternion &qa){
  q.s = qa.s;
  for(int i = 0; i < DIM; i++){
    q.v[i] = qa.v[i];
  }
}

/*!
  \brief Scale quaternion: q = alpha*qa
 */
inline void qtn_scale(quaternion &q, const quaternion &qa, 
		      const double scale){
  q.s = qa.s*scale;
  for(int i = 0; i < DIM; i++){
    q.v[i] = qa.v[i]*scale;
  }
}

/*!
  \brief Scale quaternion in place: q = alpha*q
 */
inline void qtn_scale(quaternion &q, const double scale){
  q.s *= scale;
  for(int i = 0; i < DIM; i++){
    q.v[i] *= scale;
  }
}

/*!
  \brief Add two quaternions: q = qa + alpha*qb
 */
inline void qtn_add(quaternion &q, const quaternion &qa, 
		    const quaternion &qb, const double alpha=1.0){
  q.s = qa.s + alpha * qb.s;
  for(int i = 0; i < DIM; i++){
    q.v[i] = qa.v[i] + alpha * qb.v[i];
  }
}

/*!
  \brief Add to quaternions in place: q = q + alpha*qb
 */
inline void qtn_add(quaternion &q, const quaternion &qb, 
		    const double alpha=1.0){
  q.s += alpha * qb.s;
  for(int i = 0; i < DIM; i++){
    q.v[i] += alpha * qb.v[i];
  }
}

/*!
  \brief Cuaternion conjugate
 */
inline void qtn_conj(quaternion &q, const quaternion &qa){
  q.s = qa.s;
  for(int i = 0; i < DIM; i++){
    q.v[i] = -qa.v[i];
  }
}

/*!
  \brief Cuaternion conjugate in place
 */
inline void qtn_conj(quaternion &q){
  for(int i = 0; i < DIM; i++){
    q.v[i] = -q.v[i];
  }
}

/*!
  \brief Multiply two quaternions: q = qa.(alpha*qb)
 */
inline void qtn_prod(quaternion &q, const quaternion &qa, 
		     const quaternion &qb,const double alpha=1.0){

  q.s = qa.s * qb.s;
  q.v[0] = qa.v[1] * qb.v[2] - qa.v[2] * qb.v[1];
  q.v[1] = qa.v[2] * qb.v[0] - qa.v[0] * qb.v[2];
  q.v[2] = qa.v[0] * qb.v[1] - qa.v[1] * qb.v[0];

  for(int i = 0; i < DIM; i++){
    q.s -= (qa.v[i] * qb.v[i]);

    q.v[i] += (qa.s * qb.v[i] + qb.s * qa.v[i]);
    q.v[i] *= alpha;
  }
  q.s *= alpha;
}

/*!
  \brief Multiply two quaternion in place: q = q.(alpha*qb)
 */
inline void qtn_prod(quaternion &q, const quaternion &qb, 
		     const double scale=1.0){
  quaternion qa;
  qtn_init(qa, q);
  qtn_prod(q, qa, qb, scale);
}

/*!
  \brief Compute quaternion norm
 */
inline double qtn_norm(const quaternion &q){
  return sqrt(q.s*q.s + q.v[0]*q.v[0] + q.v[1]*q.v[1] + q.v[2]*q.v[2]);
}

/*!
  \brief Compute norm of scalar part
 */
inline double qtn_norm_s(const quaternion &q){
  return ABS(q.s);
}

/*!
  \brief Compute norm of imaginary part
 */
inline double qtn_norm_v(const quaternion &q){
  return sqrt(q.v[0]*q.v[0] + q.v[1]*q.v[1] + q.v[2]*q.v[2]);
}

/*!
  \brief Compute quaternion square norm
 */
inline double qtn_sqnorm(const quaternion &q){
  return q.s*q.s + q.v[0]*q.v[0] + q.v[1]*q.v[1] + q.v[2]*q.v[2];
}

/*!
  \brief Check normal quaternion
 */
inline void qtn_isnormal(const quaternion &q, const double tol=QTOL){
  assert(ABS(qtn_norm(q) - 1.0) <= tol);
}

/*!
  \brief Normalize quaternion
 */
inline void qtn_normalize(quaternion &q){
  double qnorm = qtn_norm(q);
  if(qnorm <= 0.0){
    fprintf(stderr, "Erorr: qnorm = %.10f\n", qnorm);
    exit_job(EXIT_FAILURE);
  }
  qnorm = 1.0 / qnorm;
  qtn_scale(q, qnorm);
}

/*!
  \brief Compute quaternion inverse: q = alpha*qa^-1
 */
inline void qtn_inv(quaternion &q, quaternion qa, const double alpha=1.0){
  double q2i = qtn_sqnorm(qa);
  assert(q2i != 0.0);

  q2i = 1.0 / q2i;
  qtn_conj(qa);
  qtn_scale(q, qa, q2i*alpha);
}

/*!
  \brief Compute quaternion inverse in place: q = alpha*q^-1
 */
inline void qtn_inv(quaternion &q, const double alpha=1.0){
  double q2i = qtn_sqnorm(q);
  assert(q2i != 0.0);

  q2i = 1.0 / q2i;
  qtn_conj(q);
  qtn_scale(q, q2i*alpha);
}

/*!
  \brief Compare two quaternions
 */
inline int qtn_cmp(const quaternion &qa, const quaternion &qb, 
		    const double tol=QTOL_LARGE){
  quaternion qd;
  double dmax;

  qtn_add(qd, qa, qb, -1.0);
  dmax = (qb.s != 0.0 ? ABS(qd.s/qb.s) : ABS(qd.s));
  for(int i = 0; i < DIM; i++){
    dmax = (qb.v[i] != 0.0 ? 
	    MAX(dmax, ABS(qd.v[i]/qb.v[i])) : MAX(dmax, ABS(qd.v[i])));
  }

  return (dmax <= tol ? 1 : 0);
}

/*!
  \brief Compute rotation quaternion given angle and (normal) axis vector
 */
inline void rv_rqtn(quaternion &q, const double &phi, const double n[DIM]){
  double ds;
  double dni;
  
  dni = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  assert(dni > 0.0);
  if(ABS(dni - 1.0) > QTOL){
    printf("Warning: n must be normal vector, |n|= %.15f\n", dni);
  }
  dni = 1.0/dni;

  ds = sin(0.5*phi);
  q.s = cos(0.5*phi);
  for(int i = 0; i < DIM; i++){
    q.v[i] = ds * n[i] * dni;
  }
}


/*!
  \brief Compute rotation quaternion given rotation vector
  \details Rotation angle is given directly by the magnitude of the vector
 */
inline void rv_rqtn(quaternion &q, const double v[DIM]){
  double phi;
  double ds;
  double dni;

  phi = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  assert(phi > 0.0);
  dni = 1.0/phi;

  ds = sin(0.5*phi);
  q.s = cos(0.5*phi);
  for(int i = 0; i < DIM; i++){
    q.v[i] = ds * v[i] * dni;
  }
}

/*!
  \brief Compute rotation angle and (normal) vector for rotation quaternion
 */
inline void rqtn_rv(double &phi, double v[DIM], const quaternion &q){
  double ds;

  qtn_isnormal(q);
  phi = 2.0*asin(qtn_norm_v(q));
  qtn_vector(v, q);

  ds =  sin(0.5*phi);
  if(ds != 0.0){
    ds = 1.0/ds;
    for(int i = 0; i < DIM; i++){
      v[i] *= ds;
    }
  }else{
    phi = 0.0;
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 1.0;
  }
}

/*!
  \brief Compute random rotation quaternion (see Numerical Recipes)
 */
inline void random_rqtn(quaternion &q){
  double u0, u1, u2, u3, scale;
  RA_circle(u0, u1);
  RA_circle(u2, u3);
  scale = sqrt((1.0 - u0*u0 - u1*u1)/(u2*u2 + u3*u3));
  qtn_init(q, u3*scale, u0, u1, u2*scale);
  qtn_isnormal(q);
}

/*!
  \brief Compute rotation matrix given rotation quaternion
 */
inline void rqtn_rm(double R[DIM][DIM], const quaternion &q){

  R[0][0] = 1.0 - 2.0 * q.v[1] * q.v[1] - 2.0 * q.v[2] * q.v[2];
  R[0][1] = 2.0 * q.v[0] * q.v[1] - 2.0 * q.s * q.v[2];
  R[0][2] = 2.0 * q.s * q.v[1] + 2.0 * q.v[0] * q.v[2];

  R[1][0] = 2.0 * q.v[0] * q.v[1] + 2.0 * q.s * q.v[2];
  R[1][1] = 1.0 - 2.0 * q.v[0] * q.v[0] - 2.0 * q.v[2] * q.v[2];
  R[1][2] = -2.0 * q.s * q.v[0] + 2.0 * q.v[1] * q.v[2];

  R[2][0] = -2.0 * q.s * q.v[1] + 2.0 * q.v[0] * q.v[2];
  R[2][1] = 2.0 * q.s * q.v[0] + 2.0 * q.v[1] * q.v[2];
  R[2][2] = 1.0 - 2.0 * q.v[0] * q.v[0] - 2.0 * q.v[1] * q.v[1];
}


/*!
  \brief Compute transpose of rotation matrix given rotation quaternion
 */
inline void rqtn_rmt(double R_T[DIM][DIM], const quaternion &q){
  qtn_isnormal(q);

  R_T[0][0] = 1.0 - 2.0 * q.v[1] * q.v[1] - 2.0 * q.v[2] * q.v[2];
  R_T[0][1] = 2.0 * q.v[0] * q.v[1] + 2.0 * q.s * q.v[2];
  R_T[0][2] = -2.0 * q.s * q.v[1] + 2.0 * q.v[0] * q.v[2];

  R_T[1][0] = 2.0 * q.v[0] * q.v[1] - 2.0 * q.s * q.v[2];
  R_T[1][1] = 1.0 - 2.0 * q.v[0] * q.v[0] - 2.0 * q.v[2] * q.v[2];
  R_T[1][2] = 2.0 * q.s * q.v[0] + 2.0 * q.v[1] * q.v[2];

  R_T[2][0] = 2.0 * q.s * q.v[1] + 2.0 * q.v[0] * q.v[2];
  R_T[2][1] = -2.0 * q.s * q.v[0] + 2.0 * q.v[1] * q.v[2];
  R_T[2][2] = 1.0 - 2.0 * q.v[0] * q.v[0] - 2.0 * q.v[1] * q.v[1];
}

/*!
  \brief Compute rotation rate matrix E given quaternion
 */
inline void rqtn_rm_e(double E[DIM][DIM+1], const quaternion &q){
  qtn_isnormal(q);

  E[0][0] = -q.v[0];
  E[0][1] = q.s;
  E[0][2] = -q.v[2];
  E[0][3] = q.v[1];

  E[1][0] = -q.v[1];
  E[1][1] = q.v[2];
  E[1][2] = q.s;
  E[1][3] = -q.v[0];
  
  E[2][0] = -q.v[2];
  E[2][1] = -q.v[1];
  E[2][2] = q.v[0];
  E[2][3] = q.s;
}
/*!
  \brief Compute transpose of rotation rate matrix E given quaternion
 */
inline void rqtn_rm_et(double E_T[DIM+1][DIM], const quaternion &q){
  qtn_isnormal(q);

  E_T[0][0] = -q.v[0];
  E_T[0][1] = -q.v[1];
  E_T[0][2] = -q.v[2];

  E_T[1][0] = q.s;
  E_T[1][1] = q.v[2];
  E_T[1][2] = -q.v[1];

  E_T[2][0] = -q.v[2];
  E_T[2][1] = q.s;
  E_T[2][2] = q.v[0];

  E_T[3][0] = q.v[1];
  E_T[3][1] = -q.v[0];
  E_T[3][2] = q.s;
}

/*!
  \brief Compute rotation matrix G given quaternion
 */
inline void rqtn_rm_g(double gmatrix[DIM][DIM+1], const quaternion &q){
  qtn_isnormal(q);

  gmatrix[0][0] = -q.v[0];
  gmatrix[0][1] = q.s;
  gmatrix[0][2] = q.v[2];
  gmatrix[0][3] = -q.v[1];

  gmatrix[1][0] = -q.v[1];
  gmatrix[1][1] = -q.v[2];
  gmatrix[1][2] = q.s;
  gmatrix[1][3] = q.v[0];
  
  gmatrix[2][0] = -q.v[2];
  gmatrix[2][1] = q.v[1];
  gmatrix[2][2] = -q.v[0];
  gmatrix[2][3] = q.s;
}

/*!
  \brief Compute transpose of rotation matrix G given quaternion
 */
inline void rqtn_rm_gt(double gmatrix[DIM+1][DIM], const quaternion &q){
  qtn_isnormal(q);

  gmatrix[0][0] = -q.v[0];
  gmatrix[0][1] = -q.v[1];
  gmatrix[0][2] = -q.v[2];

  gmatrix[1][0] = q.s;
  gmatrix[1][1] = -q.v[2];
  gmatrix[1][2] = q.v[1];

  gmatrix[2][0] = q.v[2];
  gmatrix[2][1] = q.s;
  gmatrix[2][2] = -q.v[0];

  gmatrix[3][0] = -q.v[1];
  gmatrix[3][1] = q.v[0];
  gmatrix[3][2] = q.s;
}


#endif
