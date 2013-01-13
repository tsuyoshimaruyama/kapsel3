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
#include "lad3.h"

enum COORD_SYSTEM {BODY_FRAME, SPACE_FRAME};
enum COORD_TRANS {BODY2SPACE, SPACE2BODY};

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
/*void wdot(double dwdt[DIM],
	  const double w[DIM],
	  const double tau[DIM],
	  const double I[DIM][DIM]);
*/
/*!
  \brief Compute time derivative of angular velocity (in place)
 */
/*inline void wdot(double dwdt[DIM],
		 const double tau[DIM],
		 const double I[DIM][DIM]);
*/
/*!
  \brief Update orientation/ angular velocity using fourth order RK scheme
  of a free rigid body (torque=0)
  \param[in,out] q orientation quaternion
  \param[in,out] omega angular velocity (BODY FRAME)
  \param[in] I intertia tensor
  \param[in] dt time step
 */
/*void propagate_wq_RK4(quaternion &q, 
		      double omega[DIM],
		      const double I[DIM][DIM],
		      const double &dt);
*/
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
/*void propagate_w_RK4_Q_Euler(double QR[DIM][DIM], 
			     double omega[DIM],
			     const double I[DIM][DIM],
			     const double &dt);
*/
#endif
