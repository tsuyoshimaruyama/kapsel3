#include "rigid_body.h"
/*!
  \file rigid_body.cxx
  \brief Routines to solve equations of motion for orientation of a 
  symmetrical sphere (in 3D space)
  \author J. Molina
  \date 2012/04/02
  \version 1.0
*/

#define ONE 1.0
#define ZERO 0.0
//
// Get skew anti-symmetric matrix
inline void skew(double ws[DIM][DIM], 
		 const double w[DIM]){
  assert(DIM == 3);
  ws[0][0] = ZERO;
  ws[0][1] = -w[2];
  ws[0][2] = w[1];
  ws[1][0] = w[2];
  ws[1][1] = ZERO;
  ws[1][2] = -w[0];
  ws[2][0] = -w[1];
  ws[2][1] = w[0];
  ws[2][2] = ZERO;
}

//
// Make sure x is in body coordinates
// (orientation matrix)
void ensure_body_coordinates(double x[DIM],
			     const double QR[DIM][DIM],
			     const COORD_SYSTEM &coord){
  if(coord == BODY_FRAME){
    return;
  }else if(coord == SPACE_FRAME){
    v_M_prod(x, QR);
  }else {
    fprintf(stderr, "Error: unknown coordinate system\n");
    exit_job(EXIT_FAILURE);
  }
}
// (orientation quaternion)
void ensure_body_coordinates(double x[DIM],
			     const quaternion &q,
			     const COORD_SYSTEM &coord){
  if(coord == BODY_FRAME){
    return;
  }else if(coord == SPACE_FRAME){
    quaternion qx;
    quaternion q_dmy;
    qtn_init(qx, ZERO, x);
    qtn_conj(q_dmy, q);

    qtn_prod(qx, q);
    qtn_prod(q_dmy, qx);
    qtn_vector(x, q_dmy);
  }else{
    fprintf(stderr, "Error: unknown coordinate system\n");
    exit_job(EXIT_FAILURE);
  }
}

//
// Make sure x is in space coordinates
// (orientation  matrix)
void ensure_space_coordinates(double x[DIM],
			      const double QR[DIM][DIM],
			      const COORD_SYSTEM &coord){
  if(coord == SPACE_FRAME){
    return;
  }else if(coord == BODY_FRAME){
    M_v_prod(x, QR);
  }else{
    fprintf(stderr, "Error: unknown coordinate system\n");
    exit_job(EXIT_FAILURE);
  }
}
// (orientation quaternion)
void ensure_space_coordinates(double x[DIM],
			      const quaternion &q,
			      const COORD_SYSTEM &coord){
  if(coord == SPACE_FRAME){
    return;
  }else if(coord == BODY_FRAME){
    quaternion qx;
    quaternion q_dmy;
    qtn_init(qx, ZERO, x);
    qtn_conj(q_dmy, q);

    qtn_prod(qx, q_dmy);
    qtn_prod(q_dmy, q, qx);
    qtn_vector(x, q_dmy);
  }else{
    fprintf(stderr, "Error: unknown coordinate system\n");
    exit_job(EXIT_FAILURE);
  }
}

// Perform rigid body rotation
// (quaternion)
void rigid_body_rotation(double rotated[DIM], 
			 const double original[DIM], 
			 const quaternion &q, 
			 const COORD_TRANS &transform){
  quaternion qx;
  quaternion q_dmy;

  qtn_init(qx, ZERO, original);
  qtn_conj(q_dmy, q);
  
  if(transform == BODY2SPACE){
    qtn_prod(qx, q_dmy);
    qtn_prod(q_dmy, q, qx);

  }else if(transform == SPACE2BODY){
    qtn_prod(qx, q);
    qtn_prod(q_dmy, qx);

  }else{
    fprintf(stderr, "Error: unknown rigid body transformation\n");
    exit_job(EXIT_FAILURE);
  }

  qtn_vector(rotated, q_dmy);
}
inline void rigid_body_rotation(double rotated[DIM], 
				const quaternion &q, 
				const COORD_TRANS &transform){
  double original[DIM];
  v_copy(original, rotated);
  rigid_body_rotation(rotated, original, q, transform);
}

// Perform rigid body rotation
// (rotation matrix)
void rigid_body_rotation(double rotated[DIM], 
			 const double original[DIM],
			 const double QR[DIM][DIM],
			 const COORD_TRANS &transform){
  if(transform == BODY2SPACE){
    M_v_prod(rotated, QR, original);
  }else if(transform == SPACE2BODY){
    v_M_prod(rotated, original, QR);
  }else{
    fprintf(stderr, "Error: unknown rigid body transformation\n");
    exit_job(EXIT_FAILURE);
  }
}
inline void rigid_body_rotation(double rotated[DIM],
				const double QR[DIM][DIM],
				const COORD_TRANS &transform){
  double original[DIM];
  v_copy(original, rotated);
  rigid_body_rotation(rotated, original, QR, transform);
}



//
// time derivative of orientation quaternion
void qdot(quaternion &dqdt, 
	  const quaternion &q,
	  const double omega[DIM], 
	  const COORD_SYSTEM &coord){
  quaternion qw;
  qtn_init(qw, ZERO, omega);
  if(coord == SPACE_FRAME){
    qtn_prod(dqdt, qw, q, 0.5);
  } else if(coord == BODY_FRAME){
    qtn_prod(dqdt, q, qw, 0.5);
  } else {
    fprintf(stderr, "Error: unknown coordinate system\n");
    exit_job(EXIT_FAILURE);
  }
}
inline void qdot(quaternion &dqdt,
		 const double omega[DIM],
		 const COORD_SYSTEM &coord){
  quaternion q;
  qtn_init(q, dqdt);
  qdot(dqdt, q, omega, coord);
}

//
// time derivative of orientation matrix
void Qdot(double dQRdt[DIM][DIM],
	  const double QR[DIM][DIM],
	  const double omega[DIM],
	  const COORD_SYSTEM &coord){
  
  double omega_body[DIM];
  double omega_skew[DIM][DIM];

  if(coord == SPACE_FRAME){
    rigid_body_rotation(omega_body, omega, QR, SPACE2BODY);
    skew(omega_skew, omega_body);
  } else if (coord == BODY_FRAME){
    skew(omega_skew, omega);
  } else{
    fprintf(stderr, "Error: unknown coordinate system\n");
    exit_job(EXIT_FAILURE);
  }
  M_prod(dQRdt, QR, omega_skew);
}
inline void Qdot(double dQRdt[DIM][DIM],
		 const double omega[DIM],
		 const COORD_SYSTEM &coord){
  double QR[DIM][DIM];
  M_copy(QR, dQRdt);
  Qdot(dQRdt, QR, omega, coord);
}

// time derivative of angular velocity
void wdot(double dwdt[DIM], 
	  const double w[DIM],
	  const double tau[DIM],
	  const double I[DIM][DIM]){
  
  double dmy[DIM];
  double Iinv[DIM][DIM];
  M_inv(Iinv, I);
  M_v_prod(dmy, I, w);
  v_cross(dmy, w);
  v_add(dmy, tau);
  M_v_prod(dwdt, Iinv, dmy);
}
inline void wdot(double dwdt[DIM],
		 const double tau[DIM],
		 const double I[DIM][DIM]){
  double w[DIM];
  v_copy(w, dwdt);
  wdot(dwdt, w, tau, I);
}



// Rotation operator using rodrigues' formula
void rodrigues_rotation_formula(double T[DIM][DIM], 
				const double phi, 
				const double n[DIM]){
  double Id[DIM][DIM] = {{ONE, ZERO, ZERO}, 
			 {ZERO, ONE, ZERO},
			 {ZERO, ZERO, ONE}};
  double V[DIM][DIM];

  double cs = ONE - cos(phi);
  double sn = sin(phi);
  skew(V, n);
  for(int i = 0; i < DIM; i++){
    for(int j = 0; j < DIM; j++){
      T[i][j] = Id[i][j] + sn*V[i][j] +
	cs * (n[i]*n[j] - Id[i][j]);
    }
  }
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
/// TEST ROUTINES
//////////////////////////////////////////////////

void propagate_wq_RK4(quaternion &q,
		      double omega[DIM],
		      const double I[DIM][DIM],
		      const double &dt){
  quaternion kq1, kq2, kq3, kq4, dmy_q;
  double kw1[DIM], kw2[DIM], kw3[DIM], kw4[DIM], dmy_w[DIM];
  double tau[DIM] = {ZERO, ZERO, ZERO};

  qtn_init(dmy_q, q);
  v_copy(dmy_w, omega);
  qdot(kq1, dmy_q, dmy_w, BODY_FRAME);
  wdot(kw1, dmy_w, tau, I);
  qtn_scale(kq1, dt);
  v_scale(kw1, dt);

  qtn_add(dmy_q, q, kq1, 0.5); 
  v_add(dmy_w, omega, kw1, 0.5); 
  qdot(kq2, dmy_q, dmy_w, BODY_FRAME);
  wdot(kw2, dmy_w, tau, I);
  qtn_scale(kq2, dt);
  v_scale(kw2, dt);

  qtn_add(dmy_q, q, kq2, 0.5);
  v_add(dmy_w, omega, kw2, 0.5);
  qdot(kq3, dmy_q, dmy_w, BODY_FRAME);
  wdot(kw3, dmy_w, tau, I);
  qtn_scale(kq3, dt);
  v_scale(kw3, dt);

  qtn_add(dmy_q, q, kq3);
  v_add(dmy_w, omega, kw3);
  qdot(kq4, dmy_q, dmy_w, BODY_FRAME);
  wdot(kw4, dmy_w, tau, I);
  qtn_scale(kq4, dt);
  v_scale(kw4, dt);

  qtn_add(dmy_q, kq1, kq2, 2.0);
  qtn_add(dmy_q, kq3, 2.0);
  qtn_add(dmy_q, kq4);
  qtn_add(q, dmy_q, 1.0/6.0);
  qtn_normalize(q);

  v_add(dmy_w, kw1, kw2, 2.0);
  v_add(dmy_w, kw3, 2.0);
  v_add(dmy_w, kw4);
  v_add(omega, dmy_w, 1.0/6.0);
}

void propagate_w_RK4_Q_Euler(double QR[DIM][DIM], 
			     double omega[DIM],
			     const double I[DIM][DIM],
			     const double &dt){
  double rot_angle;
  double rot_vector[DIM];
  double rot_operator[DIM][DIM];
  double dmy_QR[DIM][DIM];
  double kw1[DIM], kw2[DIM], kw3[DIM], kw4[DIM], dmy_w[DIM];
  double tau[DIM] = {ZERO, ZERO, ZERO};

  // orientation update
  v_copy(rot_vector, omega);
  M_copy(dmy_QR, QR);


  rot_angle = v_norm(rot_vector);
  v_scale(rot_vector, 1.0/rot_angle);
  rot_angle *= dt;


  rodrigues_rotation_formula(rot_operator, rot_angle, rot_vector);
  M_prod(QR, dmy_QR, rot_operator);


  // omega update
  v_copy(dmy_w, omega);
  wdot(kw1, dmy_w, tau, I);
  v_scale(kw1, dt);


  v_add(dmy_w, omega, kw1, 0.5);
  wdot(kw2, dmy_w, tau, I);
  v_scale(kw2, dt);


  v_add(dmy_w, omega, kw2, 0.5);
  wdot(kw3, dmy_w, tau, I);
  v_scale(kw3, dt);


  v_add(dmy_w, omega, kw3);
  wdot(kw4, dmy_w, tau, I);
  v_scale(kw4, dt);


  v_add(dmy_w, kw1, kw2, 2.0);
  v_add(dmy_w, kw3, 2.0);
  v_add(dmy_w, kw4);
  v_add(omega, dmy_w, 1.0/6.0);

}

