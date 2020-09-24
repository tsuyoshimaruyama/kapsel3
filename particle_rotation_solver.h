/*!
  \file particle_rotation_solver.h
  \brief Solver routines for particle orientations
  \author J. Molina
  \date 2013/05/15
  \version 1.0
 */
#ifndef PARTICLE_ROTATION_SOLVER_H
#define PARTICLE_ROTATION_SOLVER_H

#include "rigid_body.h"
#include "variable.h"
/*!
  \brief Update particle orientations using the Euler
  forward method
  \details \f{align*}{
    \qtn{q}_i^{n+1} &= \qtn{q}_i^{n} + h \dot{\qtn{q}}^{n}\\
    \dot{q} &= \frac{1}{2}(0,\vec{\omega})\circ\qtn{q} = \frac{1}{2}\qtn{q}\circ(0, \vec{\omega}^\prime)
  \f}
  \param[in,out] p particle data
  \param[in] jikan time data
 */
inline void MD_solver_orientation_Euler(Particle &p, const double &dt) {
    quaternion dqdt;
    qtn_init(p.q_old, p.q);
    qdot(dqdt, p.q, p.omega, SPACE_FRAME);
    qtn_add(p.q, dqdt, dt);
    qtn_normalize(p.q);
}

/*!
  \brief Update particle orientations using Simo & Wang's second order method
  \details \f{align*}{
  \qtn{q}_i^{n+1} &= \qtn{q}_i^n +
  \frac{h}{4}\qtn{q}_i^{n}\circ(0,3\vec{\omega}_i^{\prime
  n}-\vec{\omega}_i^{\prime n-1})
  \f}
  Note that addition of angular velocity vectors only makes sense in
  body (primed) coordinates.
  \param[in,out] p particle data
  \param[in] jikan time data
  */

inline void MD_solver_orientation_SW2(Particle &p, const double &hdt) {
    double wb[DIM];
    double wb_old[DIM];
    // only add angular velocity vectors in body coordinates !
    rigid_body_rotation(wb, p.omega, p.q, SPACE2BODY);
    rigid_body_rotation(wb_old, p.omega_old, p.q_old, SPACE2BODY);
    for (int d = 0; d < DIM; d++) {
        wb[d] = 3.0 * wb[d] - wb_old[d];
    }

    quaternion dqdt;
    qtn_init(p.q_old, p.q);
    qdot(dqdt, p.q, wb, BODY_FRAME);
    qtn_add(p.q, dqdt, hdt);
    qtn_normalize(p.q);
}
/*!
  \brief Update particle orientations using a
  second-order Adams-Bashforth scheme
  \details \f{align*}{
  \qtn{q}_i^{n+1} &= \qtn{q}_i^n +
  \frac{h}{2}\left[ 3\dot{\qtn{q}}^{n} - \dot{\qtn{q}}^{n-1}\right]\\
  \dot{\qtn{q}} &= \frac{1}{2} \qtn{q}\circ (0, \vec{\omega}^\prime) = \frac{1}{2}(0,\omega)\circ\qtn{q}
  \f}
  where primes refer to the body coordinates.
  \param[in,out] p particle data
  \param[in] jikan time data
  */
inline void MD_solver_orientation_AB2(Particle &p, const double &hdt) {
    quaternion dqdt, dqdt_old;
    qdot(dqdt, p.q, p.omega, SPACE_FRAME);
    qdot(dqdt_old, p.q_old, p.omega_old, SPACE_FRAME);
    qtn_init(p.q_old, p.q);

    qtn_add(p.q, dqdt, 3.0 * hdt);
    qtn_add(p.q, dqdt_old, -hdt);
    qtn_normalize(p.q);
}

// Samuel Buss' second-order scheme
// untested
inline void MD_solver_orientation_SB2(Particle &p, const double &dt) {
    double wb[DIM];
    for (int d = 0; d < DIM; d++) {
        wb[d] = p.omega[d] + dt / 2.0 * IMOI[p.spec] * p.torque_hydro[d];
    }

    quaternion dqdt;
    qtn_init(p.q_old, p.q);
    qdot(dqdt, p.q, wb, SPACE_FRAME);
    qtn_add(p.q, dqdt, dt);
    qtn_normalize(p.q);
}

/*!
  \brief Solve Euler equations in body frame
  \f{align*}
  \dot{\omega}^x &= \left[\tau^x + \omega^y\omega^z\left(I^{yy}- I^{zz}\right)\right]/I^{xx}\\
  \dot{\omega}^y &= \left[\tau^y + \omega^z\omega^x\left(I^{zz}- I^{xx}\right)\right]/I^{yy}\\
  \dot{\omega}^z &= \left[\tau^z + \omega^x\omega^y\left(I^{xx}- I^{yy}\right)\right]/I^{zz}
  \f}
 */
inline void MD_solver_omega_Euler(double            omega[DIM],
                                  const double      torque[DIM],
                                  const double      Ib[DIM],
                                  const quaternion &q,
                                  const int         free_omega[DIM],
                                  const double &    dt) {
    double omega_b[DIM];
    double torque_b[DIM];
    rigid_body_rotation(omega_b, omega, q, SPACE2BODY);
    rigid_body_rotation(torque_b, torque, q, SPACE2BODY);

    double new_omega_b[DIM] = {omega_b[0], omega_b[1], omega_b[2]};

    // Warning: Fixed_omega flag now fixes angular velocity about rigid body frame
    //          Inconsistent with Fixed_vel flag which used lab frame
    if (free_omega[0]) new_omega_b[0] += (dt / Ib[0]) * (torque_b[0] + (Ib[1] - Ib[2]) * omega_b[1] * omega_b[2]);

    if (free_omega[1]) new_omega_b[1] += (dt / Ib[1]) * (torque_b[1] + (Ib[2] - Ib[0]) * omega_b[2] * omega_b[0]);

    if (free_omega[2]) new_omega_b[2] += (dt / Ib[2]) * (torque_b[2] + (Ib[0] - Ib[1]) * omega_b[0] * omega_b[1]);

    rigid_body_rotation(omega, new_omega_b, q, BODY2SPACE);
}
inline void MD_solver_omega_Euler_update(double            delta_omega[DIM],
                                         const double      omega[DIM],
                                         const double      torque[DIM],
                                         const double      Ib[DIM],
                                         const quaternion &q,
                                         const double &    dt) {
    double omega_b[DIM];
    double torque_b[DIM];
    rigid_body_rotation(omega_b, omega, q, SPACE2BODY);
    rigid_body_rotation(torque_b, torque, q, SPACE2BODY);

    delta_omega[0] = (dt / Ib[0]) * (torque_b[0] + (Ib[1] - Ib[2]) * omega_b[1] * omega_b[2]);
    delta_omega[1] = (dt / Ib[1]) * (torque_b[1] + (Ib[2] - Ib[0]) * omega_b[2] * omega_b[0]);
    delta_omega[2] = (dt / Ib[2]) * (torque_b[2] + (Ib[0] - Ib[1]) * omega_b[0] * omega_b[1]);

    rigid_body_rotation(delta_omega, q, BODY2SPACE);
}
/*!
  \brief Solve euler equations in body frame
*/
inline void MD_solver_omega_AB2(double            omega[DIM],
                                const double      torque[DIM],
                                const double      torque_old[DIM],
                                const double      Ib[DIM],
                                const quaternion &q,
                                const quaternion &q_old,
                                const int         free_omega[DIM],
                                const double &    hdt) {
    double omega_b[DIM];
    double torque_b[DIM];
    double torque_b_old[DIM];

    rigid_body_rotation(omega_b, omega, q, SPACE2BODY);
    rigid_body_rotation(torque_b, torque, q, SPACE2BODY);
    rigid_body_rotation(torque_b_old, torque_old, q_old, SPACE2BODY);

    double new_omega_b[DIM] = {omega_b[0], omega_b[1], omega_b[2]};

    // Warning: Fixed_omega flag now fixes angular velocity about rigid body frame
    //          Inconsistent with Fixed_vel flag which used lab frame
    if (free_omega[0])
        new_omega_b[0] +=
            (hdt / Ib[0]) * ((torque_b[0] + torque_b_old[0]) + 2.0 * (Ib[1] - Ib[2]) * omega_b[1] * omega_b[2]);

    if (free_omega[1])
        new_omega_b[1] +=
            (hdt / Ib[1]) * ((torque_b[1] + torque_b_old[1]) + 2.0 * (Ib[2] - Ib[0]) * omega_b[2] * omega_b[0]);

    if (free_omega[2])
        new_omega_b[2] +=
            (hdt / Ib[2]) * ((torque_b[2] + torque_b_old[2]) + 2.0 * (Ib[0] - Ib[1]) * omega_b[0] * omega_b[1]);

    rigid_body_rotation(omega, new_omega_b, q, BODY2SPACE);
}

#endif
