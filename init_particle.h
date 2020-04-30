/*!
  \file init_particle.h
  \brief Initialize particle properties (header file)
  \details Initializes positions, velocities, forces, torques on particles
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */

#ifndef INIT_PARTICLE_H
#define INIT_PARTICLE_H

#include <assert.h>

#include "avs_output.h"
#include "fluct.h"
#include "input.h"
#include "macro.h"
#include "md_force.h"
#include "rigid.h"
#include "variable.h"

void Init_Particle(Particle *p);
void Init_Chain(Particle *p);
void Init_Rigid(Particle *p);
void Show_parameter(Particle *p);

inline void Show_particle(Particle *p) {
    for (int n = 0; n < Particle_Number; n++) {
        fprintf(stderr,
                "%g %g %g %g %g %g\n",
                p[n].v[0],
                p[n].v[1],
                p[n].v[2],
                p[n].omega[0],
                p[n].omega[1],
                p[n].omega[2]);
    }
    fprintf(stderr, "\n\n");
}

/*!
    \brief Compute particle dipole in space frame
*/
extern void compute_particle_dipole_standard(double *mu_space, const double *mu_body, quaternion &q);

/*!
    \brief Compute particle dipole in space frame for the special case of a Quincke roller
    \details Assume that the particle dipole is always perpendicular to the plane formed by the anqular velocity vector
   \f$\vec{e}_{\omega}\f$ and the electric field direction \f$\vec{n}\f$

   \f{align*}{
     \vec{p} &=   p_0 \vec{n}\times\vec{e}_{\omega}
   \f}

   with \f$p_0\f$ the magnitude. It is thus parallel to the quincke torque \f$\tau_{\textrm{Quincke}}\f$.

    \warning We are assuming that the magnitude of the dipole is passed in the first component of mu_body
    \param[out] mu_space particle dipole in the lab frame
    \param[in]  mu_body the magnitude of the dipole (p_0, 0, 0)
    \param[in]  q orientation quaternion
*/
extern void compute_particle_dipole_quincke(double *mu_space, const double *mu_body, quaternion &q);

/*!
    \brief Generic funtion pointer used to compute particle dipole
*/
extern void (*compute_particle_dipole)(double *mu_space, const double *mu_body, quaternion &q);
#endif
