/*!
  \file operate_surface.h
  \brief Routines to control slip velocity at particle fluid boundaries
  \author J. Molina
  \date 22/06/2012
  \version 1.0
 */
#ifndef OPERATE_SURFACE_H
#define OPERATE_SURFACE_H

#include "variable.h"
#include "make_phi.h"
#include "particle_solver.h"

/*!
  \brief Enforce fluid velocity slip at particle boundary
  \details Velocity profile is enforced at the particle boundary with particle
  velocity modified to ensure momentum conservation
  \warning Particle velocites are modified!
  \param[in] u fluid velocity
  \param[out] up fluid velocity update to enforce slip boundaries
  \param[in,out] p particle data (updated to ensure momentum conservation)
  \param[in] jikan time step data (needed to update particle velocities)
 */
void Make_f_slip_particle(double const* const *u,
			  double **up,
			  Particle *p);
#endif
