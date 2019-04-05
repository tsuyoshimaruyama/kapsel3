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
#ifdef _MPI
#include <mpi.h>
#endif
#include <assert.h> 
#include "macro.h"
#include "variable.h"
#include "md_force.h"
#include "avs_output.h"
#include "input.h"
#include "fluct.h"
#include "rigid.h"
#ifdef _MPI
#include "operate_mpi_particle.h"
#endif
#include "dc.h"

void MT_Init_Particle (Particle *p);
void Init_Particle(Particle *p);
void Init_Chain(Particle *p);
void Init_Rigid(Particle *p);
void Show_parameter(Particle *p);

/* change: ouput data only in RANK ZERO */
inline void Show_particle(Particle *p){
#ifdef _MPI
    Particle_Gather (p, p_tmp, SW_OFF);
    if (procid == root) {
        Particle_qsort (p_tmp, Particle_Number);
    }
#endif
    if (procid == root) {
        for (int n = 0; n < Particle_Number; n++) {
            fprintf_single (stderr, "%g %g %g %g %g %g\n",
                     p_tmp[n].v[0], p_tmp[n].v[1], p_tmp[n].v[2],
                     p_tmp[n].omega[0], p_tmp[n].omega[1], p_tmp[n].omega[2]);
        }
    }
    fprintf_single (stderr, "\n\n");
}
#endif

