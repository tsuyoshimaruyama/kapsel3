#ifndef EWALD_MULTIPOLE_H
#define EWALD_MULTIPOLE_H

#include "rigid_body_rotation.h"
#include "lad3.h"

void init_ewald(double *a, double *b, double *c,
		double &rcut, int &kcut, double &ewald_alpha,
		int &Nparticle);

void ewald_multipole_interaction(double *E_ewald,  
				 Particle *p);
#endif
