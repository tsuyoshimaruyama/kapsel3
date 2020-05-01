#ifndef WALL_H
#define WALL_H
#include <assert.h>

#include "input.h"
#include "profile.h"

void Init_Wall(double* phi);
void Add_f_wall(Particle* p);
#endif
