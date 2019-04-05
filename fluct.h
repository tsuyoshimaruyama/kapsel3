/*!
  \file fluct.h
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Routines to compute random thermal fluctuation forces (header file)
  \todo documentation
 */
#ifndef FLUCT_H
#define FLUCT_H
#ifdef _MPI
#include <mpi.h>
#endif
#include <stdio.h>
#include <math.h>
#include "input.h"
#include "variable.h"
#include "operate_omega.h"
#include "make_phi.h"
#include "particle_solver.h"
#ifdef _MPI
#include "operate_mpi_particle.h"
#endif
#if !defined (NDEBUG)
#include "dc.h"
#else
//#include "SFMT.h"
#endif
#ifdef __cplusplus
extern "C" {
#endif
#include "mt19937ar.h"

  //void init_genrand(unsigned long s);
  //double genrand_real3(void);

#ifdef __cplusplus
}
#endif

void Add_random_force_thermostat(Particle *p, const CTime &jikan);

inline void MT_seed(const int &SW_seed, const unsigned long &seed){
    if(SW_seed == RANDOM_SEED){
        init_genrand(time(NULL));
        fprintf_single(stderr, "# MT_seed= time(NULL)\n");
    }else if(SW_seed == GIVEN_SEED){
        init_genrand(seed);
        fprintf_single(stderr, "# MT_seed= %lu\n",seed);
    }else {
        fprintf_single(stderr, "MT_seed(): invalid SW_seed.\n");
        exit_job(EXIT_FAILURE);
    }
}

inline void Gauss2(double random[2]){
    static double x1,x2;
    static double rsq, factor;
    do{
        x1 = 2.*genrand_real3()-1.;
        x2 = 2.*genrand_real3()-1.;
        rsq = x1 * x1 + x2 * x2;
    }while( rsq >= 1.0 || rsq == 0.0);
    factor = sqrt(-2.*log(rsq)/rsq);
    random[0] = factor * x1;
    random[1] = factor * x2;
}

void generate_random_number_MT(double* val, int num);

void Force_random_walk(Particle *p);

void Random_Walk(Particle *p, const double dr_factor = 0.5e-1*.5);

void Steepest_descent(Particle *p, const double dr_factor = 0.5e-1);

#endif
