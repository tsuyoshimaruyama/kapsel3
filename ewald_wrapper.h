#ifndef EWALD_WRAPPER_H
#define EWALD_WRAPPER_H

#include <stdio.h>
#include <time.h>
#include "ewald.h"

typedef struct EwaldParams {
  // input parameters
  double alpha = 8.0;
  double delta = 1.0e-16;
  double conf = 0.51;
  double epsilon = -0.01;
  void init(double _alpha, double _delta, double _conf, double _eps) {
    alpha = _alpha;
    delta = _delta;
    conf = _conf;
    epsilon = _eps;
  }
} EwaldParams;
typedef struct EwaldMem {
  int num;
  int *group_id;
  double *dval;

  double **r;
  double *q;
  double **mu;

  double **force;
  double **torque;
  double **efield;
  double ***efield_grad;
  double energy[5];
  bool charge;
  bool dipole;
  void init(const int &_num, const bool &_charge, const bool &_dipole) {
    num = _num;
    charge = _charge;
    dipole = _dipole;
    group_id = (int *)alloc_1d_int(num);
    r = (double **)alloc_2d_double(num, DIM);
    dval = (double *)alloc_1d_double(num);
    q = (charge ? (double *)alloc_1d_double(num) : NULL);
    mu = (dipole ? (double **)alloc_2d_double(num, DIM) : NULL);
    {
      for (int i = 0; i < num; i++)
        group_id[i] = -(i + 1);
      if (charge) {
        for (int i = 0; i < num; i++)
          q[i] = 0.0;
      }
      double *rr = r[0];
      for (int i = 0; i < num * DIM; i++)
        rr[i] = 0.0;
      if (dipole) {
        double *pp = mu[0];
        for (int i = 0; i < num * DIM; i++)
          pp[i] = 0.0;
      }
    }

    force = (double **)alloc_2d_double(num, DIM);
    torque = (double **)alloc_2d_double(num, DIM);
    efield = (double **)alloc_2d_double(num, DIM);
    efield_grad = (double ***)alloc_3d_double(num, DIM, DIM);
  }
  void free() {
    free_1d_int(group_id);

    if (charge) free_1d_double(q);

    free_2d_double(r);
    if (dipole) free_2d_double(mu);

    free_2d_double(force);
    free_2d_double(torque);
    free_2d_double(efield);
    free_3d_double(efield_grad);
  }
} EwaldMem;

extern EwaldParams ewald_param;
extern EwaldMem ewald_mem;
extern parallelepiped *ewald_cell;
extern ewald *ewald_sum;

void free_ewald_sum();

void init_ewald_sum(const double &lx, 
                    const double &ly, 
                    const double &lz, 
                    const int &num, 
                    const bool &charge, 
                    const bool &dipole);
void compute_ewald_sum();

#endif
