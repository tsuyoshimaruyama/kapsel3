/*!
  \file ewald_wrapper.cxx
  \author John J. Molina
  \date 2014/08/18
  \version 1.0
  \brief Wrapper to compute ewald forces, fields, field gradients, etc.
 */
#include "ewald_wrapper.h"

EwaldParams     ewald_param;
EwaldMem        ewald_mem;
parallelepiped *ewald_cell;
ewald *         ewald_sum;

void free_ewald_sum() {
    ewald_mem.free();
    delete ewald_cell;
}

void init_ewald_sum(const double &lx, const double &ly, const double &lz, const int &num) {
    double a[DIM] = {lx, 0.0, 0.0};
    double b[DIM] = {0.0, ly, 0.0};
    double c[DIM] = {0.0, 0.0, lz};
    ewald_cell    = new parallelepiped(a, b, c);
    ewald_mem.init(num);
    ewald_sum = new ewald(ewald_cell,
                          ewald_param.alpha,
                          ewald_param.epsilon,
                          ewald_param.delta,
                          ewald_param.conv,
                          ewald_mem.num,
                          ewald_param.charge,
                          ewald_param.dipole);
    ewald_sum->define_groups(ewald_mem.group_id);
    for (int i = 0; i < ewald_mem.num; i++) {
        ewald_mem.efield[i][0] = ewald_mem.efield[i][1] = ewald_mem.efield[i][2] = 0.0;
    }
}

// Before running : make sure r <- positions
//                            q <- charges
//                            mu<- dipoles
void compute_ewald_sum() {
    if (!ewald_param.enabled) return;
    double *dmy_q  = (ewald_param.charge ? ewald_mem.q : nullptr);
    double *dmy_mu = (ewald_param.dipole ? ewald_mem.mu[0] : nullptr);
    // Calculation
    ewald_sum->reset_boundary(ewald_param.epsilon);
    ewald_sum->compute(ewald_mem.energy,
                       ewald_mem.force[0],
                       ewald_mem.torque[0],
                       ewald_mem.efield[0],
                       ewald_mem.efield_grad[0][0],
                       ewald_mem.r[0],
                       dmy_q,
                       dmy_mu,
                       ewald_mem.num);
}
