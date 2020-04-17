#include "ewald_wrapper.h"

EwaldParams     ewald_param;
EwaldMem        ewald_mem;
parallelepiped *ewald_cell;
ewald *         ewald_sum;

void free_ewald_sum() {
    ewald_mem.free();
    delete ewald_cell;
}

void init_ewald_sum(const double &lx,
                    const double &ly,
                    const double &lz,
                    const int &   num,
                    const bool &  charge,
                    const bool &  dipole) {
    double a[DIM] = {lx, 0.0, 0.0};
    double b[DIM] = {0.0, ly, 0.0};
    double c[DIM] = {0.0, 0.0, lz};
    ewald_cell    = new parallelepiped(a, b, c);
    ewald_mem.init(num, charge, dipole);
    ewald_sum = new ewald(ewald_cell,
                          ewald_param.alpha,
                          ewald_param.epsilon,
                          ewald_param.delta,
                          ewald_param.conf,
                          ewald_mem.num,
                          ewald_mem.charge,
                          ewald_mem.dipole);
    ewald_sum->define_groups(ewald_mem.group_id);
    for (int rigidID = 0; rigidID < num; rigidID++) {
        // ewald_mem.q[rigidID]    = 0.0;   //charge
        ewald_mem.dval[rigidID]      = 1.0;  // dipole->1.0 charge->0.0
        ewald_mem.efield[rigidID][0] = ewald_mem.efield[rigidID][1] = ewald_mem.efield[rigidID][2] = 0.0;
    }
}

// Before running : make sure r <- positions
//                            q <- charges
//                            mu<- dipoles
void compute_ewald_sum() {
    double *dmy_q  = (ewald_mem.charge ? ewald_mem.q : NULL);
    double *dmy_mu = (ewald_mem.dipole ? ewald_mem.mu[0] : NULL);
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
