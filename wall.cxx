#include "wall.h"

void Init_Wall(double* phi_wall) {
  {  // print wall params
    char              dmy[128];
    static const char axis[] = {'x', 'y', 'z'};
    sprintf(dmy,
            "phi_flatwall_%dx%dx%d_h%c%d.dat",
            Ns[0],
            Ns[1],
            Ns[2],
            axis[wall.axis],
            static_cast<int>((wall.hi - wall.lo) / DX));

    FILE*  fwall = filecheckopen(dmy, "w");
    double dr    = DX / 3.0;
    double l     = L[wall.axis];
    double hl    = l / 2.0;
    double r     = 0.0;
    double phiw  = 0.0;
    while (r <= l) {
      phiw = (r < hl ? Phi(r, wall.lo) : 1.0 - Phi(r, wall.hi));
      fprintf(fwall, "%.5f %.5f\n", r, phiw);
      r += dr;
    }
    fclose(fwall);
  }

  {  // compute wall phi field
    double  rijk[DIM] = {0.0, 0.0, 0.0};
    double& r         = rijk[wall.axis];
    double  hl        = HL[wall.axis];
    for (int i = 0; i < NX; i++) {
      rijk[0] = static_cast<double>(i) * DX;

      for (int j = 0; j < NY; j++) {
        rijk[1] = static_cast<double>(j) * DX;

        for (int k = 0; k < NZ; k++) {
          rijk[2]      = static_cast<double>(k) * DX;
          int im       = (i * NY * NZ_) + (j * NZ_) + k;
          phi_wall[im] = (r < hl ? Phi(r, wall.lo) : 1.0 - Phi(r, wall.hi));
        }
      }
    }
  }
}
