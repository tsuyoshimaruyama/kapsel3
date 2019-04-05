#ifndef __OPERATE_MPI_PARTICLE_H__
#define __OPERATE_MPI_PARTICLE_H__
#ifdef _MPI
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "variable.h"
#include "alloc.h"
#include "operate_mpi.h"
#ifdef __cplusplus
extern "C" {
#endif
    void Build_Particle_Datatypes (void);
    void Particle_Bcast (Particle *p);
    void Particle_Gather (Particle *p, Particle *new_p, enum SW_SWITCH sw);
    void Particle_Group_Communication (Particle *p, enum SW_COMM sw_comm);
    void Particle_qsort (Particle *p, int n);
    void Particle_Sendrecv (Particle *p, Particle *new_p, enum SW_COMM sw_comm);
#ifdef __cplusplus
}
#endif
#endif /* _MPI */
#endif /* __OPERATE_MPI_H__ */
