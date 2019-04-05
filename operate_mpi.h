#ifndef __OPERATE_MPI_H__
#define __OPERATE_MPI_H__
#ifdef _MPI
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "variable.h"
#include "alloc.h"
#ifndef NDEBUG
#define MPI_ERRORCHECK(errno)                                                                  \
    do {                                                                                       \
        if ((errno) && !MPI_Error_string ((errno), errstr, &errstrlen)) {                      \
            fprintf (stderr, "MPI Error:%s (file:%s, line:%d)\n", errstr, __FILE__, __LINE__); \
            fflush (stderr);                                                                   \
            MPI_Barrier (MPI_COMM_WORLD);                                                      \
            exit (EXIT_FAILURE);                                                               \
        }                                                                                      \
    } while(0)
#else
#define MPI_ERRORCHECK(errno) 
#endif
#ifdef __cplusplus
extern "C" {
#endif
    int Set_MPI_initalize (void);
    int Set_MPI_finalize (void);
    void Set_ID_Table (void);
    void Get_mesh_array (double *src, double *dist, enum SW_SWITCH sw);
    void Set_mesh_array (double *dist, double *src);
#ifdef __cplusplus
}
#endif
#endif /* _MPI */
#endif /* __OPERATE_MPI_H__ */
