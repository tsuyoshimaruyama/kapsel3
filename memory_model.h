#ifndef __MEMORY_MODEL_H__
#define __MEMORY_MODEL_H__
#include <assert.h>
#include <limits.h>
#include "variable.h"
#include "alloc.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __cplusplus
extern "C" {
#endif
    void Range_division (void);
    int Range_check (const Index_range *i_range, Index_range *o_range);
    int Range_coord (const int *i_mesh, int *o_mesh);
    void Mem_alloc_var (double **zeta);
    void Mem_free_var (double **zeta);
#ifdef __cplusplus
}
#endif
#endif /*__MEMORY_MODEL_H__*/
