/*!
  \file avs_output.h
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Output routines for field data in AVS/Express format (header file)
 */
#ifndef AVS_OUTPUT_H
#define AVS_OUTPUT_H
#ifdef _MPI
#include <mpi.h>
#endif
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "variable.h"
#include "macro.h"
#include "input.h"
#include "fluid_solver.h"

#ifdef _MPI
#include "operate_mpi_particle.h"
#endif

enum AVS_Field {
    irregular, uniform, rectilinear
};
//////////////////////////////


/*!
  \brief Print avs parameters to stderr
 */
void Show_avs_parameter();

/*!
  \brief Write AVS .fld definition files
 */
void Init_avs(const AVS_parameters &Avs_parameters);

/*!
  \brief Set global AVS parameters
 */
void Set_avs_parameters(AVS_parameters &Avs_parameters);

/*!
  \brief Output field data for current configuration in AVS format
 */
void Output_avs(AVS_parameters &Avs_parameters, 
		double **u, 
		double *phi, 
		double *Pressure,
		double **strain,
		const CTime &time);


/*!
  \brief Output field data of charged system for current configuration
  in AVS format
 */
void Output_avs_charge(AVS_parameters &Avs_parameters, 
		       double** u, 
		       double* phi, 
		       double* colloid_charge, 
		       double* solute_charge_total,
		       double* potential,
		       const CTime &time);

#endif
