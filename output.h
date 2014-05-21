/*!
  \file resume.h
  \brief Routines to write output file (header file)
  \author J. Molina
  \date 2014/05/21
  \version 1.0
 */
#ifndef OUTPUT_H
#define OUTPUT_H

#include "input.h"
#include "variable.h"
#include "rigid.h"
#include "fft_wrapper.h"

//UDF OUTPUT WRITING
//

/*!
  \brief Output particle data for current configuration in UDF format
 */
void Output_udf(UDFManager *ufout
		,double **zeta
		,double *uk_dc
		,const Particle *p
		,const CTime &time);

#endif
