/*!
  \file resume.h
  \brief Routines to read/write restart file (header file)
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \todo documentation
 */

#ifndef RESUME_H
#define RESUME_H

#include "avs_output.h"

void Save_Restart_udf(
		      double **zeta
		      ,double *uk_dc
		      ,const Particle *p
		      ,const CTime &time
		      ,double **conc_k
		      );
void Force_restore_parameters(
			      double **zeta
			      ,double *uk_dc
			      ,Particle *p
			      ,CTime &time
			      ,double **conc_k
			      );


#endif
