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
#include "avs_output.h"
#include "avs_output_p.h"
#include "output_writer.h"

extern output_writer *writer;

void Init_output();
void Free_output();
void Show_output_parameter();

void Output_open_frame();
void Output_close_frame();
void Output_field_data(double** zeta,
		       double* uk_dc,
		       Particle* p,
		       const CTime &time);

void Output_charge_field_data(double** zeta,
			      double* uk_dc,
			      double** Concentration,
			      Particle* p,
			      const CTime &time);

void Output_particle_data(Particle*p, 
			 const CTime &time);


/*!
  \brief Output particle data for current configuration in UDF format
*/
void Output_udf(UDFManager *ufout,
		  double **zeta,
		  double *uk_dc,
		  const Particle *p,
		  const CTime &time);

#endif
