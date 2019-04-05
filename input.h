/*!
  \file input.h
  \brief Read udf input file to start simulation (header file)
  \author Y. Nakayama
  \date 2006/11/14
  \version 1.3
 */
#ifndef INPUT_H
#define INPUT_H
#ifdef _MPI
#include <mpi.h>
#endif
#include <stdio.h>
#include <iostream>
#include <string>
#include <cstring>
#include <cfloat>
#include <sys/types.h>
#include <assert.h>
#include "macro.h"
#include "alloc.h"
#include "udfmanager.h"
#include "parameter_define.h"
#include "variable.h"
extern UDFManager *ufin, *ufout, *ufres;
//extern UDFManager *ufsum;
/*!
  \brief Read udf files from command line prompt
 */
void file_get(const int argc, char *argv[]);

/*!
  \brief Read udf input file
 */
void Gourmet_file_io(const char *infile
		     ,const char *outfile
		     ,const char *sumfile
		     ,const char *deffile
		     ,const char *ctrlfile
		     ,const char *resfile
		     );

#endif
