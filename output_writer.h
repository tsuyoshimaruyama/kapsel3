#ifndef OUTPUT_WRITER_H
#define OUTPUT_WRITER_H

#include <hdf5.h>
#include "output_writer_base.h"

class hdf5_writer : public output_writer {
 public:
  hdf5_writer(const int&_nx, 
	      const int&_ny,
	      const int&_nz,
	      const int&_nz_,
	      const int&_nump
	      );

  ~hdf5_writer();

  //initialize writing one timestep
  void write_start(const CTime &time);

  //finilize writing one timestep
  void write_end();

  //write field data
  void write_field_data(double const* phi);

  //write multi-dim field data
  void write_field_data(double const* const* u, const int& dim);

  //write particle data
  void write_particle_data(Particle const* p);

  //write output parameters to stderr
  void show_parameter();

 private:
};
#endif
