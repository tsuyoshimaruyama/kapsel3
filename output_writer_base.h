#ifndef OUTPUT_WRITER_BASE_H
#define OUTPUT_WRITER_BASE_H
#include "variable.h"
class output_writer {
 public:

  //destructor of derived classes should free all allocated memory
  virtual ~output_writer(){};

  //virtual function to initialize writing one timestep
  virtual void write_start(const CTime &time) = 0;

  //virtual function to finilize writing one timestep
  virtual void write_end() = 0;

  //virtual function to write field data
  virtual void write_field_data(double const* phi) = 0;

  //virtual function to write multi-dim field data
  virtual void write_field_data(double const* const* u, const int& dim) = 0;

  //virtual function to write particle data
  virtual void write_particle_data(const Particle* p) = 0;

  //virtual function to write output parameters to stderr
  virtual void show_parameter() = 0;
};
#endif
