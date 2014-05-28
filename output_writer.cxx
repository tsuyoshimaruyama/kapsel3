#include "output_writer.h"

hdf5_writer::hdf5_writer(const int& _nx,
			 const int& _ny,
			 const int& _nz,
			 const int& _nz_,
			 const int& _nump
			 )

{
  fprintf(stderr, "#HDF5: Creating HDF5\n");
}
hdf5_writer::~hdf5_writer(){
  fprintf(stderr, "#HDF5: Closing writer\n");
}

void hdf5_writer::write_start(const CTime &time){
  fprintf(stderr, "#HDF5: Open file %d\n", time.ts);
}
void hdf5_writer::write_end(){
  fprintf(stderr, "#HDF5: Close file\n");
}
void hdf5_writer::write_field_data(double const* phi){
  fprintf(stderr, "#HDF5: Writing 1-D field data\n");
}
void hdf5_writer::write_field_data(double const* const* u, const int&dim){
  fprintf(stderr, "#HDF5: Writing N-D field data\n");
}
void hdf5_writer::write_particle_data(Particle const* p){
  fprintf(stderr, "#HDF5: Writing Particle data\n");
}
void hdf5_writer::show_parameter(){
  fprintf(stderr, "#HDF5: Writing Parameters\n");
}
