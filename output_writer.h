#ifndef OUTPUT_WRITER_H
#define OUTPUT_WRITER_H

#include <stdio.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <assert.h>
#include "macro.h"
#include "parameter_define.h"
#include "variable.h"
#include "output_writer_base.h"

class hdf5_writer : public output_writer {
 public:
  /*!
    \brief Class constructor for hdf5 writer class
    \param[in] _NX points in x
    \param[in] _NY points in y
    \param[in] _NZ points in z
    \param[in] _NZ_ total points in z, including ghost points used in fft computations
    \param[in] _DX sclae for grid spacing
    \param[in] _Nump total number of particles
    \param[in] _dt time step between frames (assumes printing at regular intervals)
    \param[in] _out_dir output directory for trajectory files
    \param[in] _out_name name of trajectory files
   */
  hdf5_writer(const int&     _NX, 
	      const int&     _NY,
	      const int&     _NZ,
	      const int&     _NZ_,
	      const double&  _DX,
	      const int&     _Nump,
	      const double&  _dt,
	      const char*    _out_dir,
	      const char*    _out_name
	      );

  /*!
    \brief Overloaded class destroyer
    \details Make sure to release all allocated memory
   */
  ~hdf5_writer();

  //
  // Functions inherited from virtual base class
  //

  /*! 
    \brief Initialization for writing one frame to disk
  */
  void write_start();

  /*!
    \brief Finish writing current frame to disk
  */
  void write_end();

  /*!
    \brief Write field data to open frame
   */
  void write_field_data(double const* phi, const char* name);

  /*!
    \brief Writes particle data to open frame
   */
  void write_particle_data(Particle *p);

  /*!
    \brief Write output parameters to stderr
  */
  void show_parameter();


  //
  // New functions
  //
  
  /*!
    \brief Print x,y,z coordiantes for ALL points 
    \details Required to obtain correct visualization using xdmf
    \param[in,out] work_v3 working memory
   */
  void write_xyz_coords(double** work_v3);

  /*!
    \brief Select hyperslab to output only partial field data
    \param[input] rank selects axis to reduce: yz slab (0), xz slab (1), xy slab (2)
    \param[input] slab_start start index for slab
    \param[input] slab_width number of slices in slab dimension
   */
  void set_hyperslab(const int &rank, const int &slab_start,  const int &slab_width);

  /*!
    \brief Select subset of particle to print
    \details if first is not a valid particle id, NO particles will be printed
   */
  void set_particle_mask(const int &first);

 private:
  /*
    Field Data
  */
  int NX;
  int NY;
  int NZ;
  int NZ_;
  float DX;
  float Origin[DIM];

  /*
    Particle Data
  */
  int nump;
  int startp;

  /*
    Time Data
   */
  int    ts;
  float  dt;
  
  /*
    Output Param
  */
  char   out_name[256];

  /*
    HDF5 Parameters
   */
  hid_t  fid;        //file id
  hid_t  gid_time;   //time frame id
  hid_t  gid_field;  //group id for field data
  hid_t  gid_part;   //group id for particle data

  //Field Output
  static const hsize_t mem_rank_field = 1;
  hsize_t mem_dims_field[mem_rank_field];
  hid_t mem_dataspace_field;

  static const hsize_t out_rank_field = DIM;
  hsize_t out_dims_field[out_rank_field];
  hid_t out_dataspace_field;

  //
  //Private Functions
  //

  /*!
    \brief check h5 error code
   */
  inline void h5_check_err(const herr_t &err){
    if(err < 0){
      fprintf(stderr, "# HDF5 Error:\n");
      H5Eprint(H5E_DEFAULT, stderr);
      exit_job(EXIT_FAILURE);
    }
  }

  /*!
    \brief Writes buffer array to disk
    \param[in] loc_id location id to write to
    \param[in] buffer array to write
    \param[in] name name of new dataset to create
    \param[in] mem_dtype data type in memory
    \param[in] mem_dataspace dataspace to use when reading buffer from memory
    \param[in] out_dtype data type on disk
    \param[in] out_dataspace dataspace to use when writing buffer to disk
   */
  inline void write_data(hid_t loc_id, const void* buffer, const char* name,
			 const hid_t mem_dtype, hid_t &mem_dataspace,
			 const hid_t out_dtype, hid_t &out_dataspace){

    hid_t status;
    hid_t out_dataset = H5Dcreate(loc_id, name, out_dtype, out_dataspace,
				  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    h5_check_err(out_dataset);

    status = H5Dwrite(out_dataset, mem_dtype, mem_dataspace, 
		      H5S_ALL, H5P_DEFAULT, buffer);
    h5_check_err(status);

    status = H5Dclose(out_dataset);
    h5_check_err(status);
  }
  /*!
    \brief Writes buffer array to disk
    \param[in] lod_id location id to write to
    \param[in] buffer array to write
    \param[in] name name of new dataset to create
    \param[in] mem_dtype data type in memory
    \param[in] out_dtype data type on disk
    \param[in] out_rank rank of array on disk
    \param[in] out_dim size of each dimension on disk
   */
  inline void write_data(hid_t loc_id, const void* buffer, const char* name,
			 const hid_t mem_dtype, const hid_t out_dtype,
			 const int &out_rank, const hsize_t* out_dim
			 ){
    hid_t status;
    hid_t out_dataspace, out_dataset;

    out_dataspace = H5Screate_simple(out_rank, out_dim, NULL);
    h5_check_err(out_dataspace);

    out_dataset   = H5Dcreate(loc_id, name, out_dtype, out_dataspace, 
			      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    h5_check_err(out_dataset);

    status = H5Dwrite(out_dataset, mem_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
    h5_check_err(status);

    status = H5Dclose(out_dataset);
    h5_check_err(status);

    status = H5Sclose(out_dataspace);
    h5_check_err(status);
  }
  
  /*!
    \brief Write scalar particle data to current frame
   */
  void write_particle_scalar(const int* data, const char* name);
  /*!
    \brief Write vector particle data to current frame
   */
  void write_particle_vector3(const double* data, const char* name);
  /*!
    \brief Write quaternion particle data to current frame
   */
  void write_particle_quaternion(const double* data, const char* name);
  /*!
    \brief Write 3x3 square matrix particle data to current frame
   */
  void write_particle_matrix3(const double* data, const char* name);
};

#endif
