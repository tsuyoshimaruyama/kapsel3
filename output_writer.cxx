#include "output_writer.h"

//hdf5 writer Constructor / Destroyer
hdf5_writer::hdf5_writer(const int&    _NX,
			 const int&    _NY,
			 const int&    _NZ,
			 const int&    _NZ_,
			 const double& _DX,
			 const int&    _nump,
			 const double& _dt,
			 const char*   _out_dir,
			 const char*   _out_name
			 )

{
  {
    //FIELD
    NX   = _NX;
    NY   = _NY;
    NZ   = _NZ;
    NZ_  = _NZ_;
    DX   = static_cast<float>(_DX);
    Origin[0] = Origin[1] = Origin[2] = 0.0;    

    //PARTICLE
    nump   = _nump;
    startp = 0;

    //TIME
    ts   = 0;
    dt   = static_cast<float>(_dt);

    //OUTPUT
    dircheckmake(_out_dir);
    sprintf(out_name, "%s/%s", _out_dir, _out_name);
  }
  
  {// HDF5
    fid = gid_field = gid_part = -1;

    //Field data read/write access specifiers
    {
      //Memory access (ignore ghost fft points)
      mem_dims_field[0] = NX*NY*NZ_;
      mem_dataspace_field = H5Screate_simple(mem_rank_field, mem_dims_field, NULL);
      h5_check_err(mem_dataspace_field);
      
      hsize_t mem_offset = 0;
      hsize_t mem_stride = NZ_;
      hsize_t mem_count  = NX*NY;
      hsize_t mem_block  = NZ;
      herr_t  status  = H5Sselect_hyperslab(mem_dataspace_field, H5S_SELECT_SET,
					    &mem_offset, &mem_stride, 
					    &mem_count, &mem_block);
      h5_check_err(status);

      //Disk access
      out_dims_field[0] = NX;
      out_dims_field[1] = NY;
      out_dims_field[2] = NZ;
      out_dataspace_field = H5Screate_simple(out_rank_field, out_dims_field, NULL);
      h5_check_err(out_dataspace_field);
    }
  }

}
hdf5_writer::~hdf5_writer(){
  herr_t status; 
  status = H5Sclose(mem_dataspace_field);
  h5_check_err(status);

  status = H5Sclose(out_dataspace_field);
  h5_check_err(status);
}

//
//Inherited functions
//
void hdf5_writer::write_start(){
  //open file & groups
  {
    char fname[256];
    sprintf(fname, "%s.%d.h5", out_name, ts);
    fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    h5_check_err(fid);
    
    //create field group
    gid_field = H5Gcreate(fid, "/field_data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    h5_check_err(gid_field);

    //create particle group
    gid_part  = H5Gcreate(fid, "/particle_data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    h5_check_err(gid_part);
  }

  //add attributes
  {
    herr_t status;
    char dmy_group[128];
    {
      sprintf(dmy_group,"/");
      status = H5LTset_attribute_int(fid, dmy_group, "ts", &ts, 1);
      h5_check_err(status);

      float dmy_float = static_cast<float>(ts)*dt;
      status = H5LTset_attribute_float(fid, dmy_group, "time", &dmy_float, 1);
      h5_check_err(status);

    }
    {
      sprintf(dmy_group, "/field_data");
      status = H5LTset_attribute_float(fid, dmy_group, "dx", &DX, 1);
      h5_check_err(status);

      int dmy_dims[DIM] = {out_dims_field[0], out_dims_field[1], out_dims_field[2]};
      status = H5LTset_attribute_int(fid, dmy_group, "nxnynz", dmy_dims, DIM);
      h5_check_err(status);

      status = H5LTset_attribute_float(fid, dmy_group, "origin", Origin, DIM);
      h5_check_err(status);
    }
  }
}
void hdf5_writer::write_end(){
  herr_t status;

  status = H5Gclose(gid_part);
  h5_check_err(status);

  status = H5Gclose(gid_field);
  h5_check_err(status);

  status = H5Fclose(fid);
  h5_check_err(status);
  ts += 1;
}
void hdf5_writer::write_field_data(const double* phi, const char* name) {
  write_data(gid_field, phi, name, 
	     H5T_NATIVE_DOUBLE, mem_dataspace_field,
	     H5T_NATIVE_FLOAT, out_dataspace_field);
}

inline void hdf5_writer::write_particle_scalar(const int * data, const char* name){
  hsize_t dmy_dims[1] = {nump - startp};
  write_data(gid_part, data, name, H5T_NATIVE_INT, H5T_NATIVE_INT, 1, dmy_dims);
}
inline void hdf5_writer::write_particle_vector3(const double* data, const char* name){
  hsize_t dmy_dims[2] = {nump - startp, 3};
  write_data(gid_part, data, name, H5T_NATIVE_DOUBLE, H5T_NATIVE_FLOAT, 2, dmy_dims);
}
inline void hdf5_writer::write_particle_quaternion(const double* data, const char* name){
  hsize_t dmy_dims[2] = {nump - startp, 4};
  write_data(gid_part, data, name, H5T_NATIVE_DOUBLE, H5T_NATIVE_FLOAT, 2, dmy_dims);
}
inline void hdf5_writer::write_particle_matrix3(const double* data, const char* name){
  hsize_t dmy_dims[3] = {nump - startp, 3, 3};
  write_data(gid_part, data, name, H5T_NATIVE_DOUBLE, H5T_NATIVE_FLOAT, 3, dmy_dims);
}
void hdf5_writer::write_particle_data(Particle *p){
  int nump_print = nump - startp;
  int nump_print3= 3*nump_print;
  int nump_print6= 6*nump_print;
  
  int*   psp    = (int*)malloc(sizeof(int)*(nump_print));
  int*   pid    = (int*)malloc(sizeof(int)*(nump_print));
  double* pdata = (double*)malloc(sizeof(double)*nump_print3*3); 
  
  int j, jj;
#pragma omp parallel for schedule(dynamic, 1) private(j, jj)
  for(int i = startp; i < nump; i++){
    j  = i - startp;
    jj = 3*j;
    
    pid[j]   = i;
    psp[j]   = p[i].spec;
    
    pdata[jj]  = p[i].x[0];
    pdata[jj+1]= p[i].x[1];
    pdata[jj+2]= p[i].x[2];
    
    pdata[nump_print3+jj]  = p[i].x_nopbc[0];
    pdata[nump_print3+jj+1]= p[i].x_nopbc[1];
    pdata[nump_print3+jj+2]= p[i].x_nopbc[2];
    
    pdata[nump_print6+jj]  = p[i].v[0];
    pdata[nump_print6+jj+1]= p[i].v[1];
    pdata[nump_print6+jj+2]= p[i].v[2];
  }
  this->write_particle_scalar(pid, "id");
  this->write_particle_scalar(psp, "spec");
  this->write_particle_vector3(&pdata[0], "r");
  this->write_particle_vector3(&pdata[nump_print3], "r_raw");
  this->write_particle_vector3(&pdata[nump_print6], "v");

#pragma opm parallel for schedule(dynamic, 1) private(j, jj)
  for(int i = startp; i < nump; i++){
    j  = i - startp;
    jj = 3*j;

    //total hydrodynamic force
    pdata[jj]  = p[i].f_hydro_previous[0] + p[i].f_slip_previous[0];
    pdata[jj+1]= p[i].f_hydro_previous[1] + p[i].f_slip_previous[1];
    pdata[jj+2]= p[i].f_hydro_previous[2] + p[i].f_slip_previous[2];

    //other forces (LJ, external, etc.)
    pdata[nump_print3+jj]  = p[i].fr_previous[0];
    pdata[nump_print3+jj+1]= p[i].fr_previous[1];
    pdata[nump_print3+jj+2]= p[i].fr_previous[2];
  }
  this->write_particle_vector3(&pdata[0], "f_hydro");
  this->write_particle_vector3(&pdata[nump_print3], "f_r");

#pragma omp parallel for schedule(dynamic, 1) private(j, jj)
  for(int i = startp; i < nump; i++){
    j = i - startp;
    jj = 3*j;

    pdata[jj]  = p[i].omega[0];
    pdata[jj+1]= p[i].omega[1];
    pdata[jj+2]= p[i].omega[2];
    
    //total hydrodynamic torque
    pdata[nump_print3+jj]  = p[i].torque_hydro_previous[0] + p[i].torque_slip_previous[0];
    pdata[nump_print3+jj+1]= p[i].torque_hydro_previous[1] + p[i].torque_slip_previous[1];
    pdata[nump_print3+jj+2]= p[i].torque_hydro_previous[2] + p[i].torque_slip_previous[2];

    /**** torque_r = 0
    pdata[nump_print6+jj]  = p[i].torque_r_previous[0];
    pdata[nump_print6+jj+1]= p[i].torque_r_previous[1];
    pdata[nump_print6+jj+2]= p[i].torque_r_previous[2];
    ****/
  }
  this->write_particle_vector3(&pdata[0], "w");
  this->write_particle_vector3(&pdata[nump_print3], "t_hydro");
  //this->write_particle_vector3(&pdata[nump_print6], "t_r);
  
#pragma omp parallel for schedule(dynamic, 1) private(j, jj)
  for(int i = startp; i < nump; i++){
    j  = i - startp;
    jj = 9*j;
    
    rqtn_rm(p[i].QR, p[i].q);
    //x axis
    pdata[jj]  = p[i].QR[0][0]; 
    pdata[jj+1]= p[i].QR[1][0]; 
    pdata[jj+2]= p[i].QR[2][0];
    //y axis
    pdata[jj+3]= p[i].QR[0][1];
    pdata[jj+4]= p[i].QR[1][1]; 
    pdata[jj+5]= p[i].QR[2][1];
    //z axis
    pdata[jj+6]= p[i].QR[0][2]; 
    pdata[jj+7]= p[i].QR[1][2];
    pdata[jj+8]= p[i].QR[2][2];
  }
  this->write_particle_matrix3(&pdata[0], "QR");
      
  free(psp);
  free(pid);
  free(pdata);
}

void hdf5_writer::show_parameter() {
  fprintf(stderr, "# *** HDF5 Parameters ***\n");
  fprintf(stderr, "# Mem Layout           = (%d, %d, %d, %d)\n", 
	  NX, NY, NZ, NZ_);
  fprintf(stderr, "# Grid spacing         = %.3g\n", DX);
  fprintf(stderr, "# Time between frames  = %.3g\n", dt);
  fprintf(stderr, "# Num Particles        = %d\n", nump);
  fprintf(stderr, "# First Print particle = %d\n", startp);
  fprintf(stderr, "# Output Name          = %s\n", out_name);
  fprintf(stderr, "# Output Dimensions    = (%lld, %lld, %lld)\n",
	  out_dims_field[0], out_dims_field[1], out_dims_field[2]);
  fprintf(stderr, "# Output Origin        = (%.3g, %.3g, %.3g)\n", 
	  Origin[0], Origin[1], Origin[2]);
  fprintf(stderr, "# ************************\n");
}

//
//New functions
//

void hdf5_writer::write_xyz_coords(double** work_v3){
  char fname[256];
  sprintf(fname, "%s_coord.h5", out_name);
  hid_t cfid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  h5_check_err(cfid);
  for(int i = 0; i < NX; i++){
    for(int j = 0; j < NY; j++){
      for(int k = 0; k < NZ; k++){
	int im = (i*NY*NZ_) + (j*NZ_) + k;
	work_v3[0][im] = static_cast<double>(i);
	work_v3[1][im] = static_cast<double>(j);
	work_v3[2][im] = static_cast<double>(k);
      }
    }
  }
  write_data(cfid, work_v3[0], "x", 
	     H5T_NATIVE_DOUBLE, mem_dataspace_field,
	     H5T_NATIVE_FLOAT, out_dataspace_field);
  write_data(cfid, work_v3[1], "y", 
	     H5T_NATIVE_DOUBLE, mem_dataspace_field,
	     H5T_NATIVE_FLOAT, out_dataspace_field);
  write_data(cfid, work_v3[2], "z", 
	     H5T_NATIVE_DOUBLE, mem_dataspace_field,
	     H5T_NATIVE_FLOAT, out_dataspace_field);
  herr_t status = H5Fclose(cfid);
  h5_check_err(status);
}

void hdf5_writer::set_hyperslab(const int &rank,
				const int &slab_start,
				const int &slab_width){
  //select hyperslab for memory access
  {
    hsize_t mem_offset;
    hsize_t mem_stride;
    hsize_t mem_count;
    hsize_t mem_block;
    if(rank == 0){ //yz slab
      mem_offset = slab_start*(NY*NZ_);
      mem_stride = NY*NZ_;
      mem_count  = slab_width;
      mem_block  = NY*NZ_;
    }else if(rank == 1){ //xz slab
      mem_offset = slab_start*NZ_;
      mem_stride = NY*NZ_;
      mem_count  = NX;
      mem_block  = slab_width;
    }else if(rank == 2){ //xy slab
      mem_offset = slab_start;
      mem_stride = NZ_;
      mem_count  = NY*NX;
      mem_block  = slab_width;
    }else{
      fprintf(stderr, "# Invalid slab rank %d\n", rank);
      exit_job(EXIT_FAILURE);
    }
    Origin[rank] = slab_start*DX;
    
    herr_t status = H5Sselect_hyperslab(mem_dataspace_field, H5S_SELECT_AND,
					&mem_offset, &mem_stride, &mem_count, &mem_block);
  }

  //modify dataspace for output
  {
    out_dims_field[rank] = slab_width;
    herr_t status = H5Sset_extent_simple(out_dataspace_field, out_rank_field, 
					 out_dims_field, out_dims_field);
  }

  //check consistency
  {
    hsize_t slab_npoints = H5Sget_select_npoints(mem_dataspace_field);
    assert(slab_npoints = out_dims_field[0]*out_dims_field[1]*out_dims_field[2]);
  }
  fprintf(stderr, "# HDF5 Fluid output: using hyperslab (rank=%d, start=%d, width=%d)\n",
	  rank, slab_start, slab_width
	  );
}

void hdf5_writer::set_particle_mask(const int &first){
  startp = (first >= 0 && first < nump ? first : nump);
  fprintf(stderr, "# HDF5 Particle output: skipping first %d particles\n", startp);
}
//end hdf5 writer
