#include "output_writer.h"

//fluid properties
const char* hdf5_writer::f_vel_name[] = {"u_x", "u_y", "u_z", "u"};
const char* hdf5_writer::f_tau_name[] = {"tau_xx", "tau_xy", "tau_xz", \
					 "tau_yy", "tau_yz", "tau_zz", "tau"};
const char* hdf5_writer::f_phi_name = "phi";
const char* hdf5_writer::f_pressure_name = "pressure";
const char* hdf5_writer::f_surface_charge_name = "surface_charge";
const char* hdf5_writer::f_solute_charge_name = "solute_charge";
const char* hdf5_writer::f_potential_charge_name= "e_potential";
//particle properties
const char* hdf5_writer::p_id_name = "id";
const char* hdf5_writer::p_spec_name = "spec";
const char* hdf5_writer::p_pos_name = "pos";
const char* hdf5_writer::p_pos_raw_name = "pos_raw";
const char* hdf5_writer::p_vel_name = "vel";
const char* hdf5_writer::p_omega_name = "omega";
const char* hdf5_writer::p_force_h_name = "force_h";
const char* hdf5_writer::p_force_r_name = "force_r";
const char* hdf5_writer::p_torque_h_name= "torque_h";
const char* hdf5_writer::p_torque_r_name= "torque_r";
const char* hdf5_writer::p_QR_name ="QR";

//hdf5 writer Constructor / Destroyer
hdf5_writer::hdf5_writer(const int&    _NX,
			 const int&    _NY,
			 const int&    _NZ,
			 const int&    _NZ_,
			 const double& _DX,
			 const int&    _nump,
			 const double& _dt,
			 const char*   _out_dir,
			 const char*   _out_name,
			 const Field_crop& _crop_field,
			 const Field_output& _print_field,
			 const Particle_output& _print_particle
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

    //TIME
    ts   = 0;
    dt   = static_cast<float>(_dt);

    //OUTPUT
    dircheckmake(_out_dir);
    sprintf(out_name, "%s", _out_name);
    sprintf(out_path, "%s/%s", _out_dir, _out_name);
    crop_field  = _crop_field;
    print_field = _print_field;
    print_particle = _print_particle;
  }
  
  {// HDF5
    fid = gid_time = gid_field = gid_part = -1;

    //Field data read/write access specifiers
    {
      //Array Dimenions
      mem_dims_field[0] = NX*NY*NZ_;
      out_dims_field[0] = NX;
      out_dims_field[1] = NY;
      out_dims_field[2] = NZ;

      //Memory access (ignore ghost fft points)
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


      //Memory acces hyperslab
      if(!crop_field.none){
	int rank = crop_field.rank;
	int slab_start = crop_field.start;
	int slab_width = crop_field.width;
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

	//Modify output parameters
	Origin[rank] = slab_start*DX;
	out_dims_field[rank] = slab_width;

	herr_t status = H5Sselect_hyperslab(mem_dataspace_field, H5S_SELECT_AND,
					    &mem_offset, &mem_stride, &mem_count, &mem_block);
	h5_check_err(status);
      }//field_crop

      //Disk Acess
      out_dataspace_field = H5Screate_simple(out_rank_field, out_dims_field, NULL);
      h5_check_err(out_dataspace_field);

      //check consistency between memory and disk spaces
      {
	hsize_t slab_npoints = H5Sget_select_npoints(mem_dataspace_field);
	assert(slab_npoints = out_dims_field[0]*out_dims_field[1]*out_dims_field[2]);
      }
    }
  }//HDF5

  {//Initialize output file
    char fname[128];
    sprintf(fname, "%s.h5", out_path);
    fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    h5_check_err(fid);
    
    //set global attributes
    {
      hid_t status;
      status = H5LTset_attribute_float(fid, "/", "dx", &DX, 1);
      h5_check_err(status);
      
      int dmy_dims[DIM] = {out_dims_field[0], out_dims_field[1], out_dims_field[2]};
      status = H5LTset_attribute_int(fid, "/", "nxnynz", dmy_dims, DIM);
      h5_check_err(status);
      
      status = H5LTset_attribute_float(fid, "/", "origin", Origin, DIM);
      h5_check_err(status);
      
      int dmy_num = nump - print_particle.first;
      status = H5LTset_attribute_int(fid, "/", "nump", &dmy_num, 1);
      h5_check_err(status);
    }
  }

}

hdf5_writer::~hdf5_writer(){
  //print grid coordinates required for xdmf format
  this -> write_xyz_coords();

  //print configure file
  this -> write_conf_file();

  herr_t status;
  status = H5Sclose(mem_dataspace_field);
  h5_check_err(status);

  status = H5Sclose(out_dataspace_field);
  h5_check_err(status);

  status = H5Fclose(fid);
  h5_check_err(status);
}

//
//Inherited functions
//
void hdf5_writer::write_start(){
  char dmy_group[128];
  //open file & groups
  {
    sprintf(dmy_group, "/frame_%d", ts);
    gid_time  = H5Gcreate(fid, dmy_group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    h5_check_err(gid_time);

    //create field group
    gid_field = H5Gcreate(gid_time, "field_data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    h5_check_err(gid_field);

    //create particle group
    gid_part  = H5Gcreate(gid_time, "particle_data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    h5_check_err(gid_part);
  }

  //add attributes
  {
    herr_t status;
    float dmy_float = static_cast<float>(ts)*dt;
    status = H5LTset_attribute_float(gid_time, dmy_group, "time", &dmy_float, 1);
    h5_check_err(status);
  }
}
void hdf5_writer::write_end(){
  herr_t status;

  status = H5Gclose(gid_part);
  h5_check_err(status);

  status = H5Gclose(gid_field);
  h5_check_err(status);

  status = H5Gclose(gid_time);
  h5_check_err(status);
  ts += 1;
}
inline void hdf5_writer::write_field_scalar(const double* phi, const char* name) {
  write_data(gid_field, phi, name, 
	     H5T_NATIVE_DOUBLE, mem_dataspace_field,
	     H5T_NATIVE_FLOAT, out_dataspace_field);
}
inline void hdf5_writer::write_particle_scalar(const int * data, const char* name){
  hsize_t dmy_dims[1] = {nump - print_particle.first};
  write_data(gid_part, data, name, H5T_NATIVE_INT, H5T_NATIVE_INT, 1, dmy_dims);
}
inline void hdf5_writer::write_particle_vectorn(const double*data, const int& dim, const char* name){
  hsize_t dmy_dims[2] = {nump - print_particle.first, dim};
  write_data(gid_part, data, name, H5T_NATIVE_DOUBLE, H5T_NATIVE_FLOAT, 2, dmy_dims);
}
inline void hdf5_writer::write_particle_matrix3(const double* data, const char* name){
  hsize_t dmy_dims[3] = {nump - print_particle.first, DIM, DIM};
  write_data(gid_part, data, name, H5T_NATIVE_DOUBLE, H5T_NATIVE_FLOAT, 3, dmy_dims);
}
void hdf5_writer::write_field_data(double** u, double* phi, double* pressure, double** tau){
  if(print_field.none) return;

  if(print_field.vel){
    this -> write_field_scalar(u[0], f_vel_name[0]);
    this -> write_field_scalar(u[1], f_vel_name[1]);
    this -> write_field_scalar(u[2], f_vel_name[2]);
  }
  if(print_field.phi){
    this -> write_field_scalar(phi, f_phi_name);
  }
  if(print_field.pressure){
    this -> write_field_scalar(pressure, f_pressure_name);
  }
  if(print_field.tau){
    this -> write_field_scalar(tau[0], f_tau_name[0]);
    this -> write_field_scalar(tau[1], f_tau_name[1]);
    this -> write_field_scalar(tau[2], f_tau_name[2]);
    this -> write_field_scalar(tau[3], f_tau_name[3]);
    this -> write_field_scalar(tau[4], f_tau_name[4]);
  }
}
void hdf5_writer::write_charge_field_data(double** u, double* phi, 
					  double* surface_charge, double* solute_charge, double* potential){
  if(print_field.none) return;

  if(print_field.vel){
    this -> write_field_scalar(u[0], f_vel_name[0]);
    this -> write_field_scalar(u[1], f_vel_name[1]);
    this -> write_field_scalar(u[2], f_vel_name[2]);
  }
  if(print_field.phi){
    this -> write_field_scalar(phi, f_phi_name);
  }
  if(print_field.charge){
    this -> write_field_scalar(surface_charge, f_surface_charge_name);
    this -> write_field_scalar(solute_charge,  f_solute_charge_name);
    this -> write_field_scalar(potential, f_potential_charge_name);
  }
}

//Todo: figure out how to print structures using hdf5 interface (avoid mem copies)!
void hdf5_writer::write_particle_data(Particle *p){
  if(print_particle.none) return;

  int &firstp    = print_particle.first;
  int nump_print = nump - firstp;
  int nump_print3= 3*nump_print;
  int nump_print6= 6*nump_print;
  
  int*   psp    = alloc_1d_int(nump_print);
  int*   pid    = alloc_1d_int(nump_print);
  double* pdata = alloc_1d_double(nump_print3*3); 
  
  int j, jj;
#pragma omp parallel for schedule(dynamic, 1) private(j, jj)
  for(int i = firstp; i < nump; i++){
    j  = i - firstp;
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
  this->write_particle_scalar(pid, p_id_name);
  this->write_particle_scalar(psp, p_spec_name);
  this->write_particle_vectorn(&pdata[0], DIM, p_pos_name);
  this->write_particle_vectorn(&pdata[nump_print3], DIM, p_pos_raw_name);
  this->write_particle_vectorn(&pdata[nump_print6], DIM, p_vel_name);

#pragma opm parallel for schedule(dynamic, 1) private(j, jj)
  for(int i = firstp; i < nump; i++){
    j  = i - firstp;
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
  this->write_particle_vectorn(&pdata[0], DIM, p_force_h_name);
  this->write_particle_vectorn(&pdata[nump_print3], DIM, p_force_r_name);

#pragma omp parallel for schedule(dynamic, 1) private(j, jj)
  for(int i = firstp; i < nump; i++){
    j = i - firstp;
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
  this->write_particle_vectorn(&pdata[0], DIM, p_omega_name);
  this->write_particle_vectorn(&pdata[nump_print3], DIM, p_torque_h_name);
  //this->write_particle_vectorn(&pdata[nump_print6], DIM, p_torque_r_name);
  
#pragma omp parallel for schedule(dynamic, 1) private(j, jj)
  for(int i = firstp; i < nump; i++){
    j  = i - firstp;
    jj = 9*j;
    
    qtn_normalize(p[i].q);
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
  this->write_particle_matrix3(&pdata[0], p_QR_name);
      
  free_1d_int(psp);
  free_1d_int(pid);
  free_1d_double(pdata);
}

void hdf5_writer::show_parameter() {
  fprintf(stderr, "# *** HDF5 Parameters ***\n");
  fprintf(stderr, "# Output Name          = %s.h5\n", out_path);
  fprintf(stderr, "# Mem Layout           = (%d, %d, %d, %d)\n", NX, NY, NZ, NZ_);
  fprintf(stderr, "# Grid spacing         = %.3g\n", DX);
  fprintf(stderr, "# Time between frames  = %.3g\n", dt);
  if(!print_particle.none){
    fprintf(stderr, "# Particles Range      = %d -> %d\n", print_particle.first, nump);
  }else{
    fprintf(stderr, "# Particle Data is being suppressed\n");
  }
  if(!print_field.none){
    if(!crop_field.none){
      fprintf(stderr, "# Field data is being cropped \n");
      fprintf(stderr, "# Slab axis  : %d (0=yz, 1=xz, 2=xy)\n", crop_field.rank);
      fprintf(stderr, "# Slab start : %d\n", crop_field.start);
      fprintf(stderr, "# Slab width : %d\n", crop_field.width);
    }
    fprintf(stderr, "# Output Dimensions    = (%lld, %lld, %lld)\n", out_dims_field[0], out_dims_field[1], out_dims_field[2]);
    fprintf(stderr, "# Output Origin        = (%.3g, %.3g, %.3g)\n", Origin[0], Origin[1], Origin[2]);
    fprintf(stderr, "# Print velocity field : %s\n", (print_field.vel ? "YES" : "NO"));
    fprintf(stderr, "# Print phi field      : %s\n", (print_field.phi ? "YES" : "NO"));
    fprintf(stderr, "# Print charge field   : %s\n", (print_field.charge  ? "YES" : "NO"));
    fprintf(stderr, "# Print Pressure field : %s\n", (print_field.pressure ? "YES" : "NO"));
    fprintf(stderr, "# Print Stress field   : %s\n", (print_field.tau ? "YES" : "NO"));
  }else{
    fprintf(stderr, "# Field Data is being suppressed\n");
  }
  fprintf(stderr, "# ***********************\n");
}

//
//New functions
//
void hdf5_writer::write_conf_file(){
  
  FILE* conf;
  char confname[128];
  sprintf(confname, "%s.conf", out_path);
  conf = filecheckopen(confname, "w");
  //Parameters
  fprintf(conf, "[Filename]\n %s\n", out_name);
  fprintf(conf, "[Grid Dimensions]\n %lld %lld %lld\n", 
	  out_dims_field[0], out_dims_field[1], out_dims_field[2]);
  fprintf(conf, "[Time Series]\n 0.0 %.6f %d\n", dt, ts);
  fprintf(conf, "[Particle Number]\n %d\n", nump - print_particle.first);

  //Particle Data
  if(!print_particle.none){
    //Particle Grid
    fprintf(conf, "[Particle Grid]\n %s\n", p_pos_name);
    fprintf(conf, "#[Particle Grid]\n#%s\n", p_pos_raw_name);
    //Scalar Particle Quantities
    fprintf(conf, "[Particle Scalar]\n %s\n", p_id_name);
    fprintf(conf, "[Particle Scalar]\n %s\n", p_spec_name);
    //Vector Particle Quantities
    fprintf(conf, "[Particle Vector]\n %s\n", p_vel_name);
    fprintf(conf, "[Particle Vector]\n %s\n", p_omega_name);
    fprintf(conf, "[Particle Vector]\n %s\n", p_force_h_name);
    fprintf(conf, "[Particle Vector]\n %s\n", p_force_r_name);
    fprintf(conf, "[Particle Vector]\n %s\n", p_torque_h_name);
    //Tensor Particle Quantities
    fprintf(conf, "[Particle Tensor]\n %s\n", p_QR_name);
  }

  // Field Data
  if(!print_field.none){
    if(print_field.vel)
      fprintf(conf, "[Field Vector]\n %s %s %s %s\n", 
	      f_vel_name[0], f_vel_name[1], f_vel_name[2], f_vel_name[3]);
    if(print_field.phi)
      fprintf(conf, "[Field Scalar]\n %s\n", f_phi_name);
    if(print_field.pressure)
      fprintf(conf, "[Field Scalar]\n %s\n", f_pressure_name);
    if(print_field.tau){
      //missing zz component to get symmetric tensor6 field...
      fprintf(conf, "[Field Scalar]\n %s\n", f_tau_name[0]);
      fprintf(conf, "[Field Scalar]\n %s\n", f_tau_name[1]);
      fprintf(conf, "[Field Scalar]\n %s\n", f_tau_name[2]);
      fprintf(conf, "[Field Scalar]\n %s\n", f_tau_name[3]);
      fprintf(conf, "[Field Scalar]\n %s\n", f_tau_name[4]);
    }
    if(print_field.charge){
      fprintf(conf, "[Field Scalar]\n %s\n", f_surface_charge_name);
      fprintf(conf, "[Field Scalar]\n %s\n", f_solute_charge_name);
      fprintf(conf, "[Field Scalar]\n %s\n", f_potential_charge_name);
    }
  }
  fclose(conf);
}

void hdf5_writer::write_xyz_coords(){
  double* work_v3[DIM];
  for(int d = 0; d < DIM; d++) 
    work_v3[d] = alloc_1d_double(NX*NY*NZ_);

  double ddx = static_cast<double>(DX);
  for(int i = 0; i < NX; i++){
    for(int j = 0; j < NY; j++){
      for(int k = 0; k < NZ; k++){
	int im = (i*NY*NZ_) + (j*NZ_) + k;
	work_v3[0][im] = static_cast<double>(i)*ddx;
	work_v3[1][im] = static_cast<double>(j)*ddx;
	work_v3[2][im] = static_cast<double>(k)*ddx;
      }
    }
  }
  hid_t gid_coord = H5Gcreate(fid, "/coord_data", 
			      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  h5_check_err(gid_coord);
  write_data(gid_coord, work_v3[0], "x", 
	     H5T_NATIVE_DOUBLE, mem_dataspace_field,
	     H5T_NATIVE_FLOAT, out_dataspace_field);
  write_data(gid_coord, work_v3[1], "y", 
	     H5T_NATIVE_DOUBLE, mem_dataspace_field,
	     H5T_NATIVE_FLOAT, out_dataspace_field);
  write_data(gid_coord, work_v3[2], "z", 
	     H5T_NATIVE_DOUBLE, mem_dataspace_field,
	     H5T_NATIVE_FLOAT, out_dataspace_field);
  herr_t status = H5Gclose(gid_coord);
  h5_check_err(status);
  for(int d = 0; d < DIM; d++)
    free_1d_double(work_v3[d]);
}
//end hdf5 writer
