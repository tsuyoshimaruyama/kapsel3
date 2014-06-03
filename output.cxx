/*!
  \file output.cxx
  \brief Wrapper rouintes to write output files
  \details Prepare all data for output, leave actual writing to specialized class
  \author J. Molina
  \date 2014/05/21
  \version 1.0
 */

#include "output.h"

output_writer *writer;
void Init_output(){

  //Init Parameters

  //field selections
  if(print_field.vel || print_field.phi || 
     (print_field.pressure && SW_EQ != Electrolyte) || 
     (print_field.tau && SW_EQ != Electrolyte) ||
     (print_field.rho && SW_EQ == Electrolyte)
     ){
    print_field.none = false;
  }else{
    print_field.none = true;
  }

  //hyper slab selection
  if(print_field_crop.rank >=0 && print_field_crop.rank < DIM && !print_field.none){
    print_field_crop.none = false;
  }else{
    print_field_crop.none = true;
  }

  //particle selection
  if(print_particle.first >= 0 && print_particle.first < Particle_Number 
     && Particle_Number > 0){
    print_particle.none = false;
  }else{
    print_particle.none = true;
  }


  //Init Writers
  if(SW_OUTFORMAT == OUT_AVS_BINARY ||  SW_OUTFORMAT == OUT_AVS_ASCII){ // LEGACY AVS
    //extended options not supported for AVS

    //no cropping
    print_field_crop.none = true;

    //print all fields
    print_field.none = false;
    print_field.vel  = true;
    print_field.phi  = true;
    print_field.rho  = true;
    print_field.pressure = true;
    print_field.tau  = true;

    //print all particle data
    print_particle.none = false;
    print_particle.first = 0;

    Set_avs_parameters(Avs_parameters);
    Init_avs(Avs_parameters);
    if(Particle_Number > 0){
      Init_avs_p(Avs_parameters);
    }
  }else if(SW_OUTFORMAT == OUT_EXT){
    if(SW_EXTFORMAT == EXT_OUT_HDF5){
      hdf5_writer* h5writer = new hdf5_writer(NX, NY, NZ, NZ_, DX,
					      Particle_Number, 
					      DT*GTS,
					      Out_dir, 
					      Out_name
					      );
      if(!print_field_crop.none)
	h5writer->set_hyperslab(print_field_crop.rank, 
				print_field_crop.start, 
				print_field_crop.width);
      if(!print_particle.none)
	h5writer->set_particle_mask(print_particle.first);
      
      //print grid coordinates required for xdmf format
      h5writer -> write_xyz_coords(work_v3);
      
      writer = h5writer;
    }
  }
}
void Free_output(){
  if(SW_OUTFORMAT == OUT_EXT){
    writer->~output_writer();
  }
}
void Show_output_parameter(){
  fprintf(stderr, "#### Output parameters #### \n");

  if(SW_UDF){
    fprintf(stderr, "# UDF output is enabled\n");
    fprintf(stderr, "# for UDF ->\t%s\n",Out_udf);
  }else {
    fprintf(stderr, "# UDF output is supressed.\n");
  }

  if(SW_OUTFORMAT == OUT_AVS_BINARY || SW_OUTFORMAT == OUT_AVS_ASCII){
    Show_avs_parameter();
  }else if(SW_OUTFORMAT == OUT_EXT){
    fprintf(stderr, "# AVS output is suppressed.\n");
    fprintf(stderr, "# Extened output enabled.\n");
    writer -> show_parameter();
  }else{
    fprintf(stderr, "# Field/Particle output is disabled.\n");
  }

  //field selection
  if(!print_field.none){
    fprintf(stderr, "# Field data is being printed\n");
    fprintf(stderr, "# Print velocity field : %s\n", 
	    (print_field.vel ? "YES" : "NO"));
    fprintf(stderr, "# Print phi field      : %s\n", 
	    (print_field.phi ? "YES" : "NO"));
    fprintf(stderr, "# Print rho field      : %s\n", 
	    ((print_field.rho && SW_EQ == Electrolyte) ? "YES" : "NO"));
    fprintf(stderr, "# Print Pressure field : %s\n",
	    ((print_field.pressure && SW_EQ != Electrolyte) ? "YES" : "NO"));
    fprintf(stderr, "# Print Stress field   : %s\n",
	    ((print_field.tau && SW_EQ != Electrolyte) ? "YES" : "NO"));
  }else{
    fprintf(stderr, "# Field data is suppressed\n");
  }

  //hyperslab
  if(!print_field_crop.none){
    fprintf(stderr, "# Field data is being cropped \n");
    fprintf(stderr, "# Slab axis  : %d (0=yz, 1=xz, 2=xy)\n", print_field_crop.rank);
    fprintf(stderr, "# Slab start : %d\n", print_field_crop.start);
    fprintf(stderr, "# Slab width : %d\n", print_field_crop.width);
  }

  //Particle selection
  if(!print_particle.none){
    fprintf(stderr, "# Particle data is being printed \n");
  }else{
    fprintf(stderr, "# Particle data is being suppressed\n");
    fprintf(stderr, "# Note: Particle data can still be found in UDF (if enabled)\n");
  }
  
  fprintf(stderr, "####                   #### \n");
}
void Output_open_frame(){
  if(SW_OUTFORMAT == OUT_EXT)
    writer -> write_start();
    
}
void Output_close_frame(){
  if(SW_OUTFORMAT == OUT_EXT)
    writer -> write_end();
}

void Output_field_data(double** zeta,
		       double* uk_dc,
		       Particle* p,
		       const CTime &time){

  double *strain[QDIM]={f_particle[0]
		      ,f_particle[1]
		      ,f_particle[2]
		      ,f_ns0[0]
		      ,f_ns0[1]
  };

  if(print_field.vel || print_field.tau){

    if(SW_EQ == Shear_Navier_Stokes_Lees_Edwards){
      Zeta_k2u_k_OBL(zeta, uk_dc, u);
   
      if(print_field.tau){
	U_k2Strain_k_OBL(u, strain);
	for(int d = 0; d < QDIM; d++){
	  A_k2a(strain[d]);
	}
	E_oblique2E(strain, false); //without mean shear flow terms
      }//print strain ?
      
      if(print_field.vel){
	U_k2u(u);
	U_oblique2u(u);             //with mean shear flow terms
      }//print u?
    }else{
      Zeta_k2u_k(zeta, uk_dc, u);
      
      if(print_field.tau){
	U_k2Strain_k(u, strain);
	for(int d = 0; d < QDIM; d++){
	  A_k2a(strain[d]);
	}
      }//print strain?
      
      if(print_field.vel){
	U_k2u(u);
      }//print u?
    }

  }// Print u / strain ?

  //! TODO: implement pressure calculation
  if(print_field.pressure){
    A_k2a(Pressure); 
  }//print pressure?

  if(print_field.phi){
    Reset_phi(phi);
    if (SW_EQ == Shear_Navier_Stokes_Lees_Edwards) {
      Make_phi_particle_OBL(phi, p);
    }else {
      Make_phi_particle(phi, p);
    }
  }//print phi?

  if(SW_OUTFORMAT == OUT_AVS_BINARY || SW_OUTFORMAT == OUT_AVS_ASCII){
    Output_avs(Avs_parameters, u, phi, Pressure, strain, time);
  }else if(SW_OUTFORMAT == OUT_EXT){
    if(print_field.vel){
      writer -> write_field_data(u[0], "ux");
      writer -> write_field_data(u[1], "uy");
      writer -> write_field_data(u[2], "uz");
    }
    if(print_field.phi)     
      writer -> write_field_data(phi, "phi");
    if(print_field.pressure)
      writer -> write_field_data(Pressure, "pressure");
    if(print_field.tau){
      writer -> write_field_data(strain[0], "tau_xx");
      writer -> write_field_data(strain[1], "tau_xy");
      writer -> write_field_data(strain[2], "tau_xz");
      writer -> write_field_data(strain[3], "tau_yy");
      writer -> write_field_data(strain[4], "tau_yz");
    }
  }
}

void Output_charge_field_data(double** zeta,
			      double* uk_dc,
			      double** Concentration,
			      Particle* p,
			      const CTime &time){
  if(print_field.vel)
    Zeta_k2u(zeta, uk_dc, u);
  
  double *potential = f_particle[0];
  double *dmy_value0 = f_particle[1];
  if(print_field.rho){
    Conc_k2charge_field(p, Concentration, potential, phi, dmy_value0);
    A2a_k(potential);
    Charge_field_k2Coulomb_potential_k_PBC(potential);
    A_k2a(potential);
    for(int n=0;n<N_spec;n++){
      A_k2a(Concentration[n]);
    }
  }//print rho?

  if(print_field.phi && !print_field.rho){
    Reset_phi(phi);
    Make_phi_particle(phi, p);
  }else if(print_field.phi && print_field.rho){
    Reset_phi(phi);
    Reset_phi(up[0]);
    Make_phi_qq_particle(phi, up[0], p);
  }//print phi/rho?

  if(print_field.rho){
    int im;
    double dmy;
    double dmy_surface_area = PI4*SQ(RADIUS);
    //compute total concentration
    for(int i = 0; i < NX; i++){
      for(int j = 0; j < NY; j++){
	for(int k = 0; k < NZ; k++){
	  im = (i*NY*NZ_) + (j*NZ_) + k;
	  dmy = 0.0;
	  for(int n = 0; n < N_spec; n++){
	    dmy += Elementary_charge*Valency[n]*Concentration[n][im];
	  }
	  up[0][im]*= dmy_surface_area;
	  up[1][im] = dmy*(1.0 - phi[im]);
	}
      }
    }
  }

  //Writers
  if(SW_OUTFORMAT == OUT_AVS_BINARY || SW_OUTFORMAT == OUT_AVS_ASCII){
    Output_avs_charge(Avs_parameters, u, phi, up[0], up[1], potential, time);
  }else if(SW_OUTFORMAT == OUT_EXT){
    if(print_field.vel){
      writer -> write_field_data(u[0], "u_x");
      writer -> write_field_data(u[1], "u_y");
      writer -> write_field_data(u[2], "u_z");
    }
    if(print_field.phi)
      writer -> write_field_data(phi, "phi");
    if(print_field.rho){
      writer -> write_field_data(up[0], "surface_charge");
      writer -> write_field_data(up[1], "rho");
      writer -> write_field_data(potential, "e_potential");
    }
  }

  //recover original state
  if(print_field.rho){
    for(int n=0; n<N_spec; n++){
      A2a_k(Concentration[n]);
    }
  }
}

//Todo: figure out how to print structures using hdf5 interface (avoid mem copies)!
void Output_particle_data(Particle* p, const CTime& time){
  if(SW_OUTFORMAT == OUT_AVS_BINARY || SW_OUTFORMAT == OUT_AVS_ASCII){
    Output_avs_p(Avs_parameters, p, time);
  }else{
    if(!print_particle.none)
      writer->write_particle_data(p);
  }
}

void Output_udf(UDFManager *ufout
                , double **zeta
                , double *uk_dc
                , const Particle *p
                , const CTime &time
		)
{
  ufout->newRecord();
  ufout->put("E", 1.0);
  ufout->put("t", time.ts);
  for(int j = 0; j < Particle_Number; j++) {
    char str[256];
    sprintf(str, "Particles[%d]", j);
    Location target(str);
    ufout->put(target.sub("R.x"), p[j].x[0]);
    ufout->put(target.sub("R.y"), p[j].x[1]);
    ufout->put(target.sub("R.z"), p[j].x[2]);
    ufout->put(target.sub("R_raw.x"), p[j].x_nopbc[0]);
    ufout->put(target.sub("R_raw.y"), p[j].x_nopbc[1]);
    ufout->put(target.sub("R_raw.z"), p[j].x_nopbc[2]);
    ufout->put(target.sub("v.x"), p[j].v[0]);
    ufout->put(target.sub("v.y"), p[j].v[1]);
    ufout->put(target.sub("v.z"), p[j].v[2]);

    qtn_isnormal(p[j].q);
    ufout->put(target.sub("q.q0"), qtn_q0(p[j].q));
    ufout->put(target.sub("q.q1"), qtn_q1(p[j].q));
    ufout->put(target.sub("q.q2"), qtn_q2(p[j].q));
    ufout->put(target.sub("q.q3"), qtn_q3(p[j].q));
    ufout->put(target.sub("omega.x"), p[j].omega[0]);
    ufout->put(target.sub("omega.y"), p[j].omega[1]);
    ufout->put(target.sub("omega.z"), p[j].omega[2]);

    ufout->put(target.sub("f_hydro.x"), p[j].f_hydro_previous[0]);
    ufout->put(target.sub("f_hydro.y"), p[j].f_hydro_previous[1]);
    ufout->put(target.sub("f_hydro.z"), p[j].f_hydro_previous[2]);
    ufout->put(target.sub("torque_hydro.x"), p[j].torque_hydro_previous[0]);
    ufout->put(target.sub("torque_hydro.y"), p[j].torque_hydro_previous[1]);
    ufout->put(target.sub("torque_hydro.z"), p[j].torque_hydro_previous[2]);

    ufout->put(target.sub("f_r.x"), p[j].fr_previous[0]);
    ufout->put(target.sub("f_r.y"), p[j].fr_previous[1]);
    ufout->put(target.sub("f_r.z"), p[j].fr_previous[2]);
    ufout->put(target.sub("torque_r.x"), 0.0);
    ufout->put(target.sub("torque_r.y"), 0.0);
    ufout->put(target.sub("torque_r.z"), 0.0);


    ufout->put(target.sub("f_slip.x"), p[j].f_slip_previous[0]);
    ufout->put(target.sub("f_slip.y"), p[j].f_slip_previous[1]);
    ufout->put(target.sub("f_slip.z"), p[j].f_slip_previous[2]);
    ufout->put(target.sub("torque_slip.x"), p[j].torque_slip_previous[0]);
    ufout->put(target.sub("torque_slip.y"), p[j].torque_slip_previous[1]);
    ufout->put(target.sub("torque_slip.z"), p[j].torque_slip_previous[2]);
  }
  if(SW_PT == rigid){
    for(int rigidID = 0; rigidID < Rigid_Number; rigidID++){
      char str[256];
      sprintf(str, "RigidParticles[%d]", rigidID);
      Location target(str);

      int rigid_first_n = Rigid_Particle_Cumul[rigidID];
      quaternion qGs;
      qtn_init(qGs, p[rigid_first_n].q);

      ufout->put(target.sub("R.x"), xGs[rigidID][0]);
      ufout->put(target.sub("R.y"), xGs[rigidID][1]);
      ufout->put(target.sub("R.z"), xGs[rigidID][2]);
      ufout->put(target.sub("R_raw.x"), xGs_nopbc[rigidID][0]);
      ufout->put(target.sub("R_raw.y"), xGs_nopbc[rigidID][1]);
      ufout->put(target.sub("R_raw.z"), xGs_nopbc[rigidID][2]);
      ufout->put(target.sub("v.x"), velocityGs[rigidID][0]);
      ufout->put(target.sub("v.y"), velocityGs[rigidID][1]);
      ufout->put(target.sub("v.z"), velocityGs[rigidID][2]);

      ufout->put(target.sub("q.q0"), qtn_q0(qGs));
      ufout->put(target.sub("q.q1"), qtn_q1(qGs));
      ufout->put(target.sub("q.q2"), qtn_q2(qGs));
      ufout->put(target.sub("q.q3"), qtn_q3(qGs));
      ufout->put(target.sub("omega.x"), omegaGs[rigidID][0]);
      ufout->put(target.sub("omega.y"), omegaGs[rigidID][1]);
      ufout->put(target.sub("omega.z"), omegaGs[rigidID][2]);

      ufout->put(target.sub("f_hydro.x"), forceGs_previous[rigidID][0]);
      ufout->put(target.sub("f_hydro.y"), forceGs_previous[rigidID][1]);
      ufout->put(target.sub("f_hydro.z"), forceGs_previous[rigidID][2]);
      ufout->put(target.sub("torque_hydro.x"), torqueGs_previous[rigidID][0]);
      ufout->put(target.sub("torque_hydro.y"), torqueGs_previous[rigidID][1]);
      ufout->put(target.sub("torque_hydro.z"), torqueGs_previous[rigidID][2]);

      ufout->put(target.sub("f_r.x"), forceGrs_previous[rigidID][0]);
      ufout->put(target.sub("f_r.y"), forceGrs_previous[rigidID][1]);
      ufout->put(target.sub("f_r.z"), forceGrs_previous[rigidID][2]);
      ufout->put(target.sub("torque_r.x"), torqueGrs_previous[rigidID][0]);
      ufout->put(target.sub("torque_r.y"), torqueGrs_previous[rigidID][1]);
      ufout->put(target.sub("torque_r.z"), torqueGrs_previous[rigidID][2]);

      ufout->put(target.sub("f_slip.x"), 0.0);
      ufout->put(target.sub("f_slip.y"), 0.0);
      ufout->put(target.sub("f_slip.z"), 0.0);
      ufout->put(target.sub("torque_slip.x"), 0.0);
      ufout->put(target.sub("torque_slip.y"), 0.0);
      ufout->put(target.sub("torque_slip.z"), 0.0);
    }
  }
}

