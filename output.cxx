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

  //legacy AVS support
  if(SW_OUTFORMAT == OUT_AVS_BINARY ||  SW_OUTFORMAT == OUT_AVS_ASCII){
    fprintf(stderr, "#AVS: set\n");
    Set_avs_parameters(Avs_parameters);
    fprintf(stderr, "#AVS: init\n");
    Init_avs(Avs_parameters);
    fprintf(stderr, "#AVS: init_p\n");
    if(Particle_Number > 0){
      Init_avs_p(Avs_parameters);
    }
    fprintf(stderr, "#AVS: finish\n");
  }else if(SW_OUTFORMAT == OUT_EXT){
    if(SW_EXTFORMAT == EXT_OUT_HDF5){
      writer = new hdf5_writer(NX, NY, NZ, NZ_,  //field dimensions
			       Particle_Number  //particle dimensions
			       );
    }
  }
}
void Free_output(){
  if(SW_OUTFORMAT == OUT_EXT){
    writer->~output_writer();
  }
}
void Show_output_parameter(){
  fprintf(stderr, "#Output parameters\n");
  if(SW_OUTFORMAT == OUT_AVS_BINARY || SW_OUTFORMAT == OUT_AVS_ASCII){
    Show_avs_parameter();
  }else if(SW_OUTFORMAT == OUT_EXT){
    fprintf(stderr, "#AVS output is suppressed.\n");
    fprintf(stderr, "#Extened output enabled.\n");
    writer -> show_parameter();
  }else{
    fprintf(stderr, "#Field/Particle output is disabled.\n");
  }
  
  if(SW_UDF){
    fprintf(stderr, "#for UDF ->\t%s\n",Out_udf);
  }else {
    fprintf(stderr, "#UDF output is supressed.\n");
  }
}

void Output_data(double** zeta,
		 double* uk_dc,
		 Particle* p,
		 const CTime &time){


  double *strain[QDIM]={f_particle[0]
		      ,f_particle[1]
		      ,f_particle[2]
		      ,f_ns0[0]
		      ,f_ns0[1]
  };

  if(SW_EQ == Shear_Navier_Stokes_Lees_Edwards){
    Zeta_k2u_k_OBL(zeta, uk_dc, u);
    U_k2Strain_k_OBL(u, strain);
    for(int d = 0; d < QDIM; d++){
      A_k2a(strain[d]);
    }
    U_k2u(u);

    E_oblique2E(strain, false); //without mean shear flow terms
    U_oblique2u(u);             //with mean shear flow terms
  }else{
    Zeta_k2u_k(zeta, uk_dc, u);
    U_k2Strain_k(u, strain);
    for(int d = 0; d < QDIM; d++){
      A_k2a(strain[d]);
    }
    U_k2u(u);
  }

  A_k2a(Pressure); //! TODO: implement pressure calculation
  {
    Reset_phi(phi);
    if (SW_EQ == Shear_Navier_Stokes_Lees_Edwards) {
      Make_phi_particle_OBL(phi, p);
    }else {
      Make_phi_particle(phi, p);
    }
  }

  if(SW_OUTFORMAT == OUT_AVS_BINARY || SW_OUTFORMAT == OUT_AVS_BINARY){
    Output_avs(Avs_parameters, u, phi, Pressure, strain, time);
    if(Particle_Number > 0){
      Output_avs_p(Avs_parameters, p, time);
    }
  }else if(SW_OUTFORMAT == OUT_EXT){
    writer -> write_start(time);
    writer -> write_field_data(u, DIM);
    writer -> write_field_data(phi);
    writer -> write_field_data(Pressure);
    writer -> write_field_data(strain, QDIM);
    if(Particle_Number > 0) 
      writer -> write_particle_data(p);
    writer -> write_end();
  }
}

void Output_charge_data(double** zeta,
			double* uk_dc,
			double** Concentration,
			Particle* p,
			const CTime &time){
  Zeta_k2u(zeta, uk_dc, u);
  
  double *potential = f_particle[0];
  double *dmy_value0 = f_particle[1];
  {
    Conc_k2charge_field(p, Concentration, potential, phi, dmy_value0);
    A2a_k(potential);
    Charge_field_k2Coulomb_potential_k_PBC(potential);
    A_k2a(potential);
    for(int n=0;n<N_spec;n++){
      A_k2a(Concentration[n]);
    }
  }
  {
    Reset_phi(phi);
    Reset_phi(up[0]);
    Make_phi_qq_particle(phi, up[0], p);
  }
  {
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
  if(SW_OUTFORMAT == OUT_AVS_BINARY || SW_OUTFORMAT == OUT_AVS_BINARY){
    Output_avs_charge(Avs_parameters, u, phi, up[0], up[1], potential, time);
    if(Particle_Number > 0){
      Output_avs_p(Avs_parameters, p, time);
    }
  }else if(SW_OUTFORMAT == OUT_EXT){
    writer -> write_start(time);
    writer -> write_end();
  }

  //recover original state
  for(int n=0; n<N_spec; n++){
    A2a_k(Concentration[n]);
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

