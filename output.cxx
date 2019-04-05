/*!
  \file output.cxx
  \brief Wrapper rouintes to write output files
  \details Prepare all data for output, leave actual writing to specialized class
  \author J. Molina
  \date 2014/05/21
  \version 1.0
 */

#include "output.h"
#ifdef WITH_EXTOUT
output_writer *writer;    //pointer to base writer
hdf5_writer   *h5writer;  //pointer to hdf5 writer
#endif 
void Init_output(Particle *p){
#ifdef WITH_EXTOUT
    writer = h5writer = NULL;
#endif
    if(SW_OUTFORMAT == OUT_AVS_BINARY ||  SW_OUTFORMAT == OUT_AVS_ASCII){ // LEGACY AVS
        print_field.none = false;
        print_field.vel = print_field.phi  = print_field.charge = print_field.pressure = print_field.tau  = true;

        for(int d = 0; d < DIM; d++){
            print_field_crop.start[d] = 0;
            print_field_crop.count[d] = Ns[d];
            print_field_crop.stride[d]= 1;
        }

        Set_avs_parameters(Avs_parameters);
        Init_avs(Avs_parameters);
        if(Particle_Number > 0) Init_avs_p(Avs_parameters);
    }else if(SW_OUTFORMAT == OUT_EXT){
#ifdef WITH_EXTOUT
    //Setup extended parameters

    //field selection
        if(SW_EQ == Electrolyte){
            print_field.pressure = false;
            print_field.tau      = false;
        }else{
            print_field.charge   = false;
        }
        print_field.none = (print_field.vel || print_field.phi || print_field.pressure || print_field.tau || print_field.charge ? false : true);

        // Initialize writers
        if(SW_EXTFORMAT == EXT_OUT_HDF5){
        //create particle and obstacle lists
            if (procid == root) {
                std::vector<int> plist, olist;
                plist.reserve(Particle_Number);
                olist.reserve(Particle_Number);
                for(int i = 0; i < Particle_Number; i++){
                    if(janus_propulsion[p[i].spec] != obstacle){
                        plist.push_back(i);
                    }else{
                        olist.push_back(i);
                    }
                }
                h5writer = new hdf5_writer(NX, NY, NZ, NZ_, DX, Particle_Number, DT*GTS, Out_dir, Out_name, print_field_crop, print_field, plist, olist, p);
                std::vector<int>().swap(plist);
                std::vector<int>().swap(olist);
                writer = static_cast<output_writer*>(h5writer);
            }
        }
#endif
    }
}

void Free_output(){
#ifdef WITH_EXTOUT
    if(SW_OUTFORMAT == OUT_EXT){
        if (procid == root) {
            writer->~output_writer();
        }
    }
#endif
}

void Show_output_parameter(){
    if(SW_UDF){
        fprintf_single(stderr, "# UDF output is enabled\n");
        fprintf_single(stderr, "# for UDF ->\t%s\n",Out_udf);
    }else {
        fprintf_single(stderr, "# UDF output is supressed.\n");
    }

    if(SW_OUTFORMAT == OUT_AVS_BINARY || SW_OUTFORMAT == OUT_AVS_ASCII){
        Show_avs_parameter();
    }else if(SW_OUTFORMAT == OUT_EXT){
#ifdef WITH_EXTOUT
        if (procid == root) {
            writer -> show_parameter();
		}
#endif
    }else{
        fprintf_single(stderr, "# Field/Particle output is disabled.\n");
    }
}

void Output_open_frame(){
#ifdef WITH_EXTOUT
    if(SW_OUTFORMAT == OUT_EXT){
        if (procid == root) {
            writer -> write_start();

        //hdf5 specific options : write data to current frame
            if(SW_EXTFORMAT == EXT_OUT_HDF5){
                if(SW_EQ == Shear_Navier_Stokes || 
                    SW_EQ == Shear_Navier_Stokes_Lees_Edwards){
                    float dmy_float = static_cast<float>(degree_oblique);
                    h5writer -> write_frame_attributes("gamma", &dmy_float, 1);
                }
            }
        }
	}
#endif
}

void Output_close_frame(){
#ifdef WITH_EXTOUT
	if(SW_OUTFORMAT == OUT_EXT){
        if (procid == root) {
            writer -> write_end();
		}
	}
#endif
}

void Output_field_data(double** zeta, double* uk_dc, Particle* p, const CTime &time){
    if(print_field.none) return;
    //double *stress[QDIM]={f_particle[0], f_particle[1], f_particle[2], f_ns0[0], f_ns0[1]};
    double **stress;
	stress = calloc_2d_double (5, mesh_size);
    Copy_v1(stress[0], f_particle[0]);
    Copy_v1(stress[1], f_particle[1]);
    Copy_v1(stress[2], f_particle[2]);
    Copy_v1(stress[3], f_ns0[0]);
    Copy_v1(stress[4], f_ns0[1]);

    if(print_field.vel || print_field.tau){
        if(SW_EQ == Shear_Navier_Stokes_Lees_Edwards){
            Zeta_k2u_k_OBL(zeta, uk_dc, u);
            if(print_field.tau){
                U_k2Stress_k_OBL(u, stress);
                for(int d = 0; d < QDIM; d++){
                    A_k2a(stress[d]);
                }
                //Stress_oblique2Stress(stress, false); //without mean shear flow terms
            }//print stress ?
            if(print_field.vel){
                U_k2u(u);
                U_oblique2u(u);      //with mean shear flow terms
            }//print u?
        }else{
            Zeta_k2u_k(zeta, uk_dc, u);
            if(print_field.tau){
                U_k2Stress_k(u, stress);
                for(int d = 0; d < QDIM; d++){
                    A_k2a(stress[d]);
                }
            }//print stress?
            if(print_field.vel){
                U_k2u(u);
            }//print u?
        }
    }// Print u / stress ?
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
        Output_avs(Avs_parameters, u, phi, Pressure, stress, time);
    }else if(SW_OUTFORMAT == OUT_EXT){
#ifdef WITH_EXTOUT
		double** u_tmp;
        double* Pressure_tmp;
        double* phi_tmp;
        double** stress_tmp;
        if(procid == root) {
            u_tmp = calloc_2d_double (3, NX * NY * NZ_);
            phi_tmp = calloc_1d_double (NX * NY * NZ_);
            Pressure_tmp = calloc_1d_double (NX * NY * NZ_);
            stress_tmp = calloc_2d_double (5, NX * NY * NZ_);
        } else {
            u_tmp = (double **) malloc ( 3 * sizeof (double *) );
            phi_tmp = (double *) malloc ( sizeof (double) );
            Pressure_tmp = (double *) malloc ( sizeof (double) );
            stress_tmp = (double **) malloc ( 5 * sizeof (double *) );
            //dmy_w = (double **) malloc ( 10 * sizeof (double *) );
        }
        Get_mesh_array (u[0], u_tmp[0], SW_OFF);
        Get_mesh_array (u[1], u_tmp[1], SW_OFF);
        Get_mesh_array (u[2], u_tmp[2], SW_OFF);
        Get_mesh_array (phi, phi_tmp, SW_OFF);
        Get_mesh_array (Pressure, Pressure_tmp, SW_OFF);
        Get_mesh_array (stress[0], stress_tmp[0], SW_OFF);
        Get_mesh_array (stress[1], stress_tmp[1], SW_OFF);
        Get_mesh_array (stress[2], stress_tmp[2], SW_OFF);
        Get_mesh_array (stress[3], stress_tmp[3], SW_OFF);
        Get_mesh_array (stress[4], stress_tmp[4], SW_OFF);
        if (procid == root) {
            //writer->write_field_data(u, phi, Pressure, stress);
            writer->write_field_data(u_tmp, phi_tmp, Pressure_tmp, stress_tmp);
            free_2d_double (stress_tmp);
            free_1d_double (Pressure_tmp);
            free_1d_double (phi_tmp);
            free_2d_double (u_tmp);
		} else {
            free (stress_tmp);
            free (Pressure_tmp);
            free (phi_tmp);
            free (u_tmp);
		}
#endif
	}
    free_2d_double (stress);
}

void Output_charge_field_data(double** zeta, double* uk_dc, double** Concentration, Particle* p, const CTime &time){
    if(print_field.none) return;
    if(print_field.vel) Zeta_k2u(zeta, uk_dc, u);

    double *potential = f_particle[0];
    double *dmy_value0 = f_particle[1];
    if(print_field.charge){
        Conc_k2charge_field(p, Concentration, potential, phi, dmy_value0);
        A2a_k(potential);
        Charge_field_k2Coulomb_potential_k_PBC(potential);
        A_k2a(potential);
        for(int n=0;n<N_spec;n++){
            A_k2a(Concentration[n]);
        }
    }//print charge ?

    if(print_field.phi && !print_field.charge){
        Reset_phi(phi);
        Make_phi_particle(phi, p);
    }else if(print_field.phi && print_field.charge){
        Reset_phi(phi);
        Reset_phi(up[0]);
        Make_phi_qq_particle(phi, up[0], p);
    }//print phi/charge?

    if(print_field.charge){
        int im;
        double dmy;
        double dmy_surface_area = PI4*SQ(RADIUS);
        //compute total concentration
        for (int i = 0; i < NPs[REAL][0]; i++) {
            for (int j = 0; j < NPs[REAL][1]; j++) {
                for (int k = 0; k < NPs[REAL][2]; k++) {
                    im = REALMODE_ARRAYINDEX(i, j, k);
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
#ifdef WITH_EXTOUT
		double** u_tmp;
        double* phi_tmp;
        double* up_0_tmp;
        double* up_1_tmp;
        double* potential_tmp;
        if(procid == root) {
            u_tmp = calloc_2d_double (3, NX * NY * NZ_);
            phi_tmp = calloc_1d_double (NX * NY * NZ_);
            up_0_tmp = calloc_1d_double (NX * NY * NZ_);
            up_1_tmp = calloc_1d_double (NX * NY * NZ_);
            potential_tmp = calloc_1d_double (NX * NY * NZ_);
        } else {
            u_tmp = (double **) malloc ( 3 * sizeof (double *) );
            phi_tmp = (double *) malloc ( sizeof (double) );
            up_0_tmp = (double *) malloc ( sizeof (double) );
            up_1_tmp = (double *) malloc ( sizeof (double) );
            potential_tmp = (double *) malloc ( sizeof (double) );
            //dmy_w = (double **) malloc ( 10 * sizeof (double *) );
        }
        Get_mesh_array (u[0], u_tmp[0], SW_OFF);
        Get_mesh_array (u[1], u_tmp[1], SW_OFF);
        Get_mesh_array (u[2], u_tmp[2], SW_OFF);
        Get_mesh_array (phi, phi_tmp, SW_OFF);
        Get_mesh_array (up[0], up_0_tmp, SW_OFF);
        Get_mesh_array (up[1], up_1_tmp, SW_OFF);
        Get_mesh_array (potential, potential_tmp, SW_OFF);
        if (procid == root) {
            //writer -> write_charge_field_data(u, phi, up[0], up[1], potential);
            writer -> write_charge_field_data(u_tmp, phi_tmp, up_0_tmp, up_1_tmp, potential_tmp);
            free_1d_double (potential_tmp);
            free_1d_double (up_1_tmp);
            free_1d_double (up_0_tmp);
            free_1d_double (phi_tmp);
            free_2d_double (u_tmp);
		} else {
            free (potential_tmp);
            free (up_1_tmp);
            free (up_0_tmp);
            free (phi_tmp);
            free (u_tmp);
		}
#endif
    }

    //recover original state
    if(print_field.charge){
        for(int n=0; n<N_spec; n++){
            A2a_k(Concentration[n]);
        }
    }
}

void Output_particle_data(Particle* p, const CTime& time){
    if(Particle_Number == 0) return;
    if(SW_OUTFORMAT == OUT_AVS_BINARY || SW_OUTFORMAT == OUT_AVS_ASCII){
        Output_avs_p(Avs_parameters, p, time);
    }else{
#ifdef WITH_EXTOUT
        Particle_Gather (p, p_tmp, SW_OFF);
        Particle_qsort (p_tmp, Particle_Number);
        if (procid == root) {
            writer->write_particle_data(p_tmp);
            writer->write_obstacle_data(p_tmp);
		}
        //writer->write_particle_data(p);
        //writer->write_obstacle_data(p);
#endif
    }
}

//Legacy udf output
void Output_udf(UDFManager *ufout, double **zeta, double *uk_dc, /*const*/ Particle *p, const CTime &time){
#if defined (_MPI)
    Particle_Gather (p, p_tmp, SW_OFF);
    Particle_qsort (p_tmp, Particle_Number);
#else
    for (int j = 0; j < Particle_Number; j++) {
        p_tmp[j] = p[j];
    }
#endif
    if (procid == root) {
        ufout->newRecord();
        ufout->put("E", 1.0);
        ufout->put("t", time.ts);
        for(int j = 0; j < Particle_Number; j++) {
            char str[256];
            sprintf(str, "Particles[%d]", j);
            Location target(str);
            ufout->put(target.sub("R.x"), p_tmp[j].x[0]);
            ufout->put(target.sub("R.y"), p_tmp[j].x[1]);
            ufout->put(target.sub("R.z"), p_tmp[j].x[2]);
            ufout->put(target.sub("R_raw.x"), p_tmp[j].x_nopbc[0]);
            ufout->put(target.sub("R_raw.y"), p_tmp[j].x_nopbc[1]);
            ufout->put(target.sub("R_raw.z"), p_tmp[j].x_nopbc[2]);
            ufout->put(target.sub("v.x"), p_tmp[j].v[0]);
            ufout->put(target.sub("v.y"), p_tmp[j].v[1]);
            ufout->put(target.sub("v.z"), p_tmp[j].v[2]);

            qtn_isnormal(p_tmp[j].q);
            ufout->put(target.sub("q.q0"), qtn_q0(p_tmp[j].q));
            ufout->put(target.sub("q.q1"), qtn_q1(p_tmp[j].q));
            ufout->put(target.sub("q.q2"), qtn_q2(p_tmp[j].q));
            ufout->put(target.sub("q.q3"), qtn_q3(p_tmp[j].q));
            ufout->put(target.sub("omega.x"), p_tmp[j].omega[0]);
            ufout->put(target.sub("omega.y"), p_tmp[j].omega[1]);
            ufout->put(target.sub("omega.z"), p_tmp[j].omega[2]);

            ufout->put(target.sub("f_hydro.x"), p_tmp[j].f_hydro_previous[0]);
            ufout->put(target.sub("f_hydro.y"), p_tmp[j].f_hydro_previous[1]);
            ufout->put(target.sub("f_hydro.z"), p_tmp[j].f_hydro_previous[2]);
            ufout->put(target.sub("torque_hydro.x"), p_tmp[j].torque_hydro_previous[0]);
            ufout->put(target.sub("torque_hydro.y"), p_tmp[j].torque_hydro_previous[1]);
            ufout->put(target.sub("torque_hydro.z"), p_tmp[j].torque_hydro_previous[2]);

            ufout->put(target.sub("f_r.x"), p_tmp[j].fr_previous[0]);
            ufout->put(target.sub("f_r.y"), p_tmp[j].fr_previous[1]);
            ufout->put(target.sub("f_r.z"), p_tmp[j].fr_previous[2]);
            ufout->put(target.sub("torque_r.x"), 0.0);
            ufout->put(target.sub("torque_r.y"), 0.0);
            ufout->put(target.sub("torque_r.z"), 0.0);

            ufout->put(target.sub("f_slip.x"), p_tmp[j].f_slip_previous[0]);
            ufout->put(target.sub("f_slip.y"), p_tmp[j].f_slip_previous[1]);
            ufout->put(target.sub("f_slip.z"), p_tmp[j].f_slip_previous[2]);
            ufout->put(target.sub("torque_slip.x"), p_tmp[j].torque_slip_previous[0]);
            ufout->put(target.sub("torque_slip.y"), p_tmp[j].torque_slip_previous[1]);
            ufout->put(target.sub("torque_slip.z"), p_tmp[j].torque_slip_previous[2]);
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
}

