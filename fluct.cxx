/*!
  \file fluct.cxx
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Routines to compute random thermal fluctuation forces
 */
#include "fluct.h"
#include "mt19937ar.h"

void Add_random_force_thermostat(Particle *p, const CTime &jikan){
    static const double Zeta_drag = 6.* M_PI * ETA * RADIUS;
    static const double Zeta_drag_rot = 8.* M_PI * ETA * POW3(RADIUS);

    const double sdv_v = sqrt(Zeta_drag * jikan.dt_md * kBT * alpha_v);
    const double sdv_omega = sqrt(Zeta_drag_rot * jikan.dt_md * kBT * alpha_o);

    static double noise_intensity_v = kT_snap_v * sdv_v;
    static double noise_intensity_o = kT_snap_o * sdv_omega;

    if(SW_PT != rigid){
        generate_random_number_MT(rand_num_particle, 3);

        for(int n=0; n<Update_Particle_Number; n++){
            if(janus_propulsion[p[n].spec] != obstacle){

                double imass = IMASS[p[n].spec];
                double imoi = IMOI[p[n].spec];
                for(int d=0; d<DIM; d++){
                    p[n].v[d] += rand_num_particle[6*(p[n].id)+d]*noise_intensity_v*imass;
                    if(ROTATION){
                        p[n].omega[d] += rand_num_particle[6*(p[n].id)+3+d]*noise_intensity_o*imoi;
                    }
                }
            }
        }
    }else{ // if rigid
        int rigidID;
        double fdt[3], ndt[3], **forceGsdt, **torqueGsdt;
        forceGsdt = alloc_2d_double(Rigid_Number, DIM);
        torqueGsdt = alloc_2d_double(Rigid_Number, DIM);
        for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
            for(int d=0; d<DIM; d++){
                forceGsdt[rigidID][d] = 0.0;
                torqueGsdt[rigidID][d] = 0.0;
            }
        }

        generate_random_number_MT(rand_num_particle, 3);
        for(int n=0; n<Update_Particle_Number; n++){

            rigidID = Particle_RigidID[n];

            for(int d=0; d<DIM; d++){
                fdt[d] = rand_num_particle[6*(p[n].id)+d]*noise_intensity_v;
                ndt[d] = rand_num_particle[6*(p[n].id)+3+d]*noise_intensity_o;
                forceGsdt[rigidID][d] += fdt[d];
                torqueGsdt[rigidID][d] += ndt[d];  // If same forces exist on a particle's surface,
                                                   // sum of torques around gravity point of constituent particles and
                                                   // the torque around gravity point of the rigid particle are same.
            }
            torqueGsdt[rigidID][0] += GRvecs[n][1] * fdt[2] - GRvecs[n][2] * fdt[1];
            torqueGsdt[rigidID][1] += GRvecs[n][2] * fdt[0] - GRvecs[n][0] * fdt[2];
            torqueGsdt[rigidID][2] += GRvecs[n][0] * fdt[1] - GRvecs[n][1] * fdt[0];
        }

        int rigid_spec;
        for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
            rigid_spec = RigidID_Components[rigidID];

            for(int d1=0; d1<DIM; d1++){
                if(Rigid_Motions_vel[rigid_spec][d1] == 0) continue; // if "fix"
                velocityGs[rigidID][d1] += forceGsdt[rigidID][d1] * Rigid_IMasses[rigidID];
            }

            for(int d1=0; d1<DIM; d1++){
                if(Rigid_Motions_omega[rigid_spec][d1] == 0) continue; // if "fix"
                for(int d2=0; d2<DIM; d2++) omegaGs[rigidID][d1] += Rigid_IMoments[rigidID][d1][d2] * torqueGsdt[rigidID][d2];
            }
        }
        for(int n=0; n<Particle_Number; n++){
            rigidID = Particle_RigidID[n];
            p[n].v[0] = velocityGs[rigidID][0] + omegaGs[rigidID][1]*GRvecs[n][2] - omegaGs[rigidID][2]*GRvecs[n][1];
            p[n].v[1] = velocityGs[rigidID][1] + omegaGs[rigidID][2]*GRvecs[n][0] - omegaGs[rigidID][0]*GRvecs[n][2];
            p[n].v[2] = velocityGs[rigidID][2] + omegaGs[rigidID][0]*GRvecs[n][1] - omegaGs[rigidID][1]*GRvecs[n][0];
        }

        free_2d_double(forceGsdt);
        free_2d_double(torqueGsdt);
    }
}

void Force_random_walk(Particle *p){
    Calc_f_Lennard_Jones(p);
}

void Random_Walk(Particle *p, const double dr_factor){
    const double dr = RADIUS * dr_factor;

    generate_random_number_MT(rand_num_particle, 2);

    for(int n=0; n<Update_Particle_Number; n++){
        double dmy[4];
        //Gauss2(dmy);
        //Gauss2(dmy+2);
        dmy[0] = rand_num_particle[4*p[n].id+0];
        dmy[1] = rand_num_particle[4*p[n].id+1];
        dmy[2] = rand_num_particle[4*p[n].id+2];
        dmy[3] = rand_num_particle[4*p[n].id+3];
        double h = 0.;
        for(int d=0; d<DIM; d++){  
            h += SQ(dmy[d]);
        }
        h = sqrt(h);
        if(h> 0.){    
            h = dr/h;
        }
        for(int d=0; d<DIM; d++){      
	        p[n].x[d] += h * dmy[d];
            p[n].x[d] -= floor(p[n].x[d]*iL_particle[d])*L_particle[d];
            assert(p[n].x[d] >= 0.);
            assert(p[n].x[d] < L_particle[d]);
        }
    }
#ifdef _MPI
    Particle_Group_Communication (p, MANY_TO_MANY);
#endif
}

void Steepest_descent(Particle *p, const double dr_factor){
    const double dr = RADIUS * dr_factor;

    Force_random_walk(p);
    for(int n=0; n<Update_Particle_Number; n++){
        double h = 0.;
        for(int d=0; d<DIM; d++){  
            h += SQ(p[n].fr[d]);
        }
        h = sqrt(h);
        if(h> 0.){    
            h = dr/h;
        }
        for(int d=0; d<DIM; d++){      
            p[n].x[d] += h * p[n].fr[d];
            p[n].x[d] -= floor(p[n].x[d]*iL_particle[d])*L_particle[d];

            assert(p[n].x[d] >= 0.);
            assert(p[n].x[d] < L_particle[d]);
        }
    }
#ifdef _MPI
    Particle_Group_Communication (p, MANY_TO_MANY);
#endif
}

void generate_random_number_MT(double* val, int num){
    double dmy[2];
    int i = 0;
	for(int n=0; n < Particle_Number; n++){
		for(int m=0; m < num; m++){
            Gauss2(dmy);
            val[i] = dmy[0];
            val[i+1] = dmy[1];
            i+=2;
		}
	}
}