/*!
  \file fluct.cxx
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Routines to compute random thermal fluctuation forces
 */
#include "fluct.h"

void Add_random_force_thermostat(Particle *p, const CTime &jikan){

  double* noise_intensity_v;
  double* noise_intensity_o;

  noise_intensity_v = alloc_1d_double(Component_Number);
  noise_intensity_o = alloc_1d_double(Component_Number);

  for(int i = 0 ; i < Component_Number ; i++){
	  double Zeta_drag = 6.* M_PI * ETA * RADII[i];
  	double Zeta_drag_rot = 8.* M_PI * ETA * POW3(RADII[i]);
	  
  	double sdv_v = sqrt(Zeta_drag * jikan.dt_md * kBT * alpha_v);
	  double sdv_omega = sqrt(Zeta_drag_rot * jikan.dt_md * kBT * alpha_o);
   
  	noise_intensity_v[i] = kT_snap_v * sdv_v;
  	noise_intensity_o[i] = kT_snap_o * sdv_omega;
  }
  if(SW_PT != rigid){
#pragma omp parallel for
    for(int n=0; n<Particle_Number; n++){
      if(janus_propulsion[p[n].spec] != obstacle){
        double dmy[6];
        Gauss2(dmy);
        Gauss2(dmy+2);
        Gauss2(dmy+4);
        
        double imass = IMASS[p[n].spec];
        double imoi = IMOI[p[n].spec];
        
        for(int d=0; d<DIM; d++){
          p[n].v[d] += dmy[d]*noise_intensity_v[p[n].spec]*imass;
          if(ROTATION){
            p[n].omega[d] += dmy[d+3]*noise_intensity_o[p[n].spec]*imoi;
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
    for(int n=0; n<Particle_Number; n++){
      double dmy[6];
      Gauss2(dmy);
      Gauss2(dmy+2);
      Gauss2(dmy+4);
      
      rigidID = Particle_RigidID[n];
      
      for(int d=0; d<DIM; d++){
        fdt[d] = dmy[d]*noise_intensity_v[p[n].spec];
        ndt[d] = dmy[d+3]*noise_intensity_o[p[n].spec];
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
