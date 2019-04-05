/*!
  \file avs_output.cxx
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Output routines for field data in AVS/Express format
 */

#include "avs_output.h"

const int Veclen = 5+5; 
const char *Label="ux uy uz phi pressure tau_xy tau_yz tau_zx tau_xx tau_yy"; // avs 出力ラベル用
const int Veclen_charge = 4+3; 
const char *Label_charge="ux uy uz phi surface_charge rho e_potential"; // avs 出力ラベル用

//const AVS_Field Field = irregular;
const AVS_Field Field = uniform;


void Show_avs_parameter(){
    if (procid == root) {
        if(SW_OUTFORMAT == OUT_AVS_BINARY){
            fprintf(stderr, "#for AVS (filetype is binary)\n");
        }else if(SW_OUTFORMAT == OUT_AVS_ASCII){
            fprintf(stderr, "#for AVS (filetype is ascii)\n");
        }else{
            fprintf(stderr, "# Uknown AVS FORMAT\n");
            exit_job(EXIT_FAILURE);
        }
        fprintf(stderr, "#directory:%s\n", Out_dir);
        fprintf(stderr, "# (mesh data)->\t{%s, %s, %s*.dat}\n"
                  ,Avs_parameters.out_fld
                  ,Avs_parameters.out_cod
                  ,Avs_parameters.out_pfx);
        if(Particle_Number > 0){
            fprintf(stderr, "# (particle data)->\t{%s, %s*.cod, %s*.dat}\n"
                  ,Avs_parameters.out_pfld
                  ,Avs_parameters.out_ppfx
                  ,Avs_parameters.out_ppfx);
        }
    }
}

/* change: Output data only in RANK ZERO */
void Init_avs(const AVS_parameters &Avs_parameters){
    if (procid == root) {

        char dmy_dir[256];
        sprintf(dmy_dir, "%s/avs", Out_dir);
        dircheckmake(dmy_dir);

        FILE *fout;
        fout=filecheckopen(Avs_parameters.fld_file,"w");
        fprintf(fout,"# AVS field file\n");
        fprintf(fout,"ndim=%d\n",DIM);
        fprintf(fout,"dim1=%d\n", Avs_parameters.nx);
        fprintf(fout,"dim2=%d\n", Avs_parameters.ny);
        fprintf(fout,"dim3=%d\n", Avs_parameters.nz);
        fprintf(fout,"nspace=%d\n", DIM);
        if(SW_EQ == Navier_Stokes ||
           SW_EQ == Shear_Navier_Stokes ||
           SW_EQ == Shear_Navier_Stokes_Lees_Edwards){
            fprintf(fout,"veclen=%d\n", Veclen);
        }else if(SW_EQ==Electrolyte){
            fprintf(fout,"veclen=%d\n", Veclen_charge);
        }
        fprintf(fout,"data=float\n");
        if(Field == irregular ){
            fprintf(fout,"field=irregular\n");
        }else if(Field == uniform ){
            fprintf(fout,"field=uniform\n");
        }else {
            fprintf(stderr, "invalid Field\n"); 
            exit_job(EXIT_FAILURE);
        }
        fprintf(fout,"nstep=%d\n",Avs_parameters.nstep);
        if(SW_EQ == Navier_Stokes ||
           SW_EQ == Shear_Navier_Stokes ||
           SW_EQ == Shear_Navier_Stokes_Lees_Edwards){
            fprintf(fout,"label = %s\n", Label);
        }else if(SW_EQ==Electrolyte){
            fprintf(fout,"label = %s\n", Label_charge);
        }

        fclose(fout);
  
        if(Field == irregular){
            fout=filecheckopen(Avs_parameters.cod_file,"wb");
            if(SW_OUTFORMAT == OUT_AVS_BINARY){
                for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
                    for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend ; j++){
                        for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend ; k++){
                            float dmy= (float)(i*DX);
                            fwrite(&dmy,sizeof(float),1,fout);
						}
                    }
                }
                for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
                    for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend ; j++){
                        for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend ; k++){
                            float dmy= (float)(j*DX);
                            fwrite(&dmy,sizeof(float),1,fout);
                        }
                    }
                }
                for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
                    for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend ; j++){
                        for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend ; k++){
                            float dmy= (float)(k*DX);
                            fwrite(&dmy,sizeof(float),1,fout);
                        }
                    }
                }
            }else if(SW_OUTFORMAT == OUT_AVS_ASCII){
                fprintf(fout,"X Y Z\n");
                for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
                    for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend ; j++){
                        for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend ; k++){
                            fprintf(fout,"%g %g %g\n", (double)i*DX, (double)j*DX, (double)k*DX);
                        }
                    }
                }
            }else{
                fprintf(stderr, "# Uknown AVS FORMAT\n");
                exit_job(EXIT_FAILURE);
            }
            fclose(fout);
        }else if(Field == uniform){
            fout=filecheckopen(Avs_parameters.cod_file,"wb");
        {
            fprintf(fout,"%g\n%g\n%g\n%g\n%g\n%g\n"
                        ,(double)(Avs_parameters.istart)*DX
                        ,(double)(Avs_parameters.iend)*DX
                        ,(double)(Avs_parameters.jstart)*DX
                        ,(double)(Avs_parameters.jend)*DX
                        ,(double)(Avs_parameters.kstart)*DX
                        ,(double)(Avs_parameters.kend)*DX);
        }
            fclose(fout);
        }else {
            fprintf(stderr, "invalid Field\n"); 
            exit_job(EXIT_FAILURE);
        }
	}
}

void Set_avs_parameters(AVS_parameters &Avs_parameters){

    Avs_parameters.nx=NX;
    Avs_parameters.ny=NY;
    Avs_parameters.nz=NZ;

    Avs_parameters.istart = 0;
    Avs_parameters.iend = NX-1;
    Avs_parameters.jstart = 0;
    Avs_parameters.jend = NY-1;
    Avs_parameters.kstart = 0;
    Avs_parameters.kend = NZ-1;

    sprintf(Avs_parameters.out_fld, "%s.fld", Out_name);
    sprintf(Avs_parameters.out_cod, "%s.cod", Out_name);
    sprintf(Avs_parameters.cod_file, "%s/%s",
    Out_dir, Avs_parameters.out_cod);

    sprintf(Avs_parameters.out_pfx, "avs/%s_", Out_name);
    sprintf(Avs_parameters.fld_file, "%s/%s", Out_dir, Avs_parameters.out_fld);

    sprintf(Avs_parameters.out_pfld, "%sp.fld", Out_name);
    //sprintf(Avs_parameters.out_pcod, "%sp.cod", Out_name);
    sprintf(Avs_parameters.out_ppfx, "avs/%sp_", Out_name);
    sprintf(Avs_parameters.pfld_file, "%s/%s", Out_dir, Avs_parameters.out_pfld);
    //sprintf(Avs_parameters.pcod_file, "%s/%s", Out_dir, Avs_parameters.out_pcod);

    Avs_parameters.nstep=(Num_snap+1);
}

inline void Binary_write(FILE *fout, AVS_parameters &Avs_parameters, const double *a){
    int im;
    for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend; k++){
        for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend; j++){
            for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
                im = (i * NY * NZ_) + (j * NZ_) + k;
                float dmy= (float)a[im];
                fwrite(&dmy,sizeof(float),1,fout);
            }
        }
    }
}

/* change: Output data only in RANK = ROOT */
inline void Add_field_description(AVS_parameters &Avs_parameters, const CTime &time, const int &veclen){
    FILE *fout;
    char line[512];
    if (procid == root) {
        fout=filecheckopen(Avs_parameters.fld_file,"a");
        fprintf(fout,"time value = \"step%dtime%g\"\n", time.ts, time.time);
        if(SW_OUTFORMAT == OUT_AVS_BINARY){
            static const int data_size=sizeof(float) * Avs_parameters.nx * Avs_parameters.ny * Avs_parameters.nz;
            if(Field == irregular){
                for(int n=0; n < DIM; n++){
                    fprintf(fout, "coord %d file = %s filetype = binary skip = %d\n", n+1, Avs_parameters.out_cod, n*data_size);
                }
            }else if(Field == uniform){
              {
                for(int n=0; n < DIM; n++){
                    fprintf(fout, "coord %d file = %s filetype = ascii skip = %d\n", n+1, Avs_parameters.out_cod, n*2);
                }
              }
            }
            for(int n=0; n < veclen; n++){
                fprintf(fout, "variable %d file = %s%d.dat filetype = binary skip = %d\n", n+1, Avs_parameters.out_pfx, time.ts, n * data_size);
            }
        }else if(SW_OUTFORMAT == OUT_AVS_ASCII){
            if(Field == irregular){
                for(int n=0; n < DIM; n++){
                    fprintf(fout, "coord %d file = %s filetype = ascii skip = 1 offset = %d stride =%d\n", n+1, Avs_parameters.out_cod, n, DIM);
                }
            }else if(Field == uniform){
                for(int n=0; n < DIM; n++){
                    fprintf(fout, "coord %d file = %s filetype = ascii skip = %d\n", n+1, Avs_parameters.out_cod, n*2);
                }
            }else {
                fprintf(stderr, "invalid Field\n"); 
                exit_job(EXIT_FAILURE);
            }
            for(int n=0; n < veclen; n++){
                fprintf(fout, "variable %d file = %s%d.dat filetype = ascii skip = 2 offset = %d stride = %d\n", n+1, Avs_parameters.out_pfx, time.ts, n, veclen);
            }
        }else{
            fprintf(stderr, "# Uknown AVS FORMAT\n");
            exit_job(EXIT_FAILURE);
        }
        fprintf(fout,"EOT\n");
        fclose(fout);
	}
}

void Output_avs(AVS_parameters &Avs_parameters, double **u, double *phi, double *Pressure, double **strain, const CTime &time){
    int im;
    double **dmy_w;

    Add_field_description(Avs_parameters,time, Veclen);

    FILE *fout;
    char line[512];
    sprintf(line,"timesteps=%d time=%f\n%s",time.ts,time.time, Label);
    sprintf(Avs_parameters.data_file,"%s/%s%d.dat", Out_dir, Avs_parameters.out_pfx, time.ts);

#if defined (_MPI)
    if(procid == root) {
        dmy_w = calloc_2d_double (10, NX * NY * NZ_);
    } else {
        dmy_w = (double **) malloc ( 10 * sizeof (double *) );
    }
    Get_mesh_array (u[0], dmy_w[0], SW_OFF);
    Get_mesh_array (u[1], dmy_w[1], SW_OFF);
    Get_mesh_array (u[2], dmy_w[2], SW_OFF);
    Get_mesh_array (phi, dmy_w[3], SW_OFF);
    Get_mesh_array (Pressure, dmy_w[4], SW_OFF);
    Get_mesh_array (strain[0], dmy_w[5], SW_OFF);
    Get_mesh_array (strain[1], dmy_w[6], SW_OFF);
    Get_mesh_array (strain[2], dmy_w[7], SW_OFF);
    Get_mesh_array (strain[3], dmy_w[8], SW_OFF);
    Get_mesh_array (strain[4], dmy_w[9], SW_OFF);
#else
    dmy_w = (double **) malloc (10 * sizeof (double *) );
    alloc_error_check (dmy_w);
    dmy_w[0] = u[0];
    dmy_w[1] = u[1];
    dmy_w[2] = u[2];
    dmy_w[3] = phi;
    dmy_w[4] = Pressure;
    dmy_w[5] = strain[0];
    dmy_w[6] = strain[1];
    dmy_w[7] = strain[2];
    dmy_w[8] = strain[3];
    dmy_w[9] = strain[4];
#endif
    if (procid == root) {
        fout=filecheckopen(Avs_parameters.data_file,"wb");

        if(SW_OUTFORMAT == OUT_AVS_BINARY){
            Binary_write (fout, Avs_parameters, dmy_w[0]);
            Binary_write (fout, Avs_parameters, dmy_w[1]);
            Binary_write (fout, Avs_parameters, dmy_w[2]);
            Binary_write (fout, Avs_parameters, dmy_w[3]);
            Binary_write (fout, Avs_parameters, dmy_w[4]);
            Binary_write (fout, Avs_parameters, dmy_w[6]); // 12
            Binary_write (fout, Avs_parameters, dmy_w[9]); // 23
            Binary_write (fout, Avs_parameters, dmy_w[7]); // 13
            Binary_write (fout, Avs_parameters, dmy_w[5]); // 11
            Binary_write (fout, Avs_parameters, dmy_w[8]); // 22
        }else if(SW_OUTFORMAT == OUT_AVS_ASCII){
            fprintf(fout,"%s\n", line);
            for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend; k++){
                for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend; j++){
	                for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
                        int im=(i*NY*NZ_)+(j*NZ_)+k;
                        fprintf(fout,"%.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g\n",
                            dmy_w[0][im], dmy_w[1][im], dmy_w[2][im],
                            dmy_w[3][im], dmy_w[4][im], dmy_w[6][im],
                            dmy_w[9][im], dmy_w[7][im], dmy_w[5][im], dmy_w[8][im]);
	                }
                }
            }
        }else{
            fprintf(stderr, "# Uknown AVS FORMAT\n");
            exit_job(EXIT_FAILURE);
        }
        fclose(fout);
#ifdef _MPI
        free_2d_double (dmy_w);
    } else {
        free (dmy_w);
#endif
    }
#ifndef _MPI
    free (dmy_w);
#endif
}

void Output_avs_charge(AVS_parameters &Avs_parameters, double** u, double* phi, double* colloid_charge
                          ,double* solute_charge_total, double* potential, const CTime &time){
    int im;
    double **dmy_w;

    Add_field_description(Avs_parameters,time, Veclen_charge);

    FILE *fout;
    char line[512];
    sprintf(line,"timesteps=%d time=%f\n%s",time.ts,time.time, Label_charge);
    sprintf(Avs_parameters.data_file,"%s/%s%d.dat", Out_dir, Avs_parameters.out_pfx, time.ts);
    double dmy_surface_area = PI4 * RADIUS * RADIUS;

#if defined (_MPI)
    if(procid == root) {
        dmy_w = calloc_2d_double (7, NX * NY * NZ_);
    } else {
        dmy_w = (double **) malloc (7 * sizeof (double *) );
    }
    Get_mesh_array (u[0], dmy_w[0], SW_OFF);
    Get_mesh_array (u[1], dmy_w[1], SW_OFF);
    Get_mesh_array (u[2], dmy_w[2], SW_OFF);
    Get_mesh_array (phi, dmy_w[3], SW_OFF);
    Get_mesh_array (colloid_charge, dmy_w[4], SW_OFF);
    Get_mesh_array (solute_charge_total, dmy_w[5], SW_OFF);
    Get_mesh_array (potential, dmy_w[6], SW_OFF);
#else
    dmy_w = (double **) malloc ( (6 + N_spec) * sizeof (double *) );
    alloc_error_check (dmy_w);
    dmy_w[0] = u[0];
    dmy_w[1] = u[1];
    dmy_w[2] = u[2];
    dmy_w[3] = phi;
    dmy_w[4] = colloid_charge;
    dmy_w[5] = solute_charge_total;
    dmy_w[6] = potential;

#endif
    if (procid == root) {	
        fout=filecheckopen(Avs_parameters.data_file,"wb");
        if(SW_OUTFORMAT == OUT_AVS_BINARY){
            Binary_write (fout, Avs_parameters, dmy_w[0]);
            Binary_write (fout, Avs_parameters, dmy_w[1]);
            Binary_write (fout, Avs_parameters, dmy_w[2]);
            Binary_write (fout, Avs_parameters, dmy_w[3]);
            Binary_write (fout, Avs_parameters, dmy_w[4]);
            Binary_write (fout, Avs_parameters, dmy_w[5]);
            Binary_write (fout, Avs_parameters, dmy_w[6]);
        }else if(SW_OUTFORMAT == OUT_AVS_ASCII){ // OUT_AVS_ASCII
            fprintf (fout, "%s\n", line);
            for (int k = Avs_parameters.kstart; k <= Avs_parameters.kend; k++) {
                for (int j = Avs_parameters.jstart; j <= Avs_parameters.jend; j++) {
                    for (int i = Avs_parameters.istart; i <= Avs_parameters.iend; i++) {
                        double dmy = 0.0;
                        im = (i * NY * NZ_) + (j * NZ_) + k;
                        fprintf(fout,"%.3g %.3g %.3g %.3g %.3g %.3g %.3g\n",
                                      u[0][im], u[1][im], u[2][im], phi[im],
                                      colloid_charge[im], solute_charge_total[im], potential[im]);
                    }
                }
			}
        }else{
            fprintf(stderr, "# Uknown AVS FORMAT\n");
            exit_job(EXIT_FAILURE);
        }
        fclose(fout);
#ifdef _MPI
        free_2d_double (dmy_w);
	} else {
        free (dmy_w);
#endif
    }
#ifndef _MPI
    free (dmy_w);
#endif
}

