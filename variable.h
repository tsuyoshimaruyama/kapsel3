/*!
  \file variable.h
  \brief Defines the global structs (CTime, Particle, Index_range)
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */
#ifndef VARIABLE_H
#define VARIABLE_H

#ifdef _MPI
#include <mpi.h>
#endif
#include <complex>
#include <stdlib.h>
#include "parameter_define.h"
#include "udfmanager.h"

#include "dc.h"

#ifdef _FFT_IMKL
#include <mkl_dfti.h>
#elif _FFT_FFTW
#include <fftw3.h>
#endif

// omega,ux,uy,phi, phi_up など場の変数を格納する
typedef double ** *Value;
typedef int ** *Value_int;

// type struct defines
typedef struct CTime {
    int ts; // time step
    double time; // time
    double dt_fluid;  // time increment
    double hdt_fluid; // 1/2 * dt
    double dt_md;  // time increment
    double hdt_md; // 1/2 * dt
} CTime;

typedef struct splineSystem{
    int      n;
    double   dx;
    double*  a;
    double*  b;
    double*  c;
    double*  d;
    double*  Q;
    double*  Aii;
    double*  Ain;
} splineSystem;

typedef struct AVS_parameters{
    int nx;
    int ny;
    int nz;
    char out_fld[128];
    char out_cod[128];
    char out_pfx[128];
    char fld_file[128];
    char cod_file[128];

    char data_file[128];

    char out_pfld[128];
    char out_ppfx[128];
    char pfld_file[128];
    int istart;
    int iend;
    int jstart;
    int jend;
    int kstart;
    int kend;
    int nstep;
} AVS_parameters;

typedef struct quaternion {
    double s; //scalar part
    double v[DIM]; //vector part
} quaternion;

typedef struct Particle {
    double mass;                       //fluid particle mass
    double surface_mass;               //surface fluid mass

    int spec;
    int id;

    double x[DIM];
    double x_previous[DIM];
    double x_nopbc[DIM];

    double v[DIM];
    double v_old[DIM];
    double v_slip[DIM];

    double f_hydro[DIM];
    double f_hydro_previous[DIM];
    double f_hydro1[DIM];
    double f_slip[DIM];
    double f_slip_previous[DIM];

    double fr[DIM];
    double fr_previous[DIM];

    double torque_r[DIM];
    double torque_r_previous[DIM];

    double omega[DIM];
    double omega_old[DIM];
    double omega_slip[DIM];

    double torque_hydro[DIM];
    double torque_hydro_previous[DIM];
    double torque_hydro1[DIM];
    double torque_slip[DIM];
    double torque_slip_previous[DIM];

    double momentum_depend_fr[DIM];

    double mass_center[DIM];           //center of mass of fluid particle


    double surface_mass_center[DIM];   //surface fluid center of mass

    double surface_dv[DIM];            //momentum change due to slip
    double surface_dw[DIM];        //ang. momentum change due to slip

    drand48_data rdata; //random routine library(Dynamic Creator) parameter

    double QR[DIM][DIM];
    double QR_old[DIM][DIM];
    double inertia[DIM][DIM];          //moment of inertia of fluid particle
    double surface_inertia[DIM][DIM];  //surface fluid moment of inertia

    quaternion q;
    quaternion q_old;
} Particle;

typedef struct Index_range {
    int istart;
    int iend;
    int jstart;
    int jend;
    int kstart;
    int kend;
} Index_range;

typedef struct Field_crop{
    int start[DIM];
    int count[DIM];
    int stride[DIM];
} Field_crop;

typedef struct Field_output{
    bool none;
    bool vel;
    bool phi;
    bool charge;
    bool pressure;
    bool tau;
} Field_output;

typedef struct Particle_output{
    bool none;
    int start;
    int count;
    int stride;
} Particle_output;

typedef std::complex <double> Complex;

//type enum defines
enum Count_SW { INIT, ADD, MEAN, SNAP_MEAN, SHOW };
enum DIVISION { divX = 0, divY = 1 };
enum EQ { Navier_Stokes, Shear_Navier_Stokes, Shear_Navier_Stokes_Lees_Edwards, Electrolyte };
enum Particle_BC { PBC_particle, Lees_Edwards, Shear_hydro };
enum Particle_IC { None, uniform_random, random_walk, FCC, BCC, user_specify };
enum Particle_IO { random_dir, space_dir, user_dir };
enum PT { spherical_particle, chain, rigid };
extern const int GIVEN_SEED;
extern const int RANDOM_SEED;
enum SW_SWITCH { SW_OFF = 0, SW_ON = 1 };
enum SW_COMM { ONE_TO_ONE_LJ, ONE_TO_ONE_SEKIBUN, ONE_TO_MANY, MANY_TO_MANY };
enum SW_time { AUTO, MANUAL };
enum SW_MODE { REAL = 0, SPECTRUM = 1 };
enum JAX {x_axis, y_axis, z_axis, no_axis};
enum JP  {motor,slip,obstacle,no_propulsion};
enum OBL_INT {linear_int, spline_int};
enum OBL_TRANSFORM {oblique2cartesian, cartesian2oblique};
enum OUTFORMAT{OUT_NONE, OUT_AVS_ASCII, OUT_AVS_BINARY, OUT_EXT};
enum EXTFORMAT{EXT_OUT_HDF5};
//////  
extern SW_time SW_TIME;
//////  
// external enumerated type
extern enum EQ SW_EQ;
extern enum Particle_IC DISTRIBUTION;
extern enum PT SW_PT;
extern enum SW_time SW_TIME;
// external structure type
extern AVS_parameters Avs_parameters;
// switch parameters
extern const char *EQ_name[];
extern const char *PT_name[];
extern int AC;
extern int BINARY;
extern int External_field;
extern int FIX_CELL;
extern int FIX_CELLxyz[DIM];
extern int Fixed_particle;
extern int G_direction;
extern int LJ_powers;
extern int LJ_truncate;
extern int N_spec;
extern int NS_source;
extern int PINNING;
extern int Poisson_Boltzmann;
extern int RESUMED;
extern int ROTATION;
extern int Shear_AC;
extern int STOKES;
extern int SW_AVS;
extern int SW_FFT;
extern int SW_UDF;
// Material parameters
extern double *MOI;
extern double *IMOI;
extern double *MASS;
extern double *MASS_RATIOS;
extern double IMASS_RATIO;
extern double *IMASS;
extern double *IMASS_RATIOS;
extern double RHO;
extern double IRHO;
extern double *RHO_particle;
extern double *S_surfaces; // pretilt scalar order
extern double *W_surfaces; // spring cst. of anchoring
extern double A;
extern double XI;
extern double IXI;
extern double HXI;
extern double A_XI;
extern double R_cutoff;
extern double A_R_cutoff;
extern double Axel;
extern double DT;
extern double DT_noise;
extern double DX;
extern double IDX;
extern double DX3;
extern double EPSILON;
extern double ETA;
extern double Delta_ETA;
extern double ikBT;
extern double kBT;
extern double NU;
extern double Nu_ratio;
extern double KMAX2;
extern double G;
extern double Inertia_stress;
extern double Ivolume;
extern double RADIUS;
extern double SIGMA;
extern double Tdump;
extern double VF;
extern double VF_LJ;
extern double T_LJ;
extern double LJ_dia;
extern double Srate_depend_LJ_cap;
extern double Min_rij;
extern double Max_force;
extern double Shear_frequency;
extern double Shear_rate;
extern double Shear_rate_eff;
extern double Shear_strain;
extern double Shear_strain_realized;
extern double dev_shear_stress[];
extern double rigid_dev_shear_stress[];
extern double &dev_shear_stress_lj;
extern double &dev_shear_stress_rot;
extern double &rigid_dev_shear_stress_lj;
extern double &rigid_dev_shear_stress_rot;
extern int Shear_strain_int;
// Step parameter(main)
extern int GTS;
extern int Num_snap;
extern int MSTEP;
extern int last_ts;
// Step parameter(sub)
extern int *Beads_Numbers;
extern int *Chain_Numbers;
extern int *Particle_Numbers;
extern int *Pinning_Numbers;
extern int *Pinning_ROT_Numbers;
extern int Component_Number;
extern int N_iteration_init_distribution;
extern int N_PIN;
extern int N_PIN_ROT;
extern int Particle_Number;
extern double Bjerrum_length;
extern double Surface_ion_number;
extern double Counterion_number;
// Mesh (Real)
extern double L[], LX, LY, LZ;
extern double HL[], HLX, HLY, HLZ;
extern double WAVE_X, WAVE_Y, WAVE_Z;
extern double L_particle[];
extern double iL_particle[];
extern double HL_particle[];
// Mesh (Real <Division>)
extern double LPs[SPACE][DIM];
extern double PREV_LPs[SPACE][DIM], NEXT_LPs[SPACE][DIM];
// Mesh (Mesh)
extern int Nmax;
extern int Nmin;
extern int Ns[], Ns_shear[], NX, NY, NZ, NZ_;
extern int HNs[], HNX, HNY, HNZ, HNZ_;
extern int N2s[], N2X, N2Y, N2Z, N2Z_;
extern int HN2s[], HN2X, HN2Y, HN2Z, HN2Z_;
extern int TRNs[], TRN_X, TRN_Y, TRN_Z;
extern int TRNs_QS[], TRN_QS_X, TRN_QS_Y, TRN_QS_Z;
extern int *KX_int, *KY_int, *KZ_int;
extern double *K2, *IK2;
// Mesh (Mesh <Division>)
extern int mesh_size;
extern int NPs[SPACE][DIM], HNQZ_;
extern int PREV_NPs[SPACE][DIM], NEXT_NPs[SPACE][DIM];
extern int *NPs_ALL;
/* Two_fluid */
extern double Mean_Bulk_concentration;
extern double Onsager_solute_coeff;
/* Electrolyte */
extern double **Concentration;
extern double **Concentration_rhs0;
extern double **Concentration_rhs1;
extern double **Surface_normal;
extern double *Onsager_coeff;
extern double *Surface_charge;
extern double *Surface_charge_e;
extern double *Total_solute;
extern double *Total_solute;
extern double *Valency;
extern double *Valency_e;
extern double Angular_Frequency;
extern double Debye_length;
extern double Dielectric_cst;
extern double E_ext[DIM];
extern double Elementary_charge;
extern double Frequency;
extern double Ion_charge;
extern double Onsager_coeff_counterion;
extern double Onsager_coeff_negative_ion;
extern double Onsager_coeff_positive_ion;
extern double Valency_counterion;
extern double Valency_negative_ion;
extern double Valency_positive_ion;
extern int **Sekibun_cell_exponential;
extern int NP_domain_exponential;
/* Temp Monitor */
extern double kT_snap_v; //fix R4A2
extern double kT_snap_o; //fix R4A2
extern double alpha_v;
extern double alpha_o;
/* AVS */
extern char Out_dir[];
extern char Out_name[];
/* UDF */
extern char *In_udf, *Sum_udf, *Out_udf, *Def_udf, *Ctrl_udf, *Res_udf;
/* FFT */
#ifdef _FFT_IMKL
extern DFTI_DESCRIPTOR_HANDLE imkl_p_fw, imkl_p_bw;
#elif _FFT_FFTW
extern fftw_plan fftw_p_fw, fftw_p_bw;
#elif _FFT_OOURA
typedef struct ooura_plan {
	int *ip;
	double *w;
	double *t;
	double ***a;
} ooura_plan;
extern ooura_plan ooura_p;
#endif
/* Sekibun cell */
extern int NP_domain;
extern int NP_domain_interface;
extern int **Sekibun_cell;
extern int **Sekibun_cell_interface;
/* Work */
extern double *Pressure, **Shear_force, **Shear_force_k, **f_particle;
extern double **ucp, *phi, *phi_sum, **up, **u, *rhop;
extern splineSystem** splineOblique;
extern double*** uspline;
extern Index_range* ijk_range_two_third_filter;
extern int n_ijk_range_two_third_filter;
extern double **work_v4, **work_v3, **work_v2, *work_v1;
extern double **f_ns0, **f_ns1, **f_ns2, **f_ns3, **f_ns4;
// OpenMP
extern int THREADNUM;
extern double *tmp_buffer1;
extern double *tmp_buffer2;
extern double **tmp_buffer_dim;
/* MPI */
extern int procs, xprocs, yprocs;
extern int procid, xid, yid;
extern int *xmesh, *ymesh;
extern int ierr;
extern int root;
extern int **tags;
extern int offset_x, offset_y;
extern int **idtbl;
extern int **lj_ref, lj_size;
extern int **sekibun_ref, sekibun_size;
extern int *sendranks, *sendcounts, *recvranks, *recvcounts;
extern int *sendbuf, *recvbuf;
extern Particle *p_tmp;
extern int *rcounts, *displs;
#ifdef _MPI
extern MPI_Datatype PT;
extern MPI_Request *ireq;
extern MPI_Status *ista;
extern MPI_Group GROUP_WORLD;
extern MPI_Group OWN_X_GROUP;
extern MPI_Comm OWN_X_COMM;
extern MPI_Group OWN_Y_GROUP;
extern MPI_Comm OWN_Y_COMM;
#endif
/* MPI Max Sekibun cell */
extern int Max_Sekibun_cells[DIM];
extern int Max_Sekibun_cell;
/* MPI Step parameter(sub) */
extern int Update_Particle_Number;
extern int Reference_Particle_Number;
extern int Local_Particle_Number;
/* MPI Debugging Error Check*/
#if defined(_MPI) && !defined(NDEBUG)
extern char errstr[MPI_MAX_ERROR_STRING + 1];
extern int errstrlen;
#endif
/* -none- */
extern long int seedval;

//////
extern OBL_INT SW_OBL_INT;
extern OUTFORMAT SW_OUTFORMAT;
extern EXTFORMAT SW_EXTFORMAT;
extern Field_crop      print_field_crop;
extern Field_output    print_field;

extern const char *OBL_INT_name[];
extern const char *OUTFORMAT_name[];
extern const char *EXTFORMAT_name[];
extern const char *JAX_name[];
extern const char *JP_name[];

//rigid variables
extern int Rigid_Number;
extern int **Rigid_Motions_vel;     // 0 (fix) or 1 (free)
extern int **Rigid_Motions_omega;   // 0 (fix) or 1 (free)
extern double **Rigid_Velocities;
extern double **Rigid_Omegas;
extern int *RigidID_Components;
extern int *Rigid_Particle_Numbers;
extern int *Rigid_Particle_Cumul;
extern double **xGs;
extern double **xGs_previous;
extern double **xGs_nopbc;
extern double *Rigid_Masses;
extern double *Rigid_IMasses;
extern double ***Rigid_Moments;
extern double ***Rigid_IMoments;
extern double ***Rigid_Moments_body;
extern double **velocityGs;
extern double **omegaGs;
extern double **forceGs;	//hydro
extern double **forceGrs;	//LJ
extern double **torqueGs;
extern double **torqueGrs;
extern double **velocityGs_old;
extern double **omegaGs_old;
extern double **forceGs_previous;
extern double **forceGrs_previous;
extern double **torqueGs_previous;
extern double **torqueGrs_previous;
extern int *Particle_RigidID;
//////

// PATCHY JANUS INTERACTIONS
extern int    SW_PATCHY;
extern int    PATCHY_POWER;
extern double PATCHY_EPSILON;
extern double PATCHY_LAMBDA;
extern double PATCHY_A_R_cutoff;
extern const double PATCHY_AXIS[];

extern int SW_JANUS;
extern int SW_JANUS_MOTOR;
extern int SW_JANUS_SLIP;
extern int MAX_SLIP_ITER;
extern double MAX_SLIP_TOL;
extern JAX *janus_axis;
extern JP *janus_propulsion;
extern double **janus_force;
extern double **janus_torque;
extern double *janus_slip_vel;
extern double *janus_slip_mode;

extern double **GRvecs;
extern double **GRvecs_body;

//////shear_degree
extern double degree_oblique;

extern double *Hydro_force;
extern double *Hydro_force_new;

extern Particle_IC DISTRIBUTION;
extern Particle_IO ORIENTATION;

extern Particle_output print_particle;

//////Hydro force for Lees-Edwards
extern double *h_obl_force;
extern double *h_obl_torque;
extern double *h_obl_force_local;
extern double *h_obl_torque_local;
extern double *h_obl_volume;
extern double *h_obl_Itrace;
extern double *h_obl_volume_local;
extern double *h_obl_Itrace_local;

typedef struct Particle_Info{
	int sign;
	double r[DIM];
	double dmy_ry;
	double dmyR;
	double dmy_phi;
	int im;
} Particle_Info;

extern Particle_Info **p_info;

//////Random numbers for particles
extern double* rand_num_particle;

#endif
