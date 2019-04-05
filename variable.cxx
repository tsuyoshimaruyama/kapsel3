#include "variable.h"
// external enumerated type
enum EQ SW_EQ;
enum Particle_IC DISTRIBUTION;
enum Particle_IO ORIENTATION;
enum PT SW_PT;
enum SW_time SW_TIME;
// external structure type
AVS_parameters Avs_parameters;
// switch parameters
const char *EQ_name[] = { "Navier_Stokes", "Shear_Navier_Stokes", "Shear_Navier_Stokes_Lees_Edwards", "Electrolyte" };
const char *PT_name[] = { "spherical_particle", "chain", "rigid" };
int AC;
int BINARY;
int External_field;
int FIX_CELL;
int FIX_CELLxyz[DIM];
int Fixed_particle = 0;
int G_direction;
int LJ_powers;
int LJ_truncate;
int N_spec;
int NS_source = 0;
int PINNING;
int Poisson_Boltzmann;
int RESUMED;
int ROTATION;
int Shear_AC;
int STOKES;
int SW_AVS;

int SW_FFT = IMKL_FFT;
int SW_UDF;
const int GIVEN_SEED = 0;
const int RANDOM_SEED = 1;
// Material parameters
double *MOI;
double *IMOI;
double *MASS;
double *MASS_RATIOS;
double IMASS_RATIO;
double *IMASS;
double *IMASS_RATIOS;
double RHO;
double IRHO;
double *RHO_particle;
double *S_surfaces; // pretilt scalar order
double *W_surfaces; // spring cst. of anchoring
double A;
double XI;
double IXI;
double HXI;
double A_XI;
double R_cutoff;
double A_R_cutoff;
double Axel;
double DT;
double DT_noise;
double DX;
double IDX;
double DX3;
double EPSILON;
double ETA;
double Delta_ETA;
double ikBT;
double kBT;
double NU;
double Nu_ratio;
double KMAX2;
double G;
double Inertia_stress;
double Ivolume;
double RADIUS;
double SIGMA;
double Tdump;
double VF;
double VF_LJ;
double T_LJ;
double LJ_dia;
double Srate_depend_LJ_cap;
double Min_rij;
double Max_force;
double Shear_frequency;
double Shear_rate;
double Shear_rate_eff;
double Shear_strain;
double Shear_strain_realized;
int Shear_strain_int;

//////
Field_crop      print_field_crop;
Field_output    print_field;

OBL_INT SW_OBL_INT;
OUTFORMAT SW_OUTFORMAT;
EXTFORMAT SW_EXTFORMAT;

//rigid variables
int Rigid_Number;
int **Rigid_Motions_vel;     // 0 (fix) or 1 (free)
int **Rigid_Motions_omega;   // 0 (fix) or 1 (free)
double **Rigid_Velocities;
double **Rigid_Omegas;
int *RigidID_Components;
int *Rigid_Particle_Numbers;
int *Rigid_Particle_Cumul;
double **xGs;
double **xGs_previous;
double **xGs_nopbc;
double *Rigid_Masses;
double *Rigid_IMasses;
double ***Rigid_Moments;
double ***Rigid_IMoments;
double ***Rigid_Moments_body;
double **velocityGs;
double **omegaGs;
double **forceGs;	//hydro
double **forceGrs;	//LJ
double **torqueGs;
double **torqueGrs;
double **velocityGs_old;
double **omegaGs_old;
double **forceGs_previous;
double **forceGrs_previous;
double **torqueGs_previous;
double **torqueGrs_previous;
int *Particle_RigidID;
//////

// PATCHY JANUS INTERACTIONS
int    SW_PATCHY;
int    PATCHY_POWER;
double PATCHY_EPSILON;
double PATCHY_LAMBDA;
double PATCHY_A_R_cutoff;

int SW_JANUS;
int SW_JANUS_MOTOR;
int SW_JANUS_SLIP;
int MAX_SLIP_ITER;
double MAX_SLIP_TOL;
JAX *janus_axis;
JP *janus_propulsion;
double **janus_force;
double **janus_torque;
double *janus_slip_vel;
double *janus_slip_mode;

double **GRvecs;
double **GRvecs_body;

//////shear_degree
double degree_oblique;


// Step parameter(main)
int GTS;
int Num_snap;
int MSTEP;
int last_ts;
// Step parameter(sub)
int *Beads_Numbers;
int *Chain_Numbers;
int *Particle_Numbers;
int *Pinning_Numbers;
int *Pinning_ROT_Numbers;
int Component_Number;
int N_iteration_init_distribution;
int N_PIN;
int N_PIN_ROT;
int Particle_Number;
double Bjerrum_length;
double Surface_ion_number;
double Counterion_number;
// Mesh (Real)
double L[DIM], LX, LY, LZ;
double HL[DIM], HLX, HLY, HLZ;
double WAVE_X, WAVE_Y, WAVE_Z;
double L_particle[DIM];
double iL_particle[DIM];
double HL_particle[DIM];
// Mesh (Real <Division>)
double LPs[SPACE][DIM];
double PREV_LPs[SPACE][DIM], NEXT_LPs[SPACE][DIM];
// Mesh (Mesh)
int Nmax;
int Nmin;
int Ns[DIM], Ns_shear[DIM], NX, NY, NZ, NZ_;
int HNs[DIM], HNX, HNY, HNZ, HNZ_;
int N2s[DIM], N2X, N2Y, N2Z, N2Z_;
int HN2s[DIM], HN2X, HN2Y, HN2Z, HN2X_, HN2Y_, HN2Z_;
int TRNs[DIM], TRN_X, TRN_Y, TRN_Z;
int TRNs_QS[DIM], TRN_QS_X, TRN_QS_Y, TRN_QS_Z;
int *KX_int, *KY_int, *KZ_int;
double *K2, *IK2;
// Mesh (Mesh <Division>)
int mesh_size;
int NPs[SPACE][DIM], HNQZ_;
int PREV_NPs[SPACE][DIM], NEXT_NPs[SPACE][DIM];
int *NPs_ALL;
/* Two_fluid */
double Mean_Bulk_concentration;
double Onsager_solute_coeff;
/* Electrolyte */
double **Concentration;
double **Concentration_rhs0;
double **Concentration_rhs1;
double **Surface_normal;
double *Onsager_coeff;
double *Surface_charge;
double *Surface_charge_e;
double *Total_solute;
double *Valency;
double *Valency_e;
double Angular_Frequency;
double Debye_length;
double Dielectric_cst;
double E_ext[DIM];
double Elementary_charge = 1.0;
double Frequency;
double Ion_charge;
double Onsager_coeff_counterion;
double Onsager_coeff_negative_ion;
double Onsager_coeff_positive_ion;
double Valency_counterion;
double Valency_negative_ion;
double Valency_positive_ion;
int **Sekibun_cell_exponential;
int NP_domain_exponential;
/* Temp Monitor */
double kT_snap_v = 1.439; //fix R4A2
double kT_snap_o = 1.22;  //fix R4A2
double alpha_v;
double alpha_o;
/* AVS */
char Out_dir[128];
char Out_name[128];
/* UDF */
char *In_udf, *Sum_udf, *Out_udf, *Def_udf, *Ctrl_udf, *Res_udf;
/* FFT */
#ifdef _FFT_IMKL
DFTI_DESCRIPTOR_HANDLE imkl_p_fw, imkl_p_bw;
#elif _FFT_FFTW
fftw_plan fftw_p_fw, fftw_p_bw;
#elif _FFT_OOURA
ooura_plan ooura_p;
#endif
/* Sekibun cell */
int NP_domain;
int NP_domain_interface;
int **Sekibun_cell;
int **Sekibun_cell_interface;
/* Work */
double *Pressure, **Shear_force, **Shear_force_k, **f_particle;
double **ucp, *phi, *phi_sum, **up, **u, *rhop;
splineSystem** splineOblique;
double*** uspline;
Index_range* ijk_range_two_third_filter;
int n_ijk_range_two_third_filter;
double **work_v4, **work_v3, **work_v2, *work_v1;
double **f_ns0, **f_ns1, **f_ns2, **f_ns3, **f_ns4;

double *Hydro_force;
double *Hydro_force_new;

// OpenMP
int THREADNUM;
double *tmp_buffer1;
double *tmp_buffer2;
double **tmp_buffer_dim;
// MPI
int procs, xprocs, yprocs;
int procid, xid, yid;
int *xmesh, *ymesh;
int ierr;
int root = 0;
int **tags;
int offset_x, offset_y;
int **idtbl;
int **lj_ref, lj_size;
int **sekibun_ref, sekibun_size;
int *sendranks, *sendcounts, *recvranks, *recvcounts;
int *sendbuf, *recvbuf;
Particle *p_tmp;
int *rcounts, *displs;
#ifdef _MPI
MPI_Datatype PT;
MPI_Request *ireq;
MPI_Status *ista;
MPI_Group GROUP_WORLD;
MPI_Group OWN_X_GROUP;
MPI_Comm OWN_X_COMM;
MPI_Group OWN_Y_GROUP;
MPI_Comm OWN_Y_COMM;
#endif
/* MPI Max Sekibun cell */
int Max_Sekibun_cells[DIM];
int Max_Sekibun_cell;
// MPI Step parameter(sub)
int Update_Particle_Number;
int Reference_Particle_Number;
int Local_Particle_Number;
/* MPI Debugging Error Check*/
#if defined(_MPI) && !defined(NDEBUG)
char errstr[MPI_MAX_ERROR_STRING + 1];
int errstrlen;
#endif

long int seedval;

const char *OBL_INT_name[]={"linear", "spline"};
const char *OUTFORMAT_name[]={"NONE", "ASCII", "BINARY", "EXTENDED"};
const char *EXTFORMAT_name[]={"HDF5"};
const char *JAX_name[]={"X", "Y", "Z", "NONE"};
const char *JP_name[]={"TUMBLER", "SQUIRMER", "OBSTACLE", "OFF"};

const double PATCHY_AXIS[DIM] = {0.0, 0.0, 1.0};

//not used
//int &TRN_QS_X=TRNs_QS[0];
//int &TRN_QS_Y=TRNs_QS[1];
//int &TRN_QS_Z=TRNs_QS[2];

double dev_shear_stress[2];
double rigid_dev_shear_stress[2];
double &dev_shear_stress_lj  = dev_shear_stress[0];
double &dev_shear_stress_rot = dev_shear_stress[1];
double &rigid_dev_shear_stress_lj  = rigid_dev_shear_stress[0];
double &rigid_dev_shear_stress_rot = rigid_dev_shear_stress[1];

UDFManager *ufin;
UDFManager *ufout;
UDFManager *ufres;

//////Hydro force for Lees-Edwards
double *h_obl_force;
double *h_obl_torque;
double *h_obl_force_local;
double *h_obl_torque_local;
double *h_obl_volume;
double *h_obl_Itrace;
double *h_obl_volume_local;
double *h_obl_Itrace_local;

//////Random numbers for particles
double* rand_num_particle;

Particle_Info **p_info;