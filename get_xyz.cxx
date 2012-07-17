#include <assert.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "alloc.h"
#include "udfmanager.h"
#include "macro.h"
#include "quaternion.h"
#include "rigid_body.h"
#include "lad3.h"
#define MAXBUFFER 100
#define NDIM 3
using namespace std;
enum ID_OPTION_P {POSITION, ORIENTATION, ALL_POS, NO_POS};
enum ID_OPTION_V {VELOCITY, OMEGA, ALL_VEL, NO_VEL};
enum ID_OPTION_F {FORCE, TORQUE, ALL_FRC, NO_FRC};
enum ID_OPTION_N {JANUS_X, JANUS_Y, JANUS_Z, NO_JANUS};
enum OPTION {in_pos = 1, in_vel = 2, in_frc = 3, in_nt=4, in_fin = 5, in_fout = 6};
double ex[NDIM] = {1.0, 0.0, 0.0};
double ey[NDIM] = {0.0, 1.0, 0.0};
double ez[NDIM] = {0.0, 0.0, 1.0};

void newframe(ofstream &outfile,
	      double &t,
	      int &ntot);

void write_xyz(ofstream &outfile, 
	       double &t,
	       const int &pid,
	       const int &sid,
	       double r[NDIM], 
	       double QR[NDIM][NDIM],
	       double v[NDIM],
	       double vs[NDIM],
	       double w[NDIM],
	       double frc[NDIM],
	       double frc_slip[NDIM],
	       double tau[NDIM],
	       double tau_slip[NDIM],
	       double ni[NDIM],
	       double n0[NDIM],
	       int &ntot,
	       ID_OPTION_P &pflag,
	       ID_OPTION_V &vflag,
	       ID_OPTION_F &fflag,
	       ID_OPTION_N &nflag
	       );

void init_xyz(ofstream &outfile, char *fname, ID_OPTION_P &pflag, ID_OPTION_V &vflag, ID_OPTION_F &fflag,
	      ID_OPTION_N &nflag);

void close_xyz(ofstream &outfile);


inline void wrong_invocation(){
  cout << "Usage: get_xyz -OPT_pos -OPT_vel -OPT_frc -OPT_janus UDFfile XYZfile" << endl;
  cout << "Options:" << endl;
  cout << "OPT_pos:" << endl;
  cout << "  -r\t positions" << endl;
  cout << "  -q\t orientations" << endl;
  cout << "  -rq\t all pos/orientation info" << endl;
  cout << endl;
  cout << "OPT_vel:" << endl;
  cout << "  -v\t linear velocities" << endl;
  cout << "  -w\t angular velocities" << endl;
  cout << "  -vw\t all velocity infro" << endl;
  cout << endl;
  cout << "OPT_frc:" << endl;
  cout << "  -f\t hydrodynamic forces" << endl;
  cout << "  -t\t hydrodynamic torques" << endl;
  cout << "  -ft\t all force/torque info" << endl;
  cout << "OPT_janus:" << endl;
  cout << "  -jx\t janus x axis" << endl;
  cout << "  -jy\t janus y axis" << endl;
  cout << "  -jz\t janus z axis" << endl;
  exit(1);
}

int main(int argc, char* argv[])
{
  
  //check invocation
  ID_OPTION_P idp = NO_POS;
  ID_OPTION_V idv = NO_VEL;
  ID_OPTION_F idf = NO_FRC;
  ID_OPTION_N idn = NO_JANUS;
  double * janus_axis = NULL;
  
  if(argc != 7){
    wrong_invocation();
  }
  if(strcmp(argv[in_pos], "-r") == 0){
    idp = POSITION;
    cout << "Print position" << endl;
  }else if(strcmp(argv[in_pos], "-q") == 0){
    idp = ORIENTATION;
    cout << "Print orientation" << endl;
  }else if(strcmp(argv[in_pos], "-rq") == 0){
    idp = ALL_POS;
    cout << "Print all position info" << endl;
  }else if(strcmp(argv[in_pos], "-n")){
    idp = NO_POS;
    cout << "Skipping position info" << endl;
  }else{
    wrong_invocation();
  }

  if(strcmp(argv[in_vel], "-v") == 0){
    idv = VELOCITY;
    cout << "Print velocity" << endl;
  }else if(strcmp(argv[in_vel], "-w") == 0){
    idv = OMEGA;
    cout << "Print angular velocity" << endl;
  }else if(strcmp(argv[in_vel], "-vw") == 0){
    idv = ALL_VEL;
    cout << "Print all velocity info" << endl;
  }else if(strcmp(argv[in_vel], "-n") == 0){
    idv = NO_VEL;
    cout << "Skipping velocity info" << endl;
  }else{
    wrong_invocation();
  }

  if(strcmp(argv[in_frc], "-f") == 0){
    idf = FORCE;
    cout << "Print hydrodynamic force" << endl;
  }else if(strcmp(argv[in_frc], "-t") == 0){
    idf = TORQUE;
    cout << "Print hydrodynamic torque" << endl;
  }else if(strcmp(argv[in_frc], "-ft") == 0){
    idf = ALL_FRC;
    cout << "Print hydrodynamic force info" << endl;
  }else if(strcmp(argv[in_frc], "-n") == 0){
    idf = NO_FRC;
    cout << "Skipping force info" << endl;
  }else{
    wrong_invocation();
  }

  if(strcmp(argv[in_nt], "-jx") == 0){
    idn = JANUS_X;
    janus_axis = ex;
    cout << "Print Janus x axis" << endl;
  }else if(strcmp(argv[in_nt], "-jy") == 0){
    idn = JANUS_Y;
    janus_axis = ey;
    cout << "Print Janus y axis" << endl;
  }else if(strcmp(argv[in_nt], "-jz") == 0){
    idn = JANUS_Z;
    janus_axis = ez;
    cout << "Print Janus z axis" << endl;
  }else if(strcmp(argv[in_nt], "-n") == 0){
    idn = NO_JANUS;
    cout << "Skipping Janus info" << endl;
  }else{
    wrong_invocation();
  }

  if(idp == NO_POS && idv == NO_VEL && idf == NO_FRC && idn == NO_JANUS){
    cout << "Nothing to do" << endl;
    exit(0);
  }

  //make sure udf file exists
  UDFManager *ufin;
  if(file_check(argv[in_fin])){
    ufin = new UDFManager(argv[in_fin]);
  }
  ofstream outfile;
  init_xyz(outfile, argv[in_fout], idp, idv, idf, idn);
  

  int nx, ny, nz;
  int dt, frames, records;
  int ntotal;
  int nspec;
  int *pnum;
  int *spec_id;
  double *bslip;
  const string sep = "\t\t";  

  {
    //mesh size
    Location target("mesh");
    int npx, npy, npz;
    ufin -> get(target.sub("NPX"), npx);
    ufin -> get(target.sub("NPY"), npy);
    ufin -> get(target.sub("NPZ"), npz);
    nx = 1 << npx;
    ny = 1 << npy;
    nz = 1 << npz;
  }
  {
    //particle species data
    Location target("object_type");
    string str;
    ufin -> get(target.sub("type"), str);
    if(str != "spherical_particle"){
      cout << "Error unknown particle type" << endl;
      exit(1);
    }
    nspec = ufin -> size("object_type.spherical_particle.Particle_spec[]");

    pnum = alloc_1d_int(nspec);
    bslip = alloc_1d_double(nspec);
    int nmax = 0;
    ntotal = 0;
    for(int i = 0; i < nspec; i++){
      char str[256];
      sprintf(str, "object_type.spherical_particle.Particle_spec[%d]",i);
      Location target(str);
      ufin -> get(target.sub("Particle_number"), pnum[i]);
      ufin -> get(target.sub("janus_slip_vel"), bslip[i]);
      fprintf(stderr, "SLIIIIP: %10.8g ", bslip[i]);
      bslip[i] = 1.0 / (2./3 * bslip[i]);
      fprintf(stderr, " (%10.8g) \n", bslip[i]);
      nmax = MAX(nmax, pnum[i]);
      ntotal += pnum[i];
    }
    spec_id = (int*) malloc(sizeof(int) * ntotal);
    int sum = 0;
    for(int i = 0; i < nspec; i++){
      for(int j = 0; j < pnum[i]; j++){
	spec_id[sum] = i;
	sum++;
      }
    }
    assert(sum == ntotal);
  }
  {
    //frame data
    Location target("output");
    ufin -> get(target.sub("GTS"), dt);
    ufin -> get(target.sub("Num_snap"), frames);
    records = ufin -> totalRecord();
  }
  {
    //read particle data
    double r[NDIM];
    double v[NDIM];
    double vs[NDIM];
    double w[NDIM];
    double QR[NDIM][NDIM];
    double frc[NDIM];
    double frc_slip[NDIM];
    double tau[NDIM];
    double tau_slip[NDIM];
    double n0[NDIM];
    double ni[NDIM];
    quaternion q;
    double t;
    double q0,q1,q2,q3;
    for(int i = 0; i < records; i++){
      char str[256];
      sprintf(str, "#%d", i);
      ufin -> jump(str);
      ufin -> get("t", t);

      newframe(outfile, t, ntotal);
      for(int j = 0; j < ntotal; j++){
	char str[256];
	sprintf(str, "Particles[%d]", j);
	Location target(str);

	if(idp == POSITION || idp == ALL_POS){
	  ufin -> get(target.sub("R.x"), r[0]);
	  ufin -> get(target.sub("R.y"), r[1]);
	  ufin -> get(target.sub("R.z"), r[2]);
	}
	
	if(idp == ORIENTATION || idp == ALL_POS){
	  ufin -> get(target.sub("q.q0"), q0);
	  ufin -> get(target.sub("q.q1"), q1);
	  ufin -> get(target.sub("q.q2"), q2);
	  ufin -> get(target.sub("q.q3"), q3);
	  qtn_init(q, q0 ,q1, q2, q3);
	  qtn_isnormal(q, QTOL_LARGE);
	  rqtn_rm(QR, q);
	}

	if(idv == VELOCITY || idv == ALL_VEL){
	  ufin -> get(target.sub("v.x"), v[0]);
	  ufin -> get(target.sub("v.y"), v[1]);
	  ufin -> get(target.sub("v.z"), v[2]);
	  for(int d = 0; d < NDIM; d++){
	    vs[d] = v[d] * bslip[spec_id[j]];
	  }
	}

	if(idv == OMEGA || idv == ALL_VEL){
	  ufin -> get(target.sub("omega.x"), w[0]);
	  ufin -> get(target.sub("omega.y"), w[1]);
	  ufin -> get(target.sub("omega.z"), w[2]);
	}

	if(idf == FORCE || idf == ALL_FRC){
	  ufin -> get(target.sub("f_hydro.x"), frc[0]);
	  ufin -> get(target.sub("f_hydro.y"), frc[1]);
	  ufin -> get(target.sub("f_hydro.z"), frc[2]);

	  ufin -> get(target.sub("f_slip.x"), frc_slip[0]);
	  ufin -> get(target.sub("f_slip.y"), frc_slip[1]);
	  ufin -> get(target.sub("f_slip.z"), frc_slip[2]);
	}

	if(idf == TORQUE || idf == ALL_FRC){
	  ufin -> get(target.sub("torque_hydro.x"), tau[0]);
	  ufin -> get(target.sub("torque_hydro.y"), tau[1]);
	  ufin -> get(target.sub("torque_hydro.z"), tau[2]);

	  ufin -> get(target.sub("torque_slip.x"), tau_slip[0]);
	  ufin -> get(target.sub("torque_slip.y"), tau_slip[1]);
	  ufin -> get(target.sub("torque_slip.z"), tau_slip[2]);
	}

	if(idn != NO_JANUS){
	  ufin -> get(target.sub("q.q0"), q0);
	  ufin -> get(target.sub("q.q1"), q1);
	  ufin -> get(target.sub("q.q2"), q2);
	  ufin -> get(target.sub("q.q3"), q3);
	  qtn_init(q, q0, q1, q2, q3);
	  qtn_isnormal(q, QTOL_LARGE);
	  qtn_normalize(q);
	  rigid_body_rotation(ni, janus_axis, q, BODY2SPACE);
	  if(i == 0){
	    v_copy(n0, ni);
	  }
	}
	write_xyz(outfile, t, j+1, spec_id[j], r, QR, v, vs, w, frc, frc_slip, tau, tau_slip, ni, n0, ntotal, idp, idv, idf, idn);
      }
    }
  }

  free_1d_int(pnum);
  close_xyz(outfile);
  return 0;
}

void newframe(ofstream &outfile,
	      double &t,
	      int &ntot){
  //  outfile << endl;
}

void write_xyz(ofstream &outfile, 
	       double &t,
	       const int &pid,
	       const int &sid,
	       double r[NDIM], 
	       double QR[NDIM][NDIM],
	       double v[NDIM],
	       double vs[NDIM],
	       double w[NDIM],
	       double frc[NDIM],
	       double frc_slip[NDIM],
	       double tau[NDIM],
	       double tau_slip[NDIM],
	       double ni[NDIM],
	       double n0[NDIM],
	       int &ntot,
	       ID_OPTION_P &pflag,
	       ID_OPTION_V &vflag,
	       ID_OPTION_F &fflag,
	       ID_OPTION_N &nflag){
  char str[4096];
  char dmy_str[256];
  double dmy_ndot = ni[0]*n0[0] + ni[1]*n0[1] + ni[2]*n0[2];
  double dmy_theta = (dmy_ndot < 1.0 ? acos(dmy_ndot) : 0.0);

  sprintf(str, "%d %d ", pid, sid);
  if(pflag == POSITION || pflag == ALL_POS){
    sprintf(dmy_str, "%.6g  %.6g  %.6g  ", r[0], r[1], r[2]);
    strcat(str, dmy_str);
  }
  if(vflag == VELOCITY || vflag == ALL_VEL){
    sprintf(dmy_str, "%.6g  %.6g  %.6g  %.6g  %.6g  ", v[0], v[1], v[2],
	    sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]),
	    sqrt(vs[0]*vs[0] + vs[1]*vs[1] + vs[2]*vs[2])
	    );
    strcat(str, dmy_str);
  }

  if(pflag == ORIENTATION || pflag == ALL_POS){
    sprintf(dmy_str, "%.6g  %.6g  %.6g  %.6g  %.6g  %.6g  %.6g  %.6g  %.6g  ", 
	    QR[0][0], QR[1][0], QR[2][0],
	    QR[0][1], QR[1][1], QR[2][1],
	    QR[0][2], QR[1][2], QR[2][2]);
    strcat(str, dmy_str);
  }
  if(vflag == OMEGA || vflag == ALL_VEL){
    sprintf(dmy_str, "%.6g  %.6g  %.6g  ", w[0], w[1], w[2]);
    strcat(str, dmy_str);
  }

  if(fflag == FORCE || fflag == ALL_FRC){
    sprintf(dmy_str, "%.6g  %.6g  %.6g  %.6g  %.6g  %.6g  ", frc[0], frc[1], frc[2], frc_slip[0], frc_slip[1], frc_slip[2]);
    strcat(str, dmy_str);
  }
  if(fflag == TORQUE || fflag == ALL_FRC){
    sprintf(dmy_str, "%.6g  %.6g  %.6g  %.6g  %.6g  %.6g  ", tau[0], tau[1], tau[2], tau_slip[0], tau_slip[1], tau_slip[2]);
    strcat(str, dmy_str);
  }
  if(nflag != NO_JANUS){
    sprintf(dmy_str, "%.6g %.6g %.6g %.15g %.6g", ni[0], ni[1], ni[2], dmy_ndot, dmy_theta * 180.0 / M_PI);
    strcat(str, dmy_str);
  }

  if(pflag != POSITION && pflag != ORIENTATION && pflag != ALL_POS && pflag != NO_POS){
    cout << "Wrong position flag!" << endl;
    exit(1);
  }
  if(vflag != VELOCITY && vflag != OMEGA && vflag != ALL_VEL && vflag != NO_VEL){
    cout << "Wrong velocity flag!" << endl;
    exit(1);
  }
  if(fflag != FORCE && fflag != TORQUE && fflag != ALL_FRC && fflag != NO_FRC){
    cout << "Wrong  force flag!" << endl;
    exit(1);
  }
  if(nflag != JANUS_X && nflag != JANUS_Y && nflag != JANUS_Z && nflag != NO_JANUS){
    cout << "Wrong janus flag!" << endl;
    exit(1);
  }

  outfile << str << endl;
}

void init_xyz(ofstream &outfile, char *fname, ID_OPTION_P &pflag, ID_OPTION_V &vflag, ID_OPTION_F &fflag,
	      ID_OPTION_N &nflag){

  outfile.open(fname);
  char str[256];    
  sprintf(str, "## 1:p_id 2:spec_id ");
  if(pflag == POSITION || pflag == ALL_POS){
    strcat(str, "3:r_x 4:r_y 5:r_z ");
  }
  if(vflag == VELOCITY || vflag == ALL_VEL){
    strcat(str, "6:v_x 7:v_y 8:v_z 9:v 10:vs ");

  }
  if(pflag == ORIENTATION || pflag == ALL_POS){
    strcat(str, "11:e1_x 12:e1_y 13:e1_z 14:e2_x 15:e2_y 16:e2_z 17:e3_x 18:e3_y 19:e3_z ");
  }
  if(vflag == OMEGA || vflag == ALL_VEL){
    strcat(str, "20:w_x 21:w_y 22:w_z ");
  }
  if(fflag == FORCE || fflag == ALL_FRC){
    strcat(str, "23:F_x 24:F_y 25:F_z 26:Fs_x 27:Fs_y 28:Fs_z ");
  }
  if(fflag == TORQUE || fflag == ALL_FRC){
    strcat(str, "29:N_x 30:N_y 31:N_z 32:Ns_x 33:Ns_y 34:Ns_z ");
  }
  if(nflag != NO_JANUS){
    strcat(str, "35:e_x 36:e_y 37:e_z 38:c_ee 39:theta");
  }
  outfile << str << endl;

}

void close_xyz(ofstream &outfile){
  outfile.close();
}
