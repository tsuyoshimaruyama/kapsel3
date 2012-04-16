#include <assert.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "alloc.h"
#include "udfmanager.h"
#include "macro.h"
#include "quaternion.h"
#define MAXBUFFER 100
#define NDIM 3
using namespace std;
enum ID_OPTION {POSITION, ORIENTATION, VELOCITY, OMEGA, POSVEL,ANGULAR_POSVEL, ALL, INVALID};

void newframe(ofstream &outfile,
	      double &t,
	      int &ntot);

void write_xyz(ofstream &outfile, 
	       double &t,
	       double r[NDIM], 
	       double QR[NDIM][NDIM],
	       double v[NDIM],
	       double w[NDIM],
	       int &ntot,
	       ID_OPTION &flag
	       );

void init_xyz(ofstream &outfile, char *fname, ID_OPTION &flag);

void close_xyz(ofstream &outfile);


inline void wrong_invocation(){
  cout << "Usage: get_xyz -OPT UDFfile XYZfile" << endl;
  cout << "Options:" << endl;
  cout << "  -r\t positions" << endl;
  cout << "  -q\t orientations" << endl;
  cout << "  -v\t velocites" << endl;
  cout << "  -w\t angular velocities" << endl;
  cout << "  -rv\t positions + velocities" << endl;
  cout << "  -qw\t orientation + angular velocity" << endl;
  cout << "  -all\t all info" << endl;
  exit(1);
}

int main(int argc, char* argv[])
{
  
  //check invocation
  ID_OPTION id = INVALID;
  if(argc != 4){
    wrong_invocation();
  }
  if(strcmp(argv[1], "-r") == 0){
    id = POSITION;
    cout << "Print position" << endl;
  }
  if(strcmp(argv[1], "-q") == 0){
    id = ORIENTATION;
    cout << "Print orientation" << endl;
  }
  if(strcmp(argv[1], "-v") == 0){
    id = VELOCITY;
    cout << "Print velocity" << endl;
  }
  if(strcmp(argv[1], "-w") == 0){
    id = OMEGA;
    cout << "Print angular velocity" << endl;
  }
  if(strcmp(argv[1], "-rv") == 0){
    id = POSVEL;
    cout << "Print linear coordinates" << endl;
  }
  if(strcmp(argv[1], "-qw") == 0){
    id = ANGULAR_POSVEL;
    cout << "Print angular coordinates" << endl;
  }
  if(strcmp(argv[1], "-all") == 0){
    id = ALL;
    cout << "Print all coordinates" << endl;
  }
  if(id == INVALID) {
    wrong_invocation();
  }

  //make sure udf file exists
  UDFManager *ufin;
  if(file_check(argv[2])){
    //    cout << "Reading UDF file: "<< argv[2] << endl;
    ufin = new UDFManager(argv[2]);
  }
  ofstream outfile;
  //  cout << "Writing XYZ file: " << argv[3] << endl;
  //outfile.open(argv[3]);
  init_xyz(outfile, argv[3], id);
  

  int nx, ny, nz;
  int dt, frames, records;
  int ntotal;
  int nspec;
  int *pnum;
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
    int nmax = 0;
    ntotal = 0;
    for(int i = 0; i < nspec; i++){
      char str[256];
      sprintf(str, "object_type.spherical_particle.Particle_spec[%d]",i);
      Location target(str);
      ufin -> get(target.sub("Particle_number"), pnum[i]);
      nmax = MAX(nmax, pnum[i]);
      ntotal += pnum[i];
    }
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
    double w[NDIM];
    double QR[NDIM][NDIM];
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

	ufin -> get(target.sub("R.x"), r[0]);
	ufin -> get(target.sub("R.y"), r[1]);
	ufin -> get(target.sub("R.z"), r[2]);
	ufin -> get(target.sub("v.x"), v[0]);
	ufin -> get(target.sub("v.y"), v[1]);
	ufin -> get(target.sub("v.z"), v[2]);

	ufin -> get(target.sub("q.q0"), q0);
	ufin -> get(target.sub("q.q1"), q1);
	ufin -> get(target.sub("q.q2"), q2);
	ufin -> get(target.sub("q.q3"), q3);
	ufin -> get(target.sub("omega.x"), w[0]);
	ufin -> get(target.sub("omega.y"), w[1]);
	ufin -> get(target.sub("omega.z"), w[2]);
	qtn_init(q, q0 ,q1, q2, q3);
	qtn_isnormal(q, QTOL_LARGE);
	qtn_rm(QR, q);

	write_xyz(outfile, t, r, QR, v, w, ntotal, id);
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
  outfile << "# " << t << "  " << ntot << endl;
}

void write_xyz(ofstream &outfile, 
	       double &t,
	       double r[NDIM], 
	       double QR[NDIM][NDIM],
	       double v[NDIM],
	       double w[NDIM],
	       int &ntot,
	       ID_OPTION &flag){
  char str[4096];
  char dmy_str[256];
  

  sprintf(str, "%f ", t);
  if(flag == POSITION || flag == POSVEL || flag == ALL){
    sprintf(dmy_str, "%.6g  %.6g  %.6g ", r[0], r[1], r[2]);
    strcat(str, dmy_str);
  }
  if(flag == VELOCITY || flag == POSVEL || flag == ALL){
    sprintf(dmy_str, "%.6g  %.6g  %.6g ", v[0], v[1], v[2]);
    strcat(str, dmy_str);
    
  }
  if(flag == ORIENTATION || flag == ANGULAR_POSVEL || flag == ALL){
    sprintf(dmy_str, "%.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g ", 
	    QR[0][0], QR[1][0], QR[2][0],
	    QR[0][1], QR[1][1], QR[2][1],
	    QR[0][2], QR[1][2], QR[2][2]);
    strcat(str, dmy_str);
  }
  if(flag == OMEGA || flag == ANGULAR_POSVEL || flag == ALL){
    sprintf(dmy_str, "%.6g %.6g %.6g ", w[0], w[1], w[2]);
    strcat(str, dmy_str);
  }
  if(flag != POSITION && flag != ORIENTATION && flag != VELOCITY
     && flag != OMEGA && flag != POSVEL && flag != ANGULAR_POSVEL 
     && flag != ALL){
    cout << "Wrong flag!" << endl;
    exit(1);
  }
  outfile << str << endl;
}

void init_xyz(ofstream &outfile, char *fname, ID_OPTION &flag){

  outfile.open(fname);
  char str[256];    
  sprintf(str, "## ");
  if(flag == POSITION || flag == POSVEL || flag == ALL){
    strcat(str, "r_x, r_y, r_z ");
    
  }
  if(flag == VELOCITY || flag == POSVEL || flag == ALL){
    strcat(str, "v_x, v_y, v_z ");

  }
  if(flag == ORIENTATION || flag == ANGULAR_POSVEL || flag == ALL){
    strcat(str, "e1_x e1_y e1_z e2_x e2_y e2_z e3_x e3_y e3_z ");

  }
  if(flag == OMEGA || flag == ANGULAR_POSVEL || flag == ALL){
    strcat(str, "w_x, w_y, w_z");

  }
  if(flag != POSITION && flag != ORIENTATION && flag != VELOCITY
     && flag != OMEGA && flag != POSVEL && flag != ANGULAR_POSVEL 
     && flag != ALL){
    cout << "Wrong flag!" << endl;
    exit(1);
  }
  outfile << str << endl;

}

void close_xyz(ofstream &outfile){
  outfile.close();
}
