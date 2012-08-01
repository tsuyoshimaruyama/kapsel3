#ifndef AVS_PROCESS_H
#define AVS_PROCESS_H

#include "avs_utils.h"
extern int Cs[DIM];
extern int &CX; //
extern int &CY; //
extern int &CZ; //
extern int HCs[DIM];
extern int &HCX; //
extern int &HCY; //
extern int &HCZ; //
extern bool error_calc;

void wrong_invocation(){
  fprintf(stderr, "usage: pavs -p PID -l Ns -i UDF [-rms]\n");
  fprintf(stderr, "       -p   PID\t ID (1...N) of centered particle\n");
  fprintf(stderr, "       -l   Ns \t Box size length\n");
  fprintf(stderr, "       -i   UDF\t Input file\n");
  fprintf(stderr, "       -err Perform error calculation [optional]\n");
  exit_job(EXIT_FAILURE);
}

void init_io(int argc, char *argv[], UDFManager* &udf_in, int &pid){
  
  char* fname;
  int pset, lset, fset, vset;
  pset = lset = fset = vset = 0;
  error_calc = false;
  for(int i = 1; i < argc; i++){
    if(strcmp(argv[i], "-p") == 0){
      if(i+1 == argc || argv[i+1][0] == '-'){
	fprintf(stderr, "Error: No particle selected\n");
	wrong_invocation();
      }
      pid = atoi(argv[++i]);
      if(--pid < 0){
	fprintf(stderr, "Incorrect id\n");
	wrong_invocation();
      }
      pset = 1;
    }else if(strcmp(argv[i],"-l") == 0){
      if(i+1 == argc || argv[i+1][0] == '-'){
	fprintf(stderr, "Error: No box size selected\n");
	wrong_invocation();
      }
      Cs[0] = atoi(argv[++i]);
      if(Cs[0] <= 0){
	fprintf(stderr, "Incorrect size\n");
	wrong_invocation();
      }      
      if(Cs[0]%2 == 0){
	Cs[0]++;
      }
      Cs[1] = Cs[2] = Cs[0];
      HCs[0] = HCs[1] = HCs[2] = Cs[0] / 2;
      lset = 1;
    }else if(strcmp(argv[i], "-i") == 0){
      if(i+1 == argc || argv[i+1][0] == '-'){
	fprintf(stderr, "Error: No udf file selected\n");
	wrong_invocation();
      }
      i++;
      fname = argv[i];
      if(file_check(argv[i])){
	udf_in = new UDFManager(argv[i]);
      }else{
	fprintf(stderr, "Error: Cannot open file\n");
      }
      fset = 1;
    }else if(strcmp(argv[i], "-err") == 0){
      i++;
      error_calc = true;
    }
  }

  if(!pset){
    fprintf(stderr, "Error: No particle selected\n");
    wrong_invocation();
  }
  if(!lset){
    fprintf(stderr, "Error: No box size selected\n");
    wrong_invocation();
  }
  if(!fset){
    fprintf(stderr, "Error: No udf file selected\n");
    wrong_invocation();
  }

  fprintf(stderr, " # Centered Particle: %8d\n", pid + 1);
  fprintf(stderr, " # Box Size         : (%d,%d,%d)\n", 
	  Cs[0], Cs[1], Cs[2]);
  fprintf(stderr, " # HBox Size        : (%d,%d,%d)\n", 
	  HCs[0], HCs[1], HCs[2]);
  fprintf(stderr, " # UDF File         : %s\n", fname);

}


#endif
