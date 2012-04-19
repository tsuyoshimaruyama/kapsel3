#include <iostream>
#include <fstream>
#include <string>
#include "rigid_body.h"

using namespace std;
int main(int argc, char**argv){
  ofstream o_qr, o_rr;
  o_qr.open("quaternion_rotation.dat");
  o_rr.open("matrix_rotation.dat");
  printf("Rigid Body Test\n");

  //Initial conditions: anular momentum and orientation
  double p0[DIM] = {0.0, 0.0, 4.0}; //angular momentum (space frame)
  double n_theta = M_PI/3.0;        //rotation vector for initial orientation
  double n_phi = 0.0;
  double n[DIM] = {sin(n_theta)*cos(n_phi), 
		   sin(n_theta)*sin(n_phi),
		   cos(n_theta)};
  double phi = M_PI/6.0;            //rotation angle
  quaternion q0;
  double R0[DIM][DIM] = {{1.0, 0.0, 0.0},
			 {0.0, 1.0, 0.0},
			 {0.0, 0.0, 1.0}};

  rv_qtn(q0, phi, n); //orientation - quaternion repr.
  qtn_rm(R0, q0);     //orientation - so(3)

  double Ip[DIM] = {2.0, 2.0, 8.0}; // principles intertia moments (body axes)
  double I[DIM][DIM] = {{Ip[0], 0.0, 0.0}, // inertia tensor (body axes)
			{0.0, Ip[1], 0.0},
			{0.0, 0.0, Ip[2]}};
  double II[DIM][DIM] = {{1.0/Ip[0], 0.0, 0.0},
			 {0.0, 1.0/Ip[1], 0.0},
			 {0.0, 0.0, 1.0/Ip[2]}};
  double w0[DIM];
  double O, T, dt;
  rigid_body_rotation(w0, p0, R0, SPACE2BODY);
  M_v_prod(w0, II);      //initial angular velocity (body axes)
  O = (Ip[2] - Ip[0])/Ip[0]*w0[2]; //angular velocity of precesion
  T = (O != 0.0) ? (2*M_PI / O) : (2*M_PI / v_norm(w0));
  dt = T/1000;

  printf("#w0 = (%.5f, %.5f, %.5f)\n", 
	  w0[0], w0[1], w0[2]);
  printf("#q0 = (%.5f, %.5f, %.5f, %.5f)\n",
	  qtn_q0(q0), qtn_q1(q0), qtn_q2(q0), qtn_q3(q0));
  printf("#|q|= %.5f\n", qtn_norm(q0));
  printf("#Omega = %.5f\n", O);
  printf("#dt = %.5f\n", dt);

  quaternion q;
  double w[DIM], w2[DIM];
  double Rt[DIM][DIM], Rt2[DIM][DIM];
  v_copy(w, w0);
  v_copy(w2, w0);
  qtn_init(q, q0);
  qtn_rm(Rt, q);
  M_copy(Rt2, Rt);
  o_qr.precision(8);
  o_rr.precision(8);
  string spc = "\t ";
  o_qr << Rt[0][0] << spc << Rt[1][0] << spc << Rt[2][0] << spc
       << Rt[0][1] << spc << Rt[1][1] << spc << Rt[2][1] << spc
       << Rt[0][2] << spc << Rt[1][2] << spc << Rt[2][2] << endl;
  o_rr << Rt2[0][0] << spc << Rt2[1][0] << spc << Rt2[2][0] << spc
       << Rt2[0][1] << spc << Rt2[1][1] << spc << Rt2[2][1] << spc
       << Rt2[0][2] << spc << Rt2[1][2] << spc << Rt2[2][2] << endl;
  for(int steps = 0; steps < 10000; steps++){
    propagate_wq_RK4(q, w, I, dt);
    propagate_w_RK4_Q_Euler(Rt2, w2, I, dt);
    M_isValidRotation(Rt2);
    o_qr << Rt[0][0] << spc << Rt[1][0] << spc << Rt[2][0] << spc
	 << Rt[0][1] << spc << Rt[1][1] << spc << Rt[2][1] << spc
	 << Rt[0][2] << spc << Rt[1][2] << spc << Rt[2][2] << endl;
    o_rr << Rt2[0][0] << spc << Rt2[1][0] << spc << Rt2[2][0] << spc
	 << Rt2[0][1] << spc << Rt2[1][1] << spc << Rt2[2][1] << spc
	 << Rt2[0][2] << spc << Rt2[1][2] << spc << Rt2[2][2] << endl;

    printf("%d\n", steps);
    qtn_rm(Rt, q);

  }
  o_qr.close();
  o_rr.close();
  return 0;
}
