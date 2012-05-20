#include "ewald_multipole.h"

static double e1 = {1, 0, 0};
static double e2 = {0, 1, 0};
static double e3 = {0, 0, 1};

int N; 
double h[DIM][DIM], h_i[DIM][DIM];

int kc, kc2, kc_i, kc2_i;
double rc, rc2, rc_i, rc2_i;

double V, Vi;
double alpha, alpha2, alpha_i, alpha2_i;

int nk_domain;
int **Ewald_k_cell;
double *Ewald_exp_cell;

double **f_k, **f_r;                // forces
double **t_k, **t_r, t_surface;        // torques
double e_k, e_r, e_self, e_surface;    // energies

inline void MI_distance(const double *r1, const double *r2,
			double *r12, double &dr12){
  double s1[DIM], s2[DIM], s12[DIM];

  M_v_prod(s1, h_i, r1);
  M_v_prod(s2, h_i, r2);
  for(int d = 0; d < DIM; d++){
    s12[d] = s1[d] - s2[d];
    s12[d] = s12[d] - (double)Nint(s12[d]);
  }
  M_v_prod(r12, h, s12);
  dr12 = v_norm(r12);
}

// 2*rc < smallest perpendicular width of the unit cell
inline void valid_rcut(const double &rcut, 
		       const double *a, const double *b, const double *c){
  double wa, wb, wc;
  double bc[DIM], ca[DIM], ab[DIM];
  v_cross(bc, b, c);
  v_cross(ca, c, a);
  v_cross(ab, a, b);
  
  v_inner_prod(wa, bc);
  v_inner_prod(wb, ca);
  v_inner_prod(wc, ab);
  
  wa = ABS(wa) / v_norm(bc);
  wb = ABS(wb) / v_norm(ca);
  wc = ABS(wc) / v_norm(ab);
  if(rcut > 0.5 * MIN(wa, MIN(wb, wc))){
    fprintf(stderr, "Error: decrese rcut for ewald summation %f (%f)\n"
	    rcut, 0.5*MIN(wa, MIN(wb, wc)));
    exit_job(EXIT_FAILURE);
  }
}

inline void n2k(const int *nv, double *kv){
  double dnv[DIM];
  dnv[0] = (double) nv[0];
  dnv[1] = (double) nv[1];
  dnv[2] = (double) nv[2];
  M_v_prod(kv, h_i, dnv, PI2);
}
// Create Ewald Sekibun cells
void Ewald_sekibun_cell(int &nk_domain, 
			int **k_cell, 
			double *exp_cell, 
			int &nmax){

  int n_mesh = 8*((nmax + 2)*(nmax+2)*(nmax+2));
  int **dmy_k = alloc_2d_int(n_mesh, DIM);
  double **dmy_exp = alloc_1d_double(n_mesh);

  int mesh = 0;
  double kv[DIM], nv[DIM];

  // x > 0 
  for(int i = 1; i <= nmax, i++){

    nv[0] = i;
    for(int j = -nmax; j <= nmax, j++){

      nv[1] = j;
      for(int k = -nmax; k <= nmax, k++){

	nv[2] = k;
	n2k(nv, kv);
      }
    }
  }

  // x = 0 && y > 0
  {
    int i = 0;
    nv[0] = i;
    for(int j = 1; j <= nmax, j++){

      nv[1] = j;
      for(int k = -nmax; k <= nmax; k++){

	nv[2] = k;
	n2k(nv, kv);
      }
    }
  }

  // x = 0 && y == 0
  {
    int i = 0, j = 0;
    nv[0] = i;
    nv[1] = j;
    for(int k = 0; k <= nmax; k++){

      nv[2] = k;
      n2k(nv, kv);
    }
  }

}

// Initialize ewald parameters 
// Should only be called once !
void init_ewald(double *a, double *b, double *c,
		double &rcut, int &kcut, double &ewald_alpha, 
		int &Nparticle){

  // box matrix
  h[0][0] = a[0];
  h[1][0] = a[1];
  h[2][0] = a[2];

  h[0][1] = b[0];
  h[1][1] = b[1];
  h[2][1] = b[2];

  h[0][2] = c[0];
  h[1][2] = c[1];
  h[2][2] = c[2];

  V = M_det(h);
  Vi = 1.0/V;
  M_inv(h_i, h);
  
  // real space cutoff
  valid_rcut(rcut, a, b, c);
  rc = rcut;
  rc2 = rc*rc;
  rc_i = 1.0/rcut;
  rc2_i = 1.0/rc2;

  // reciprocal space cutoff
  kc = kcut;
  kc2 = kc*kc;
  kc_i = 1.0/kc;
  kc2_i = 1.0/kc2;

  // ewald damping parameter
  alpha = ewald_alpha;
  alpha2 = alpha*alpha;
  alpha_i = 1.0/alpha;
  alpha2_i = 1.0/alpha2;

  // particle forces & torques
  N = Nparticle;
  f_k = alloc_2d_double(N, DIM);
  f_r = alloc_2d_double(N, DIM);
  t_k = alloc_2d_double(N, DIM);
  t_r = alloc_2d_double(N, DIM);
  t_surface = alloc_2d_double(N, DIM);
}

// Reset ewald energy, force, torque
inline void reset_ewald(){
  e_k = 0.0;
  e_r = 0.0;
  e_self = 0.0;
  e_surface = 0.0;

  for(int n = 0; n < N; n++){
    for(int d = 0; d < DIM; d++){
      f_k[n][d] = 0.0;
      f_r[n][d] = 0.0;

      t_k[n][d] = 0.0;
      t_r[n][d] = 0.0;
      t_surface[n][d] = 0.0;
    }
  }
}

void ewald_multipole_real(const Particle *p){
}

void ewald_multipole_recip(const Particle *p){
  for(int n = 0; )
}

void ewald_multipole_surface(const Particle *p){
}

void ewald_multipole_self(const Particle *p){
}

void ewald_multipole_interaction(double *E_ewald,
				 Particle *p){
  reset_ewald();

  ewald_multipole_real(p);

  ewald_multipole_recip(p);

  ewald_multipole_surface(p);
  
  ewald_multipole_self(p);

  E_ewald[0] = e_r + e_k + e_self + e_surface;
  E_ewald[1] = e_r;
  E_ewald[2] = e_k;
  E_ewald[3] = e_self;
  E_ewald[4] = e_surface;
  for(int n = 0; n < N; n++){
    for(int d = 0; d < DIM; d+=){
      p[n].ewald_force[d] = f_r[n][d] + f_k[n][d];
      p[n].ewald_torque[d] = t_r[n][d] + t_k[n][d] + t_surface[n][d];
    }
  }
}
