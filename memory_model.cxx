#include "memory_model.h"
/*領域分割パラメタ設定*/
void Range_division (void) {

    int *npx, *npy;
    int *nqy, *nqz_;
    int nps_max[SPACE][DIM];
    int dmy[DIM];
    int dmy_mod;
    int id = 0;
#ifndef _MPI
    procs = 1;
    xprocs = 1;
    yprocs = 1;
    procid = 0;
    xid = 0;
    yid = 0;
#endif
#ifdef _OPENMP
    THREADNUM = omp_get_max_threads();
#else
    THREADNUM = 1;
#endif
    NPs_ALL = alloc_1d_int(procs * SPACE * DIM);
    npx = alloc_1d_int(xprocs);
    npy = alloc_1d_int(xprocs);
    for(int i = 0; i < xprocs; i++) {
        npx[i] = (int) (NX / xprocs);
        npy[i] = (int) (NY / xprocs);
    }
    dmy_mod = NX % xprocs;
    if(dmy_mod) {
/*
        for(int i = 0; i < dmy_mod; i++) {
            npx[i]++;
        }
*/
        npx[xprocs - 1] += dmy_mod;
    }
    dmy_mod = NY % xprocs;
    if(dmy_mod) {
/*
        for(int i = 0; i < dmy_mod; i++) {
            npy[i]++;
        }
*/
        npy[xprocs - 1] += dmy_mod;
    }
    nqy = alloc_1d_int(yprocs);
    nqz_ = alloc_1d_int(yprocs);
    for(int i = 0; i < yprocs; i++) {
        nqy[i] = (int)(NY / yprocs);
        nqz_[i] = (int)(HNZ_ / yprocs);
    }
    dmy_mod = NY % yprocs;
    if(dmy_mod) {
/*
        for(int i = 0; i < dmy_mod; i++) {
            nqy[i]++;
        }
*/
        nqy[yprocs - 1] += dmy_mod;
    }
    dmy_mod = HNZ_ % yprocs;
    if(dmy_mod) {
/*
        for(int i = 0; i < dmy_mod; i++) {
            nqz_[i]++;
        }
*/
        nqz_[yprocs - 1] += dmy_mod;
    }
    HNQZ_ = nqz_[yid];
    for(int i = 0; i < yprocs; i++) {
        nqz_[i] *= 2;
    }
// Mesh
    NPs[REAL][0] = npx[xid];
    NPs[REAL][1] = nqy[yid];
    NPs[REAL][2] = NZ;
    NPs[SPECTRUM][0] = NX;
    NPs[SPECTRUM][1] = npy[xid];
    NPs[SPECTRUM][2] = nqz_[yid];
    for (int i = 0; i < xprocs; i++){
        for (int j = 0; j < yprocs; j++){
            NPs_ALL[((i * yprocs + j) * SPACE * DIM) + (REAL * DIM) + 0] = npx[i];
            NPs_ALL[((i * yprocs + j) * SPACE * DIM) + (REAL * DIM) + 1] = nqy[j];
            NPs_ALL[((i * yprocs + j) * SPACE * DIM) + (REAL * DIM) + 2] = NZ;
            NPs_ALL[((i * yprocs + j) * SPACE * DIM) + (SPECTRUM * DIM) + 0] = NX;
            NPs_ALL[((i * yprocs + j) * SPACE * DIM) + (SPECTRUM * DIM) + 1] = npy[i];
            NPs_ALL[((i * yprocs + j) * SPACE * DIM) + (SPECTRUM * DIM) + 2] = nqz_[j];
        }
    }
    // Mesh MAX
    for (int i = 0; i < SPACE; i++){
        for (int j = 0; j < DIM; j++){
            nps_max[i][j] = INT_MIN;
        }
    }
    for (int i = 0; i < procs; i++){
        for (int j = 0; j < SPACE; j++){
            for (int k = 0; k < DIM; k++){
                if (nps_max[j][k] < NPs_ALL[(i * SPACE * DIM) + (j * DIM) + k]) {
                    nps_max[j][k] = NPs_ALL[(i * SPACE * DIM) + (j * DIM) + k];
                }
            }
        }
    }
    dmy[0] = nps_max[REAL][0] * nps_max[REAL][1] * NZ_;
    dmy[1] = nps_max[REAL][0] * NY * nps_max[SPECTRUM][2];
    dmy[2] = NX * nps_max[SPECTRUM][1] * nps_max[SPECTRUM][2];
    mesh_size = 0;
    for(int i = 0; i < DIM; i++){
        if (mesh_size <= dmy[i]) mesh_size = dmy[i];
    }
// Process MIN - MAX (mesh)
    PREV_NPs[REAL][0] = 0;
    PREV_NPs[REAL][1] = 0;
    PREV_NPs[REAL][2] = 0;
    PREV_NPs[SPECTRUM][0] = 0;
    PREV_NPs[SPECTRUM][1] = 0;
    PREV_NPs[SPECTRUM][2] = 0;
    for(int i = 0; i < xid; i++) {
        PREV_NPs[REAL][0] += npx[i];
        PREV_NPs[SPECTRUM][1] += npy[i];
    }
    for(int i = 0; i < yid; i++) {
        PREV_NPs[REAL][1] += nqy[i];
        PREV_NPs[SPECTRUM][2] += nqz_[i];
    }
    NEXT_NPs[REAL][0] = 0;
    NEXT_NPs[REAL][1] = 0;
    NEXT_NPs[REAL][2] = NZ;
    NEXT_NPs[SPECTRUM][0] = NX;
    NEXT_NPs[SPECTRUM][1] = 0;
    NEXT_NPs[SPECTRUM][2] = 0;
    for(int i = 0; i <= xid; i++) {
        NEXT_NPs[REAL][0] += npx[i];
        NEXT_NPs[SPECTRUM][1] += npy[i];
    }
    for(int i = 0; i <= yid; i++) {
        NEXT_NPs[REAL][1] += nqy[i];
        NEXT_NPs[SPECTRUM][2] += nqz_[i];
    }
#ifdef _MPI
    rcounts = alloc_1d_int(procs);
    displs = alloc_1d_int(procs + 1);
    displs[0] = 0;
    for(int i = 0; i < xprocs; i++) {
        for(int j = 0; j < yprocs; j++) {
            id = i * yprocs + j;
            rcounts[id] = npx[i] * nqy[j] * NZ_;
            displs[id + 1] = displs[id] + rcounts[id];
        }
    }
#endif
// Process MIN - MAX (real)
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < DIM; j++){
            LPs[i][j] = (NPs[i][j] * DX);
            PREV_LPs[i][j] = (PREV_NPs[i][j] * DX);
            NEXT_LPs[i][j] = (NEXT_NPs[i][j] * DX);
        }
    }
    free_1d_int(npx);
    free_1d_int(npy);
    free_1d_int(nqy);
    free_1d_int(nqz_);
}
/*分割した領域の範囲チェック*/
int Range_check (const Index_range *i_range, Index_range *o_range) {
    int res = ANS_FALSE;

    // X
    if (i_range->iend < PREV_NPs[SPECTRUM][0] || NEXT_NPs[SPECTRUM][0] <= i_range->istart ||
        i_range->jend < PREV_NPs[SPECTRUM][1] || NEXT_NPs[SPECTRUM][1] <= i_range->jstart ||
        i_range->kend < PREV_NPs[SPECTRUM][2] || NEXT_NPs[SPECTRUM][2] <= i_range->kstart ) {
        res = ANS_FALSE;
        o_range->istart = -1;
        o_range->iend = -2;
        o_range->jstart = -1;
        o_range->jend = -2;
        o_range->kstart = -1;
        o_range->kend = -2;
    } else {
        res = ANS_TRUE;
        if (PREV_NPs[SPECTRUM][0] <= i_range->istart && i_range->istart < NEXT_NPs[SPECTRUM][0]) {
            o_range->istart = i_range->istart - PREV_NPs[SPECTRUM][0];
        } else {
            o_range->istart = 0;
        }
        if (PREV_NPs[SPECTRUM][0] <= i_range->iend && i_range->iend < NEXT_NPs[SPECTRUM][0]) {
            o_range->iend = i_range->iend - PREV_NPs[SPECTRUM][0];
        } else {
            o_range->iend = NPs[SPECTRUM][0] - 1;
        }
        if (PREV_NPs[SPECTRUM][1] <= i_range->jstart && i_range->jstart < NEXT_NPs[SPECTRUM][1]) {
            o_range->jstart = i_range->jstart - PREV_NPs[SPECTRUM][1];
        } else {
            o_range->jstart = 0;
        }
        if (PREV_NPs[SPECTRUM][1] <= i_range->jend && i_range->jend < NEXT_NPs[SPECTRUM][1]) {
            o_range->jend = i_range->jend - PREV_NPs[SPECTRUM][1];
        } else {
            o_range->jend = NPs[SPECTRUM][1] - 1;
        }
        if (PREV_NPs[SPECTRUM][2] <= i_range->kstart && i_range->kstart < NEXT_NPs[SPECTRUM][2]) {
            o_range->kstart = i_range->kstart - PREV_NPs[SPECTRUM][2];
        } else {
            o_range->kstart = 0;
        }
        if (PREV_NPs[SPECTRUM][2] <= i_range->kend && i_range->kend < NEXT_NPs[SPECTRUM][2]) {
            o_range->kend = i_range->kend - PREV_NPs[SPECTRUM][2];
        } else {
            o_range->kend = NPs[SPECTRUM][2] - 1;
        }
    }
    return res;
}
/*分割した粒子の範囲チェック*/
int Range_coord (const int *i_mesh, int *o_mesh) {
    int res = ANS_FALSE;

    if ((i_mesh[0] >= PREV_NPs[REAL][0] && NEXT_NPs[REAL][0] > i_mesh[0])
     && (i_mesh[1] >= PREV_NPs[REAL][1] && NEXT_NPs[REAL][1] > i_mesh[1])
     && (i_mesh[2] >= PREV_NPs[REAL][2] && NEXT_NPs[REAL][2] > i_mesh[2])) {
        res = ANS_TRUE;
        o_mesh[0] = i_mesh[0] - PREV_NPs[REAL][0];
        o_mesh[1] = i_mesh[1] - PREV_NPs[REAL][1];
        o_mesh[2] = i_mesh[2] - PREV_NPs[REAL][2];
    } else {
        res = ANS_FALSE;
        o_mesh[0] = -1;
        o_mesh[1] = -1;
        o_mesh[2] = -1;
    }
    return res;
}
/*領域管理関数群*/
void Mem_alloc_var (double **zeta) {
    p_tmp = (Particle *) malloc (Particle_Number * sizeof (Particle) );
    alloc_error_check (p_tmp);
    if (SW_EQ == Navier_Stokes || SW_EQ == Shear_Navier_Stokes || SW_EQ == Shear_Navier_Stokes_Lees_Edwards) {
        ucp = (double **) malloc (sizeof (double *) * DIM);
        alloc_error_check (ucp);
        for (int d = 0; d < DIM; d++) {
            ucp[d] = alloc_1d_double (mesh_size);
        }
    } else if (SW_EQ == Electrolyte) {
        Valency = alloc_1d_double (N_spec);
        Onsager_coeff = alloc_1d_double (N_spec);
        Valency_e = alloc_1d_double (N_spec);
        Total_solute = alloc_1d_double (N_spec);
        Surface_normal = (double **) malloc (sizeof (double *) * DIM);
        alloc_error_check (Surface_normal);
        Concentration = (double **) malloc (sizeof (double *) * N_spec);
        alloc_error_check (Concentration);
        Concentration_rhs0 = (double **) malloc (sizeof (double *) * N_spec);
        alloc_error_check (Concentration_rhs0);
        Concentration_rhs1 = (double **) malloc (sizeof (double *) * N_spec);
        alloc_error_check (Concentration_rhs1);
        for (int n = 0; n < N_spec; n++) {
            Concentration[n] = alloc_1d_double (mesh_size);
            Concentration_rhs0[n] = alloc_1d_double (mesh_size);
            Concentration_rhs1[n] = alloc_1d_double (mesh_size);
        }
        for (int d = 0; d < DIM; d++) {
            Surface_normal[d] = alloc_1d_double (mesh_size);
        }
    }
    Pressure = alloc_1d_double (mesh_size);
    f_ns0 = (double **) malloc (sizeof (double *) * DIM - 1);
    alloc_error_check (f_ns0);
    f_ns1 = (double **) malloc (sizeof (double *) * DIM - 1);
    alloc_error_check (f_ns1);
    f_ns2 = (double **) malloc (sizeof (double *) * DIM - 1);
    alloc_error_check (f_ns2);
    f_ns3 = (double **) malloc (sizeof (double *) * DIM - 1);
    alloc_error_check (f_ns3);
    f_ns4 = (double **) malloc (sizeof (double *) * DIM - 1);
    alloc_error_check (f_ns4);
    Shear_force = (double **) malloc (sizeof (double *) * DIM);
    alloc_error_check (Shear_force);
    Shear_force_k = (double **) malloc (sizeof (double *) * DIM);
    alloc_error_check (Shear_force_k);
    f_particle = (double **) malloc (sizeof (double *) * DIM);
    alloc_error_check (f_particle);
    u = (double **) malloc (sizeof (double *) * DIM);
    alloc_error_check (u);
    up = (double **) malloc (sizeof (double *) * DIM);
    alloc_error_check (up);
    work_v3 = (double **) malloc(sizeof(double *) * DIM);
    alloc_error_check (work_v3);
    work_v4 = (double **) malloc(sizeof(double *) * DIM);
    alloc_error_check (work_v4);
    for (int d = 0; d < DIM - 1; d++) {
        f_ns0[d] = alloc_1d_double (mesh_size);
        f_ns1[d] = alloc_1d_double (mesh_size);
        f_ns2[d] = alloc_1d_double (mesh_size);
        f_ns3[d] = alloc_1d_double (mesh_size);
        f_ns4[d] = alloc_1d_double (mesh_size);
        zeta[d] = alloc_1d_double (mesh_size);
    }
    for (int d = 0; d < DIM; d++) {
        Shear_force[d] = alloc_1d_double (mesh_size);
        Shear_force_k[d] = alloc_1d_double (mesh_size);
        f_particle[d] = alloc_1d_double (mesh_size);
        u[d] = alloc_1d_double (mesh_size);
        up[d] = alloc_1d_double (mesh_size);
        work_v3[d] = alloc_1d_double(NX*NY*NZ_);
        work_v4[d] = alloc_1d_double(NX*NY*NZ_);
    }
    work_v2 = (double **) malloc(sizeof(double *) * (DIM - 1));
    for(int d=0; d < DIM - 1; d++){
        work_v2[d] = alloc_1d_double(mesh_size);
    }

    phi = alloc_1d_double (mesh_size);
    rhop = alloc_1d_double (mesh_size);
    IK2 = alloc_1d_double (mesh_size);
    K2 = alloc_1d_double (mesh_size);
    KX_int = alloc_1d_int (mesh_size);
    KY_int = alloc_1d_int (mesh_size);
    KZ_int = alloc_1d_int (mesh_size);

    phi_sum = alloc_1d_double(mesh_size);
    work_v1 = alloc_1d_double(mesh_size);

    if (SW_EQ == Shear_Navier_Stokes_Lees_Edwards) {
        Hydro_force = alloc_1d_double(mesh_size);
        Hydro_force_new = alloc_1d_double(mesh_size);
        h_obl_force = alloc_1d_double (DIM*Particle_Number);
        h_obl_force_local = alloc_1d_double (DIM*Particle_Number);
        if(ROTATION){
            h_obl_torque = alloc_1d_double (DIM*Particle_Number);
            h_obl_torque_local = alloc_1d_double (DIM*Particle_Number);
        }
        h_obl_volume = alloc_1d_double (Particle_Number);
        h_obl_Itrace = alloc_1d_double (Particle_Number);
        h_obl_volume_local = alloc_1d_double (Particle_Number);
        h_obl_Itrace_local = alloc_1d_double (Particle_Number);
	}

    rand_num_particle = alloc_1d_double (Particle_Number*6);

#ifdef _OPENMP
    tmp_buffer1 = alloc_1d_double (THREADNUM * mesh_size);
    tmp_buffer2 = alloc_1d_double (THREADNUM * mesh_size);
    tmp_buffer_dim = alloc_2d_double (DIM, THREADNUM * mesh_size);
#endif
}


void Mem_free_var (double **zeta) {
    free (p_tmp);
    free_1d_int (NPs_ALL);
#ifdef _MPI
    free_1d_int (rcounts);
    free_1d_int (displs);
#endif
    free_1d_double (rand_num_particle);
    if (SW_EQ == Shear_Navier_Stokes_Lees_Edwards) {
        free_1d_double (h_obl_Itrace_local);
        free_1d_double (h_obl_volume_local);
        free_1d_double (h_obl_Itrace);
        free_1d_double (h_obl_volume);
        if(ROTATION){
            free_1d_double (h_obl_torque_local);
            free_1d_double (h_obl_torque);
		}
        free_1d_double (h_obl_force_local);
        free_1d_double (h_obl_force);
        free_1d_double (Hydro_force_new);
        free_1d_double (Hydro_force);
	}
    for (int d = 0; d < DIM; d++) {
        free_1d_double (work_v4[d]);
        free_1d_double (work_v3[d]);
    }
    for (int d = 0; d < DIM - 1; d++) {
        free_1d_double (work_v2[d]);
    }
    free(work_v4);
    free(work_v3);
    free(work_v2);
#ifdef _OPENMP
    free_1d_double (tmp_buffer1);
    free_1d_double (tmp_buffer2);
    free_2d_double (tmp_buffer_dim);
#endif
    free_1d_int (KZ_int);
    free_1d_int (KY_int);
    free_1d_int (KX_int);
    free_1d_double (K2);
    free_1d_double (IK2);
    free_1d_double (rhop);
    free_1d_double (phi);
    for (int d = 0; d < DIM; d++) {
        free_1d_double (up[d]);
        free_1d_double (u[d]);
        free_1d_double (f_particle[d]);
        free_1d_double (Shear_force_k[d]);
        free_1d_double (Shear_force[d]);
    }
    free (up);
    free (u);
    free (f_particle);
    free (Shear_force);
    free (Shear_force_k);
    for (int d = 0; d < DIM - 1; d++) {
        free_1d_double (zeta[d]);
        free_1d_double (f_ns4[d]);
        free_1d_double (f_ns3[d]);
        free_1d_double (f_ns2[d]);
        free_1d_double (f_ns1[d]);
        free_1d_double (f_ns0[d]);
    }
    free (zeta);
    free (f_ns4);
    free (f_ns3);
    free (f_ns2);
    free (f_ns1);
    free (f_ns0);
    free_1d_double (Pressure);
    if (SW_EQ == Navier_Stokes || SW_EQ == Shear_Navier_Stokes) {
        for (int d = 0; d < DIM; d++) {
            free_1d_double (ucp[d]);
        }
        free (ucp);
    } else if (SW_EQ == Electrolyte) {
        for (int d = 0; d < DIM; d++) {
            free_1d_double (Surface_normal[d]);
        }
        for (int n = 0; n < N_spec; n++) {
            free_1d_double (Concentration_rhs1[n]);
            free_1d_double (Concentration_rhs0[n]);
            free_1d_double (Concentration[n]);
        }
        free (Concentration_rhs1);
        free (Concentration_rhs0);
        free (Concentration);
        free (Surface_normal);
        free_1d_double (Total_solute);
        free_1d_double (Valency_e);
        free_1d_double (Onsager_coeff);
        free_1d_double (Valency);
    }
    switch (SW_PT) {
        case spherical_particle:
            free_1d_int (Particle_Numbers);
            free_1d_double (MASS_RATIOS);
            free_1d_double (RHO_particle);
            free_1d_double (MASS);
            free_1d_double (IMASS);
            free_1d_double (IMASS_RATIOS);
            free_1d_double (MOI);
            free_1d_double (IMOI);
            free_1d_double (S_surfaces);
            free_1d_double (W_surfaces);
            free_1d_double (Surface_charge);
            free_1d_double (Surface_charge_e);
            break;
        case chain:
            free_1d_int (Particle_Numbers);
            free_1d_int (Beads_Numbers);
            free_1d_int (Chain_Numbers);
            free_1d_double (MASS_RATIOS);
            free_1d_double (RHO_particle);
            free_1d_double (MASS);
            free_1d_double (IMASS);
            free_1d_double (IMASS_RATIOS);
            free_1d_double (MOI);
            free_1d_double (IMOI);
            free_1d_double (S_surfaces);
            free_1d_double (W_surfaces);
            free_1d_double (Surface_charge);
            free_1d_double (Surface_charge_e);
            break;
        default:
            break;
    }
    if (PINNING) {
        free_1d_int (Pinning_ROT_Numbers);
        free_1d_int (Pinning_Numbers);
    }
    //free_1d_int (ip);
    //free_1d_double (t);
    //free_1d_double (w);
    delete[] ijk_range_two_third_filter;
    free_1d_double(phi_sum);
    free_1d_double(work_v1);
}
