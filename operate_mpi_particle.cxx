#include "operate_mpi_particle.h"
#ifdef _MPI
/*粒子領域内外のチェック*/
int chk_coord (const double *x) {
    int ans = ANS_FALSE;

    if ((x[0] >= PREV_LPs[REAL][0] && x[0] < NEXT_LPs[REAL][0]) &&
        (x[1] >= PREV_LPs[REAL][1] && x[1] < NEXT_LPs[REAL][1]) &&
        (x[2] >= PREV_LPs[REAL][2] && x[2] < NEXT_LPs[REAL][2])) {
        ans = ANS_TRUE;
    } else {
        ans = ANS_FALSE;
    }
    return ans;
}
/*ID重複のチェック(重複防止)*/
int chk_id (int idmax, int *id, int pid) {
    int i, ans = ANS_TRUE;

    for (i = 0; i <= idmax; i++) {
        if (id[i] == pid) {
            ans = ANS_FALSE;
            break;
        }
    }
    return ans;
}
/*粒子ソート(比較関数)*/
int particle_comp (const void *p0, const void *p1) {
    const Particle *p0_ = (Particle *) p0;
    const Particle *p1_ = (Particle *) p1;
    if (p0_->id < p1_->id) {
        return -1;
    } else if (p0_->id > p1_->id) {
        return 1;
    } else {
        return 0;
    }
}
/*粒子ソート(ソート関数)*/
void Particle_qsort (Particle *p, int n) {
    qsort (p, n, sizeof (Particle), particle_comp);
}
/*Datatype definition for communication*/
void Build_Particle_Datatypes (void) {
    Particle dummy_particle;
#if !defined (NDEBUG)
    mt_struct dummy_mt;
    int cnt_mt[7], cnt_particle[5];
    MPI_Aint disp_mt[7], disp_particle[5];
    MPI_Datatype type_mt[7], type_particle[5];
    MPI_Datatype VECTOR_3D, MT;

    // 3d-vector array of double
    ierr = MPI_Type_contiguous (DIM, MPI_DOUBLE, &VECTOR_3D);
    MPI_ERRORCHECK (ierr);
    // setting random number of particle structure
    cnt_mt[0] = 1;
    cnt_mt[1] = 4;
    cnt_mt[2] = 3;
    cnt_mt[3] = 4;
    cnt_mt[4] = 2;
    cnt_mt[5] = 2;
    cnt_mt[6] = (int) sizeof (uint32_t *);
    ierr = MPI_Address (&dummy_mt.aaa, &disp_mt[0]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_mt.mm, &disp_mt[1]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_mt.wmask, &disp_mt[2]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_mt.shift0, &disp_mt[3]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_mt.maskB, &disp_mt[4]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_mt.i, &disp_mt[5]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_mt.state, &disp_mt[6]);
    MPI_ERRORCHECK (ierr);
    for (int i = 6; i >= 0; --i) {
        disp_mt[i] -= disp_mt[0];
    }
    type_mt[0] = MPI_UNSIGNED;
    type_mt[1] = MPI_INT;
    type_mt[2] = MPI_UNSIGNED;
    type_mt[3] = MPI_INT;
    type_mt[4] = MPI_UNSIGNED;
    type_mt[5] = MPI_INT;
    type_mt[6] = MPI_CHAR;
    ierr = MPI_Type_struct (7, cnt_mt, disp_mt, type_mt, &MT);
    MPI_ERRORCHECK (ierr);
    // particle structure
    cnt_particle[0] = 1;
    cnt_particle[1] = 2;
    cnt_particle[2] = MT_ARRAY;
    cnt_particle[3] = 21;
    cnt_particle[4] = 1;
    ierr = MPI_Address (&dummy_particle.eff_mass_ratio, &disp_particle[0]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_particle.spec, &disp_particle[1]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_particle.state, &disp_particle[2]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_particle.x, &disp_particle[3]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_particle.rdata, &disp_particle[4]);
    MPI_ERRORCHECK (ierr);
    for (int i = 4; i >= 0; --i) {
        disp_particle[i] -= disp_particle[0];
    }
    type_particle[0] = MPI_DOUBLE;
    type_particle[1] = MPI_INT;
    type_particle[2] = MPI_UNSIGNED;
    type_particle[3] = VECTOR_3D;
    type_particle[4] = MT;
    ierr = MPI_Type_struct (5, cnt_particle, disp_particle, type_particle, &PT);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Type_commit (&PT);
    MPI_ERRORCHECK (ierr);
#else
    drand48_data dummy_data;
	quaternion dummy_q;
    int cnt_data[2], cnt_particle[6];
    MPI_Aint disp_data[2], disp_particle[6];
    MPI_Datatype type_data[2], type_particle[6];
    MPI_Datatype VECTOR_3D, DRAND, QUATERNION, MATRIX_3x3;

    // 3d-vector array of double
    ierr = MPI_Type_contiguous (DIM, MPI_DOUBLE, &VECTOR_3D);
    MPI_ERRORCHECK (ierr);
    // matrix 3x3 array of double
    ierr = MPI_Type_contiguous (9, MPI_DOUBLE, &MATRIX_3x3);
    MPI_ERRORCHECK (ierr);
    // setting random number of particle structure
    cnt_data[0] = 8;
    cnt_data[1] = 1;
    ierr = MPI_Address (&dummy_data.__x, &disp_data[0]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_data.__a, &disp_data[1]);
    MPI_ERRORCHECK (ierr);
    for (int i = 1; i >= 0; --i) {
        disp_data[i] -= disp_data[0];
    }
    type_data[0] = MPI_SHORT;
    type_data[1] = MPI_LONG_LONG_INT;
    ierr = MPI_Type_struct (2, cnt_data, disp_data, type_data, &DRAND);
    MPI_ERRORCHECK (ierr);
    // setting quanternion
    cnt_data[0] = 1;
    cnt_data[1] = 1;
	ierr = MPI_Address (&dummy_q.s, &disp_data[0]);
    MPI_ERRORCHECK (ierr);
	ierr = MPI_Address (&dummy_q.v, &disp_data[1]);
    MPI_ERRORCHECK (ierr);
    for (int i = 1; i >= 0; --i) {
        disp_data[i] -= disp_data[0];
    }
    type_data[0] = MPI_DOUBLE;
    type_data[1] = VECTOR_3D;
    ierr = MPI_Type_struct (2, cnt_data, disp_data, type_data, &QUATERNION);
    MPI_ERRORCHECK (ierr);
    // particle structure
    cnt_particle[0] = 2;   //number of "double" type data in "Particle"
    cnt_particle[1] = 2;   //int
    cnt_particle[2] = 28;  //double[DIM]
    cnt_particle[3] = 1;   //rand48
    cnt_particle[4] = 4;   //double[DIM][DIM]
    cnt_particle[5] = 2;   //quanternion
    ierr = MPI_Address (&dummy_particle.mass, &disp_particle[0]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_particle.spec, &disp_particle[1]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_particle.x, &disp_particle[2]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_particle.rdata, &disp_particle[3]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_particle.QR, &disp_particle[4]);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Address (&dummy_particle.q, &disp_particle[5]);
    MPI_ERRORCHECK (ierr);
    for (int i = 5; i >= 0; --i) {
        disp_particle[i] -= disp_particle[0];
    }
    type_particle[0] = MPI_DOUBLE;
    type_particle[1] = MPI_INT;
    type_particle[2] = VECTOR_3D;
    type_particle[3] = DRAND;
    type_particle[4] = MATRIX_3x3;
    type_particle[5] = QUATERNION;
    ierr = MPI_Type_struct (6, cnt_particle, disp_particle, type_particle, &PT);
    MPI_ERRORCHECK (ierr);
    ierr = MPI_Type_commit (&PT);
    MPI_ERRORCHECK (ierr);
#endif

}
/* Particles send - recv (from any rank, to any rank)*/
void Particle_Sendrecv (Particle *p, Particle *new_p, enum SW_COMM sw_comm) {
    int id_base, id, x, y, num = 0;

    if (sw_comm == ONE_TO_MANY || sw_comm == MANY_TO_MANY) {
        if (procid == root) {
            fprintf (stderr, "SEND/RECV SW_COMM Error");
            fflush (stderr);
        }
        MPI_Barrier (MPI_COMM_WORLD);
        exit (EXIT_FAILURE);
    }
    int *sendcounts = alloc_1d_int (procs);
    int *recvcounts = alloc_1d_int (procs);
    int *sendranks = alloc_1d_int (procs);
    int *recvranks = alloc_1d_int (procs);
    int *sendbuf = alloc_1d_int (2 * procs);
    int *recvbuf = alloc_1d_int (2 * procs);
    Particle **send_p = (Particle **) malloc (procs * sizeof (Particle *));
    alloc_error_check (send_p);
    for (int i = 0; i < procs; i++) {
        send_p[i] = (Particle *) calloc (Update_Particle_Number, sizeof (Particle) );
        alloc_error_check (send_p[i]);
    }
    if (procs > 1) {
        //初期化
        for (int i = 0; i < procs; i++) {
            sendranks[i] = MPI_PROC_NULL;
            sendcounts[i] = 0;
            recvranks[i] = MPI_PROC_NULL;
            recvcounts[i] = 0;
        }
        //送信準備
        if (sw_comm == ONE_TO_ONE_LJ) {
            for (int n = 0; n < Update_Particle_Number; n++) {
                x = (int)(p[n].x[0] * IDX);
                y = (int)(p[n].x[1] * IDX);
                //粒子が本来いるべきプロセスを求める
                id_base = ID(xmesh[x], ymesh[y]);
                for (int i = 0; i < lj_size; i++) {
                    //求めたプロセスの周囲に存在するプロセスにも送信するため
                    //送信用バッファに詰め込む(LJ Cutoff用)
                    id = LJ_CUTOFF_REFID(id_base, i);
                    sendranks[id] = id;
                    recvranks[id] = procid;
                    send_p[id][sendcounts[id]] = p[n];
                    sendcounts[id]++;
                    recvcounts[id]++;
                }
            }
        } else if (sw_comm == ONE_TO_ONE_SEKIBUN) {
            for (int n = 0; n < Update_Particle_Number; n++) {
                x = (int)(p[n].x[0] * IDX);
                y = (int)(p[n].x[1] * IDX);
                //粒子が本来いるべきプロセスを求める
                id_base = ID(xmesh[x], ymesh[y]);
                for (int i = 0; i < sekibun_size; i++) {
                    //求めたプロセスの周囲に存在するプロセスにも送信するため
                    //送信用バッファに詰め込む(Sekibun_cell用)
                    id = SEKIBUN_REFID(id_base, i);
                    sendranks[id] = id;
                    recvranks[id] = procid;
                    send_p[id][sendcounts[id]] = p[n];
                    sendcounts[id]++;
                    recvcounts[id]++;
                }
            }
        }
        //粒子の送信ランク、個数の通知
        for (int i = 0; i < procs; i++) {
            sendbuf[2 * i] = recvranks[i];
            sendbuf[2 * i + 1] = recvcounts[i];
        }

        ierr = MPI_Alltoall (sendbuf, 2, MPI_INT, recvbuf, 2, MPI_INT, MPI_COMM_WORLD);

        MPI_ERRORCHECK (ierr);
        for (int i = 0; i < procs; i++) {
            recvranks[i] = recvbuf[2 * i];
            recvcounts[i] = recvbuf[2 * i + 1];
        }
        //粒子の送受信
        //場から粒子の作用を計算する際、粒子格納順序を保障する必要があるため
        //プロセスIDの昇順に通信を行い粒子を格納する必要がある
        num = 0;
        for (int i = 0; i < procs; i++) {

            ierr = MPI_Irecv ( (new_p + num), recvcounts[i], PT, recvranks[i],
                              TAG (procid, recvranks[i]), MPI_COMM_WORLD,
                              (ireq + i) );

            MPI_ERRORCHECK (ierr);
            num += recvcounts[i];
        }
        for (int i = 0; i < procs; i++) {

            ierr = MPI_Isend (send_p[i], sendcounts[i], PT, sendranks[i],
                              TAG (sendranks[i], procid), MPI_COMM_WORLD,
                              (ireq + procs + i) );

            MPI_ERRORCHECK (ierr);
        }

        ierr = MPI_Waitall ( (2 * procs), ireq, ista);

        MPI_ERRORCHECK (ierr);
    } else {
        for (int n = 0; n < Update_Particle_Number; n++) {
            new_p[n] = p[n];
        }
        num = Particle_Number;
    }
    Update_Particle_Number = num;
    Reference_Particle_Number = 0;
    Local_Particle_Number = num;
#if !defined (NDEBUG)
    for (int i = 0; i < Local_Particle_Number; i++) {
        new_p[i].rdata.state = &(new_p[i].state[0]);
    }
#endif
    for (int i = 0; i < procs; i++) {
        free (send_p[i]);
    }
    free (send_p);
    free_1d_int (sendcounts);
    free_1d_int (recvcounts);
    free_1d_int (sendranks);
    free_1d_int (recvranks);
    free_1d_int (sendbuf);
    free_1d_int (recvbuf);
}
/* Particles gather (from all, to (sw = SW_ON) ? all : (sw = SW_OFF) ? root) */
void Particle_Gather (Particle *p, Particle *new_p, enum SW_SWITCH sw) {
    int *recvcounts, *displs;
    if (procs > 1) {
        recvcounts = alloc_1d_int (procs);
        displs = alloc_1d_int (procs + 1);
        // gather particles number of each process
        // SW_ON : All, SW_OFF : single(root)
        if (sw == SW_OFF) {

            ierr = MPI_Gather (&Update_Particle_Number, 1, MPI_INT, recvcounts,
                               1, MPI_INT, root, MPI_COMM_WORLD);

            MPI_ERRORCHECK (ierr);
        } else {

            ierr = MPI_Allgather (&Update_Particle_Number, 1, MPI_INT,
                                  recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

            MPI_ERRORCHECK (ierr);
        }
        displs[0] = 0;
        for (int prc = 0; prc < procs; prc++) {
            displs[prc + 1] = displs[prc] + recvcounts[prc];
        }
        // gather particles of each process
        // SW_ON : All, SW_OFF : single(root)
        if (sw == SW_OFF) {

            ierr = MPI_Gatherv (p, Update_Particle_Number, PT, new_p, recvcounts,
                                displs, PT, root, MPI_COMM_WORLD);

            MPI_ERRORCHECK (ierr);
        } else {

            ierr = MPI_Allgatherv (p, Update_Particle_Number, PT, new_p,
                                   recvcounts, displs, PT, MPI_COMM_WORLD);

            MPI_ERRORCHECK (ierr);
        }
        free_1d_int (recvcounts);
        free_1d_int (displs);
    } else {
        for (int n = 0; n < Update_Particle_Number; n++) {
            new_p[n] = p[n];
        }
    }
#if !defined (NDEBUG)
    for (int n = 0; n < Particle_Number; n++) {
        new_p[n].rdata.state = &(new_p[n].state[0]);
    }
#endif

}
/* Particles broadcast (from root, to all)*/
void Particle_Bcast (Particle *p) {

    if (procs > 1) {

        ierr = MPI_Bcast (p, Particle_Number, PT, root, MPI_COMM_WORLD);

        MPI_ERRORCHECK (ierr);
    }
#if !defined (NDEBUG)
    for (int i = 0; i < Particle_Number; i++) {
        p[i].rdata.state = &(p[i].state[0]);
    }
#endif

}
/* Particles collective communication */
void Particle_Group_Communication (Particle *p, enum SW_COMM sw_comm) {
    static Particle zero = {0.0};
    int num, pnum;
    int *id;

    id = alloc_1d_int (Particle_Number);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (Particle_Number, id, p_tmp, zero)
#endif
    for (int i = 0; i < Particle_Number; i++){
        id[i] = 0;
        p_tmp[i] = zero;
    }
    // Particle communication switches
    if (sw_comm == ONE_TO_ONE_LJ || sw_comm == ONE_TO_ONE_SEKIBUN ) {
        Particle_Sendrecv (p, p_tmp, sw_comm);
        pnum = Local_Particle_Number;
    } else if (sw_comm == ONE_TO_MANY) {
        Particle_Bcast (p);
        pnum = Particle_Number;
        for (int n = 0; n < pnum; n++) {
            p_tmp[n] = p[n];
        }
    } else if (sw_comm == MANY_TO_MANY) {
        Particle_Gather (p, p_tmp, SW_ON);
        pnum = Particle_Number;
    } else {
        ;
    }
#ifdef _OPENMP
#pragma omp parallel for default(none) shared (Particle_Number, p, zero)
#endif
    for (int i = 0; i < Particle_Number; i++){
        p[i] = zero;
    }
    num = 0;
    // Selection of update particle
    for (int n = 0; n < pnum; n++) {
        if (chk_coord (p_tmp[n].x) ) {
            p[num] = p_tmp[n];
#if !defined (NDEBUG)
            p[n].rdata.state = &(p[n].state[0]);
#endif
            id[n]++;
            num++;
        }
    }
    Update_Particle_Number = num;
    // Selection of reference particle
    for (int n = 0; n < pnum; n++) {
        if (!id[n]) {
            p[num] = p_tmp[n];
#if !defined (NDEBUG)
            p[n].rdata.state = &(p[n].state[0]);
#endif
            num++;
        }
    }
    Reference_Particle_Number = num - Update_Particle_Number;
    Local_Particle_Number = Update_Particle_Number + Reference_Particle_Number;
    free_1d_int (id);

}
#endif /* _MPI */
