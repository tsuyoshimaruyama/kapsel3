/*フーリエ変換*/
#ifdef __INTEL_COMPILER
#include <mathimf.h>
#include <tgmath.h>
#endif /* __INTEL_COMPILER */
#include <string.h>
#include <complex>
#include "fft_wrapper_base.h"

#define REALMODE_ARRAYINDEX_DIM(d,i,j,k) (((d) * NPs[REAL][0] * NPs[REAL][1] * NZ_) + ((i) * NPs[REAL][1] * NZ_) + ((j) * NZ_) + (k))
#define SPECTRUMMODE_ARRAYINDEX_DIM(d,i,j,k) (((d) * NPs[SPECTRUM][0] * NPs[SPECTRUM][1] * NPs[SPECTRUM][2]) + ((i) * NPs[SPECTRUM][1] * NPs[SPECTRUM][2]) + ((j) * NPs[SPECTRUM][2]) + (k))
static void Dfti_Parameter_Set();
static void Dfti_SetValue();
#if (defined (_OPENMP) || defined(_MPI))
static MKL_LONG status;
#endif
#if (defined (_OPENMP) && !defined(_MPI))
#ifdef _FFT_IMKL
static DFTI_DESCRIPTOR_HANDLE forward_dfti;
static DFTI_DESCRIPTOR_HANDLE backward_dfti;
#endif
#endif
#if defined(_MPI)
static DFTI_DESCRIPTOR_HANDLE forward_dfti1;
static DFTI_DESCRIPTOR_HANDLE forward_dfti2;
static DFTI_DESCRIPTOR_HANDLE forward_dfti3;
static DFTI_DESCRIPTOR_HANDLE backward_dfti1;
static DFTI_DESCRIPTOR_HANDLE backward_dfti2;
static DFTI_DESCRIPTOR_HANDLE backward_dfti3;
static int npx, npy, nqy, nqz_, hnqz_;
static int max_npx, max_npy, max_nqy, max_nqz_, max_hnqz_;
int *npx_all, *npy_all, *nqy_all, *nqz_all;
int *npx_proc, *npy_proc, *nqy_proc, *nqz_proc;
static int cnt1, cnt2, size;
double *fft, *sbuf, *rbuf;
double *fft_free, *sbuf_free, *rbuf_free;
#endif
#ifndef NDEBUG
#define DFTI_ERRORCHECK(errno) \
    do { \
        if((errno) && !DftiErrorClass((errno), DFTI_NO_ERROR)) { \
            fprintf (stderr, "FFT Error:%s (line:%d)\n", DftiErrorMessage((errno)), __LINE__); \
            exit (EXIT_FAILURE); \
        } \
    } while(0)
#else
#define DFTI_ERRORCHECK(errno)
#endif
//******************************************************************************
// hi-speed memory copy
inline double* calloc_1d_double_alignment(double **p1, double **p2, int size, int alignment) {
	*p2 = calloc_1d_double(size, alignment);
	*p1 = *p2;
	return *p2;
}
//******************************************************************************

//#if defined(_FFT_MPI_IMKL) || defined(_FFT_OMP_MPI_IMKL)
void Init_fft_omp_mpi_imkl(void) {
	Dfti_Parameter_Set();
	Dfti_SetValue();
}
void Dfti_Parameter_Set() {
#if defined(_MPI)
	int nps_max[SPACE][DIM];
	npx_all = alloc_1d_int(xprocs);
	npy_all = alloc_1d_int(xprocs);
	nqy_all = alloc_1d_int(yprocs);
	nqz_all = alloc_1d_int(yprocs);
	npx_proc = calloc_1d_int(xprocs);
	npy_proc = calloc_1d_int(xprocs);
	nqy_proc = calloc_1d_int(yprocs);
	nqz_proc = calloc_1d_int(yprocs);
#endif
#if (defined (_OPENMP) && !defined(_MPI))
#ifdef _FFT_IMKL
	forward_dfti = NULL;
	backward_dfti = NULL;
#endif
#elif defined(_MPI)
	for (int i = 0; i < SPACE; i++) {
		for (int j = 0; j < DIM; j++) {
			nps_max[i][j] = INT_MIN;
		}
	}
	for (int i = 0; i < procs; i++) {
		for (int j = 0; j < SPACE; j++) {
			for (int k = 0; k < DIM; k++) {
				if (nps_max[j][k] < NPs_ALL[(i * SPACE * DIM) + (j * DIM) + k]) {
					nps_max[j][k] = NPs_ALL[(i * SPACE * DIM) + (j * DIM) + k];
				}
			}
		}
	}
	for (int xp = 0; xp < xprocs; xp++) {
		npx_all[xp] = NPs_ALL[((xp * yprocs + 0) * SPACE * DIM) + (REAL * DIM) + 0];
		npy_all[xp] = NPs_ALL[((xp * yprocs + 0) * SPACE * DIM) + (SPECTRUM * DIM) + 1];
	}
	for (int yp = 0; yp < yprocs; yp++) {
		nqy_all[yp] = NPs_ALL[((0 * yprocs + yp) * SPACE * DIM) + (REAL * DIM) + 1];
		nqz_all[yp] = NPs_ALL[((0 * yprocs + yp) * SPACE * DIM) + (SPECTRUM * DIM) + 2];
	}
	for (int xp = 1; xp < xprocs; xp++) {
		npx_proc[xp] = npx_proc[xp - 1] + npx_all[xp - 1];
		npy_proc[xp] = npy_proc[xp - 1] + npy_all[xp - 1];
	}
	for (int yp = 1; yp < yprocs; yp++) {
		nqy_proc[yp] = nqy_proc[yp - 1] + nqy_all[yp - 1];
		nqz_proc[yp] = nqz_proc[yp - 1] + nqz_all[yp - 1];
	}
	max_npx = nps_max[REAL][0];
	max_npy = nps_max[SPECTRUM][1];
	max_nqy = nps_max[REAL][1];
	max_nqz_ = nps_max[SPECTRUM][2];
	npx = NPs[REAL][0];
	npy = NPs[SPECTRUM][1];
	nqy = NPs[REAL][1];
	nqz_ = NPs[SPECTRUM][2];
	hnqz_ = HNQZ_;
	cnt1 = max_npx * max_nqy * max_nqz_;
	cnt2 = max_npx * max_npy * max_nqz_;
	size = (yprocs * cnt1 < xprocs * cnt2) ? xprocs * cnt2 : yprocs * cnt1;
	fft = calloc_1d_double_alignment(&fft_free, &fft, mesh_size, 16);
	sbuf = calloc_1d_double_alignment(&sbuf_free, &sbuf, size, 16);
	rbuf = calloc_1d_double_alignment(&rbuf_free, &rbuf, size, 16);
	forward_dfti1 = NULL;
	forward_dfti2 = NULL;
	forward_dfti3 = NULL;
	backward_dfti1 = NULL;
	backward_dfti2 = NULL;
	backward_dfti3 = NULL;
#endif
}
void Dfti_SetValue() {
#if (defined (_OPENMP) && !defined(_MPI))
#ifdef _FFT_IMKL
	const MKL_LONG length[DIM] = { NX, NY, NZ };
	const MKL_LONG istride[DIM + 1] = { 0, NY * NZ, NZ, 1 };
	const MKL_LONG ostride[DIM + 1] = { 0, NY * HNZ_, HNZ_, 1 };
#endif
#elif defined(_MPI)
	const MKL_LONG length_x = NX;
	const MKL_LONG length_y = NY;
	const MKL_LONG length_z = NZ;
	const MKL_LONG dist_x = 1;
	const MKL_LONG dist_y = 1;
	const MKL_LONG idist_z = NZ_;
	const MKL_LONG odist_z = HNZ_;
	const MKL_LONG stride_x[2] = { 0, npy * hnqz_ };
	const MKL_LONG stride_y[2] = { 0, hnqz_ };
	const MKL_LONG stride_z[2] = { 0, 1 };
	const MKL_LONG ntrans_x = npy * hnqz_;
	const MKL_LONG ntrans_y = hnqz_;
	const MKL_LONG ntrans_z = npx * nqy;
#endif

#if (defined (_OPENMP) && !defined(_MPI))
#ifdef _FFT_IMKL
	//**** A2a_k **************************************************************************************/
	DftiCreateDescriptor(&forward_dfti, DFTI_DOUBLE, DFTI_REAL, 3, length);
	DftiSetValue(forward_dfti, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	DftiSetValue(forward_dfti, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
	DftiSetValue(forward_dfti, DFTI_INPUT_STRIDES, istride);
	DftiSetValue(forward_dfti, DFTI_OUTPUT_STRIDES, ostride);
	DftiCommitDescriptor(forward_dfti);
	//**** A_k2a **************************************************************************************/
	DftiCreateDescriptor(&backward_dfti, DFTI_DOUBLE, DFTI_REAL, 3, length);
	DftiSetValue(backward_dfti, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	DftiSetValue(backward_dfti, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
	DftiSetValue(backward_dfti, DFTI_INPUT_STRIDES, ostride);
	DftiSetValue(backward_dfti, DFTI_OUTPUT_STRIDES, istride);
	DftiCommitDescriptor(backward_dfti);
#endif
#elif defined(_MPI)
	//**** A2a_k First ********************************************************************************/
	DftiCreateDescriptor(&forward_dfti1, DFTI_DOUBLE, DFTI_REAL, 1, length_z);
	DftiSetValue(forward_dfti1, DFTI_PLACEMENT, DFTI_INPLACE);
	DftiSetValue(forward_dfti1, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
	DftiSetValue(forward_dfti1, DFTI_INPUT_STRIDES, stride_z);
	DftiSetValue(forward_dfti1, DFTI_OUTPUT_STRIDES, stride_z);
	DftiSetValue(forward_dfti1, DFTI_INPUT_DISTANCE, idist_z);
	DftiSetValue(forward_dfti1, DFTI_OUTPUT_DISTANCE, odist_z);
	DftiSetValue(forward_dfti1, DFTI_NUMBER_OF_TRANSFORMS, ntrans_z);
	DftiSetValue(forward_dfti1, DFTI_THREAD_LIMIT, 1);
	DftiCommitDescriptor(forward_dfti1);
	//**** A2a_k Second *******************************************************************************/
	DftiCreateDescriptor(&forward_dfti2, DFTI_DOUBLE, DFTI_COMPLEX, 1, length_y);
	DftiSetValue(forward_dfti2, DFTI_PLACEMENT, DFTI_INPLACE);
	DftiSetValue(forward_dfti2, DFTI_COMPLEX_STORAGE, DFTI_COMPLEX_COMPLEX);
	DftiSetValue(forward_dfti2, DFTI_INPUT_STRIDES, stride_y);
	DftiSetValue(forward_dfti2, DFTI_OUTPUT_STRIDES, stride_y);
	DftiSetValue(forward_dfti2, DFTI_INPUT_DISTANCE, dist_y);
	DftiSetValue(forward_dfti2, DFTI_OUTPUT_DISTANCE, dist_y);
	DftiSetValue(forward_dfti2, DFTI_NUMBER_OF_TRANSFORMS, ntrans_y);
	DftiSetValue(forward_dfti2, DFTI_THREAD_LIMIT, 1);
	DftiCommitDescriptor(forward_dfti2);
	//**** A2a_k Third ********************************************************************************/
	DftiCreateDescriptor(&forward_dfti3, DFTI_DOUBLE, DFTI_COMPLEX, 1, length_x);
	DftiSetValue(forward_dfti3, DFTI_PLACEMENT, DFTI_INPLACE);
	DftiSetValue(forward_dfti3, DFTI_COMPLEX_STORAGE, DFTI_COMPLEX_COMPLEX);
	DftiSetValue(forward_dfti3, DFTI_INPUT_STRIDES, stride_x);
	DftiSetValue(forward_dfti3, DFTI_OUTPUT_STRIDES, stride_x);
	DftiSetValue(forward_dfti3, DFTI_INPUT_DISTANCE, dist_x);
	DftiSetValue(forward_dfti3, DFTI_OUTPUT_DISTANCE, dist_x);
	DftiSetValue(forward_dfti3, DFTI_NUMBER_OF_TRANSFORMS, ntrans_x);
	DftiSetValue(forward_dfti3, DFTI_THREAD_LIMIT, 1);
	DftiCommitDescriptor(forward_dfti3);
	//**** A_k2a First ********************************************************************************/
	DftiCreateDescriptor(&backward_dfti1, DFTI_DOUBLE, DFTI_COMPLEX, 1, length_x);
	DftiSetValue(backward_dfti1, DFTI_PLACEMENT, DFTI_INPLACE);
	DftiSetValue(backward_dfti1, DFTI_COMPLEX_STORAGE, DFTI_COMPLEX_COMPLEX);
	DftiSetValue(backward_dfti1, DFTI_INPUT_STRIDES, stride_x);
	DftiSetValue(backward_dfti1, DFTI_OUTPUT_STRIDES, stride_x);
	DftiSetValue(backward_dfti1, DFTI_INPUT_DISTANCE, dist_x);
	DftiSetValue(backward_dfti1, DFTI_OUTPUT_DISTANCE, dist_x);
	DftiSetValue(backward_dfti1, DFTI_NUMBER_OF_TRANSFORMS, ntrans_x);
	DftiSetValue(backward_dfti1, DFTI_THREAD_LIMIT, 1);
	DftiCommitDescriptor(backward_dfti1);
	//**** A_k2a Second *******************************************************************************/
	DftiCreateDescriptor(&backward_dfti2, DFTI_DOUBLE, DFTI_COMPLEX, 1, length_y);
	DftiSetValue(backward_dfti2, DFTI_PLACEMENT, DFTI_INPLACE);
	DftiSetValue(backward_dfti2, DFTI_COMPLEX_STORAGE, DFTI_COMPLEX_COMPLEX);
	DftiSetValue(backward_dfti2, DFTI_INPUT_STRIDES, stride_y);
	DftiSetValue(backward_dfti2, DFTI_OUTPUT_STRIDES, stride_y);
	DftiSetValue(backward_dfti2, DFTI_INPUT_DISTANCE, dist_y);
	DftiSetValue(backward_dfti2, DFTI_OUTPUT_DISTANCE, dist_y);
	DftiSetValue(backward_dfti2, DFTI_NUMBER_OF_TRANSFORMS, ntrans_y);
	DftiSetValue(backward_dfti2, DFTI_THREAD_LIMIT, 1);
	DftiCommitDescriptor(backward_dfti2);
	//**** A_k2a Third ********************************************************************************/
	DftiCreateDescriptor(&backward_dfti3, DFTI_DOUBLE, DFTI_REAL, 1, length_z);
	DftiSetValue(backward_dfti3, DFTI_PLACEMENT, DFTI_INPLACE);
	DftiSetValue(backward_dfti3, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
	DftiSetValue(backward_dfti3, DFTI_INPUT_STRIDES, stride_z);
	DftiSetValue(backward_dfti3, DFTI_OUTPUT_STRIDES, stride_z);
	DftiSetValue(backward_dfti3, DFTI_INPUT_DISTANCE, odist_z);
	DftiSetValue(backward_dfti3, DFTI_OUTPUT_DISTANCE, idist_z);
	DftiSetValue(backward_dfti3, DFTI_NUMBER_OF_TRANSFORMS, ntrans_z);
	DftiSetValue(backward_dfti3, DFTI_THREAD_LIMIT, 1);
	DftiCommitDescriptor(backward_dfti3);
#endif
}
void Dfti_finalize() {

#if (defined (_OPENMP) && !defined(_MPI))
#ifdef _FFT_IMKL
	DftiFreeDescriptor(&forward_dfti);
	DftiFreeDescriptor(&backward_dfti);
#endif
#elif defined(_MPI)
	DftiFreeDescriptor(&forward_dfti1);
	DftiFreeDescriptor(&forward_dfti2);
	DftiFreeDescriptor(&forward_dfti3);
	DftiFreeDescriptor(&backward_dfti1);
	DftiFreeDescriptor(&backward_dfti2);
	DftiFreeDescriptor(&backward_dfti3);
	free_1d_double(sbuf_free);
	free_1d_double(rbuf_free);
	free_1d_double(fft_free);
	free_1d_int(npx_all);
	free_1d_int(npy_all);
	free_1d_int(nqy_all);
	free_1d_int(nqz_all);
	free_1d_int(npx_proc);
	free_1d_int(npy_proc);
	free_1d_int(nqy_proc);
	free_1d_int(nqz_proc);
#endif
}
//#endif

void A2a_k_1D(double *a) {
#ifndef _MPI
#ifdef _FFT_IMKL
	long status = DftiComputeForward(imkl_p_fw, a);
#elif _FFT_FFTW
	fftw_execute_dft_r2c(fftw_p_fw, a, reinterpret_cast<fftw_complex *>(a));
#elif _FFT_OOURA
	initview_3d_double(NX, NY, NZ_, a, ooura_p.a);
	rdft3d(NX, NY, NZ, 1, ooura_p.a, ooura_p.t, ooura_p.ip, ooura_p.w);
	rdft3dsort(NX, NY, NZ, 1, ooura_p.a);
	Complex *ak = reinterpret_cast<Complex *>(a);
#pragma omp parallel for
	for (int i = 0; i < NX * NY * HNZ_; i++)
		ak[i] = std::conj(ak[i]);
#endif
#else

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(NPs, NZ_, a, fft)
#endif

	for (int i = 0; i < NPs[REAL][0]; i++) {
		for (int j = 0; j < NPs[REAL][1]; j++) {
			for (int k = 0; k < NPs[REAL][2]; k++) {
				fft[REALMODE_ARRAYINDEX(i, j, k)] = a[REALMODE_ARRAYINDEX(i, j, k)];
			}
		}
	}

	// A2a_k First FFT(NZ)
	status = DftiComputeForward(forward_dfti1, fft);
	DFTI_ERRORCHECK(status);
	// packing data a to sbuf
	if (yprocs > 1) {
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(max_npx, max_nqy, max_nqz_, npx, nqy, nqz_all, nqz_proc, NZ_,       \
           yprocs, fft, sbuf)
#endif
		for (int i = 0; i < npx; i++) {
			for (int j = 0; j < nqy; j++) {
				for (int yp = 0; yp < yprocs; yp++) {
					for (int k = 0; k < nqz_all[yp]; k++) {
						sbuf[(yp * max_npx * max_nqy * max_nqz_) +
							(i * max_nqy * max_nqz_) + (j * max_nqz_) + k] =
							fft[(i * nqy * NZ_) + (j * NZ_) + (k + nqz_proc[yp])];
					}
				}
			}
		}
		// Transpose

		ierr = MPI_Alltoall(sbuf, cnt1, MPI_DOUBLE, rbuf, cnt1, MPI_DOUBLE,
			OWN_Y_COMM);

		MPI_ERRORCHECK(ierr);
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(max_npx, max_nqy, max_nqz_, npx, NY, nqy_all, nqy_proc, nqz_,       \
           yprocs, fft, rbuf)
#endif
		for (int i = 0; i < npx; i++) {
			for (int yp = 0; yp < yprocs; yp++) {
				for (int j = 0; j < nqy_all[yp]; j++) {
					for (int k = 0; k < nqz_; k++) {
						fft[(i * NY * nqz_) + ((j + nqy_proc[yp]) * nqz_) + k] =
							rbuf[(yp * max_npx * max_nqy * max_nqz_) +
							(i * max_nqy * max_nqz_) + (j * max_nqz_) + k];
					}
				}
			}
		}
	}
	// A2a_k Second FFT(NY)
	for (int i = 0; i < npx; i++) {

		status = DftiComputeForward(forward_dfti2, (fft + (i * NY * nqz_)));

		DFTI_ERRORCHECK(status);
	}
	// packing data a to sbuf
	if (xprocs > 1) {
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(max_npx, max_npy, max_nqz_, npx, npy_all, npy_proc, nqz_, NY,       \
           xprocs, fft, sbuf)
#endif
		for (int i = 0; i < npx; i++) {
			for (int xp = 0; xp < xprocs; xp++) {
				for (int j = 0; j < npy_all[xp]; j++) {
					for (int k = 0; k < nqz_; k++) {
						sbuf[(xp * max_npx * max_npy * max_nqz_) +
							(i * max_npy * max_nqz_) + (j * max_nqz_) + k] =
							fft[(i * NY * nqz_) + ((j + npy_proc[xp]) * nqz_) + k];
					}
				}
			}
		}
		// Transpose

		ierr = MPI_Alltoall(sbuf, cnt2, MPI_DOUBLE, rbuf, cnt2, MPI_DOUBLE, OWN_X_COMM);

		MPI_ERRORCHECK(ierr);
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(max_npx, max_npy, max_nqz_, npx_all, npx_proc, npy, nqz_, xprocs,   \
           fft, rbuf)
#endif
		for (int xp = 0; xp < xprocs; xp++) {
			for (int i = 0; i < npx_all[xp]; i++) {
				for (int j = 0; j < npy; j++) {
					for (int k = 0; k < nqz_; k++) {
						fft[((i + npx_proc[xp]) * npy * nqz_) + (j * nqz_) + k] =
							rbuf[(xp * max_npx * max_npy * max_nqz_) +
							(i * max_npy * max_nqz_) + (j * max_nqz_) + k];
					}
				}
			}
		}
	}
	// A2a_k Third FFT(NX)

	status = DftiComputeForward(forward_dfti3, fft);

	DFTI_ERRORCHECK(status);
	/*
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(NPs, HNQZ_, a, fft)
#endif
	for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
		for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
			for (int k = 0; k < HNQZ_; k++) {
				a[SPECTRUMMODE_ARRAYINDEX(i, j, (2 * k))] = fft[SPECTRUMMODE_ARRAYINDEX(i, j, (2 * k))];
				a[SPECTRUMMODE_ARRAYINDEX(i, j, ((2 * k) + 1))] = -1.0 * fft[SPECTRUMMODE_ARRAYINDEX(i, j, ((2 * k) + 1))];
			}
		}
	}
	*/
	for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
		for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
			for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
				a[SPECTRUMMODE_ARRAYINDEX(i, j, k)] = fft[SPECTRUMMODE_ARRAYINDEX(i, j, k)];
			}
		}
	}
#endif
}
void A2a_k_nD(double **a, int dim) {
#ifndef _MPI
	for (int d = 0; d < dim; d++) {
		A2a_k_1D(a[d]);
	}
#else
	const int cnt1_dim = dim * cnt1;
	const int cnt2_dim = dim * cnt2;
	double *fft_dim, *sbuf_dim, *rbuf_dim;
	double *fft_dim_free, *sbuf_dim_free, *rbuf_dim_free;
	calloc_1d_double_alignment(&fft_dim_free, &fft_dim, dim * mesh_size, 16);
	calloc_1d_double_alignment(&sbuf_dim_free, &sbuf_dim, dim * size, 16);
	calloc_1d_double_alignment(&rbuf_dim_free, &rbuf_dim, dim * size, 16);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(dim, NPs, NZ_, a, fft_dim)
#endif
	for (int i = 0; i < NPs[REAL][0]; i++) {
		for (int d = 0; d < dim; d++) {
			for (int j = 0; j < NPs[REAL][1]; j++) {
				for (int k = 0; k < NPs[REAL][2]; k++) {
					fft_dim[REALMODE_ARRAYINDEX_DIM(d, i, j, k)] =
						a[d][REALMODE_ARRAYINDEX(i, j, k)];
				}
			}
		}
	}
	// A2a_k First FFT(NZ)
	for (int d = 0; d < dim; d++) {

		status =
			DftiComputeForward(forward_dfti1, (fft_dim + (d * npx * nqy * NZ_)));

		DFTI_ERRORCHECK(status);
	}
	// packing data a to sbuf_dim
	if (yprocs > 1) {
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(dim, max_npx, max_nqy, max_nqz_, npx, nqy, nqz_all, nqz_proc, NZ_,  \
           yprocs, fft_dim, sbuf_dim)
#endif
		for (int i = 0; i < npx; i++) {
			for (int d = 0; d < dim; d++) {
				for (int j = 0; j < nqy; j++) {
					for (int yp = 0; yp < yprocs; yp++) {
						for (int k = 0; k < nqz_all[yp]; k++) {
							sbuf_dim[(yp * dim * max_npx * max_nqy * max_nqz_) +
								(d * max_npx * max_nqy * max_nqz_) +
								(i * max_nqy * max_nqz_) + (j * max_nqz_) + k] =
								fft_dim[(d * npx * nqy * NZ_) + (i * nqy * NZ_) + (j * NZ_) +
								(k + nqz_proc[yp])];
						}
					}
				}
			}
		}
		// Transpose

		ierr = MPI_Alltoall(sbuf_dim, cnt1_dim, MPI_DOUBLE, rbuf_dim, cnt1_dim,
			MPI_DOUBLE, OWN_Y_COMM);

		MPI_ERRORCHECK(ierr);
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(dim, max_npx, max_nqy, max_nqz_, npx, nqy_all, nqy_proc, nqz_, NY,  \
           yprocs, fft_dim, rbuf_dim)
#endif
		for (int i = 0; i < npx; i++) {
			for (int d = 0; d < dim; d++) {
				for (int yp = 0; yp < yprocs; yp++) {
					for (int j = 0; j < nqy_all[yp]; j++) {
						for (int k = 0; k < nqz_; k++) {
							fft_dim[(d * npx * NY * nqz_) + (i * NY * nqz_) +
								((j + nqy_proc[yp]) * nqz_) + k] =
								rbuf_dim[(yp * dim * max_npx * max_nqy * max_nqz_) +
								(d * max_npx * max_nqy * max_nqz_) +
								(i * max_nqy * max_nqz_) + (j * max_nqz_) + k];
						}
					}
				}
			}
		}
	}
	// A2a_k Second FFT(NY)
	for (int d = 0; d < dim; d++) {
		for (int i = 0; i < npx; i++) {

			status = DftiComputeForward(
				forward_dfti2, (fft_dim + (d * npx * NY * nqz_) + (i * NY * nqz_)));

			DFTI_ERRORCHECK(status);
		}
	}
	// packing data a to sbuf_dim
	if (xprocs > 1) {
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(dim, max_npx, max_npy, max_nqz_, npx, npy_all, npy_proc, nqz_, NY,  \
           xprocs, fft_dim, sbuf_dim)
#endif
		for (int i = 0; i < npx; i++) {
			for (int d = 0; d < dim; d++) {
				for (int xp = 0; xp < xprocs; xp++) {
					for (int j = 0; j < npy_all[xp]; j++) {
						for (int k = 0; k < nqz_; k++) {
							sbuf_dim[(xp * dim * max_npx * max_npy * max_nqz_) +
								(d * max_npx * max_npy * max_nqz_) +
								(i * max_npy * max_nqz_) + (j * max_nqz_) + k] =
								fft_dim[(d * npx * NY * nqz_) + (i * NY * nqz_) +
								((j + npy_proc[xp]) * nqz_) + k];
						}
					}
				}
			}
		}
		// Transpose

		ierr = MPI_Alltoall(sbuf_dim, cnt2_dim, MPI_DOUBLE, rbuf_dim, cnt2_dim,
			MPI_DOUBLE, OWN_X_COMM);

		MPI_ERRORCHECK(ierr);
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(dim, max_npx, max_npy, max_nqz_, npx_all, npx_proc, npy, nqz_, NX,  \
           xprocs, fft_dim, rbuf_dim)
#endif
		for (int xp = 0; xp < xprocs; xp++) {
			for (int d = 0; d < dim; d++) {
				for (int i = 0; i < npx_all[xp]; i++) {
					for (int j = 0; j < npy; j++) {
						for (int k = 0; k < nqz_; k++) {
							fft_dim[(d * NX * npy * nqz_) +
								((i + npx_proc[xp]) * npy * nqz_) + (j * nqz_) + k] =
								rbuf_dim[(xp * dim * max_npx * max_npy * max_nqz_) +
								(d * max_npx * max_npy * max_nqz_) +
								(i * max_npy * max_nqz_) + (j * max_nqz_) + k];
						}
					}
				}
			}
		}
	}
	// A2a_k Third FFT(NX)
	for (int d = 0; d < dim; d++) {

		status =
			DftiComputeForward(forward_dfti3, (fft_dim + (d * NX * npy * nqz_)));

		DFTI_ERRORCHECK(status);
	}
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(dim, NPs, HNQZ_, a, fft_dim)
#endif
	/*
	for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
		for (int d = 0; d < dim; d++) {
			for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
				for (int k = 0; k < HNQZ_; k++) {
					a[d][SPECTRUMMODE_ARRAYINDEX(i, j, (2 * k))] =
						fft_dim[SPECTRUMMODE_ARRAYINDEX_DIM(d, i, j, (2 * k))];
					a[d][SPECTRUMMODE_ARRAYINDEX(i, j, ((2 * k) + 1))] =
						-1.0 *
						fft_dim[SPECTRUMMODE_ARRAYINDEX_DIM(d, i, j, ((2 * k) + 1))];
				}
			}
		}
	}
	*/
	for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
		for (int d = 0; d < dim; d++) {
			for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
				for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
					a[d][SPECTRUMMODE_ARRAYINDEX(i, j, k)] = fft_dim[SPECTRUMMODE_ARRAYINDEX_DIM(d, i, j, k)];
				}
			}
		}
	}
	free_1d_double(sbuf_dim_free);
	free_1d_double(rbuf_dim_free);
	free_1d_double(fft_dim_free);
#endif
}
void A_k2a_1D(double *a) {
#ifndef _MPI
#ifdef _FFT_IMKL
	long status = DftiComputeBackward(imkl_p_bw, a);
#elif _FFT_FFTW
	static const double scale = 1.0 / (NX * NY * NZ);
	Complex *ak = reinterpret_cast<Complex *>(a);
#pragma omp parallel for
	for (int i = 0; i < NX * NY * HNZ_; i++)
		ak[i] *= scale;
	fftw_execute_dft_c2r(fftw_p_bw, reinterpret_cast<fftw_complex *>(a), a);
#elif _FFT_OOURA
	static const double scale = 2.0 / (NX * NY * NZ);
	Complex *ak = reinterpret_cast<Complex *>(a);
#pragma omp parallel for
	for (int i = 0; i < NX * NY * HNZ_; i++)
		ak[i] = scale * std::conj(ak[i]);

	initview_3d_double(NX, NY, NZ_, a, ooura_p.a);
	rdft3dsort(NX, NY, NZ, -1, ooura_p.a);
	rdft3d(NX, NY, NZ, -1, ooura_p.a, ooura_p.t, ooura_p.ip, ooura_p.w);
#endif
#else
	double scale = 1.0 / (NX * NY * NZ);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(NPs, HNQZ_, a, fft)
#endif
	for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
		for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
			for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
				fft[SPECTRUMMODE_ARRAYINDEX(i, j, k)] = a[SPECTRUMMODE_ARRAYINDEX(i, j, k)];
			}
		}
	}
	// A2a_k First FFT(NX)

	status = DftiComputeBackward(backward_dfti1, fft);

	DFTI_ERRORCHECK(status);
	// packing data a to sbuf
	if (xprocs > 1) {
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(max_npx, max_npy, max_nqz_, npx_all, npx_proc, npy, nqz_, xprocs,   \
           fft, sbuf)
#endif
		for (int xp = 0; xp < xprocs; xp++) {
			for (int i = 0; i < npx_all[xp]; i++) {
				for (int j = 0; j < npy; j++) {
					for (int k = 0; k < nqz_; k++) {
						sbuf[(xp * max_npx * max_npy * max_nqz_) +
							(i * max_npy * max_nqz_) + (j * max_nqz_) + k] =
							fft[((i + npx_proc[xp]) * npy * nqz_) + (j * nqz_) + k];
					}
				}
			}
		}
		// Transpose

		ierr = MPI_Alltoall(sbuf, cnt2, MPI_DOUBLE, rbuf, cnt2, MPI_DOUBLE,
			OWN_X_COMM);

		MPI_ERRORCHECK(ierr);
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(max_npx, max_npy, max_nqz_, npx, npy_all, npy_proc, nqz_, NY,       \
           xprocs, fft, rbuf)
#endif
		for (int i = 0; i < npx; i++) {
			for (int xp = 0; xp < xprocs; xp++) {
				for (int j = 0; j < npy_all[xp]; j++) {
					for (int k = 0; k < nqz_; k++) {
						fft[(i * NY * nqz_) + ((j + npy_proc[xp]) * nqz_) + k] =
							rbuf[(xp * max_npx * max_npy * max_nqz_) +
							(i * max_npy * max_nqz_) + (j * max_nqz_) + k];
					}
				}
			}
		}
	}
	// A_k2a Second FFT(NY)
	for (int i = 0; i < npx; i++) {

		status = DftiComputeBackward(backward_dfti2, (fft + (i * NY * nqz_)));

		DFTI_ERRORCHECK(status);
	}
	// packing data a to sbuf
	if (yprocs > 1) {
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(max_npx, max_nqy, max_nqz_, npx, nqy_all, nqy_proc, nqz_, NY,       \
           yprocs, fft, sbuf)
#endif
		for (int i = 0; i < npx; i++) {
			for (int yp = 0; yp < yprocs; yp++) {
				for (int j = 0; j < nqy_all[yp]; j++) {
					for (int k = 0; k < nqz_; k++) {
						sbuf[(yp * max_npx * max_nqy * max_nqz_) +
							(i * max_nqy * max_nqz_) + (j * max_nqz_) + k] =
							fft[(i * NY * nqz_) + ((j + nqy_proc[yp]) * nqz_) + k];
					}
				}
			}
		}
		// Transpose

		ierr = MPI_Alltoall(sbuf, cnt1, MPI_DOUBLE, rbuf, cnt1, MPI_DOUBLE,
			OWN_Y_COMM);

		MPI_ERRORCHECK(ierr);
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(max_npx, max_nqy, max_nqz_, npx, nqy, nqz_all, nqz_proc, NZ_,       \
           yprocs, fft, rbuf)
#endif
		for (int i = 0; i < npx; i++) {
			for (int j = 0; j < nqy; j++) {
				for (int yp = 0; yp < yprocs; yp++) {
					for (int k = 0; k < nqz_all[yp]; k++) {
						fft[(i * nqy * NZ_) + (j * NZ_) + (k + nqz_proc[yp])] =
							rbuf[(yp * max_npx * max_nqy * max_nqz_) +
							(i * max_nqy * max_nqz_) + (j * max_nqz_) + k];
					}
				}
			}
		}
	}
	// A2a_k First FFT(NZ)

	status = DftiComputeBackward(backward_dfti3, fft);

	DFTI_ERRORCHECK(status);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(NPs, NZ_, a, fft, scale)
#endif
	for (int i = 0; i < NPs[REAL][0]; i++) {
		for (int j = 0; j < NPs[REAL][1]; j++) {
			for (int k = 0; k < NPs[REAL][2]; k++) {
				a[REALMODE_ARRAYINDEX(i, j, k)] =
					fft[REALMODE_ARRAYINDEX(i, j, k)] * scale;
			}
		}
	}
#endif
}
void A_k2a_nD(double **a, int dim) {
#if !defined(_MPI)
	for (int d = 0; d < dim; d++) {
		A_k2a_1D(a[d]);
	}
#else
	const int cnt1_dim = dim * cnt1;
	const int cnt2_dim = dim * cnt2;
	double scale = 1.0 / (NX * NY * NZ);
	double *fft_dim, *sbuf_dim, *rbuf_dim;
	double *fft_dim_free, *sbuf_dim_free, *rbuf_dim_free;
	calloc_1d_double_alignment(&fft_dim_free, &fft_dim, dim * mesh_size, 16);
	calloc_1d_double_alignment(&sbuf_dim_free, &sbuf_dim, dim * size, 16);
	calloc_1d_double_alignment(&rbuf_dim_free, &rbuf_dim, dim * size, 16);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(dim, NPs, HNQZ_, a, fft_dim)
#endif
	/*
	for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
		for (int d = 0; d < dim; d++) {
			for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
				for (int k = 0; k < HNQZ_; k++) {
					fft_dim[SPECTRUMMODE_ARRAYINDEX_DIM(d, i, j, (2 * k))] =
						a[d][SPECTRUMMODE_ARRAYINDEX(i, j, (2 * k))];
					fft_dim[SPECTRUMMODE_ARRAYINDEX_DIM(d, i, j, (2 * k) + 1)] =
						-1.0 * a[d][SPECTRUMMODE_ARRAYINDEX(i, j, ((2 * k) + 1))];
				}
			}
		}
	}
	*/
	for (int i = 0; i < NPs[SPECTRUM][0]; i++) {
		for (int d = 0; d < dim; d++) {
			for (int j = 0; j < NPs[SPECTRUM][1]; j++) {
				for (int k = 0; k < NPs[SPECTRUM][2]; k++) {
					fft_dim[SPECTRUMMODE_ARRAYINDEX_DIM(d, i, j, k)] = a[d][SPECTRUMMODE_ARRAYINDEX(i, j, k)];
				}
			}
		}
	}
	// A2a_k First FFT(NX)
	for (int d = 0; d < dim; d++) {

		status =
			DftiComputeBackward(backward_dfti1, (fft_dim + (d * NX * npy * nqz_)));

		DFTI_ERRORCHECK(status);
	}
	// packing data a to sbuf_dim
	if (xprocs > 1) {
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(dim, max_npx, max_npy, max_nqz_, npx_all, npx_proc, npy, nqz_, NX,  \
           xprocs, fft_dim, sbuf_dim)
#endif
		for (int xp = 0; xp < xprocs; xp++) {
			for (int d = 0; d < dim; d++) {
				for (int i = 0; i < npx_all[xp]; i++) {
					for (int j = 0; j < npy; j++) {
						for (int k = 0; k < nqz_; k++) {
							sbuf_dim[(xp * dim * max_npx * max_npy * max_nqz_) +
								(d * max_npx * max_npy * max_nqz_) +
								(i * max_npy * max_nqz_) + (j * max_nqz_) + k] =
								fft_dim[(d * NX * npy * nqz_) +
								((i + npx_proc[xp]) * npy * nqz_) + (j * nqz_) + k];
						}
					}
				}
			}
		}
		// Transpose

		ierr = MPI_Alltoall(sbuf_dim, cnt2_dim, MPI_DOUBLE, rbuf_dim, cnt2_dim,
			MPI_DOUBLE, OWN_X_COMM);

		MPI_ERRORCHECK(ierr);
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(dim, max_npx, max_npy, max_nqz_, npx, npy_all, npy_proc, nqz_, NY,  \
           xprocs, fft_dim, rbuf_dim)
#endif
		for (int i = 0; i < npx; i++) {
			for (int d = 0; d < dim; d++) {
				for (int xp = 0; xp < xprocs; xp++) {
					for (int j = 0; j < npy_all[xp]; j++) {
						for (int k = 0; k < nqz_; k++) {
							fft_dim[(d * npx * NY * nqz_) + (i * NY * nqz_) +
								((j + npy_proc[xp]) * nqz_) + k] =
								rbuf_dim[(xp * dim * max_npx * max_npy * max_nqz_) +
								(d * max_npx * max_npy * max_nqz_) +
								(i * max_npy * max_nqz_) + (j * max_nqz_) + k];
						}
					}
				}
			}
		}
	}
	// A_k2a Second FFT(NY)
	for (int d = 0; d < dim; d++) {
		for (int i = 0; i < npx; i++) {

			status = DftiComputeBackward(
				backward_dfti2, (fft_dim + (d * npx * NY * nqz_) + (i * NY * nqz_)));

			DFTI_ERRORCHECK(status);
		}
	}
	// packing data a to sbuf_dim
	if (yprocs > 1) {
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(dim, max_npx, max_nqy, max_nqz_, npx, nqy_all, nqy_proc, nqz_, NY,  \
           yprocs, fft_dim, sbuf_dim)
#endif
		for (int i = 0; i < npx; i++) {
			for (int d = 0; d < dim; d++) {
				for (int yp = 0; yp < yprocs; yp++) {
					for (int j = 0; j < nqy_all[yp]; j++) {
						for (int k = 0; k < nqz_; k++) {
							sbuf_dim[(yp * dim * max_npx * max_nqy * max_nqz_) +
								(d * max_npx * max_nqy * max_nqz_) +
								(i * max_nqy * max_nqz_) + (j * max_nqz_) + k] =
								fft_dim[(d * npx * NY * nqz_) + (i * NY * nqz_) +
								((j + nqy_proc[yp]) * nqz_) + k];
						}
					}
				}
			}
		}
		// Transpose

		ierr = MPI_Alltoall(sbuf_dim, cnt1_dim, MPI_DOUBLE, rbuf_dim, cnt1_dim,
			MPI_DOUBLE, OWN_Y_COMM);

		MPI_ERRORCHECK(ierr);
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(dim, max_npx, max_nqy, max_nqz_, npx, nqy, nqz_all, nqz_proc, NZ_,  \
           yprocs, fft_dim, rbuf_dim)
#endif
		for (int i = 0; i < npx; i++) {
			for (int d = 0; d < dim; d++) {
				for (int j = 0; j < nqy; j++) {
					for (int yp = 0; yp < yprocs; yp++) {
						for (int k = 0; k < nqz_all[yp]; k++) {
							fft_dim[(d * npx * nqy * NZ_) + (i * nqy * NZ_) + (j * NZ_) +
								(k + nqz_proc[yp])] =
								rbuf_dim[(yp * dim * max_npx * max_nqy * max_nqz_) +
								(d * max_npx * max_nqy * max_nqz_) +
								(i * max_nqy * max_nqz_) + (j * max_nqz_) + k];
						}
					}
				}
			}
		}
	}
	// A2a_k First FFT(NZ)
	for (int d = 0; d < dim; d++) {

		status =
			DftiComputeBackward(backward_dfti3, (fft_dim + (d * npx * nqy * NZ_)));

		DFTI_ERRORCHECK(status);
	}
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(dim, NPs, NZ_, a, fft_dim, scale)
#endif
	for (int i = 0; i < NPs[REAL][0]; i++) {
		for (int d = 0; d < dim; d++) {
			for (int j = 0; j < NPs[REAL][1]; j++) {
				for (int k = 0; k < NPs[REAL][2]; k++) {
					a[d][REALMODE_ARRAYINDEX(i, j, k)] =
						fft_dim[REALMODE_ARRAYINDEX_DIM(d, i, j, k)] * scale;
				}
			}
		}
}
	free_1d_double(sbuf_dim_free);
	free_1d_double(rbuf_dim_free);
	free_1d_double(fft_dim_free);
#endif
}
