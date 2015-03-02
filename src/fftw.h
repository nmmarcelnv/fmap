#ifndef _FFTW_H_
#define _FFTW_H_

#include <fftw3.h>
#include <stdbool.h>
#define fftw_real double
#define NUM_THREADS 1

#ifdef __cplusplus
extern "C" {
#endif

void fft3d_r2c(int L, int M, int N, fftw_real a[L][M][2*(N/2+1)]);
void fft3d_c2r(int L, int M, int N, fftw_real c[L][M][2*(N/2+1)]);
void fft3d_add_inplace(int L, int M, int N,fftw_real a[L][M][2*(N/2+1)], fftw_real b[L][M][2*(N/2+1)]);

#ifdef __cplusplus
}
#endif

#endif  /* _FFTW_H_ */
