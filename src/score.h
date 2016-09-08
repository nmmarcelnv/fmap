#ifndef _SCORE_H_
#define _SCORE_H_

#include "common.h"

#ifdef __cplusplus
extern "C" {
#endif

void SetBool(const int l, bool vol[l][l][l], bool label);
void SetZero(const int l, fftw_real grd[l][l][2*(l/2+1)]);
void vRep(const int l,bool vol[l][l][l], double cutoff, fftw_real grd[l][l][2*(l/2+1)], const double scl, const char* str, const double kBT);
void AddTo1(const int l,fftw_real grd1[l][l][2*(l/2+1)], const double scl1, fftw_real grd2[l][l][2*(l/2+1)], const double scl2);
void SelErn(const int l,bool vol[l][l][l],fftw_real grd1[l][l][2*(l/2+1)],const double ang[3],const double erncut);
void BltzSumFilt(const int l,bool vol[l][l][l],fftw_real sav[l][l][2*(l/2+1)], fftw_real local[l][l][2*(l/2+1)], const double kBT);
void Bltz2Ern(const int l,fftw_real sav[l][l][2*(l/2+1)],const double kBT);
//spot.c
double GetSpotCut(int l,bool vol[l][l][l], fftw_real grd[l][l][2*(l/2+1)],int n,int ntop, const double kBT);
void VEInd(const int l, bool vol[l][l][l], fftw_real ele[l][l][2*(l/2+1)], fftw_real vdw[l][l][2*(l/2+1)], int nbin, double mat[nbin][nbin]);
#ifdef __cplusplus
}
#endif

#endif  /* _SCORE_H_ */
