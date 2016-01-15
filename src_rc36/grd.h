#ifndef _GRD_H_
#define _GRD_H_

#include "common.h"

#ifdef __cplusplus
extern "C" {
#endif

int PBC(const int x, const int mx);
void val2grd(const int n, const double v[n], const int px[n], const int py[n], const int pz[n], const int l, fftw_real* grid);
void pbc1(const int n, const int x[n], int px[n], const int lc, const int shift);
void RGrdR12(const int l, fftw_real* grid, PRO* pro, PARM sys, double (*rdf)(double,double));
void RGrdR6(const int l, fftw_real* grid, PRO* pro, PARM sys, double (*rdf)(double,double));
void RGrdEle(const int l, fftw_real* grid, PRO* pro, PARM sys, double (*rdf)(double,double));
void RGrdVol(const int l, fftw_real* grid, PRO* pro, PARM sys, double (*rdf)(double,double));
void LGrdR12(const int l, fftw_real* grid, PRO *pro);
void LGrdR6(const int l, fftw_real* grid, PRO *pro);
void LGrd2ndEle(const int l, fftw_real* grid, PRO *pro);

#ifdef __cplusplus
}
#endif

#endif  /* _GRD_H_ */
