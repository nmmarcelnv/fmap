#ifndef _CROSS_H_
#define _CROSS_H_

#include "common.h"

#ifdef __cplusplus
extern "C" {
#endif

void Cross(PRO *rec, PRO *lig, PARM sys, const int l,fftw_real sav[l][l][2*(l/2+1)], const double ang[3]);

#ifdef __cplusplus
}
#endif

#endif  /* _CROSS_H_ */
