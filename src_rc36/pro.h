#ifndef _PRO_H_
#define _PRO_H_

#include <malloc.h>
#include <math.h>
#include <assert.h>

#include "pdb.h"

#ifdef __cplusplus
extern "C" {
#endif

//struct protein;
struct protein{
        int n;
        int *xi,*yi,*zi;
        double *xf,*yf,*zf;
        double *q,*r,*Asq,*Bsq,*kap;
};
typedef struct protein PRO;
PRO* ProCreate(const int n, ATOM atoms[n]);
void ProUpdate(PRO* pro, const int n, double xyz[n][3]);
void ProFree(PRO* pro);

#ifdef __cplusplus
}
#endif

#endif  /* _PRO_H_ */
