#ifndef _DRT_H_
#define _DRT_H_

#include "common.h"

#ifdef __cplusplus
extern "C" {
#endif

double vollnk(void* data1, void* data2, double d2,double dx2, void* sysv);
double elelnk(void* data1, void* data2, double d2,double dx2, void* sysv);
double vdwlnk(void* data1, void* data2, double d2,double dx2, void* sysv);
double volmul(void* data1, void* data2, double offset[],double dx2, void* sysv);
double elemul(void* data1, void* data2, double offset[],double dx2, void* sysv);
double vdwmul(void* data1, void* data2, double offset[],double dx2, void* sysv);
double drtf(int n, ATOM atoms1[], int m, ATOM atoms2[], double offset[3],double dx2,void* sysv,
	double (*mul)(void*,void*,double [], double, void*));
double drtlnkf(int n, ATOM atoms1[], int m, ATOM atoms2[], double offset[3],double dx2,void* sysv, CLINKED *lnk,
	double (*mul)(void*,void*, double, double, void*));
void softlnkijk1(int nCrd, ATOM Crds[], int nPro, ATOM Pros[], int l, double dx, PARM* sys,CLINKED* lnk,
	int nv, double ijk[nv][3], double ern[4][nv], double ernCrd[4][nCrd], double ernPro[4][nPro], int pos);
void softijk1(int nCrd, ATOM Crds[], int nPro, ATOM Pros[], int l, double dx, PARM* sys,CLINKED* lnk,
	int nv, double ijk[nv][3], double ern[4][nv], double ernCrd[4][nCrd], double ernPro[4][nPro], int pos);
void ijkrep(int nCrd, ATOM Crds[], int nPro, ATOM Pros[], PARM* sys, int nv, double angs[nv][3], double ijk[nv][3], double score[nv], double ern[4][nv],double ernCrd[4][nCrd], double ernPro[4][nPro], FILE *fp);

#ifdef __cplusplus
}
#endif

#endif  /* _DRT_H_ */
