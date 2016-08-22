#ifndef _PDB_H_
#define _PDB_H_

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAXLENLINE 255
#define INFOLEN 31
#define INFOLEN_1 30

typedef struct {
	char info[INFOLEN];
        double xyz[3];
        double q;
        double r;
	double Asq; //sqrt(A)
	double Bsq; //sqrt(B)
	double kap; //kappa in debye huckle
	double ern; //to record ernegy in drt calculation
}ATOM;

int FileExist(char *pdbfn);
int CountAtoms(char *pdbfn);
int ReadPdb(char *pdbfn, int n, ATOM atoms[]);
int ReadPqr(char *pdbfn, int n, ATOM atoms[]);
void CalCtd(int nAtom, ATOM atoms[], double cen[3]);
void ToCtd(int nAtom, ATOM atoms[], double cen[3]);
void Unit2dx(int nAtom, ATOM atoms[],const double dx);
void ShowPqr(FILE *fp, int nAtom, ATOM atoms[]);
void ShowPqrdx(FILE *fp, int nAtom, ATOM atoms[], double xyz[][3], double dx);
void SclRad(int nAtom, ATOM atoms[],const double scl);
void SclChg(int nAtom, ATOM atoms[],const double scl);
void SetKap(int nAtom, ATOM atoms[],const double kap);

#ifdef __cplusplus
}
#endif

#endif  /* _PDB_H_ */
