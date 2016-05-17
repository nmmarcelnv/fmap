#ifndef _ROT_H_
#define _ROT_H_

#include "pdb.h"
#include <math.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void Euler2Rot(double psi, double theta, double phi, double rot[9]);
void RotXYZ(int nAtom, ATOM atoms[], double xyznew[][3], double rot[9]);
void RotPro(int nAtom, ATOM atoms[], ATOM pro[], double rot[9]);
int CountAng(char *pdbfn);
int SetAng(char *pdbfn, double Angs[][3]);

#ifdef __cplusplus
}
#endif

#endif  /* _ROT_H_ */
