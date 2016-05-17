#ifndef _PDBNMR_H_
#define _PDBNMR_H_

#include "pdb.h"

#ifdef __cplusplus
extern "C" {
#endif


int MDLCountAtoms(char *pdbfn);
int MDLReadxyz(FILE *file, int n, double xyz[][3]);

#endif /*_PDBNMR_H_*/
