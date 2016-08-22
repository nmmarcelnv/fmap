#ifndef _SDAMMIO_H_
#define _SDAMMIO_H_

#include "pdb.h"

#ifdef __cplusplus
extern "C" {
#endif


int SDACountAtoms(char *pdbfn);
int SDAReadxyz(FILE *file, int n, double tr[][3], double xr[][3], double yr[][3]);

#endif /*_SDAMMIO_H_*/
