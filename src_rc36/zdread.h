/* cc create_lig.c -lm -o create_lig -O3 */
/* Rong Chen 1/15/2002 */
#ifndef _ZDREAD_H_
#define _ZDREAD_H_

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "pdb.h"

#ifdef __cplusplus
extern "C" {
#endif
	
typedef struct
{
        int N;
        double spacing;
        double rand[3];
        double r[3];
        double l[3];
        char rec[MAXLENLINE];
        char lig[MAXLENLINE];
} ZDCOM;

int ZDread(char *zdockfile, ZDCOM *zdcom, double angl[][3], double tran[][3], double zdscore[], bool load);

#ifdef __cplusplus
}
#endif

#endif  /* _ZDREAD_H_ */
