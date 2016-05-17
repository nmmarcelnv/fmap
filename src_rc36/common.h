#ifndef _COMMON_H_
#define _COMMON_H_

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <assert.h>

#ifdef _OPENMP
	#include <omp.h>
#endif /* _OPENMP */
#ifdef USE_MPI
        #include <mpi.h>
#endif /* USE_MPI */

#include "pdb.h"
#include "pro.h"
#include "rdf.h"
#include "rot.h"
#include "parm.h"
#include "fftw.h"
#include "grd.h"
#include "cross.h"
#include "score.h"
#include "cputime.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

#endif  /* _COMMON_H_ */
