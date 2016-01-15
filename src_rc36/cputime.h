#ifndef _CPUTIME_H_
#define _CPUTIEM_H_

#include <time.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif

void Times(bool start, int n,...);
void TimeRep(FILE *fp);

#ifdef __cplusplus
}
#endif

#endif  /* _CPUTIME_H_ */
