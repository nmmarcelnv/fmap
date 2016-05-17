#ifndef _RDF_H_
#define _RDF_H_

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

double RdfR6m(const double r2, const double skip);
double RdfR12m(const double r2, const double skip);
double RdfDebye(const double r2, const double kap);
double RdfVol(const double r2, const double rad);

#ifdef __cplusplus
}
#endif

#endif  /* _RDF_H_ */
