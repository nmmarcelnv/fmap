#include "rdf.h"

double RdfR6m(const double r2, const double skip){
	double r6=r2*r2*r2;
        return -1.0/r6;
}

double RdfR12m(const double r2, const double skip){
	double r6=r2*r2*r2;
	double r12=r6*r6;
	return 1.0/r12;
}

double RdfDebye(const double r2, const double kap){
	double r=sqrt(r2);
	return exp(-1.0*r*kap)/r;
}

double RdfVol(const double r2, const double rad){
	double rad2=rad*rad;
	return (r2>rad2)?0.0:1.0;
}
