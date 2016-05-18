#ifndef _PARAMETER_H_
#define _PARAMETER_H_

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	double rup;
	double rEup;
	double rlow;
	double dx;
	double sdie;
	double kap;
	double kBT;
	double escl;
	double vscl;
	double erncut;
}PARM;

double GetKappa(double IonicStrength, double SDie, double Temp);
double GetkBT(double Temp);
double GetWaterDie(double Temp);

#ifdef __cplusplus
}
#endif

#endif  /* _PARAMETER_H_ */
