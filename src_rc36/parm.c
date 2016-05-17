#include "parm.h"

      /* Compute parameters: 
       *
       * kappa^2 = (8 pi N_A e_c^2) I_s / (1000 eps_w k_B T)
       * kappa   = 0.325567 * I_s^{1/2}   angstroms^{-1}
       * deblen  = 1 / kappa
       *         = 3.071564378 * I_s^{1/2}   angstroms
       * \bar{kappa}^2 = eps_w * kappa^2 
       * zmagic  = (4 * pi * e_c^2) / (k_B T)   (we scale the diagonal later)
       *         = 7046.528838
       * From APBS/generic/vpbe.c
       */

double GetKappa(double IonicStrength, double SDie, double Temp){
          const double N_A = 6.022045000e+23;
          const double e_c = 4.803242384e-10;
          const double k_B = 1.380662000e-16;
          const double Pi = 3.14159265358979323846;
          return sqrt(IonicStrength*1.0e-16*(8.0*Pi*N_A*e_c*e_c)/(1000.0*SDie*k_B*Temp));
}

double GetkBT(double Temp){
	const double N_A = 6.022045000e+23;
	//SI units, 2010 CODATA value, J/K = m2·kg/(s2·K) in SI base units
	//const double k_B = 1.380648813e-23;
	//cal/K 	1 steam table calorie = 4.1868 J
	const double k_B = 3.297623030e-24;
	return N_A*k_B*Temp/1000.0;
}

//Malmberg, C. & Maryott, A. Dielectric constant of water from 0 to 100~C. J. RES. NAT. BUR. STAN., NIST, 1956 ,56 ,1
double GetWaterDie(double Temp){
	const double c2k=273.15;
        double t=Temp-c2k;
        //sdie=87.740-0.40008*t + 9.398*(1e-4)*t^2 -1.410*(1e-6)*t^3
	double sdie = 87.740-0.40008*t+9.398e-4*t*t+1.410e-6*t*t*t;
	return sdie;
}
