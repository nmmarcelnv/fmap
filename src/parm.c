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
