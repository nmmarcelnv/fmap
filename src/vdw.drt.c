#include "common.h"

double vdwmul(void* data1, void* data2, double offset[3],double dx2, void* sysv){
        ATOM* atm1=data1;
        ATOM* atm2=data2;
        PARM* sys=sysv;
        double x=atm1->xyz[0]-atm2->xyz[0]-offset[0];
        double y=atm1->xyz[1]-atm2->xyz[1]-offset[1];
        double z=atm1->xyz[2]-atm2->xyz[2]-offset[2];
        double d2=x*x+y*y+z*z;
        //dx2=(atm1->r+1.2)*(atm1->r+1.2);
        //if (d2<dx2) return 0.0;
        d2=(d2>dx2)?d2:dx2;
        double d6=d2*d2*d2;
        double d12=d6*d6;
	double Cut2=sys->rup*sys->rup;
        if (d2>Cut2){
                return 0.0f;
        }else{
                return atm1->Asq*atm2->Asq/d12-atm1->Bsq*atm2->Bsq/d6;
        }
}

double vdwlnk(void* data1, void* data2, double d2,double dx2,void* sysv){
        ATOM* atm1=data1;
        ATOM* atm2=data2;
        PARM* sys=sysv;
        //dx2=(atm1->r+1.2)*(atm1->r+1.2);
        //if (d2<dx2) return 0.0;
        d2=(d2>dx2)?d2:dx2;
        double d6=d2*d2*d2;
        double d12=d6*d6;
	double Cut2=sys->rup*sys->rup;
        if (d2>Cut2){
                return 0.0f;
        }else{
                return atm1->Asq*atm2->Asq/d12-atm1->Bsq*atm2->Bsq/d6;
        }
}
