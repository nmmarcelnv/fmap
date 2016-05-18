#include "common.h"

double elemul(void* data1, void* data2, double offset[3],double dx2, void* sysv){
        ATOM* atm1=data1;
        ATOM* atm2=data2;
        PARM* sys=sysv;
        //showatom(atm1);
        //showatom(atm2);
        double x=atm1->xyz[0]-atm2->xyz[0]-offset[0];
        double y=atm1->xyz[1]-atm2->xyz[1]-offset[1];
        double z=atm1->xyz[2]-atm2->xyz[2]-offset[2];
        double d2=x*x+y*y+z*z;
        //printf("%16e%16e%16e%16e\n",d2,atm1->q,atm2->q,atm1->q*atm2->q*exp(-1.0*dd*kap)/dd);
        double Cut2=sys->rEup*sys->rEup;
        if (d2>Cut2){
                return 0.0f;
        }else{
                double dd=sqrt(d2>dx2?d2:dx2);
                return atm1->q*atm2->q*exp(-1.0*dd*sys->kap)/dd;
        }
}

double elelnk(void* data1, void* data2, double d2, double dx2, void* sysv){
        ATOM* atm1=data1;
        ATOM* atm2=data2;
       	PARM* sys=sysv;
        double Cut2=sys->rEup*sys->rEup;
        if (d2>Cut2){
                return 0.0f;
        }else{
                double dd=sqrt(d2>dx2?d2:dx2);
                return atm1->q*atm2->q*exp(-1.0*dd*sys->kap)/dd;
        }
}
