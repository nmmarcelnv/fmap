#include "common.h"

double volmul(void* data1, void* data2, double offset[3],double dx2, void* sysv){
        ATOM* atm1=data1;
        ATOM* atm2=data2;
        double x=atm1->xyz[0]-atm2->xyz[0]-offset[0];
        double y=atm1->xyz[1]-atm2->xyz[1]-offset[1];
        double z=atm1->xyz[2]-atm2->xyz[2]-offset[2];
        double d2=x*x+y*y+z*z;
        double r=atm1->r+atm2->r;
        double r2=r*r;
        if (d2>r2){
                return 0.0;
        }else{
                return 1.0;
        }
}

double vollnk(void* data1, void* data2, double d2,double dx2, void* sysv){
        ATOM* atm1=data1;
        ATOM* atm2=data2;
        double r=atm1->r+atm2->r;
        double r2=r*r;
        if (d2>r2){
                return 0.0;
        }else{
                return 1.0;
        }
}
