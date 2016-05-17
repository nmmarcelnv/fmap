#include<ctype.h>
#include<stdlib.h>
#include<stdio.h>
#include<stddef.h>
#include<float.h>
#include<string.h>
#include<math.h>
#include "rand.h"

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3set(int *inext, int *inextp, long ma[56]){
        long mj/*,mk*/;
        //int i,ii,k;

        if (++(*inext) == 56) *inext=1;
        if (++(*inextp) == 56) *inextp=1;
        mj=ma[*inext]-ma[*inextp];
        if (mj < MZ) mj += MBIG;
        ma[*inext]=mj;
        return mj*FAC;
}

double ran1(RANPARM *ranp){
	long mj,mk;
	int i,ii,k;
	
	if (ranp->id <0 ){
		mj=MSEED-(ranp->id < 0 ? -ranp->id : ranp->id);
		mj %= MBIG;
		ranp->ma[55]=mj;
		mk=1;
		for (i=1;i <= 54;i++) {
                        ii=(21*i) % 55;
                        ranp->ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ranp->ma[ii];
                }
                for (k=1;k <= 4;k++)
                        for (i=1;i <= 55;i++) {
                                ranp->ma[i] -= ranp->ma[1+(i+30) % 55];
                                if (ranp->ma[i] < MZ) ranp->ma[i] += MBIG;
                        }
                ranp->in=0;
                ranp->inp=31;
                ranp->id=1;
	}
	if (++ranp->in == 56) ranp->in=1;
        if (++ranp->inp == 56) ranp->inp=1;
        mj=ranp->ma[ranp->in]-ranp->ma[ranp->inp];
        if (mj < MZ) mj += MBIG;
        ranp->ma[ranp->in]=mj;
        return mj*FAC;
}

double ran3parm(long *idum,int *inext, int *inextp, long ma[56])
{
        long mj,mk;
        int i,ii,k;

	if (*idum < 0){
                mj=MSEED-(*idum < 0 ? -*idum : *idum);
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i <= 54;i++) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ma[ii];
                }
                for (k=1;k <= 4;k++)
                        for (i=1;i <= 55;i++) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                        }
                *inext=0;
                *inextp=31;
                *idum=1;
        }
        if (++(*inext) == 56) *inext=1;
        if (++(*inextp) == 56) *inextp=1;
        mj=ma[*inext]-ma[*inextp];
        if (mj < MZ) mj += MBIG;
        ma[*inext]=mj;
        return mj*FAC;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */
double ran3(long *idum)
{
        static int inext,inextp;
        static long ma[56];
        static int iff=0;
        long mj,mk;
        int i,ii,k;

        if (*idum < 0 || iff == 0) {
	//if (*idum < 0){
                iff=1;
                mj=MSEED-(*idum < 0 ? -*idum : *idum);
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i <= 54;i++) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ma[ii];
                }
                for (k=1;k <= 4;k++)
                        for (i=1;i <= 55;i++) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                        }
                inext=0;
                inextp=31;
                *idum=1;
        }
        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG;
        ma[inext]=mj;
        return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

int readconf(FILE *fp0,RANPARM *ranp){
	int i;
	if(feof(fp0)){
		return 0;
	}
	else {
		if (fscanf(fp0,"%10ld%10d%10d",&ranp->id,&ranp->in,&ranp->inp)!=3){
			return 0;
		}
		for (i=1; i<56; i++){
			if(fscanf(fp0,"%10ld",&ranp->ma[i])!=1){
			return 0;
			}
		}
		if (fscanf(fp0,"\n%[^\n]",(char *)&ranp->cf)!=1){
			return 0;
		}
		return 1;
	}
}

void writeconf(FILE *fp0, RANPARM *ranp){
	int i;
	fprintf(fp0,"%10ld%10d%10d",ranp->id,ranp->in,ranp->inp);
	for (i=1; i<56; i++){
		fprintf(fp0,"%10ld",ranp->ma[i]);
		if ((i+5)%10==0){
		fprintf(fp0,"\n");
		}
	}
	fprintf(fp0,"%s\n",ranp->cf);
}

void ransave(RANPARM *ranp1,RANPARM *ranp2){
	int i;
	ranp2->id=ranp1->id;
	ranp2->in=ranp1->in;
	ranp2->inp=ranp1->inp;
	for (i=0;i<56;i++){
		ranp2->ma[i]=ranp1->ma[i];
	}
}

#define TWOPI 6.2831853
void ranp2tr(RANPARM *ranp,double tr[6], double radius){
	double rad,costh,sinth,phi;
	//printf("%f\n",ran1(ranp));
        rad=radius*ran1(ranp);
        costh=-1.0+2.0*ran1(ranp);
        sinth=sqrt(1.0-costh*costh);
        phi=TWOPI*ran1(ranp);
        tr[0]=rad*sinth*cos(phi); //x0B
        tr[1]=rad*sinth*sin(phi); //y0B
        tr[2]=rad*costh; //z0B
        tr[3]=-1.0+2.0*ran1(ranp); //costhB
        tr[4]=TWOPI*ran1(ranp); //phiB 
        tr[5]=TWOPI*ran1(ranp)-TWOPI/2.0; //chiB
	//printf("conf: %8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n",rad,costh,phi,tr[3],tr[4],tr[5]);//costhB,phiB,chiB)
}

void ranp2cf(RANPARM *ranp,double cf[6], double radius){
	cf[0]=radius*ran1(ranp); //r
	cf[1]=-1.0+2.0*ran1(ranp); //costh
	cf[2]=TWOPI*ran1(ranp); //phi
	cf[3]=-1.0+2.0*ran1(ranp); //costhB
	cf[4]=TWOPI*ran1(ranp); //chiB 
	cf[5]=TWOPI*ran1(ranp)-TWOPI/2.0; //chiB
}
#undef TWOPI 

void transf(double costhB, double phiB, double chiB,double rmat[3][3]){
	double sinthB;
	double z[3],x1[3],y1[3];
	sinthB=sqrt(1.0-costhB*costhB);
        z[0]=sinthB*cos(phiB);
        z[1]=sinthB*sin(phiB);
        z[2]=costhB;
        x1[0]=z[2]/sqrt(z[0]*z[0]+z[2]*z[2]);
        x1[1]=0.0;
        x1[2]=-z[0]/sqrt(z[0]*z[0]+z[2]*z[2]);
        y1[0]=-z[0]*z[1]/sqrt(z[0]*z[0]+z[2]*z[2]);
        y1[1]=(z[0]*z[0]+z[2]*z[2])/sqrt(z[0]*z[0]+z[2]*z[2]);
        y1[2]=-z[1]*z[2]/sqrt(z[0]*z[0]+z[2]*z[2]);
        rmat[0][0]=cos(chiB)*x1[0]+sin(chiB)*y1[0];
        rmat[1][0]=cos(chiB)*x1[1]+sin(chiB)*y1[1];
        rmat[2][0]=cos(chiB)*x1[2]+sin(chiB)*y1[2];
        rmat[0][2]=z[0];
        rmat[1][2]=z[1];
        rmat[2][2]=z[2];
        rmat[0][1]=z[1]*rmat[2][0]-z[2]*rmat[1][0];
        rmat[1][1]=z[2]*rmat[0][0]-z[0]*rmat[2][0];
        rmat[2][1]=z[0]*rmat[1][0]-z[1]*rmat[0][0];
	/*
	printf("%8.3f%8.3f%8.3f\n",z[0],z[1],z[2]);
	printf("%8.3f%8.3f%8.3f\n",x1[0],x1[1],x1[2]);
	printf("%8.3f%8.3f%8.3f\n",y1[0],y1[1],y1[2]);
	printf("%8.3f%8.3f%8.3f\n",rmat[0][0],rmat[1][0],rmat[2][0]);
	printf("%8.3f%8.3f%8.3f\n",rmat[0][1],rmat[1][1],rmat[2][1]);
	printf("%8.3f%8.3f%8.3f\n",rmat[0][2],rmat[1][2],rmat[2][2]);
	*/
}

void genconf(double lxyz[][3], int nums, double tranmat[6],double ltrxyz[][3]){
	int i;
	double rmat[3][3];
	transf(tranmat[3],tranmat[4],tranmat[5],rmat);
	for (i=0;i<nums;i++){
		ltrxyz[i][0]=tranmat[0]+lxyz[i][0]*rmat[0][0]+lxyz[i][1]*rmat[0][1]+lxyz[i][2]*rmat[0][2];
		ltrxyz[i][1]=tranmat[1]+lxyz[i][0]*rmat[1][0]+lxyz[i][1]*rmat[1][1]+lxyz[i][2]*rmat[1][2];
		ltrxyz[i][2]=tranmat[2]+lxyz[i][0]*rmat[2][0]+lxyz[i][1]*rmat[2][1]+lxyz[i][2]*rmat[2][2];
	}
}


