#include "common.h"
#include "linked.h"

double drtlnkf(int n, ATOM atoms1[], int m, ATOM atoms2[], double offset[3],double dx2,void* sys, CLINKED *lnk,
	double (*mul)(void*,void*, double, double, void*)){
	int i,j;
	double sum=0.0;
	double mx,my,mz;
	int idx[n];
	double r2[n];
	double cur;
	for (j=0;j<m;j++){
		mx=offset[0]+atoms2[j].xyz[0];
		my=offset[1]+atoms2[j].xyz[1];
		mz=offset[2]+atoms2[j].xyz[2];
		n=lnk_nearest(lnk,mx,my,mz,idx,r2);
		//if (n>0){
		//	fprintf(stderr,"lnk: %d offset: %8.3f%8.3f%8.3f %8.3f%8.3f%8.3f\n",n,offset[0],offset[1],offset[2],mx,my,mz);
		//}
		for (i=0;i<n;i++){
			//printf("AA: %16e\n",mul(&atoms1[i],&atoms2[j],offset));
			//fprintf(stderr,"r2 %d %f\n",idx[i],sqrt(r2[i]));
			cur=mul(&atoms1[idx[i]],&atoms2[j],r2[i],dx2,sys);
			sum+=cur;
			//atoms1[idx[i]].ern+=cur;
			//atoms2[j].ern+=cur;
		}
	}
	//printf("pp: %16e\n",sum);
	return sum;
}

void drtlnk(int n, ATOM atoms1[], int m, ATOM atoms2[], int l,double dx, fftw_complex* grd0, int loc,void* sys,CLINKED *lnk,
	double (*mul)(void*,void*,double, double ,void*)){
	double dx2=1.0;
	double offset[3];
	fftw_complex (*grd)[l][l]=(void*)grd0;
	int ii,jj,kk;
	for (ii=0;ii<l;ii++){
		offset[0]=ii*dx;
                for (jj=0;jj<l;jj++){
			offset[1]=jj*dx;
                        for (kk=0;kk<l;kk++){
				offset[2]=kk*dx;
				double s=drtlnkf(n,atoms1,m,atoms2,offset,dx2,sys,lnk,mul);
				grd[ii][jj][kk][loc]=s;
				//printf("%d %d %d %16e\n",i,j,k,s);
			}
		}
	}
}

void drtlnkijk(int n, ATOM atoms1[], int m, ATOM atoms2[], int l,double dx, int ne, int ijk[ne][3], double ern[ne], 
	void* sys, CLINKED *lnk,
        double (*mul)(void*,void*,double,double, void*)){
        double dx2=1.0f;
        int i;
        for (i=0;i<ne;i++){
        	double offset[3];
        	int ii,jj,kk;
                ii=ijk[i][0];
                offset[0]=ii*dx;
                jj=ijk[i][1];
                offset[1]=jj*dx;
                kk=ijk[i][2];
                offset[2]=kk*dx;
                double s=drtlnkf(n,atoms1,m,atoms2,offset,dx2,sys,lnk,mul);
                ern[i]=s;
        }
}
