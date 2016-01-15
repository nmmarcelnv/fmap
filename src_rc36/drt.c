#include "common.h"

double drtf(int n, ATOM atoms1[], int m, ATOM atoms2[], double offset[3],double dx2,void* sys,
	double (*mul)(void*,void*,double [], double, void*)){
	int i,j;
	double sum=0.0;
	double cur;
	for (i=0;i<n;i++){
		for (j=0;j<m;j++){
			//printf("AA: %16e\n",mul(&atoms1[i],&atoms2[j],offset));
			
			cur=mul(&atoms1[i],&atoms2[j],offset,dx2,sys);
			sum+=cur;
			//atoms1[i].ern+=cur;
			//atoms2[j].ern+=cur;
		}
	}
	//printf("pp: %16e\n",sum);
	return sum;
}

void drt(int n, ATOM atoms1[], int m, ATOM atoms2[], int l,double dx, fftw_complex* grd0, int loc,void* sys,
	double (*mul)(void*,void*,double [],double, void*)){
	double dx2=1.0f;
	double offset[3];
	fftw_complex (*grd)[l][l]=(void*)grd0;
	int ii,jj,kk;
	for (ii=0;ii<l;ii++){
		offset[0]=ii*dx;
                for (jj=0;jj<l;jj++){
			offset[1]=jj*dx;
                        for (kk=0;kk<l;kk++){
				offset[2]=kk*dx;	
				double s=drtf(n,atoms1,m,atoms2,offset,dx2,sys,mul);
				grd[ii][jj][kk][loc]=s;
				//printf("%d %d %d %16e\n",i,j,k,s);
			}
		}
	}
}

void drtijk(int n, ATOM atoms1[], int m, ATOM atoms2[], int l,double dx, int ne, int ijk[ne][3], double ern[ne], void* sys,
        double (*mul)(void*,void*,double [],double, void*)){
	double dx2=1.0f;
	double offset[3];
	int i;
	int ii,jj,kk;
	for (i=0;i<ne;i++){
		ii=ijk[i][0];
		offset[0]=ii*dx;
		jj=ijk[i][1];
		offset[1]=jj*dx;
		kk=ijk[i][2];
		offset[2]=kk*dx;
		double s=drtf(n,atoms1,m,atoms2,offset,dx2,sys,mul);
		ern[i]=s;
	}
}
