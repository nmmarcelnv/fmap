#include "common.h"

void getminmax(int l,bool vol[l][l][l], fftw_real grd[l][l][2*(l/2+1)],fftw_real minmax[2]){
	int i,j,k;
	minmax[0]=FLT_MAX; //min
	minmax[1]=-FLT_MAX; //max
	for (i=0;i<l;i++){
        	for (j=0;j<l;j++){
                        for (k=0;k<l;k++){
                                if (vol[i][j][k]==false){
					if (grd[i][j][k]<minmax[0]){
						minmax[0]=grd[i][j][k];
					}
					if (grd[i][j][k]>minmax[1]){
						minmax[1]=grd[i][j][k];
					}
				}
			}
		}
	}
}

int Index(fftw_real minmax[2],int n,fftw_real dx,fftw_real v){
	//printf("dx:%16.8f\n",dx);
	return (int)((v-minmax[0])/dx);
}

int SpotRep(fftw_real minmax[2],int n,fftw_real dx,long cnts[n+1], fftw_real en[n+1], int ntop){
	int i;
	fftw_real val;
	fftw_real sum=0.0;
	long sumc=0;
	int itop=0;
	for (i=0;i<n+1;i++){
		val=minmax[0]+(0.5+i)*dx;
		sum+=en[i];
		sumc+=cnts[i];
		if (ntop>sumc){
			itop=i;
		}
		//printf("Spot: %16.8f%16ld%16ld%16.8e%16.8e\n",val,cnts[i],sumc,en[i],sum);
	}
	return itop;
}

double GetSpotCut(int l,bool vol[l][l][l], fftw_real grd[l][l][2*(l/2+1)],int n,int ntop, const double kBT){
	int i,j,k;
	int vi;
	fftw_real en[n+1];
	long cnts[n+1];
	fftw_real mm[2];
	getminmax(l,vol,grd,mm);
	double dx=(mm[1]-mm[0])/(double)(n);
	printf("mm:%16.8f%16.8f\n",mm[0],mm[1]);
	for (i=0;i<n+1;i++){
		cnts[i]=0;
		en[i]=0.0;
	}
	for (i=0;i<l;i++){
                for (j=0;j<l;j++){
                        for (k=0;k<l;k++){
                                if (vol[i][j][k]==false){
					vi=Index(mm,n,dx,grd[i][j][k]);
					cnts[vi]++;
					en[vi]+=myEXP(grd[i][j][k]/-kBT);

				}
			}
		}
	}
	int itop=SpotRep(mm,n,dx,cnts,en,ntop);
	/*
	int hl=l/2;
	int ii,jj,kk;
	for (i=0;i<l;i++){
		ii=(i+hl)%l;
                for (j=0;j<l;j++){
			jj=(j+hl)%l;
                        for (k=0;k<l;k++){
				kk=(k+hl)%l;
                                if (vol[i][j][k]==false){
					vi=Index(mm,n,dx,grd[i][j][k]);
					if (vi<=itop){
						printf("Vi:%8d%8d%8d%16.8f%8d\n",ii,jj,kk,grd[i][j][k],vi);
					}
				}
			}
		}
	}
	*/
	return en[itop];
}
