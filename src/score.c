#include "score.h"

double myEXP(const double x){
        if (x>FLT_MAX_EXP){
                //fprintf(stderr,"%f > FLT_MAX_EXP: %d\n",x,FLT_MAX_EXP);
                return 0.0;
        }else{
                return exp(x);
        }
}

void showscorenet(const char* func, double sum, double frac, double u,double u1){
        fprintf(stderr,"%s\tSum: %16e Occ:%16.12f Ave:%16.8f Net:%16.8f\n",func,sum,frac,u,u1);
}

void toBool(int l, fftw_real *grid, bool *gridbool, fftw_real cutoff){
	fftw_real (*grd)[l][2*(l/2+1)]=(void*)grid;
	bool (*gbl)[l][l]=(void*)gridbool;
	int i,j,k;
	#pragma omp parallel for private(j,k)
	for (i=0;i<l;i++){
		for (j=0;j<l;j++){
			for (k=0;k<l;k++){
				if (grd[i][j][k]>cutoff){
					gbl[i][j][k]=true;
				}else{
					gbl[i][j][k]=false;
				}
			}
		}
	}
}

void SetBool(const int l, bool vol[l][l][l], bool label){
	int i,j,k;
	#pragma omp parallel for private(j,k)
	for (i=0;i<l;i++){
		for (j=0;j<l;j++){
			for (k=0;k<l;k++){
				vol[i][j][k]=label;
			}
		}
	}
}

void SetZero(const int l, fftw_real grd[l][l][2*(l/2+1)]){
	int i,j,k;
	#pragma omp parallel for private(j,k)
	for (i=0;i<l;i++){
                for (j=0;j<l;j++){
                        for (k=0;k<2*(l/2+1);k++){
                                grd[i][j][k]=0.0;
                        }
                }
        }
}

void zeroarr(int n, double *v){
        int i;
        for (i=0;i<n;i++){
                v[i]=0;
        }
}

void vRep(const int l, bool vol[l][l][l], const double cutoff,\
	  fftw_real grd[l][l][2*(l/2+1)], const double scl, const char* str, const double kBT){
        int i,j,k;
        double sum=0.0;
        double count=0.0;
	double dl3=(double)l*(double)l*(double)l;
	//inital vol if str[:3]==vol
	if (!strncmp(str,"vol",3)){
		toBool(l,&grd[0][0][0],&vol[0][0][0],cutoff);
	}
	#pragma omp parallel for private(j,k) reduction (+ : sum,count)
        for (i=0;i<l;i++){
                for (j=0;j<l;j++){
                        for (k=0;k<l;k++){
                                if (vol[i][j][k]){
                                        count+=1.0;
                                }else{
                                        sum+=myEXP(grd[i][j][k]*scl/-kBT);
                                }
                        }
                }
        }
        double u,un;
        if (!strncmp(str,"vol",3)){
                u=-kBT*log(1.0-count/dl3);
                un=0.0;
        }else{
                u=-kBT*log(sum/(dl3));
                un=-kBT*log(sum/(dl3-count));
        }
        showscorenet(str,sum,count/(dl3),u,un);
}

void AddTo1(int l, fftw_real grd1[l][l][2*(l/2+1)], const double scl1, fftw_real grd2[l][l][2*(l/2+1)], const double scl2){
        int i,j,k;
	#pragma omp parallel for private(j,k)
        for (i=0;i<l;i++){
                for (j=0;j<l;j++){
                        for (k=0;k<l;k++){
                                        grd1[i][j][k]=grd1[i][j][k]*scl1+grd2[i][j][k]*scl2;
                        }
                }
        }
}

void SelErn(const int l, bool vol[l][l][l], fftw_real grd1[l][l][2*(l/2+1)],const double ang[3],const double erncut){
	int i,j,k;
	for (i=0;i<l;i++){
		for (j=0;j<l;j++){
			for (k=0;k<l;k++){
				if (vol[i][j][k]){
					continue;
				}else{
					if (grd1[i][j][k]<erncut){
						printf("%10.6f%10.6f%10.6f%8d%8d%8d%10.3f\n",ang[0],ang[1],ang[2],i,j,k,grd1[i][j][k]);
					}
				}
			}
		}
	}
}

void BltzSumFilt(const int l, bool vol[l][l][l], fftw_real sav[l][l][2*(l/2+1)], fftw_real local[l][l][2*(l/2+1)], const double kBT){
	int i,j,k;
	#pragma omp parallel for private(j,k)
	for (i=0;i<l;i++){
                for (j=0;j<l;j++){
                        for (k=0;k<l;k++){
				if (vol[i][j][k]){
					continue;
				}else{
					sav[i][j][k]+=myEXP(local[i][j][k]/-kBT);
				}
			}
		}
	}
}

void getVi(fftw_real val, int *vi){
	const double cutoff=40.0;
	const double scl=2000.0/cutoff;
	if (val>-cutoff&& val <cutoff){
		*vi=(int)(floor(val*scl+0.5))+2000;
	}else{
		printf("%s:%d val:out of range  |%f|< %f\n", __FILE__, __LINE__,val,cutoff);
	}
}

void VEInd(const int l, bool vol[l][l][l], fftw_real ele[l][l][2*(l/2+1)], fftw_real vdw[l][l][2*(l/2+1)], int nbin, double mat[nbin][nbin]){
	int i,j,k;
	#pragma omp parallel for private(j,k)
	for (i=0;i<l;i++){
               for (j=0;j<l;j++){
                        for (k=0;k<l;k++){
                                if (vol[i][j][k]){
                                        continue;
				}else{
					int ei,vi;
					getVi(ele[i][j][k],&ei);
					getVi(vdw[i][j][k],&vi);
					#pragma omp atomic
					mat[ei][vi]+=1.0;
				}
			}
		}
	}
}

void Bltz2Ern(const int l,fftw_real sav[l][l][2*(l/2+1)],const double kBT){
	int i,j,k;
	#pragma omp parallel for private(j,k)
	for (i=0;i<l;i++){
		for (j=0;j<l;j++){
			for (k=0;k<l;k++){
				if (sav[i][j][k]>=1.0){
					sav[i][j][k]=-kBT*log(sav[i][j][k]);
				}
			}
		}
	}
}

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
	//printf("mm:%16.8f%16.8f\n",mm[0],mm[1]);
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
	return mm[0]+(itop+1)*dx;
}
