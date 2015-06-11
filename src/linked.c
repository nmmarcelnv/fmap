#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "linked.h"

#define EMPTY -1

void Grdxyz(double x, double y, double z,double bound[6], double spacing, int* gx, int* gy, int* gz);
void getbound(double xyz[][3],int num, double bound[6]);
int calcsize(double bound[6], double spacing, int box[3]);
int Indxyz(double x, double y, double z, CLINKED lnk);
void buildlinked(double xyz[][3], int num, CLINKED lnk);

struct CLINKED{
	int num;
	int box[3];
	double boxlen[3];
	int lcxyz,lcyz,lcz;
	int* head;
	int* lscl;
	double spacing;
	double bound[6];
	double outer[6];
	int pbc;
	int nb;
	double (*xyz)[3];
};
 
int pbc1shf(int x, int l, int* shift){
	int px;
	if (x>=0 && x<l){
		*shift=0;
		px=x;
	}else{
		px=x%l;
		if (px<0){
			px += l;
			*shift=-1;
		}else{
			*shift=1;
		}
	}
	return px;
}

void Grdxyz(double x, double y, double z,double bound[6], double spacing, int* gx, int* gy, int* gz){
        *gx=(int)round(((x-bound[0])/spacing));
        *gy=(int)round(((y-bound[1])/spacing));
        *gz=(int)round(((z-bound[2])/spacing));
	//fprintf(stderr,"Grdxyz: %8.3f%8.3f%8.3f%8d%8d%8d\n",x,y,z,*gx,*gy,*gz);
}

void getbound(double xyz[][3],int num, double bound[6]){
        int i,j;
        for (j=0;j<3;j++){
                bound[j]=FLT_MAX;
                bound[j+3]=-FLT_MAX;
        }
        for (i=0;i<num;i++){
                for (j=0;j<3;j++){
                        if (xyz[i][j]<bound[j]) bound[j]=xyz[i][j];
                        if (xyz[i][j]>bound[j+3]) bound[j+3]=xyz[i][j];
                }
        }
}

int calcsize(double bound[6], double spacing, int box[3]){
        int j;
        for (j=0;j<3;j++){
                box[j]=(int)ceil((bound[j+3]-bound[j])/spacing);
        }
        return box[0]*box[1]*box[2];
}

int Indxyz(double x, double y, double z, CLINKED lnk){
        int cx,cy,cz;
	int shift;
        Grdxyz(x,y,z,lnk->bound,lnk->spacing,&cx,&cy,&cz);
	//fprintf(stderr,"%d %d %d\n",cx,cy,cz);
	if (lnk->pbc){
		cx=pbc1shf(cx,lnk->box[0],&shift);
		cy=pbc1shf(cy,lnk->box[1],&shift);
		cz=pbc1shf(cz,lnk->box[2],&shift);
	}
        return cx*lnk->lcyz+cy*lnk->lcz+cz;
}

void buildlinked(double xyz[][3], int num, CLINKED lnk){
        int i,c;
        for (i=0;i<lnk->lcxyz;i++){
                lnk->head[i]=EMPTY;
        }
        for (i=0;i<num;i++){
                c=Indxyz(xyz[i][0],xyz[i][1],xyz[i][2],lnk);
                lnk->lscl[i]=lnk->head[c];
                lnk->head[c]=i;
        }
}
//box[i]*spacing == PBC boundary
//nb*spacing >=cutoff
CLINKED lnk_create(int num, double xyz[][3], int nb, double spacing, int pbc, int box[3]){
	CLINKED lnk;
	int i;
	//struct CLINKED instead of CLINKED, otherwise malloc(): memory corruption
	if(!(lnk = malloc(sizeof(struct CLINKED)))) {
                return 0;
        }
	lnk->num=num;
	lnk->nb=nb;
	lnk->spacing=spacing;
	lnk->pbc=pbc;
	lnk->xyz=&xyz[0];
	if (lnk->pbc){
		for (i=0;i<3;i++){
			
			lnk->box[i]=box[i];
			lnk->bound[i]=-0.5*lnk->box[i]*lnk->spacing;
			lnk->bound[i+3]=0.5*lnk->box[i]*lnk->spacing;
			lnk->boxlen[i]=lnk->box[i]*lnk->spacing;
		}
		lnk->lcxyz=lnk->box[0]*lnk->box[1]*lnk->box[2];
	}else{
		getbound(xyz,num,lnk->bound);
		for (i=0;i<3;i++){
			lnk->outer[i]=lnk->bound[i]-lnk->nb*lnk->spacing;
			lnk->outer[i+3]=lnk->bound[i+3]+lnk->nb*lnk->spacing;
		}
		lnk->lcxyz=calcsize(lnk->bound,lnk->spacing,lnk->box);
	}
	lnk->lcyz=lnk->box[1]*lnk->box[2];
	lnk->lcz=lnk->box[2];
	if(!(lnk->head=malloc(sizeof(int)*lnk->lcxyz))){
		return 0;
	}
	if(!(lnk->lscl=malloc(sizeof(int)*num))){
                return 0;
        }
	buildlinked(xyz,num,lnk);
	return lnk;
}

void lnk_free(CLINKED *lnk)
{
	assert(lnk && *lnk);
	assert((*lnk)->head);
        free((*lnk)->head);
	(*lnk)->head=NULL;
	assert((*lnk)->lscl);
	free((*lnk)->lscl);
	(*lnk)->lscl=NULL;
	free((*lnk)->lscl);
	assert(*lnk);
        free(*lnk);
	(*lnk)=NULL;
}

double pbcdst(double x, double boxl){
	double frac,fold;
	frac=modf(fabs(x)/boxl,&fold);
	frac=(frac<0.5)?frac:1-frac;
	return frac*boxl;
}

int lnk_nearest(CLINKED lnk, double x, double y, double z, int maxnb, int idx[], double r2[]){
	int cnt=0;
	int cx,cy,cz;
	int ii,jj,kk;
	int c1,vj;
	int pi,pj,pk;
	int shx,shy,shz;
	double mx,my,mz;
	int nb;
	vj=EMPTY; //avoid vj undifined when xyz is out: 2014/01/08 qsb
	nb=lnk->nb;
	pi=0;pj=0;pk=0;
	shx=0;shy=0;shz=0;
	if (!lnk->pbc){
		if (x<lnk->outer[0]||x>lnk->outer[3]||y<lnk->outer[1]||y>lnk->outer[4]
			||z<lnk->outer[2]||z>lnk->outer[5]){
			return 0;
		}
	}
	Grdxyz(x,y,z,lnk->bound,lnk->spacing,&cx,&cy,&cz);
	for (ii=cx-nb;ii<=cx+nb;ii++){
	  if (lnk->pbc) pi=pbc1shf(ii,lnk->box[0],&shx);
       	  for (jj=cy-nb;jj<=cy+nb;jj++){
            if (lnk->pbc) pj=pbc1shf(jj,lnk->box[1],&shy);
            for (kk=cz-nb;kk<=cz+nb;kk++){
		  if (lnk->pbc){
		    pk=pbc1shf(kk,lnk->box[2],&shz);
                    c1=pi*lnk->lcyz+pj*lnk->lcz+pk;
                    vj=lnk->head[c1];
		  }else{  
		    if (ii>=0&&ii<lnk->box[0]&&jj>=0&&jj<lnk->box[1]&&kk>=0&&kk<lnk->box[2]){
			c1=ii*lnk->lcyz+jj*lnk->lcz+kk;
                        vj=lnk->head[c1];
		    }
                  }
                  while (vj!=EMPTY){
			idx[cnt]=vj;
			mx=x-lnk->xyz[vj][0];
			my=y-lnk->xyz[vj][1];
			mz=z-lnk->xyz[vj][2];
			if (lnk->pbc){
				mx=pbcdst(mx,lnk->boxlen[0]);
				my=pbcdst(my,lnk->boxlen[1]);
				mz=pbcdst(mz,lnk->boxlen[2]);
			}
			r2[cnt]=mx*mx+my*my+mz*mz;
			if (cnt<maxnb-1){
				cnt++;
			}else{
				fprintf(stderr,"%s:%d INCREASE NB, %d IS TOO SMALL\n",__FILE__,__LINE__,maxnb);
			}
			vj=lnk->lscl[vj];
		  }
	    }
          }
	}
	return cnt;
}
