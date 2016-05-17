#include "pro.h"
/*
struct protein{
	int n;
	int *xi,*yi,*zi;
	double *xf,*yf,*zf;
	double *q,*r,*Asq,*Bsq,*kap;
};
*/
//xi=floor(x); xf=x-xi 
void xyzif(double xyz[3], int *xi, int *yi, int *zi, double *xf, double *yf, double *zf){
	//*xi=(int)(ceil(xyz[0]));
	*xi=(int)(floor(xyz[0]));
	//*xi=(int)(xyz[0]);
	//*yi=(int)(ceil(xyz[1]));
	*yi=(int)(floor(xyz[1]));
	//*yi=(int)(xyz[1]);
	//*zi=(int)(ceil(xyz[2]));
	*zi=(int)(floor(xyz[2]));
	//*zi=(int)(xyz[2]);
	*xf=xyz[0]-(*xi);
	*yf=xyz[1]-(*yi);
	*zf=xyz[2]-(*zi);
}

PRO* ProCreate(const int n, ATOM atoms[n]){
	PRO* pro;
	int i;
	if(!(pro = malloc(sizeof(PRO)))) return NULL;

	pro->n=n;
	if(!(pro->xi =malloc(n*sizeof(int)))) return NULL;
	if(!(pro->yi =malloc(n*sizeof(int)))) return NULL;
	if(!(pro->zi =malloc(n*sizeof(int)))) return NULL;
	if(!(pro->xf =malloc(n*sizeof(double)))) return NULL;
	if(!(pro->yf =malloc(n*sizeof(double)))) return NULL;
	if(!(pro->zf =malloc(n*sizeof(double)))) return NULL;
	if(!(pro->q =malloc(n*sizeof(double)))) return NULL;
	if(!(pro->r =malloc(n*sizeof(double)))) return NULL;
	if(!(pro->Asq =malloc(n*sizeof(double)))) return NULL;
	if(!(pro->Bsq =malloc(n*sizeof(double)))) return NULL;
	if(!(pro->kap =malloc(n*sizeof(double)))) return NULL;

	#pragma omp parallel for
	for (i=0;i<n;i++){
		xyzif(atoms[i].xyz,&(pro->xi[i]),&(pro->yi[i]),&(pro->zi[i]),&(pro->xf[i]),&(pro->yf[i]),&(pro->zf[i]));
		pro->q[i]=atoms[i].q;
		pro->r[i]=atoms[i].r;
		pro->Asq[i]=atoms[i].Asq;
		pro->Bsq[i]=atoms[i].Bsq;
		pro->kap[i]=atoms[i].kap;
	}

#ifdef DEBUG
	assert(pro->n>0);
	assert(pro!=NULL);
	assert(pro->xi!=NULL);
	assert(pro->yi!=NULL);
	assert(pro->zi!=NULL);	
	assert(pro->xf!=NULL);
	assert(pro->yf!=NULL);
	assert(pro->zf!=NULL);
	assert(pro->q!=NULL);
        assert(pro->r!=NULL);
        assert(pro->Asq!=NULL);
	assert(pro->Bsq!=NULL);
        assert(pro->kap!=NULL);
	printf("%s\t:%d\t:Number of atoms:\t%d\n",__FILE__,__LINE__,pro->n);
#endif
	return pro;
}
void ProUpdate(PRO* pro, const int n, double xyz[n][3]){
	int i;
#ifdef DEBUG
	assert(n==pro->n);
#endif
	#pragma omp parallel for
	for (i=0;i<n;i++){
		xyzif(xyz[i],&(pro->xi[i]),&(pro->yi[i]),&(pro->zi[i]),&(pro->xf[i]),&(pro->yf[i]),&(pro->zf[i]));
	}
}

void ProFree(PRO* pro){
	if (pro){
		free(pro->xi);
		free(pro->xf);
		free(pro->yi);
		free(pro->zi);
		free(pro->q);
		free(pro->r);
		free(pro->Asq);
		free(pro->Bsq);
		free(pro->kap);
		free(pro);
	}
}
