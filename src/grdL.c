#include "grd.h"

int setIndWght(char* fname, int n, int Ind[][10], double Wght[][10]){
                FILE *file;
        int pos;
        //int ns;
        if (FileExist(fname)==0){
                exit(EXIT_FAILURE);
        }
        file=fopen(fname, "r");
        pos=0;
        while(20==fscanf(file,"%d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf",&Ind[pos][0],&Wght[pos][0],&Ind[pos][1],&Wght[pos][1],&Ind[pos][2],&Wght[pos][2],&Ind[pos][3],&Wght[pos][3],&Ind[pos][4],&Wght[pos][4],&Ind[pos][5],&Wght[pos][5],&Ind[pos][6],&Wght[pos][6],&Ind[pos][7],&Wght[pos][7],&Ind[pos][8],&Wght[pos][8],&Ind[pos][9],&Wght[pos][9])){
                /*
                if (ns!=20){
                        fprintf(stderr,"%s pos: %d ns: %d != 20\n",__FUNCTION__,pos,ns);
                        fprintf(stderr,"%d %20.16lf %d %20.16lf\n",Ind[pos][0],Wght[pos][0],Ind[pos][9],Wght[pos][9]);
                }
                assert(ns==20);
                */
                pos++;
        }
        fclose(file);
#ifdef DEBUG
        assert(pos==n);
	assert(pos>0);
#endif
        return pos;
}

//calcualted index(i) for 1d[100] from 3d[10][10][10] with finer grid at 1/10th
void xyz2pos(const int n, const double xf[n], const double yf[n], const double zf[n], int pos[n]){
	#pragma omp parallel for schedule(static)	
	for (int i=0;i<n;i++){
		pos[i]=100*((int)(xf[i]*10.0))+10*((int)(yf[i]*10.0))+(int)(zf[i]*10.0);
#ifdef DEBUG
		assert(pos[i]<1000);
#endif
	}
}

//calcualte index(i,j,k) for 3d[4][4][4] from 1d[16]
void ind2ijk(const int n, const int loc[n], int xoff[n], int yoff[n], int zoff[n]){
	#pragma omp for schedule(static) nowait
	for (int i=0;i<n;i++){
		xoff[i]=loc[i]/4/4;
                yoff[i]=(loc[i]-xoff[i]*4*4)/4;
		zoff[i]=loc[i]-xoff[i]*4*4-yoff[i]*4;
#ifdef DEBUG
		assert(xoff[i]<4);
		assert(yoff[i]<4);
		assert(zoff[i]<4);
#endif
	}
}

void setlocwgt(const int n, int pos[n], const int m, int loc[n], double w[n], int Ind[1000][10], double Wght[1000][10]){
	#pragma omp for schedule(static) nowait
	for (int i=0;i<n;i++){
		loc[i]=Ind[pos[i]][m];
		w[i]=Wght[pos[i]][m];
	}
}

void mult2(const int n, const double w[n], const double coef[n], double wn[n]){
	#pragma omp for schedule(static) nowait
	for (int i=0;i<n;i++){
		wn[i]=w[i]*coef[i];
	}
}

void pbc1n(const int n, const int x[n], int px[n], const int lc, const int off[n], const int shift){
        #pragma omp for schedule(static) nowait
        for (int i=0;i<n;i++){
                px[i]=PBC(x[i]+off[i]+shift,lc);
	}
}

void grdL2nd(const int l, fftw_real* grid,\
             const int n, const int xi[n], const int yi[n], const int zi[n] ,\
	     const double xf[n], const double yf[n], const double zf[n], const double coef1[n]){
	int pos[n],xoff[n],yoff[n],zoff[n],loc[n],px[n],py[n],pz[n];
	double w[n],wn[n];
	static int nlattice;
        static int Ind[1000][10];
        static double Wght[1000][10];
	int m;
#ifdef DEBUG
	printf("%i %i %i %f %f %f\n",xi[n-1],yi[n-1],zi[n-1],xf[n-1],yf[n-1],zf[n-1]);
#endif
        if (nlattice==0){
                nlattice=setIndWght("lattice.txt",1000,Ind,Wght);
        }
	xyz2pos(n,xf,yf,zf,pos);
	#pragma omp parallel private(m)
	for (m=0;m<10;m++){
		setlocwgt(n,pos,m,loc,w,Ind,Wght);
		ind2ijk(n,loc,xoff,yoff,zoff);
		mult2(n,w,coef1,wn);
		pbc1n(n,xi,px,l,xoff,-1);
		pbc1n(n,yi,py,l,yoff,-1);
		pbc1n(n,zi,pz,l,zoff,-1);
		val2grd(n,wn,px,py,pz,l,grid);
	}
}

//xc=1-xf
void setcomp(const int n, const double xf[n], double xc[n]){
	#pragma omp for schedule(static) nowait
	for (int i=0;i<n;i++){
		xc[i]=1.0-xf[i];
#ifdef DEBUG
		assert(xc[i]<=1.0);
		assert(xc[i]>=0.0);
#endif
	}
}

void mult3(const int n, const double xc[n], const double yc[n], const double zc[n], double w[n]){
	#pragma omp for schedule(static) nowait
	for (int i=0;i<n;i++){
		w[i]=xc[i]*yc[i]*zc[i];
	}
}

void grdL(const int l, fftw_real* grid,\
        const int n, const int xi[], const int yi[], const int zi[], \
	const double xf[], const double yf[], const double zf[], const double coef1[]){
	int px[n],py[n],pz[n],pxn[n],pyn[n],pzn[n];
	double xc[n], yc[n], zc[n], w[n], wn[n];
#ifdef DEBUG
	printf("%i %i %i %f %f %f\n",xi[n-1],yi[n-1],zi[n-1],xf[n-1],yf[n-1],zf[n-1]);
#endif
	#pragma omp parallel
	{
	setcomp(n,xf,xc);
	setcomp(n,yf,yc);
	setcomp(n,zf,zc);
	pbc1(n,xi,px,l,0);
	pbc1(n,yi,py,l,0);
	pbc1(n,zi,pz,l,0);
	pbc1(n,xi,pxn,l,1);
        pbc1(n,yi,pyn,l,1);
        pbc1(n,zi,pzn,l,1);
	mult3(n,xc,yc,zc,w);mult2(n,w,coef1,wn);val2grd(n,wn,px,py,pz,l,grid);
	mult3(n,xc,yc,zf,w);mult2(n,w,coef1,wn);val2grd(n,wn,px,py,pzn,l,grid);
	mult3(n,xc,yf,zc,w);mult2(n,w,coef1,wn);val2grd(n,wn,px,pyn,pz,l,grid);
	mult3(n,xc,yf,zf,w);mult2(n,w,coef1,wn);val2grd(n,wn,px,pyn,pzn,l,grid);
	mult3(n,xf,yc,zc,w);mult2(n,w,coef1,wn);val2grd(n,wn,pxn,py,pz,l,grid);
        mult3(n,xf,yc,zf,w);mult2(n,w,coef1,wn);val2grd(n,wn,pxn,py,pzn,l,grid);
        mult3(n,xf,yf,zc,w);mult2(n,w,coef1,wn);val2grd(n,wn,pxn,pyn,pz,l,grid);
        mult3(n,xf,yf,zf,w);mult2(n,w,coef1,wn);val2grd(n,wn,pxn,pyn,pzn,l,grid);
	}
}

void LGrdR12(const int l, fftw_real* grid, PRO *pro){
        grdL(l,grid,pro->n,pro->xi,pro->yi,pro->zi,pro->xf,pro->yf,pro->zf,pro->Asq);
}

void LGrdR6(const int l, fftw_real* grid, PRO *pro){
        grdL(l,grid,pro->n,pro->xi,pro->yi,pro->zi,pro->xf,pro->yf,pro->zf,pro->Bsq);
}

void LGrd2ndEle(const int l, fftw_real* grid, PRO *pro){
	grdL2nd(l,grid,pro->n,pro->xi,pro->yi,pro->zi,pro->xf,pro->yf,pro->zf,pro->q);
}
