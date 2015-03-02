#include "common.h"
#include "linked.h"
#include "drt.h"

void ernzero(int n, ATOM atoms[n]){
	int i;
	for (i=0;i<n;i++){
		atoms[i].ern=0.0;
	}
}

void ernaddPro(int n, ATOM atoms[n], double ern[n], double scl){
	int i;
	for (i=0;i<n;i++){
		ern[i]+=atoms[i].ern*scl;
	}
}

void ernadd(int n, double e1[n],double e2[n],double e3[n]){
	int i;
	for (i=0;i<n;i++){
		e3[i]=e1[i]+e2[i];
	}
}

void ernscl(int n, double e1[n], double scl){
	int i;
	for (i=0;i<n;i++){
		e1[i]=e1[i]*scl;
	}
}

void ViIJK(FILE *fp,int nv, double angs[nv][3], double ijk[nv][3], double score[nv], double ern[4][nv]){
	int i;
        for (i=0;i<nv;i++){
		fprintf(fp,"%10.6f%10.6f%10.6f%8.0f%8.0f%8.0f%10.3lf%16.8le%16.8le%16.8le%16.8le\n",angs[i][0],angs[i][1],angs[i][2],ijk[i][0],ijk[i][1],ijk[i][2],score[i],ern[3][i],ern[0][i],ern[1][i],ern[2][i]);
	}		
}

double VolChk(double val){
	return (val>0.001)?1.0:0.0;
}

void ShowErn(FILE *fp, int nAtom, ATOM atoms[nAtom], double ern[4][nAtom]){
        int i;
	double bf;
        for (i=0;i<nAtom;i++){
		if (ern[3][i]>100){
			bf=-99.99;
		}else{
			bf=-ern[3][i];
		}
                fprintf(fp,"%30s%8.3lf%8.3lf%8.3lf%6.2lf%6.2lf%16.8le%16.8le%16.8le%16.8le\n",
                        atoms[i].info,atoms[i].xyz[0],atoms[i].xyz[1],atoms[i].xyz[2],
                        VolChk(ern[0][i]),bf,ern[0][i],ern[1][i],ern[2][i],ern[3][i]);
        }
}

void softlnkijk1(int nCrd, ATOM Crds[], int nPro, ATOM Pros[], int l, double dx, PARM* sys,CLINKED* lnk,
	int nv, double ijk[nv][3], double ern[4][nv], double ernCrd[4][nCrd], double ernPro[4][nPro], int pos){
	const double escl=sys->escl;
	const double vscl=sys->vscl;
	double offset[3];
	double rlow2=sys->rlow*sys->rlow;
	int j;
	for (j=0;j<3;j++){
		offset[j]=ijk[pos][j]*dx;
	}
	ernzero(nCrd, Crds);
	ernzero(nPro, Pros);
	ern[0][pos]=drtlnkf(nCrd, Crds, nPro, Pros,offset,rlow2,sys,lnk,vollnk);
	ernaddPro(nCrd, Crds,ernCrd[0],1.0);
	ernaddPro(nPro, Pros,ernPro[0],1.0);
	
	ernzero(nCrd, Crds);
        ernzero(nPro, Pros);
        ern[1][pos]=drtlnkf(nCrd, Crds, nPro, Pros,offset,rlow2,sys,lnk,elelnk);
        ernaddPro(nCrd, Crds,ernCrd[1],escl);
        ernaddPro(nPro, Pros,ernPro[1],escl);

	ernzero(nCrd, Crds);
        ernzero(nPro, Pros);
        ern[2][pos]=drtlnkf(nCrd, Crds, nPro, Pros,offset,rlow2,sys,lnk,vdwlnk);
        ernaddPro(nCrd, Crds,ernCrd[2],vscl);
        ernaddPro(nPro, Pros,ernPro[2],vscl);
}

void softijk1(int nCrd, ATOM Crds[], int nPro, ATOM Pros[], int l, double dx, PARM* sys,CLINKED* lnk,
	int nv, double ijk[nv][3], double ern[4][nv], double ernCrd[4][nCrd], double ernPro[4][nPro], int pos){
	const double escl=sys->escl;
	const double vscl=sys->vscl;
	double offset[3];
	double rlow2=sys->rlow*sys->rlow;
	int j;
	for (j=0;j<3;j++){
		offset[j]=ijk[pos][j]*dx;
	}
	ernzero(nCrd, Crds);
	ernzero(nPro, Pros);
	ern[0][pos]=drtf(nCrd, Crds, nPro, Pros,offset,rlow2,sys,volmul);
	ernaddPro(nCrd, Crds,ernCrd[0],1.0);
	ernaddPro(nPro, Pros,ernPro[0],1.0);
	
	ernzero(nCrd, Crds);
        ernzero(nPro, Pros);
        ern[1][pos]=drtf(nCrd, Crds, nPro, Pros,offset,rlow2,sys,elemul);
        ernaddPro(nCrd, Crds,ernCrd[1],escl);
        ernaddPro(nPro, Pros,ernPro[1],escl);

	ernzero(nCrd, Crds);
        ernzero(nPro, Pros);
        ern[2][pos]=drtf(nCrd, Crds, nPro, Pros,offset,rlow2,sys,vdwmul);
        ernaddPro(nCrd, Crds,ernCrd[2],vscl);
        ernaddPro(nPro, Pros,ernPro[2],vscl);
}

void ijkrep(int nCrd, ATOM Crds[], int nPro, ATOM Pros[], PARM* sys, int nv, double angs[nv][3], double ijk[nv][3], double score[nv], double ern[4][nv],double ernCrd[4][nCrd], double ernPro[4][nPro], FILE *fp){
	const double escl=sys->escl;
        const double vscl=sys->vscl;
	ernscl(nv,ern[1],escl);	
	ernscl(nv,ern[2],vscl);	
	ernadd(nv,ern[1],ern[2],ern[3]);
        ernadd(nCrd,ernCrd[1],ernCrd[2],ernCrd[3]);
        ernadd(nPro,ernPro[1],ernPro[2],ernPro[3]);
        ViIJK(fp,nv,angs,ijk,score,ern);
        ShowErn(fp,nCrd,Crds,ernCrd);
        ShowErn(fp,nPro,Pros,ernPro);
}
