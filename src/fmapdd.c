#include "common.h"
#include "zdread.h"
#include "linked.h"
#include "drt.h"

void usage(char *prog){
        printf("Usage: %s fmap.out ion tempK\n",prog);
}

void zeroarr(int n, double *v){
	int i;
	for (i=0;i<n;i++){
		v[i]=0;
	}
}

int main(int argc, char **argv){
#ifdef DEBUG
    printf("RUNNING DEBUG BUILD\n");
#endif
	 if (argc<6){
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }
	ZDCOM zd;
	int nv=ZDread(argv[1],&zd,(void*)(NULL),(void*)(NULL),(void*)(NULL),false);
	double Angs[nv][3],tran[nv][3],score[nv];
	ZDread(argv[1],&zd,Angs,tran,score,true);
        char* CrdFn=zd.rec;
        char* ProFn=zd.lig;
        int l=zd.N;
	double dx=zd.spacing;
	double lnkc2=4.0*dx;
        double ion=atof(argv[2]);
	double tempK=atof(argv[3]);
	double esclrel=atof(argv[4]);
	double vsclrel=atof(argv[5]);
        PARM sys;
        sys.rup=12.0;
	sys.rlow=1.0;
	sys.dx=dx;
	sys.sdie=78.5;
	sys.kap=(double)GetKappa(ion,(double)sys.sdie,tempK);
	sys.escl=esclrel*332.0/sys.sdie;
	sys.vscl=vsclrel;
	sys.kBT=GetkBT(tempK);
	//sys.kBT=0.591;
	int nCrd=CountAtoms(CrdFn);
        ATOM Crds[nCrd];
        ReadPqr(CrdFn,nCrd,Crds);
	SetKap(nCrd,Crds,sys.kap);
        ToCtd(nCrd,Crds,zd.r);
	CLINKED* lnk;
	int bl=(int)(round)(2*l*dx/lnkc2);
        int box[3]={bl,bl,bl};
	double xyz[nCrd][3];
        int i;
        for (i=0;i<nCrd;i++){
                xyz[i][0]=Crds[i].xyz[0];
                xyz[i][1]=Crds[i].xyz[1];
                xyz[i][2]=Crds[i].xyz[2];
        }
        lnk=lnk_create(xyz,nCrd,lnkc2,1,box);
        int nb=(int)ceil(2.0*sys.rup/lnkc2);
        lnk_setnb(lnk,nb);
	int nPro=CountAtoms(ProFn);
        ATOM Pros[nPro];
        ReadPqr(ProFn,nPro,Pros);
	SetKap(nPro,Pros,sys.kap);
        ToCtd(nPro,Pros,zd.l);
	printf("%d\t%f\n",l,dx);
        printf("%f\t%f\t%f\n",0.0,0.0,0.0);
        printf("%s\t%8.3f%8.3f%8.3f\n",CrdFn,zd.r[0],zd.r[1],zd.r[2]);
        printf("%s\t%8.3f%8.3f%8.3f\n",ProFn,zd.l[0],zd.l[1],zd.l[2]);
	double ern[4][nv];
	double ernCrd[4][nCrd];
	double ernPro[4][nPro];
	zeroarr(nCrd*4,&(ernCrd[0][0]));
	zeroarr(nPro*4,&(ernPro[0][0]));
	#pragma omp parallel 
	{
	#pragma omp for schedule(dynamic,1) nowait
	for (i=0;i<nv;i++){
		double rot[9];
		ATOM ProCur[nPro];
		Euler2Rot(Angs[i][0],Angs[i][1],Angs[i][2],rot);
		RotPro(nPro,Pros,ProCur,rot);
		softlnkijk1(nCrd,Crds,nPro,ProCur,l,dx,&sys,lnk,nv,tran,ern,ernCrd,ernPro,i);
	}
	}
	ijkrep(nCrd,Crds,nPro,Pros,&sys,nv,Angs,tran,score,ern,ernCrd,ernPro,stdout);
	lnk_free(lnk);
	exit(EXIT_SUCCESS);
}	
