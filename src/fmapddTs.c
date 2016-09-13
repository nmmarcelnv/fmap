#include "common.h"
#include "zdread.h"
#include "rand.h"
#include "linked.h"
#include "drt.h"

void usage(char *prog){
        printf("Usage: %s subA.vdw subB.vdw maxrad conf.out nconf ion tempK eScl vScl [rEcut rVcut]\n",prog);
}

int main(int argc, char **argv){
#ifdef DEBUG
	fprintf(stderr,"RUNNING DEBUG BUILD\n");
#endif
	if (argc<10){
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }
	ZDCOM zd;
	//int nv=ZDread(argv[1],&zd,(void*)(NULL),(void*)(NULL),(void*)(NULL),false);
	double maxrad=atof(argv[3]);
	long nv=atol(argv[5]);
	double Angs[nv][3],tran[nv][3],score[nv];
	RANPARM ranp1;
	FILE* conffp=fopen(argv[4],"r");
	double conf[nv][6];	
        int i;
	for (i=0;i<nv;i++){
		if (readconf(conffp,&ranp1)==1){
			ranp2tr(&ranp1,conf[i],maxrad);
			tran[i][0]=0.0;tran[i][1]=0.0;tran[i][2]=0.0;
		}
	}
	fclose(conffp);
	//ZDread(argv[1],&zd,Angs,tran,score,true);
	strcpy(zd.rec,argv[1]);
	strcpy(zd.lig,argv[2]);
	zd.N=200;
	zd.spacing=0.6;
	zd.r[0]=0.0;zd.r[1]=0.0;zd.r[2]=0.0;
	zd.l[0]=0.0;zd.l[1]=0.0;zd.l[2]=0.0;
        char* CrdFn=zd.rec;
        char* ProFn=zd.lig;
        int l=zd.N;
	double dx=zd.spacing;
	double lnkc2=4.0*dx;
        double ion=atof(argv[6]);
	double tempK=atof(argv[7]);
	double esclrel=atof(argv[8]);
	double vsclrel=atof(argv[9]);
	double rEcut=12.0;
	double rVcut=12.0;
	if (argc>=11){
		rEcut=atof(argv[10]);
	}
	if (argc>=12){
		rVcut=atof(argv[11]);
	}
        PARM sys;
        sys.rup=rVcut;
	sys.rEup=rEcut;
	sys.rlow=1.0;
	sys.sdie=GetWaterDie(tempK);
	sys.kap=GetKappa(ion,sys.sdie,tempK);
	sys.escl=esclrel*332.0/sys.sdie;
	sys.vscl=vsclrel;
	sys.kBT=GetkBT(tempK);
	//sys.kBT=0.591;
	int nCrd=CountAtoms(CrdFn);
        ATOM Crds[nCrd];
        ReadPqr(CrdFn,nCrd,Crds);
	SetKap(nCrd,Crds,sys.kap);
        ToCtd(nCrd,Crds,zd.r);
	CLINKED lnk;
	int bl=(int)(round)(2*l*dx/lnkc2);
        int box[3]={bl,bl,bl};
	double xyz[nCrd][3];
        for (i=0;i<nCrd;i++){
                xyz[i][0]=Crds[i].xyz[0];
                xyz[i][1]=Crds[i].xyz[1];
                xyz[i][2]=Crds[i].xyz[2];
        }
	double rcut=(sys.rEup>sys.rup)?sys.rEup:sys.rup;
        int nb=(int)ceil(2.0*rcut/lnkc2);
        lnk=lnk_create(nCrd,xyz,nb,lnkc2/2.0,1,box);
	int nPro=CountAtoms(ProFn);
        ATOM Pros[nPro];
        ReadPqr(ProFn,nPro,Pros);
	SetKap(nPro,Pros,sys.kap);
        ToCtd(nPro,Pros,zd.l);
	//printf("%d\t%f\n",l,dx);
        //printf("%f\t%f\t%f\n",0.0,0.0,0.0);
        //printf("%s\t%8.3f%8.3f%8.3f\n",CrdFn,zd.r[0],zd.r[1],zd.r[2]);
        //printf("%s\t%8.3f%8.3f%8.3f\n",ProFn,zd.l[0],zd.l[1],zd.l[2]);
	double ern[4][nv];
	double ernCrd[4][nCrd];
	double ernPro[4][nPro];
	zeroarr(nCrd*4,&(ernCrd[0][0]));
	zeroarr(nPro*4,&(ernPro[0][0]));
	#pragma omp parallel
	{
	ATOM ProCur[nPro];
        int j;
        for (j=0;j<nPro;j++){
                ProCur[j]=Pros[j];
        }
	#pragma omp for schedule(static)
	for (i=0;i<nv;i++){
		double rot[9];
		double xyz[nPro][3],xyztr[nPro][3];
		for (j=0;j<nPro;j++){
			xyz[j][0]=Pros[j].xyz[0];
			xyz[j][1]=Pros[j].xyz[1];
			xyz[j][2]=Pros[j].xyz[2];
		}
		genconf(xyz,nPro,conf[i],xyztr);
		for (j=0;j<nPro;j++){
                	ProCur[j].xyz[0]=xyztr[j][0];
                	ProCur[j].xyz[1]=xyztr[j][1];
                	ProCur[j].xyz[2]=xyztr[j][2];
                }
		
		//Euler2Rot(Angs[i][0],Angs[i][1],Angs[i][2],rot);
		//RotPro(nPro,Pros,ProCur,rot);
		//ShowPqr(stderr,nPro,ProCur);
		softlnkijk1(nCrd,Crds,nPro,ProCur,l,dx,&sys,lnk,nv,tran,ern,ernCrd,ernPro,i);
	}
	}
	for (i=0;i<nv;i++){
		tran[i][0]=conf[i][0];tran[i][1]=conf[i][1];tran[i][2]=conf[i][2];
		Angs[i][0]=conf[i][3];Angs[i][1]=conf[i][4];Angs[i][2]=conf[i][5];
	}
	ijkrep(nCrd,Crds,nPro,Pros,&sys,nv,Angs,tran,score,ern,ernCrd,ernPro,stdout);
	lnk_free(&lnk);
	exit(EXIT_SUCCESS);
}	
