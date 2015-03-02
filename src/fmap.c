#include "common.h"

#ifdef _OPENMP
	fftw_real *recVol, *recR12,*recR6, *recEle;
#endif

void usage(char *prog){
        printf("Usage: %s box.pdb pro.pdb ngrid gridsize ion SclRad SclQ tempK ang.dat eScl vScl Ecut\n",prog);
}

int main(int argc, char **argv){
#ifdef DEBUG
    printf("RUNNING DEBUG BUILD\n");
#endif
	 if (argc<13){
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }
        char* CrdFn=argv[1];
        char* ProFn=argv[2];
        int l=atoi(argv[3]);
        double dx=atof(argv[4]);
        double ion=atof(argv[5]);
	double rscl=atof(argv[6]);
	double qscl=atof(argv[7]);
	double tempK=atof(argv[8]);
	char* ang=argv[9];
	double esclrel=atof(argv[10]);
	double vsclrel=atof(argv[11]);
	double erncut=atof(argv[12]);
	fftw_init_threads();
	fftw_plan_with_nthreads(omp_get_max_threads());
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
	sys.erncut=erncut; //-1.0*sys.kBT*(3.0*log((double)(l))-log(10.0));
	double cenRec[3],cenLig[3];
	int nCrd=CountAtoms(CrdFn);
        ATOM Crds[nCrd];
        ReadPqr(CrdFn,nCrd,Crds);
	SetKap(nCrd,Crds,sys.kap);
	SclRad(nCrd,Crds,rscl);
	SclChg(nCrd,Crds,qscl);
        //CalCtd(nCrd,Crds,cenRec);
        //ToCtd(nCrd,Crds,cenRec);
        cenRec[0]=0.0;cenRec[1]=0.0;cenRec[2]=0.0;
	Unit2dx(nCrd,Crds,dx);
	int nPro=CountAtoms(ProFn);
        ATOM Pros[nPro];
        ReadPqr(ProFn,nPro,Pros);
	SetKap(nPro,Pros,sys.kap);
	SclRad(nPro,Pros,rscl);
        //CalCtd(nPro,Pros,cenLig);
        //ToCtd(nPro,Pros,cenLig);
        cenLig[0]=0.0;cenLig[1]=0.0;cenLig[2]=0.0;
	Unit2dx(nPro,Pros,dx);
	printf("%d\t%f\n",l,dx);
        printf("%f\t%f\t%f\n",0.0,0.0,0.0);
        printf("%s\t%8.3f%8.3f%8.3f\n",CrdFn,cenRec[0],cenRec[1],cenRec[2]);
        printf("%s\t%8.3f%8.3f%8.3f\n",ProFn,cenLig[0],cenLig[1],cenLig[2]);
	PRO *rec,*lig;
	rec=ProCreate(nCrd,Crds);
	lig=ProCreate(nPro,Pros);
	fftw_real sav[l][l][2*(l/2+1)];
	SetZero(l,sav);
        int nAng=CountAng(ang);
        double Angs[nAng][3];
	SetAng(ang,Angs);
	int i;
	double xyzLig[nPro][3];
	double rot[9];
#ifdef _OPENMP
	fftw_real extVol[l][l][2*(l/2+1)],extEle[l][l][2*(l/2+1)];
	fftw_real extR12[l][l][2*(l/2+1)],extR6[l][l][2*(l/2+1)];
	recVol=&(extVol[0][0][0]);
	recEle=&(extEle[0][0][0]);
	recR12=&(extR12[0][0][0]);
	recR6=&(extR6[0][0][0]);
#endif
	for (i=0;i<nAng;i++){
		Euler2Rot(Angs[i][0],Angs[i][1],Angs[i][2],rot);
		RotXYZ(nPro,Pros,xyzLig,rot);
		ProUpdate(lig,nPro,xyzLig);
		Cross(rec,lig,sys,l,sav,Angs[i]);
	}
	double ang0[3]={0.0,0.0,0.0};
	bool vol[l][l][l];
	SetBool(l,vol,false);
	Bltz2Ern(l,sav,sys.kBT);
	vRep(l,vol,0.001,sav,1.0,"sav",sys.kBT);
	SelErn(l,vol,sav,ang0,sys.erncut);	
	ProFree(rec);
	ProFree(lig);	
	fftw_cleanup_threads();
	exit(EXIT_SUCCESS);
}	
