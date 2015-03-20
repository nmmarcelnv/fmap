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
	Times(true,3,0,8,9);

#ifdef USE_MPI
	int rank,size;
        MPI_Init(&argc,&argv);
        MPI_Comm_size(MPI_COMM_WORLD,&size);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif /* USE_MPI */
	Times(false,1,8);
        
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
//#ifdef _OPENMP
//	printf("TACC Nthread:%d\n",omp_get_max_threads());
//	printf("TACC Size:%d\n",size);
//	printf("TACC rank:%d\n",rank);
//#endif

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

	int nCrd,nPro,nAng;
#ifdef USE_MPI
	if (rank==0){
#endif /* USE_MPI */
		nCrd=CountAtoms(CrdFn);
		nPro=CountAtoms(ProFn);
                nAng=CountAng(ang);
#ifdef USE_MPI
        }
	Times(true,1,8);
	MPI_Bcast(&nCrd,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nPro,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nAng,1,MPI_INT,0,MPI_COMM_WORLD);
	
	//https://www.msi.umn.edu/workshops/mpi/hands-on/derived-datatypes/struct/assign
	//define MPI type for ATOM modified for MPI-2
	MPI_Datatype MPI_ATOM;
	int MPI_ATOM_len[2]={INFOLEN,9};
	MPI_Datatype MPI_ATOM_old[2]={MPI_CHAR,MPI_DOUBLE};
	MPI_Aint  MPI_ATOM_idx[2];
	//http://stackoverflow.com/questions/9864510/struct-serialization-in-c-and-transfer-over-mpi
	MPI_ATOM_idx[0]=offsetof(ATOM,info);
	MPI_ATOM_idx[1]=offsetof(ATOM,xyz);
        MPI_Type_create_struct(2,MPI_ATOM_len,MPI_ATOM_idx,MPI_ATOM_old,&MPI_ATOM);
	MPI_Type_commit(&MPI_ATOM);
	Times(false,1,8);
#endif /* USE_MPI */

        ATOM Crds[nCrd];
        ATOM Pros[nPro];
        double Angs[nAng][3];
#ifdef USE_MPI
	if (rank==0){
#endif /* USE_MPI */
		double cenRec[3]= { 0.0 };
		double cenLig[3]= { 0.0 };
        	ReadPqr(CrdFn,nCrd,Crds);
		SetKap(nCrd,Crds,sys.kap);
		SclRad(nCrd,Crds,rscl);
		SclChg(nCrd,Crds,qscl);
		Unit2dx(nCrd,Crds,dx);

        	ReadPqr(ProFn,nPro,Pros);
		SetKap(nPro,Pros,sys.kap);
		SclRad(nPro,Pros,rscl);
		Unit2dx(nPro,Pros,dx);

		printf("%d\t%f\n",l,dx);
        	printf("%f\t%f\t%f\n",0.0,0.0,0.0);
        	printf("%s\t%8.3f%8.3f%8.3f\n",CrdFn,cenRec[0],cenRec[1],cenRec[2]);
        	printf("%s\t%8.3f%8.3f%8.3f\n",ProFn,cenLig[0],cenLig[1],cenLig[2]);

		SetAng(ang,Angs);
#ifdef USE_MPI
	}
	Times(false,1,8);
	MPI_Bcast(Crds,nCrd,MPI_ATOM,0,MPI_COMM_WORLD);
	MPI_Bcast(Pros,nPro,MPI_ATOM,0,MPI_COMM_WORLD);
	MPI_Bcast(Angs,nAng*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
	Times(false,1,8);
#endif /* USE_MPI */

	PRO *rec,*lig;
	rec=ProCreate(nCrd,Crds);
	lig=ProCreate(nPro,Pros);
	fftw_real sav[l][l][2*(l/2+1)];
	SetZero(l,sav);
	
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
	Times(false,1,0);
	Times(true,1,7);
	for (i=0;i<nAng;i++){
#ifdef USE_MPI
		if (i%size==rank){
#endif /* USE_MPI */
			Euler2Rot(Angs[i][0],Angs[i][1],Angs[i][2],rot);
			RotXYZ(nPro,Pros,xyzLig,rot);
			ProUpdate(lig,nPro,xyzLig);
			Cross(rec,lig,sys,l,sav,Angs[i]);
#ifdef USE_MPI
		}
#endif /* USE_MPI */
	}
	Times(false,1,7);
	bool vol[l][l][l];
	SetBool(l,vol,false);
	Bltz2Ern(l,sav,sys.kBT);
	vRep(l,vol,0.001,sav,1.0,"sav",sys.kBT);
#ifdef USE_MPI
	if (size==1){
#endif /* USE_MPI */
		double ang0[3]={0.0,0.0,0.0};
		SelErn(l,vol,sav,ang0,sys.erncut);
#ifdef USE_MPI
	}
#endif /* USE_MPI */
	ProFree(rec);
	ProFree(lig);	
	fftw_cleanup_threads();
#ifdef USE_MPI
	MPI_Type_free(&MPI_ATOM);
	Times(false,1,9);
	TimeRep(stderr);
	MPI_Finalize();
#endif /* USE_MPI */
	exit(EXIT_SUCCESS);
}	
