#include "common.h"

#ifdef _OPENMP
	fftw_real *recVol, *recR12,*recR6, *recEle;
#endif

void usage(char *prog){
        printf("Usage: %s box.pdb profns ngrid gridsize ion SclRad SclQ tempK eScl vScl Ecut [rEcut rVcut]\n",prog);
}

int main(int argc, char **argv){
#ifdef DEBUG
    fprintf(stderr,"RUNNING DEBUG BUILD\n");
#endif
	if (argc<12){
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
        char* ProFns=argv[2];
        int l=atoi(argv[3]);
        double dx=atof(argv[4]);
        double ion=atof(argv[5]);
	double rscl=atof(argv[6]);
	double qscl=atof(argv[7]);
	double tempK=atof(argv[8]);
	//char* ang=argv[9];
	double esclrel=atof(argv[9]);
	double vsclrel=atof(argv[10]);
	double erncut=atof(argv[11]);
	double rEcut=12.0;
	double rVcut=12.0;
	if (argc>=13){
		rEcut=atof(argv[12]);
	}
	if (argc>=14){
		rVcut=atof(argv[13]);
	}

	fftw_init_threads();
	fftw_plan_with_nthreads(omp_get_max_threads());
//#ifdef _OPENMP
//	printf("TACC Nthread:%d\n",omp_get_max_threads());
//	printf("TACC Size:%d\n",size);
//	printf("TACC rank:%d\n",rank);
//#endif

	PARM sys;	
	sys.rup=rVcut;
	sys.rEup=rEcut;
	sys.rlow=1.0;
	sys.dx=dx;
	sys.sdie=GetWaterDie(tempK);
	sys.kap=GetKappa(ion,sys.sdie,tempK);
	sys.escl=esclrel*332.0/sys.sdie;
	sys.vscl=vsclrel;
	sys.kBT=GetkBT(tempK);
	//sys.kBT=0.591;
	sys.erncut=erncut; //-1.0*sys.kBT*(3.0*log((double)(l))-log(10.0));

	int nCrd,nFns;
#ifdef USE_MPI
	if (rank==0){
#endif /* USE_MPI */
		nCrd=CountAtoms(CrdFn);
                nFns=CountAng(ProFns);
#ifdef USE_MPI
        }
	Times(true,1,8);
	MPI_Bcast(&nCrd,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nFns,1,MPI_INT,0,MPI_COMM_WORLD);
	
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
        char Fns[nFns][MAXLENLINE];
#ifdef USE_MPI
	if (rank==0){
#endif /* USE_MPI */
		double cenRec[3]= { 0.0 };
		double cenLig[3]= { 0.0 };
        	ReadPqr(CrdFn,nCrd,Crds);
		SetKap(nCrd,Crds,sys.kap);
		SclRad(nCrd,Crds,rscl);
		SclChg(nCrd,Crds,qscl);
        	CalCtd(nCrd,Crds,cenRec);
        	ToCtd(nCrd,Crds,cenRec);
		Unit2dx(nCrd,Crds,dx);

		printf("%d\t%f\n",l,dx);
        	printf("%f\t%f\t%f\n",0.0,0.0,0.0);
        	printf("%s\t%8.3f%8.3f%8.3f\n",CrdFn,cenRec[0],cenRec[1],cenRec[2]);
        	printf("%s\t%8.3f%8.3f%8.3f\n",ProFns,cenLig[0],cenLig[1],cenLig[2]);

		SetFns(ProFns,Fns);
#ifdef USE_MPI
	}
	Times(false,1,8);
	MPI_Bcast(Crds,nCrd,MPI_ATOM,0,MPI_COMM_WORLD);
	MPI_Bcast(Fns,nFns*MAXLENLINE,MPI_CHAR,0,MPI_COMM_WORLD);
	Times(false,1,8);
#endif /* USE_MPI */

	PRO *rec;
	rec=ProCreate(nCrd,Crds);
	fftw_real sav[l][l][2*(l/2+1)];
	SetZero(l,sav);
	
	int i;
#ifdef _OPENMP
	fftw_real extVol[l][l][2*(l/2+1)],extEle[l][l][2*(l/2+1)];
	fftw_real extR12[l][l][2*(l/2+1)],extR6[l][l][2*(l/2+1)];
	recVol=&(extVol[0][0][0]);
	recEle=&(extEle[0][0][0]);
	recR12=&(extR12[0][0][0]);
	recR6=&(extR6[0][0][0]);
	int nb=4001;
	double mat[nb][nb];
#endif
	Times(false,1,0);
	Times(true,1,7);
	for (i=0;i<nFns;i++){
#ifdef USE_MPI
		if (i%size==rank){
#endif /* USE_MPI */
			int nPro;
			double cenLig[3]= { 0.0 };
			//printf("Fns[i]:%s\n",Fns[i]);
			nPro=CountAtoms(Fns[i]);
        		ATOM Pros[nPro];
        		ReadPqr(Fns[i],nPro,Pros);
			SetKap(nPro,Pros,sys.kap);
			SclRad(nPro,Pros,rscl);
        		CalCtd(nPro,Pros,cenLig);
        		ToCtd(nPro,Pros,cenLig);
			Unit2dx(nPro,Pros,dx);
			PRO *lig;
			lig=ProCreate(nPro,Pros);
			double ang[3]={0.0};
			Cross(rec,lig,sys,l,sav,ang,nb,mat);
			ProFree(lig);	
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
		FILE* mfp=fopen("mat.bin","wb");
		fwrite(mat,sizeof(double),nb*nb,mfp);
		fclose(mfp);
#ifdef USE_MPI
	}
#endif /* USE_MPI */
	ProFree(rec);
	fftw_cleanup_threads();
	Times(false,1,9);
	TimeRep(stderr);
#ifdef USE_MPI
	MPI_Type_free(&MPI_ATOM);
	MPI_Finalize();
#endif /* USE_MPI */
	exit(EXIT_SUCCESS);
}	
