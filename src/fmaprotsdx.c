#include "common.h"
#include "rand.h"
#include "pdbnmr.h"

#ifdef _OPENMP
	fftw_real *recVol, *recR12,*recR6, *recEle;
#endif

void usage(char *prog){
        printf("Usage: %s box.pdb pro.pdb ngrid gridsize ion SclRad SclQ tempK ang.dat eScl vScl Ecut seed PBCScl [rEcut rVcut]\n",prog);
}

int main(int argc, char **argv){
#ifdef DEBUG
    fprintf(stderr,"RUNNING DEBUG BUILD\n");
#endif
	if (argc<15){
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
	char* CenFn=argv[9];
	double esclrel=atof(argv[10]);
	double vsclrel=atof(argv[11]);
	double erncut=atof(argv[12]);
	int seed=atoi(argv[13]);
	double pbcscl=atof(argv[14]);
	double rEcut=12.0;
	double rVcut=12.0;
	if (argc>=16){
		rEcut=atof(argv[15]);
	}
	if (argc>=17){
		rVcut=atof(argv[16]);
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

	int nCrd,nPro,nCen;
#ifdef USE_MPI
	if (rank==0){
#endif /* USE_MPI */
		nCrd=CountAtoms(CrdFn);
		nPro=CountAtoms(ProFn);
                nCen=MDLCountAtoms(CenFn);
#ifdef USE_MPI
        }
	Times(true,1,8);
	MPI_Bcast(&nCrd,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nPro,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nCen,1,MPI_INT,0,MPI_COMM_WORLD);
	
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
        //double Angs[nAng][3];
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

        	ReadPqr(ProFn,nPro,Pros);
		SetKap(nPro,Pros,sys.kap);
		SclRad(nPro,Pros,rscl);
        	CalCtd(nPro,Pros,cenLig);
        	ToCtd(nPro,Pros,cenLig);
		Unit2dx(nPro,Pros,dx);

		printf("%d\t%f\n",l,dx);
        	printf("%f\t%f\t%f\n",0.0,0.0,0.0);
        	printf("%s\t%8.3f%8.3f%8.3f\n",CrdFn,cenRec[0],cenRec[1],cenRec[2]);
        	printf("%s\t%8.3f%8.3f%8.3f\n",ProFn,cenLig[0],cenLig[1],cenLig[2]);

		//SetAng(ang,Angs);
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
	int nb=4001;
	double mat[nb][nb];
#endif
	Times(false,1,0);
	Times(true,1,7);
	FILE *nmrfp;
	nmrfp=fopen(CenFn,"r");
	double cenxyz[nCen][3];
	RANPARM ranp1;
	ranp1.id=seed;
        ran1(&ranp1);
	int j,m;
	double conf[6];
	double ang[3]={0.0};
	for(;;){
#ifdef USE_MPI
		if (i%size==rank){
#endif /* USE_MPI */
			int nc=MDLReadxyz(nmrfp,nCen,cenxyz);
			if (nc!=nCen){
				break;
			}
			int nstep=nPro/nc;
			int start=0;
			double m1xyz[nstep][3],xyz[nstep][3];
			for (j=0;j<nc;j++){
				for (m=0;m<nstep;m++){
					m1xyz[m][0]=Pros[start+m].xyz[0];
					m1xyz[m][1]=Pros[start+m].xyz[1];
					m1xyz[m][2]=Pros[start+m].xyz[2];
				}
				ranp2tr(&ranp1,conf,0.0);
                		genconf(m1xyz,nstep,conf,xyz);
				for (m=0;m<nstep;m++){
                                        xyzLig[start+m][0]=xyz[m][0]+cenxyz[j][0]*pbcscl/dx;
                                        xyzLig[start+m][1]=xyz[m][1]+cenxyz[j][1]*pbcscl/dx;
                                        xyzLig[start+m][2]=xyz[m][2]+cenxyz[j][2]*pbcscl/dx;
                                }
				start+=nstep;
			}
			//ShowPqrdx(stdout,nPro,Pros, xyzLig, dx);
			ProUpdate(lig,nPro,xyzLig);
			Cross(rec,lig,sys,l,sav,ang,nb,mat);
#ifdef USE_MPI
		}
#endif /* USE_MPI */
	}
	fclose(nmrfp);
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
	ProFree(lig);	
	fftw_cleanup_threads();
	Times(false,1,9);
	TimeRep(stderr);
#ifdef USE_MPI
	MPI_Type_free(&MPI_ATOM);
	MPI_Finalize();
#endif /* USE_MPI */
	exit(EXIT_SUCCESS);
}	
