#include "common.h"
#include "zdread.h"
#include "linked.h"
#include "drt.h"

void usage(char *prog){
        printf("Usage: %s fmap.out ion tempK eScl vScl\n",prog);
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
	Times(true,3,0,8,9);

#ifdef USE_MPI
	int rank,size;
        MPI_Init(&argc,&argv);
        MPI_Comm_size(MPI_COMM_WORLD,&size);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	Times(false,1,8);
#endif
		
        double ion=atof(argv[2]);
	double tempK=atof(argv[3]);
	double esclrel=atof(argv[4]);
	double vsclrel=atof(argv[5]);
        PARM sys;
        sys.rup=12.0;
	sys.rlow=1.0;
	sys.sdie=78.5;
	sys.kap=(double)GetKappa(ion,(double)sys.sdie,tempK);
	sys.escl=esclrel*332.0/sys.sdie;
	sys.vscl=vsclrel;
	sys.kBT=GetkBT(tempK);
	//sys.kBT=0.591;

	ZDCOM zd;
	int nv,nCrd,nPro;
	int l;
	double dx;
	double lnkc2;
#ifdef USE_MPI
	if (rank==0){
#endif
		nv=ZDread(argv[1],&zd,(void*)(NULL),(void*)(NULL),(void*)(NULL),false);
#ifdef USE_MPI
	}
	Times(true,1,8);
	MPI_Bcast(&nv,1,MPI_INT,0,MPI_COMM_WORLD);
	Times(false,1,8);
#endif
	
	double Angs[nv][3],tran[nv][3],score[nv];
#ifdef USE_MPI
	if (rank==0){
#endif
		ZDread(argv[1],&zd,Angs,tran,score,true);
        	l=zd.N;
		dx=zd.spacing;
		lnkc2=4.0*dx;
		nCrd=CountAtoms(zd.rec);
		nPro=CountAtoms(zd.lig);
#ifdef USE_MPI
	}
	Times(true,1,8);
	MPI_Bcast(&l,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&dx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&lnkc2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&nCrd,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&nPro,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(Angs,nv*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(tran,nv*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
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
#endif
		
	sys.dx=dx;
        ATOM Crds[nCrd];
        ATOM Pros[nPro];
#ifdef USE_MPI
	if (rank==0){
#endif
        	ReadPqr(zd.rec,nCrd,Crds);
        	ToCtd(nCrd,Crds,zd.r);
		SetKap(nCrd,Crds,sys.kap);
        	ReadPqr(zd.lig,nPro,Pros);
		SetKap(nPro,Pros,sys.kap);
        	ToCtd(nPro,Pros,zd.l);
		printf("%d\t%f\n",l,dx);
        	printf("%f\t%f\t%f\n",0.0,0.0,0.0);
        	printf("%s\t%8.3f%8.3f%8.3f\n",zd.rec,zd.r[0],zd.r[1],zd.r[2]);
        	printf("%s\t%8.3f%8.3f%8.3f\n",zd.lig,zd.l[0],zd.l[1],zd.l[2]);
#ifdef USE_MPI
	}
	Times(true,1,8);
	MPI_Bcast(Crds,nCrd,MPI_ATOM,0,MPI_COMM_WORLD);
        MPI_Bcast(Pros,nPro,MPI_ATOM,0,MPI_COMM_WORLD);
	Times(false,1,8);
#endif
	
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

	double ern[4][nv];
	double ernCrd[4][nCrd];
	double ernPro[4][nPro];
	zeroarr(nv*4,&(ern[0][0]));
	zeroarr(nCrd*4,&(ernCrd[0][0]));
	zeroarr(nPro*4,&(ernPro[0][0]));
	Times(false,1,0);
        Times(true,1,7);
	#pragma omp parallel 
	{
	#pragma omp for schedule(dynamic,1) nowait
	for (i=0;i<nv;i++){
		double rot[9];
		ATOM ProCur[nPro];
#ifdef USE_MPI
		if (i%size==rank){
#endif	
			Euler2Rot(Angs[i][0],Angs[i][1],Angs[i][2],rot);
			RotPro(nPro,Pros,ProCur,rot);
			softlnkijk1(nCrd,Crds,nPro,ProCur,l,dx,&sys,lnk,nv,tran,ern,ernCrd,ernPro,i);
#ifdef USE_MPI
		}
#endif	
	}
	}
        Times(false,1,7);
#ifdef USE_MPI
	Times(true,1,8);
	double ernglob[4][nv];
	MPI_Reduce(ern,ernglob,4*nv,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	if (rank==0){
		ijkrep(nCrd,Crds,nPro,Pros,&sys,nv,Angs,tran,score,ernglob,ernCrd,ernPro,stdout);
	}
	MPI_Type_free(&MPI_ATOM);
	Times(false,1,8);
#else
	ijkrep(nCrd,Crds,nPro,Pros,&sys,nv,Angs,tran,score,ern,ernCrd,ernPro,stdout);
#endif
	lnk_free(lnk);
	Times(false,1,9);
        TimeRep(stderr);
#ifdef USE_MPI
	MPI_Finalize();
#endif	
	exit(EXIT_SUCCESS);
}	
