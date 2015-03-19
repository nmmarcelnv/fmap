#include "cputime.h"

#ifdef _OPENMP
        #include <omp.h>
#endif /* _OPENMP */
#ifdef USE_MPI
        #include <mpi.h>
#endif /* USE_MPI */

#define NTIMER 10
static char* time_name[NTIMER]={
	"Init",         //0
	"InitRec1",     //1
	"InitRecs",     //2
	"InitLigs",     //3
	"FFTForw",      //4
	"FFTBack",      //5
	"FFTAdd",       //6
	"MAJOR",        //7
	"MPIComm",      //8
	"Total"};       //9

static double clock_start[NTIMER];
static double clock_sum[NTIMER];
#ifdef _OPENMP
static double time_start[NTIMER];
static double time_sum[NTIMER];
#endif

void timestart(int i){
	double val;
	val = (double) clock() / (double) CLOCKS_PER_SEC;
	clock_start[i] = val;
#ifdef _OPENMP
	val = omp_get_wtime();
	time_start[i] = val;
#endif
}

void timeaccum(int i){
	double val;
	val = (double) clock() / (double) CLOCKS_PER_SEC;
	clock_sum[i] += val - clock_start[i];
#ifdef _OPENMP
        val = omp_get_wtime();
        time_sum[i] += val - time_start[i];
#endif
}

void Times(bool start, int n, ...){
	void (*func)(int);
	func=start?timestart:timeaccum;
	va_list arguments;
	va_start ( arguments, n );
	int i;
	for (i = 0; i < n; i++ ){
        	func(va_arg ( arguments, int )); 
    	}
}

void TimeRep(FILE *fp){
	int i;
	for (i=0;i<NTIMER;i++){
		fprintf(fp,"Clck local %s:\t%16.3f s. (%12.8f %)\n",time_name[i],clock_sum[i],clock_sum[i]/clock_sum[NTIMER-1]*100);
#ifdef _OPENMP
		fprintf(fp,"Time local %s:\t%16.3f s. (%12.8f %)\n",time_name[i],time_sum[i],time_sum[i]/time_sum[NTIMER-1]*100);
#endif
	}
#ifdef USE_MPI
	double clock_glb[NTIMER];
	MPI_Reduce(clock_sum,clock_glb,NTIMER,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#ifdef _OPENMP
	double time_glb[NTIMER];
	MPI_Reduce(time_sum,time_glb,NTIMER,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#endif
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if (rank==0){
		for (i=0;i<NTIMER;i++){
		fprintf(fp,"Clck global %s:\t%16.3f s. (%12.8f %)\n",time_name[i],clock_glb[i],clock_glb[i]/clock_glb[NTIMER-1]*100);
#ifdef _OPENMP
		fprintf(fp,"Time global %s:\t%16.3f s. (%12.8f %)\n",time_name[i],time_glb[i],time_glb[i]/time_glb[NTIMER-1]*100);
#endif
		}
	}
#endif
}
