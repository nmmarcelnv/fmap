#include "common.h"

void usage(char *prog){
        printf("Usage: echo escl vscl|%s mat.bin nbin dx lowboundary tempK\n",prog);
}

void readmat(char*fn, int l, double mat[l][l]){
	if (FileExist(fn)==0){
                exit(EXIT_FAILURE);
        }
	FILE *fp;
	fp=fopen(fn, "rb");
	fread(mat,sizeof(double),l*l,fp);
	fclose(fp);
}

double matscl(int l, double mat[l][l], double dx, double lowb, double kBT, double escl, double vscl){
	double bzs=0.0;
	int i,j;
	for (i=0;i<l;i++){
		double ei=dx*i+lowb;
		for (j=0;j<l;j++){
			double vi=dx*j+lowb;
			if (mat[i][j]>0.5){
				double val=exp((ei*escl+vi*vscl)/-kBT)*mat[i][j];
				bzs+=val;
#ifdef DEBUG
			printf("%8.3f%8.3f%8.3f%16d%16e\n",ei,vi,ei+vi,mat[i][j],val);
#endif
			}
		}
	}
	return bzs;
}

int main(int argc, char **argv){
#ifdef DEBUG
    fprintf(stderr,"RUNNING DEBUG BUILD\n");
#endif
	if (argc<6){
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }
 
        char* Fn=argv[1];
        int l=atoi(argv[2]);
        double dx=atof(argv[3]);
        double lowb=atof(argv[4]);
	double tempK=atof(argv[5]);
	double mat[l][l];
	double escl=1.0;
	double vscl=1.0;

	double kBT=GetkBT(tempK);
	readmat(Fn,l,mat);
	while (2==scanf("%lf %lf",&escl,&vscl)){
		double s=matscl(l,mat,dx,lowb,kBT,escl,vscl);
		printf("#%8.3f%8.3f%16e\n",escl,vscl,s);
	}
	exit(EXIT_SUCCESS);
}	
