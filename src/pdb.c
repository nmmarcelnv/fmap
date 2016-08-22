#include "pdb.h"

/*checking if a file exists*/
int FileExist(char *pdbfn){
        FILE *file;
        if ((file=fopen(pdbfn, "r")) == NULL){
                fprintf(stderr,"FILE NOT FOUND: %s\n",pdbfn);
                return 0;
        }else{
                fclose(file);
                return 1;
        };
}

/*get number of atoms*/
int CountAtoms(char *pdbfn){
        FILE *file;
        char line[MAXLENLINE];
        int pos;
        if (FileExist(pdbfn)==0){
                exit(EXIT_FAILURE);
        }
        file=fopen(pdbfn, "r");
        pos=0;
        while (fgets(line, MAXLENLINE, file) != NULL){
                if (!strncmp(line, "ATOM",4) || !strncmp(line, "HETATM",6)){
                        pos=pos+1;
                }
        }
        fclose(file);
        return pos;
}

int ReadPdb(char *pdbfn, int n, ATOM atoms[]){
        FILE *file;
        char line[MAXLENLINE];
        if (FileExist(pdbfn)==0){
                exit(EXIT_FAILURE);
        }
        file=fopen(pdbfn, "r");
        int i=0;
        while (fgets(line, MAXLENLINE, file) != NULL){
                if (!strncmp(line, "ATOM",4) || !strncmp(line, "HETATM",6)){
                        strncpy(atoms[i].info,line,INFOLEN_1);
                        atoms[i].info[INFOLEN_1]=0;
                        sscanf(line+30,"%8lf%8lf%8lf",
                                &atoms[i].xyz[0],&atoms[i].xyz[1],&atoms[i].xyz[2]);
                        i++;
                }
        }
        fclose(file);
        i=(i==n)?i:0;
#ifdef DEBUG
        fprintf(stderr,"%s:%d:Number of atoms:%d\n",__FILE__,__LINE__,i);
#endif
        return i;
}

int ReadPqr(char *pdbfn, int n, ATOM atoms[]){
	FILE *file;
        char line[MAXLENLINE];
        if (FileExist(pdbfn)==0){
                exit(EXIT_FAILURE);
        }
        file=fopen(pdbfn, "r");
        int i=0;
        while (fgets(line, MAXLENLINE, file) != NULL){
                if (!strncmp(line, "ATOM",4) || !strncmp(line, "HETATM",6)){
			strncpy(atoms[i].info,line,INFOLEN_1);
                        atoms[i].info[INFOLEN_1]=0;
			sscanf(line+30,"%8lf%8lf%8lf%8lf%8lf%*8s%*16f%*16f%16lf%16lf",
				&atoms[i].xyz[0],&atoms[i].xyz[1],&atoms[i].xyz[2],&atoms[i].q,&atoms[i].r,&atoms[i].Asq,&atoms[i].Bsq);
			i++;
	        }
        }
        fclose(file);
	i=(i==n)?i:0;
#ifdef DEBUG
        fprintf(stderr,"%s\t:%d\t:Number of atoms:\t%d\n",__FILE__,__LINE__,i);
#endif
	return i;
}

void CalCtd(int nAtom, ATOM atoms[], double cen[3]){
        int i,j;
        for (j=0;j<3;j++){
                cen[j]=0.0;
        }
        for (i=0;i<nAtom;i++){
		for (j=0;j<3;j++){
        		cen[j]+=atoms[i].xyz[j];
		}	
        }
        for (j=0;j<3;j++){
                cen[j]/=nAtom;
        }
#ifdef DEBUG
	fprintf(stderr,"%s\t:%d\t:Center:\t%8.3f%8.3f%8.3f\n",__FILE__,__LINE__,cen[0],cen[1],cen[2]);
#endif
}

void ToCtd(int nAtom, ATOM atoms[], double cen[3]){
        int i,j;
	#pragma omp parallel for private(j)
        for (i=0;i<nAtom;i++){
		for (j=0;j<3;j++){
        		atoms[i].xyz[j]-=cen[j];
		}
        }
}

void Unit2dx(int nAtom, ATOM atoms[],const double dx){
	int i,j;
	#pragma omp parallel for private(j)
        for (i=0;i<nAtom;i++){
		for (j=0;j<3;j++){
			atoms[i].xyz[j]/=dx;
		}
	}
}

void ShowPqr(FILE *fp, int nAtom, ATOM atoms[]){
        int i;
        for (i=0;i<nAtom;i++){
                fprintf(fp,"%30s%8.3lf%8.3lf%8.3lf%8.4lf%8.4lf%16.8le%16.8le\n",
                        atoms[i].info,atoms[i].xyz[0],atoms[i].xyz[1],atoms[i].xyz[2],
                        atoms[i].q,atoms[i].r,atoms[i].Asq,atoms[i].Bsq);
        }
}

void ShowPqrdx(FILE *fp, int nAtom, ATOM atoms[], double xyz[][3], double dx){
        int i;
        for (i=0;i<nAtom;i++){
                fprintf(fp,"%30s%8.3lf%8.3lf%8.3lf%8.4lf%8.4lf%16.8le%16.8le\n",
                        atoms[i].info,xyz[i][0]*dx,xyz[i][1]*dx,xyz[i][2]*dx,
                        atoms[i].q,atoms[i].r,atoms[i].Asq,atoms[i].Bsq);
        }
}

void SclRad(int nAtom, ATOM atoms[],const double scl){
        int i;
	#pragma omp parallel for
        for (i=0;i<nAtom;i++){
                        atoms[i].r*=scl;
        }
}
void SclChg(int nAtom, ATOM atoms[],const double scl){
        int i;
	#pragma omp parallel for
        for (i=0;i<nAtom;i++){
                        atoms[i].q*=scl;
        }
}

void SetKap(int nAtom, ATOM atoms[],const double kap){
	int i;
	#pragma omp parallel for
        for (i=0;i<nAtom;i++){
                        atoms[i].kap=kap;
        }
}
