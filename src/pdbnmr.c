#include "pdbnmr.h"

/*get number of atoms in a MODEL*/
int MDLCountAtoms(char *pdbfn){
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
		if (!strncmp(line, "ENDMDL",5)){
			break;
		}
        }
        fclose(file);
        return pos;
}

int MDLReadxyz(FILE *file, int n, double xyz[][3]){
        char line[MAXLENLINE];
        int i=0;
        while (fgets(line, MAXLENLINE, file) != NULL){
                if (!strncmp(line, "ATOM",4) || !strncmp(line, "HETATM",6)){
                        sscanf(line+30,"%8lf%8lf%8lf",
                                &xyz[i][0],&xyz[i][1],&xyz[i][2]);
                        i++;
                }
		if (!strncmp(line, "ENDMDL",5)){
                        break;
                }
        }
        i=(i==n)?i:0;
#ifdef DEBUG
        fprintf(stderr,"%s:%d:Number of atoms:%d\n",__FILE__,__LINE__,i);
#endif
        return i;
}
