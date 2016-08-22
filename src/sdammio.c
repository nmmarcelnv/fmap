#include <stdlib.h>
#include <string.h>
#include "pdb.h"
#include "rand.h"
#include "sdammio.h"

int SDAReadline(char* line,double tr[3], double xr[3], double yr[3]){
	int id=0;
#ifdef DEBUG
	fprintf(stderr,"%s",line);
#endif
	sscanf(line,"%*d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf ",&id,&tr[0],&tr[1],&tr[2],&xr[0],&xr[1],&xr[2],&yr[0],&yr[1],&yr[2]);
	return id;
}

/*get number of atoms in a MODEL*/
int SDACountAtoms(char *pdbfn){
        FILE *file;
        char line[MAXLENLINE];
        int pos,num,id;
        if (FileExist(pdbfn)==0){
                exit(EXIT_FAILURE);
        }
        file=fopen(pdbfn, "r");
        pos=0;
	num=0;
	int preid=0;
	double tr[3],xr[3],yr[3];
        while (fgets(line, MAXLENLINE, file) != NULL){
                if (!strncmp(line, "#",1)){
                        pos=pos+1;
			//##version0.2 sdamm           trajectory  30   1  6 15  1  0  1
			if (pos==2){
				sscanf(line,"%*s %*s %*s %d %*d %*d %*d %*d %*d",&num);	
			}
                }else{
			id=SDAReadline(line,tr,xr,yr);
			if (id>preid){
				preid=id;
			}else{
				break;
			}
		}
        }
	if (num!=preid){
		fprintf(stderr,"%s:%d:Number of atoms:%d != id %d \n",__FILE__,__LINE__,num,id);
	}
        fclose(file);
        return num;
}

int SDAReadxyz(FILE *file, int n, double tr[][3], double xr[][3], double yr[][3]){
        char line[MAXLENLINE];
        int i=0;
	int id=0;
        while (fgets(line, MAXLENLINE, file) != NULL){
                if (strncmp(line, "#",1)){
			id=SDAReadline(line,tr[i],xr[i],yr[i]);
#ifdef DEBUG
			fprintf(stderr,"%f %f %f %f %f %f %f %f %f\n",tr[i][0],tr[i][1],tr[i][2],xr[i][0],xr[i][1],xr[i][2],yr[i][0],yr[i][1],yr[i][2]);
#endif
                        i++;
                }
		if (id==n){
                        break;
                }
        }
        i=(i==n)?i:0;
#ifdef DEBUG
        fprintf(stderr,"%s:%d:Number of atoms:%d\n",__FILE__,__LINE__,i);
#endif
        return i;
}
