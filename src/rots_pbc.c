#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include "rand.h"
#define DIM 3
#define MAXCORD 9999
#define INFOLEN 31
#define INFOLEN_1 30
#define MAXGRID 500
#define PI 3.1415926
#define sqr(x) ((x)*(x))
#define cub(x) ((x)*(x)*(x))
int set_coord(char * pdbfn, double xyz[][3], char info[][INFOLEN])
{
        FILE *file;
        char line[MAXLENLINE];
        char word[7];
        if ((file=fopen(pdbfn, "r")) == NULL)
        {
                fprintf(stderr,"FILE NOT FOUND: %s",pdbfn);
                exit(1);
        };
        int pos;
        int i;
        char numstr[9];
        //char testchar[INFOLEN];
        pos=0;

        while (fgets(line, MAXLENLINE, file) != NULL)
        {
                memset(word, 0, 7 * sizeof(char));
                strncpy(word, line, 6);
                if (!strcmp(word, "ATOM  ") || !strcmp(word, "HETATM"))
                {
                        for (i=0;i<3;i++)
                        {
                                memset(numstr, 0, 9 * sizeof(char));
                                strncpy(numstr, line + 30 + 8*i, 8);
                                xyz[pos][i]=atof(numstr);
                        }
                        strncpy(info[pos], line, INFOLEN_1);
                        /*
                        strncpy(testchar, line + 6, INFOLEN_1);
                        printf("%s\n%s\n",testchar,line);
                        */
                        pos=pos+1;
                }
        }
        fclose(file);
        return pos;
}

void calCtd(double xyz[][DIM], int nums, double center[DIM])
{
        int i;
        int j;
        for (j=0; j<DIM; j++)
        {
                center[j]=0;
        }
        for (i=0; i<nums; i++)
        {
                for (j=0; j<DIM; j++)
                {
                        center[j] += xyz[i][j];
                }
        }
        for (j=0; j<DIM; j++)
        {
                center[j]=center[j]/nums;
        }
}

void toCtd(double xyz[][DIM], int nums, double center[DIM], double cenxyz[][DIM]){
        int i,j;
        for (i=0; i<nums; i++)
        {
                for (j=0; j<DIM; j++)
                {
                        cenxyz[i][j]= xyz[i][j]-center[j];
                }
        }
}

double getmaxrad(double xyz[][DIM],int nums){
	int i;
	double maxrad2,rad2i;
	maxrad2=0.0f;
	for (i=0; i<nums; i++){
		rad2i=xyz[i][0]*xyz[i][0]+xyz[i][1]*xyz[i][1]+xyz[i][2]*xyz[i][2];
		if (rad2i>maxrad2){
			maxrad2=rad2i;
		}
	}
	return sqrt(maxrad2);
}
double getminrad(double xyz[][DIM],int nums){
	int i;
	double minrad2,rad2i;
	minrad2=FLT_MAX;
	for (i=0; i<nums; i++){
		rad2i=xyz[i][0]*xyz[i][0]+xyz[i][1]*xyz[i][1]+xyz[i][2]*xyz[i][2];
		if (rad2i<minrad2){
			minrad2=rad2i;
		}
	}
	return sqrt(minrad2);
}

double dist2(double x1,double y1, double z1, double x2,double y2, double z2){
	return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
}

void setbyid(double xyz[][DIM], double curxyz_1[DIM], int curid){
        int j;
        for (j=0; j<DIM; j++){
                curxyz_1[j]=xyz[curid][j];
        }
}

double vol_sphere(double rad){
	return 4.0*PI*rad*rad*rad/3.0;
}

//phi volnume fraction //radiu ratio
double allowed_frac(double phi,double zeta){
	double frac;
	frac=-log(1-phi)\
	     + (3*zeta+3*sqr(zeta)+cub(zeta))*phi/(1-phi)\
	     + (9*sqr(zeta)/2.0+3*cub(zeta))*sqr(phi/(1-phi))\
	     + 3*cub(zeta)*cub(phi/(1-phi));
	return exp(-frac);
}

void zerogrid(int grid1[MAXGRID][MAXGRID][MAXGRID]){
        int i,j,k;
        for (i=1; i<MAXGRID; i++)
                for (j=1; j<MAXGRID; j++)
                        for (k=1; k<MAXGRID; k++)
                                grid1[i][j][k]=0;
}

int pc_map1(int x, int mx){
	if (x<0){
		return mx+x;
	}
	if (x>=mx){
		return x-mx;
	}
	return x;
}

/*Get grid xyz from real xyz, return 1 if out of range of grid*/
int real2grid(double curxyz[3], double gridxyz[3], int nxyz[3], double dxyz[3]){
        int i;
        for (i=0; i<3; i++){
                gridxyz[i] = (int)(curxyz[i]/dxyz[i]);
                if (gridxyz[i]<0 || gridxyz[i]>nxyz[i]){
                        return 1;
                }
        }
        return 0;
}

int pc_fill2grid(int fgrid[MAXGRID][MAXGRID][MAXGRID],double curxyz[3],\
             int nxyz[3],double dxyz[3], double rad, int label){
        double gridxyz[3];
        int i,j,k,xyzmin[3],xyzmax[3];
	int pci,pcj,pck;
        double distx,disty,distz;
        double rad2;
	int nfilled;
	nfilled=0;
	rad2=rad*rad;
        if (real2grid(curxyz,gridxyz,nxyz,dxyz)==1){
                fprintf(stderr,"Out of boundy in function real2grid %8.3f%8.3f%8.3f\n",curxyz[0],curxyz[1],curxyz[2]);
                exit(1);
        }
        for (i=0; i<3; i++){
                xyzmin[i]=floor(gridxyz[i]-rad/dxyz[i]);
                xyzmax[i]=ceil(gridxyz[i]+rad/dxyz[i]);
        }
	//fprintf(stderr,"%8d%8d%8d%8d%8d%8d\n",xyzmin[0],xyzmin[1],xyzmin[2],xyzmax[0],xyzmax[1],xyzmax[2]);
        for (i=xyzmin[0]; i<=xyzmax[0]; i++){
                for (j=xyzmin[1]; j<=xyzmax[1]; j++){
                        for (k=xyzmin[2]; k<=xyzmax[2]; k++){
				//fprintf(stderr,"%8d%8d%8d\n",i,j,k);
                                distx=i*dxyz[0]-curxyz[0];
                                disty=j*dxyz[1]-curxyz[1];
                                distz=k*dxyz[2]-curxyz[2];
                                if ((distx*distx+disty*disty+distz*distz)\
                                                <= rad2){
					pci=pc_map1(i,nxyz[0]);
					pcj=pc_map1(j,nxyz[1]);
					pck=pc_map1(k,nxyz[2]);
					if (fgrid[pci][pcj][pck] == label) continue;	
                                        fgrid[pci][pcj][pck]=label;
					nfilled++;
				}
			}
                }
	}
	return nfilled;
}

long pc_fullgrid(int fgrid[MAXGRID][MAXGRID][MAXGRID],\
                double xyz[MAXCORD][3], int nums,\
                int nxyz[3],double dxyz[3], double radius[MAXCORD],int label){
        int i;
        double curxyz[3];
	long filled;
	filled=0;
        for (i=0; i<nums; i++){
                setbyid(xyz,curxyz,i);
		//fprintf(stderr,"%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",xyz[i][0],xyz[i][1],xyz[i][2],curxyz[0],curxyz[1],curxyz[2]);
                filled += pc_fill2grid(fgrid,curxyz,nxyz,dxyz,radius[i],label);
        }
	return filled;
}

int grid_crash(int grid[MAXGRID][MAXGRID][MAXGRID],int txyz[][3],int num, int mx, int my, int mz){
	int i,ix,iy,iz;
	for (i=0;i<num;i++){
		ix=pc_map1(txyz[i][0]+mx,MAXGRID);
		iy=pc_map1(txyz[i][1]+my,MAXGRID);
		iz=pc_map1(txyz[i][2]+mz,MAXGRID);
		if (grid[ix][iy][iz]==1){
			return 1;
		}
	}
	return 0;
}

long test_grid_nr(int grid[MAXGRID][MAXGRID][MAXGRID], int txyz[][3], int num){
	int i,j,k;
	long nfree=0;
	long nburied=0;
	long ntest=0;
	long npass=0;
	for (i=0;i<MAXGRID;i++){
		for (j=0;j<MAXGRID;j++){
			for (k=0;k<MAXGRID;k++){
				if (grid[i][j][k]==0){
					nfree++;
					continue;
				}
				if (grid[i][j][k]>1){ //1 and 2
					nburied++;
					continue;
				}
				ntest++;
				if (grid_crash(grid,txyz,num,i,j,k)){
					 continue;
				}
				npass++;
			}
		}
	}
	//fprintf(stderr,"%16ld%16ld%16ld%16ld\n",nfree,nburied,ntest,npass);
	return nfree+npass;
}

void showxyz(double xyz[][3], int nums, FILE* fout){
        int i;
        for (i=0;i<nums;i++){
                fprintf(fout,"ATOM  %5d  Na  Na   %5d   %8.3f%8.3f%8.3f\n",\
                                i,i,xyz[i][0],xyz[i][1],xyz[i][2]);
        }
}

void gridize(double xyz[][3], int txyz[][3], int num, double dxyz[3]){
	int i,j;
	for (i=0;i<num;i++){
		for (j=0;j<3;j++){
			txyz[i][j]=(int)(xyz[i][j]/dxyz[j]);
		}	
	}
}

void write_coord(char info[][INFOLEN],double xyz0[][3],\
                int nums, FILE* fp0){
        int i;
        for (i=0; i<nums; i++){
                fprintf(fp0,"%30s%8.3f%8.3f%8.3f\n",info[i],xyz0[i][0],xyz0[i][1],xyz0[i][2]);
        }
}

void trans(double xyz[][3], int nums, double mx, double my,double mz){
	int i;
	for (i=0;i<nums;i++){
		xyz[i][0]+=mx;
		xyz[i][1]+=my;
		xyz[i][2]+=mz;
	}
}

int main(int argc, char **argv){
	double xyz[MAXCORD][3];
	int i;
	double xyz1[MAXCORD][3],xyz2[MAXCORD][3];
	char info1[MAXCORD][INFOLEN];
	int nAtom1;
	double ctd[3];
	double conf[6];
	RANPARM ranp1;
	nAtom1 = set_coord(argv[1], xyz1, info1);
	long seed=atol(argv[2]);
	int starti=atoi(argv[3]);
	ranp1.id=seed;
	ran1(&ranp1);	
	calCtd(xyz1,nAtom1,ctd);
	toCtd(xyz1,nAtom1,ctd,xyz2);
	int nAtom2;
	double xyzc[MAXCORD][3];
	char infoc[MAXCORD][INFOLEN];
	nAtom2 = set_coord(argv[4], xyzc, infoc);
	double scl=atof(argv[5]);
	for (i=0;i<starti;i++){
		ranp2tr(&ranp1,conf,0.0);
	}
	for (i=0;i<nAtom2;i++){
		ranp2tr(&ranp1,conf,0.0);
		genconf(xyz2,nAtom1,conf,xyz);
		trans(xyz,nAtom1,xyzc[i][0]*scl,xyzc[i][1]*scl,xyzc[i][2]*scl);
		write_coord(info1,xyz,nAtom1,stdout);
	}
	exit(0);
}
