#include "rot.h"

/**********************************************/
/* Function rotateAtom  */
/*  rotates around 3 euler angles  */
/**********************************************/
void rotateAtom (double oldX, double oldY, double oldZ,
                 double *newX, double *newY, double *newZ,
                 double psi, double theta, double phi )
{
  double r11, r21, r31, r12, r22, r32, r13, r23, r33;
  r11 = cos(psi)*cos(phi)  -  sin(psi)*cos(theta)*sin(phi);
  r21 = sin(psi)*cos(phi)  +  cos(psi)*cos(theta)*sin(phi);
  r31 = sin(theta)*sin(phi);

  r12 = -cos(psi)*sin(phi)  -  sin(psi)*cos(theta)*cos(phi);
  r22 = -sin(psi)*sin(phi)  +  cos(psi)*cos(theta)*cos(phi);
  r32 = sin(theta)*cos(phi);

  r13 = sin(psi)*sin(theta);
  r23 = -cos(psi)*sin(theta);
  r33 = cos(theta);

  *newX = r11 * oldX + r12 * oldY + r13 * oldZ;
  *newY = r21 * oldX + r22 * oldY + r23 * oldZ;
  *newZ = r31 * oldX + r32 * oldY + r33 * oldZ;

} /* rotateAtom */

void Euler2Rot(double psi, double theta, double phi, double rot[9]){
  double r11, r21, r31, r12, r22, r32, r13, r23, r33;
  r11 = cos(psi)*cos(phi)  -  sin(psi)*cos(theta)*sin(phi);
  r21 = sin(psi)*cos(phi)  +  cos(psi)*cos(theta)*sin(phi);
  r31 = sin(theta)*sin(phi);

  r12 = -cos(psi)*sin(phi)  -  sin(psi)*cos(theta)*cos(phi);
  r22 = -sin(psi)*sin(phi)  +  cos(psi)*cos(theta)*cos(phi);
  r32 = sin(theta)*cos(phi);

  r13 = sin(psi)*sin(theta);
  r23 = -cos(psi)*sin(theta);
  r33 = cos(theta);
  rot[0]=r11;rot[1]=r12;rot[2]=r13;
  rot[3]=r21;rot[4]=r22;rot[5]=r23;
  rot[6]=r31;rot[7]=r32;rot[8]=r33;
}

void RotXYZ(int nAtom, ATOM atoms[], double xyznew[][3], double rot[9]){
	int i;
	#pragma omp parallel for	
	for (i=0;i<nAtom;i++){
		xyznew[i][0]=atoms[i].xyz[0]*rot[0]+atoms[i].xyz[1]*rot[1]+atoms[i].xyz[2]*rot[2];
		xyznew[i][1]=atoms[i].xyz[0]*rot[3]+atoms[i].xyz[1]*rot[4]+atoms[i].xyz[2]*rot[5];
		xyznew[i][2]=atoms[i].xyz[0]*rot[6]+atoms[i].xyz[1]*rot[7]+atoms[i].xyz[2]*rot[8];
	}
}

void RotPro(int nAtom, ATOM atoms[], ATOM pro[], double rot[9]){
        int i;
        //#pragma omp parallel for        
        for (i=0;i<nAtom;i++){
		//pro[i]=atoms[i];
                pro[i].xyz[0]=atoms[i].xyz[0]*rot[0]+atoms[i].xyz[1]*rot[1]+atoms[i].xyz[2]*rot[2];
                pro[i].xyz[1]=atoms[i].xyz[0]*rot[3]+atoms[i].xyz[1]*rot[4]+atoms[i].xyz[2]*rot[5];
                pro[i].xyz[2]=atoms[i].xyz[0]*rot[6]+atoms[i].xyz[1]*rot[7]+atoms[i].xyz[2]*rot[8];
        }
}

int CountAng(char *pdbfn){
        FILE *file;
        char line[MAXLENLINE];
        int pos;
        if (FileExist(pdbfn)==0){
                exit(EXIT_FAILURE);
        }
        file=fopen(pdbfn, "r");
        pos=0;
        while (fgets(line, MAXLENLINE, file) != NULL){
                pos=pos+1;
        }
        return pos;
}

int SetAng(char *pdbfn, double Angs[][3]){
        FILE *file;
        char line[MAXLENLINE];
        int pos;
        if (FileExist(pdbfn)==0){
                exit(EXIT_FAILURE);
        }
        file=fopen(pdbfn, "r");
        pos=0;
        while (fgets(line, MAXLENLINE, file) != NULL){
                sscanf(line,"%lf %lf %lf",&Angs[pos][0],&Angs[pos][1],&Angs[pos][2]);
                pos=pos+1;
        }
        return pos;
}
