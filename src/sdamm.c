#include <math.h>
#include <stdio.h>
#include "sdamm.h"
#include "rand.h"

//maths.f:44: subroutine cross(ab,a,b)	
void cross(double ab[3], double a[3], double b[3]){
	ab[0] = a[1]*b[2] - a[2]*b[1];
	ab[1] = a[2]*b[0] - a[0]*b[2];
	ab[2] = a[0]*b[1] - a[1]*b[0];	
}

void norm(double a[3], double b[3]){
	double s=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	b[0]=a[0]/s;
	b[1]=a[1]/s;
	b[2]=a[2]/s;
}

void xryr2rot(double xr[3], double yr[3], double rotm[3][3]){
	norm(xr,rotm[0]);
	norm(yr,rotm[1]);
	double zr[3];
	//same as in restart.f90:138
	cross(rotm[2],rotm[0],rotm[1]);
	//cross(zr,rotm[0],rotm[1]);
	//norm(zr,rotm[2]);
#ifdef DEBUG
	fprintf(stderr,"rotx: %f %f %f\n",rotm[0][0],rotm[0][1],rotm[0][2]);
	fprintf(stderr,"roty: %f %f %f\n",rotm[1][0],rotm[1][1],rotm[1][2]);
	fprintf(stderr,"rotz: %f %f %f\n",rotm[2][0],rotm[2][1],rotm[2][2]);
#endif
}

void sdammTR(double tr[3], double xr[3], double yr[3], int nAtom, double xyzold[][3], double xyznew[][3]){
	double rotm[3][3];
	xryr2rot(xr,yr,rotm);
	xyzs_t_r(xyzold, nAtom, tr, rotm, xyznew);
}
