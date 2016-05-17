#ifndef _RAND_H_
#define _RAND_H_

#ifdef __cplusplus
extern "C" {
#endif

#define MAXLENLINE 255

typedef struct {
	long id;
	int in;
	int inp;
	long ma[56];
	char cf[MAXLENLINE];
}RANPARM;

double ran1(RANPARM *ranp);
int readconf(FILE *fp0,RANPARM *ranp);
void writeconf(FILE *fp0, RANPARM *ranp);
void ransave(RANPARM *ranp1,RANPARM *ranp2);
void ranp2tr(RANPARM *ranp, double tr[6], double radius);
void ranp2cf(RANPARM *ranp, double cf[6], double radius);
void transf(double costhB, double phiB, double chiB,double rmat[3][3]);
void genconf(double lxyz[][3], int nums, double tranmat[6],double ltrxyz[][3]);

#ifdef __cplusplus
}
#endif

#endif	/* _RAND_H_ */
