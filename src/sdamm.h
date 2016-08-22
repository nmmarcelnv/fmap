#ifndef _SDAMM_H_
#define _SDAMM_H_

#ifdef __cplusplus
extern "C" {
#endif

void cross(double ab[3], double a[3], double b[3]);
void norm(double a[3], double b[3]);
void xryr2rot(double xr[3], double yr[3], double rotm[3][3]);
void sdammTR(double tr[3], double xr[3], double yr[3], int nAtom, double xyzold[][3], double xyznew[][3]);

#ifdef __cplusplus
}
#endif

#endif  /* _SDAMM_H_ */
