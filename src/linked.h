#ifndef _LINKED_H_
#define _LINKED_H_

#ifdef __cplusplus
extern "C" {
#endif

struct linked;
typedef struct linked CLINKED;
CLINKED* lnk_create(double xyz[][3], int num, double cutoff, int pbc, int box[3]);
void lnk_free(CLINKED* lnk);
int lnk_nearest(CLINKED* lnk, double x, double y, double z, int idx[],double r2[]);
void lnk_setnb(CLINKED* lnk,int nb);
#ifdef __cplusplus
}
#endif

#endif  /* _LINKED_H_ */
