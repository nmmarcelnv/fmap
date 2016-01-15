#ifndef _LINKED_H_
#define _LINKED_H_

//follow idiom from typedef struct Stack *Stack_T;
#define T CLINKED
typedef struct T* T;

#ifdef __cplusplus
extern "C" {
#endif

T lnk_create(int num, double xyz[num][3], int nb, double spacing, int pbc, int box[3]);
int lnk_nearest(T lnk, double x, double y, double z, int maxnb, int idx[maxnb], double r2[maxnb]);
void lnk_free(T *lnk);

#ifdef __cplusplus
}
#endif

#undef T

#endif  /* _LINKED_H_ */
