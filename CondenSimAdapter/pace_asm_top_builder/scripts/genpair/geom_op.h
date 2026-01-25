#ifndef GEOM_OP_H
#define GEOM_OP_H



#define XX 0
#define YY 1
#define ZZ 2
#define PI 3.1415926


#define VECADD(c,a,b) c[XX]=a[XX]+b[XX];c[YY]=a[YY]+b[YY];c[ZZ]=a[ZZ]+b[ZZ]
#define VECSUB(c,a,b) c[XX]=a[XX]-b[XX];c[YY]=a[YY]-b[YY];c[ZZ]=a[ZZ]-b[ZZ]
#define VECLEN(a)     sqrt(a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ])
#define VECDOT(a,b)   (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ])
#define VECOUT(c,a,b)  c[XX]=a[YY]*b[ZZ]-a[ZZ]*b[YY];c[YY]=a[ZZ]*b[XX]-a[XX]*b[ZZ];c[ZZ]=a[XX]*b[YY]-a[YY]*b[XX]
#define VECELONG(c,a,v)  c[XX]=a[XX]*v;c[YY]=a[YY]*v;c[ZZ]=a[ZZ]*v
#define PRVEC(a)     printf("%f %f %f\n", (&a)[XX], (&a)[YY], (&a)[ZZ])
#define VECCPY(c,a)    c[XX]=a[XX];c[YY]=a[YY];c[ZZ]=a[ZZ]



float fsgn(float a);
float bond(float *ci,float *c1);
float angle(float *ci,float *c1,float *c2);
float dihedral(float *ci,float *c1,float *c2,float *c3);
void calxyz(float *ci, float *c1, float *c2, float *c3, float b, float a, float d);
void  addH_typ0(int *list_A,float *cord, float *hyg);
void addH_typ1(int *list_A, float *cord, float *hyg);
void addH_typ2(int *list_A, float *cord, float *hyg);
void addH_typ4(int *list_A, float *cord, float *hyg);
void addH_typ3(int *list_A, float *cord, float *hyg, float dih);
double cal_NOE(int *list_A, int *list_B,float *cord);


#endif
