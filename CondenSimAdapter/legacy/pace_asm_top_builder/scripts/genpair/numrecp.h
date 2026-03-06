#ifndef MATRIX_H
#define MATRIX_H




float **init_matrix(int m, int n);
double **init_matrix_d (int m, int n);
void free_matrix (float **ip, int m );
void free_matrix_d (double **ip, int m );

void jacobi( float *a, int n, float *d,float *v, int *nrot);

void eigsrt(float *d,float *v,int n);
void   rotvec (int n,float *v, float *t);

float *vector(int m);
void free_vector(float *ip);
float pythag(float a, float b);
float max (float a, float b);
void svdcmp(float **a, int m,int n,float *w,float **v);
void matrix_mulp(float **a, int ma, int na, float **b, int mb, int nb, float **ip);
//ip=a*b
void matrix_transp(float **a, int m, int n, float **ip);
//at-->ip
void matrix_copy (float **a, int m , int n, float **ip);
//a-->ip

void agaus( double **a, double *b, int n, double *x);

void fit_plane (float *x, int numx,double *nx, double *ny, double *nz, double *p );

double random_gen ();
void fit_line_simp (float *pos_pl, int ncomp, double *a, double *b);
void dist2line(double *a, double *b, double *c, double *sl);
void line_slope_point (double *a, double *b, double *kk);


#endif






