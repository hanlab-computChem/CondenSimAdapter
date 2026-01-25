#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "numrecp.h"


#define ROTATE(a,i,j,k,l) g=a[i*n+j];h=a[k*n+l];a[i*n+j]=g-s*(h+g*tau);\
	a[k*n+l]=h+s*(g-h*tau);



float **init_matrix(int m, int n)
//build matrix of a[1..m][1..n]

{

    float **ip;
    int i;
     ip=malloc((m+1)*sizeof(float *));
     for (i=0;i<=m;i++)
     {
       ip[i]=malloc((n+1)*sizeof (float));

      }


     return ip;

}


double **init_matrix_d(int m, int n)
//build matrix of a[1..m][1..n]

{

    double  **ip;
    int i;
     ip=malloc((m+1)*sizeof(double *));
     for (i=0;i<=m;i++)
     {
       ip[i]=malloc((n+1)*sizeof (double));

      }


     return ip;

}




void free_matrix (float **ip, int m )
//free a[1..m][1..n]
{
   int i;

    for (i=0;i<=m; i++)
      {

        free(ip[i]);
      }

     free(ip);

}


void free_matrix_d (double **ip, int m )
//free a[1..m][1..n]
{
   int i;

    for (i=0;i<=m; i++)
      {

        free(ip[i]);
      }

     free(ip);

}



void jacobi( float *a, int n, float *d,float *v, int *nrot)
//matrix is for [0..n-1][0..n-1]
{
	int j,iq,ip,i;
	float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
        float zero=1e-8; 	

	b=malloc(sizeof(float)*n);
	z=malloc(sizeof(float)*n);

	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip*n+iq]=0.0;
		v[ip*n+ip]=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip*n+ip];
		z[ip]=0.0;
	}
	*nrot=0;

//        for (ip=0;ip<n;ip++){
 //          for (iq=0;iq<n;iq++){ printf("%10.3f ",a[ip*n+iq]);}
  //        printf("\n");
   //      }


	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=0;ip<=n-2;ip++) {
			for (iq=ip+1;iq<=n-1;iq++)
				sm += fabs(a[ip*n+iq]);
		}
               
		if (sm < zero) {
			free(z);
			free(b);
     //                   printf("ok");
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<=n-2;ip++) {
			for (iq=ip+1;iq<=n-1;iq++) {
				g=100.0*fabs(a[ip*n+iq]);
				if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
					&& fabs(d[iq])+g == fabs(d[iq]))
					a[ip*n+iq]=0.0;
				else if (fabs(a[ip*n+iq]) > tresh) {
					h=d[iq]-d[ip];
					if (fabs(h)+g == fabs(h))
						t=(a[ip*n+iq])/h;
					else {
						theta=0.5*h/(a[ip*n+iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip*n+iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip*n+iq]=0.0;
					for (j=0;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n-1;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=0;j<=n-1;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=0;ip<=n-1;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
        free(z);
        free(b);
	printf("Too many iterations in routine JACOBI\n");
}

#undef ROTATE


void eigsrt(float *d,float *v,int n)


{
	int k,j,i;
	float p;

	for (i=0;i<n-1;i++) {
                k=i;
		p=d[k];
		for (j=i+1;j<=n-1;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=0;j<=n-1;j++) {
				p=v[j*n+i];
				v[j*n+i]=v[j*n+k];
				v[j*n+k]=p;
			}
		}
	}
}




void   rotvec (int n,float *v, float *t)                                        
{
      float *s;
      s=malloc(n*sizeof(float));  
                            
      int i,j;
                                                               
       for (i=0;i<n;i++){                                                   
         s[i]=v[i];                                                   
         v[i]=0.0;                                                      
      }                                                         

      for ( i=0; i<n;i++){                                                        
         for ( j=0;j<n;j++){                                                    
            v[i]=v[i]+s[j]*t[i*n+j];                                       
         }                                                       
      }                                                          


     free(s);                                                          
}                                                              



float *vector(int m)
//build vec for a[1..m]
{
   float *ip;
    ip=malloc((m+1)*sizeof(float));
    return ip;
}

void free_vector(float *ip)
{
  free(ip);
}




float pythag(float a, float b)
{
  float absa, absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa>absb) return absa*sqrt( 1.0+(absb/absa)*(absb/absa));
  else return (absb==0.0 ? 0.0 : absb* sqrt ( 1.0 + (absa/absb)*(absa/absb))); 

}

float max (float a, float b)
{
   if (a>=b) return a;
   else return b;

}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void svdcmp(float **a, int m,int n,float *w,float **v)
{
	int flag,i,its,j,jj,k,l,nm;
	float c,f,h,s,x,y,z;
	float anorm=0.0,g=0.0,scale=0.0;
        float *rv1;

	if (m < n) {printf("SVDCMP: You must augment A with extra zero rows");exit(1);}
	rv1=malloc((n+1)*sizeof(float));
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i]= a[k][i]/scale;
					s = s+a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				if (i != n) {
					for (j=l;j<=n;j++) {
						for (s=0.0,k=i;k<=m;k++) s =s+ a[k][i]*a[k][j];
						f=s/h;
						for (k=i;k<=m;k++) a[k][j] = a[k][j]+ f*a[k][i];
					}
				}
				for (k=i;k<=m;k++) a[k][i] = a[k][i] *scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale = scale + fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] = a[i][k] / scale;
					s = s+ a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				if (i != m) {
					for (j=l;j<=m;j++) {
						for (s=0.0,k=l;k<=n;k++) s = s + a[j][k]*a[i][k];
						for (k=l;k<=n;k++) a[j][k] = a[j][k] + s*rv1[k];
					}
				}
				for (k=l;k<=n;k++) a[i][k] = a[i][k] * scale;
			}
		}
		anorm=max(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s = s + a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] = v[k][j] + s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=n;i>=1;i--) {
		l=i+1;
		g=w[i];
		if (i < n)
			for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			if (i != n) {
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=m;k++) s = s+ a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
					for (k=i;k<=m;k++) a[k][j] =a[k][j] + f*a[k][i];
				}
			}
			for (j=i;j<=m;j++) a[j][i] = a[j][i] *g;
		} else {
			for (j=i;j<=m;j++) a[j][i]=0.0;
		}
		++(a[i][i]);
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if (fabs(rv1[l])+anorm == anorm) {
					flag=0;
					break;
				}
				if (fabs(w[nm])+anorm == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					if (fabs(f)+anorm != anorm) {
						g=w[i];
						h=pythag(f,g);
						w[i]=h;
						h=1.0/h;
						c=g*h;
						s=(-f*h);
						for (j=1;j<=m;j++) {
							y=a[j][nm];
							z=a[j][i];
							a[j][nm]=y*c+z*s;
							a[j][i]=z*c-y*s;
						}
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
				}
				break;
			}
			if (its == 30) {printf("No convergence in 30 SVDCMP iterations");exit(1);}
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vector(rv1);
}

#undef SIGN





void matrix_mulp(float **a, int ma, int na, float **b, int mb, int nb, float **ip)

//[1..m][1..n]
{
     if (na!=mb) {printf("2nd Dim of A is not equal to 1st Dim of B\n");exit(1);}
 

      int i,j,k;
      float summ;

       for (i=1;i<=ma;i++)
        {
          for (k=1;k<=nb;k++)
           {
            summ=0;
             for (j=1;j<=na;j++)
              {
               summ+=a[i][j]*b[j][k];
              }
             ip[i][k]=summ;
           
            }

         }



}




void  matrix_transp(float **a, int m, int n, float **ip)
//[1..m][1..n]
{
    int i,j;
     for (i=1;i<=m;i++)
      for (j=1;j<=n;j++)
          ip[j][i]=a[i][j];


} 





void matrix_copy (float **a, int m , int n, float **ip)
//[1..m][1..n]
{
    int i,j;
     for (i=1;i<=m;i++)
      for (j=1;j<=n;j++)
       ip[i][j]=a[i][j];


}





void agaus( double **a, double *b, int n, double *x)
//solve linear eq by Gaussian Elimination
// program will terminate if singular value happens
{
   int *js;
   js=malloc((n+1)*sizeof(int));
    int l;
    int i,j,k,is;
    double d,t;

    for(k=1;k<=n-1;k++)
    {
       d=0.0;
       for (i=k;i<=n;i++)
        for (j=k;j<=n;j++)
           if (fabs(a[i][j])>d) 
             {
               d=fabs(a[i][j]);
               js[k]=j;
               is=i;
             }
       if ((d+1.0)==1.0)
         {printf("Singular Problem!\n");exit(1);}
       else
        {
         if (js[k]!=k)
          {
            for (i=1;i<=n;i++)
              {
                t=a[i][k];
                a[i][k]=a[i][js[k]];
                a[i][js[k]]=t; 
              }
          }

         if (is!=k)
          {
            for (j=k;j<=n;j++)
             {
               t=a[k][j];
               a[k][j]=a[is][j];
               a[is][j]=t;
 
             }
            t=b[k];
            b[k]=b[is];
            b[is]=t;

          }


        }

       for(j=k+1;j<=n;j++)
           a[k][j]=a[k][j]/a[k][k];

       b[k]=b[k]/a[k][k];
 
       for (i=k+1;i<=n;i++)
        {
        for (j=k+1;j<=n;j++)
           a[i][j]=a[i][j]-a[i][k]*a[k][j];

         b[i]=b[i]-a[i][k]*b[k];

        }


    }

       if ((fabs(a[n][n])+1.0)==1.0)
         {printf("Singular Problem!\n");exit(1);}

    //backsubstitution
    x[n]=b[n]/a[n][n];

    for (i=n-1;i>=1;i--)
     {
//      if ((fabs(a[i][i])+1.0)==1.0)
  //      {printf("Singular Problem!\n");exit(1);}
      t=0.0;
       for (j=i+1;j<=n;j++)
         t+=a[i][j]*x[j];
       x[i]=(b[i]-t);

     }

    
    for (k=n-1;k>=1;k--)
      if (js[k]!=k)
        {
         t=x[k];
         x[k]=x[js[k]];
         x[js[k]]=t;

        }




  free(js);


}



void fit_plane (float *x, int numx,double *nx, double *ny, double *nz, double *p )
//[0..n-1]
{
   double sxx, syx, szx, sx;
   double sxy, syy, szy, sy;
   double sxz, syz, szz, sz;

      int i;

     sxx=0;    
     for (i=0;i<numx;i++)
       sxx+=x[3*i+0]*x[3*i+0];


     syx=0;    
     for (i=0;i<numx;i++)
       syx+=x[3*i+1]*x[3*i+0];
   

     szx=0;    
     for (i=0;i<numx;i++)
       szx+=x[3*i+2]*x[3*i+0];


     sxy=0;    
     for (i=0;i<numx;i++)
       sxy+=x[3*i+0]*x[3*i+1];

     syy=0;    
     for (i=0;i<numx;i++)
       syy+=x[3*i+1]*x[3*i+1];

     szy=0;    
     for (i=0;i<numx;i++)
       szy+=x[3*i+2]*x[3*i+1];

     sxz=0;    
     for (i=0;i<numx;i++)
       sxz+=x[3*i+0]*x[3*i+2];


     syz=0;    
     for (i=0;i<numx;i++)
       syz+=x[3*i+1]*x[3*i+2];


     szz=0;    
     for (i=0;i<numx;i++)
       szz+=x[3*i+2]*x[3*i+2];


     sx=0;    
     for (i=0;i<numx;i++)
       sx+=x[3*i+0];


     sy=0;    
     for (i=0;i<numx;i++)
       sy+=x[3*i+1];


     sz=0;    
     for (i=0;i<numx;i++)
       sz+=x[3*i+2];


    double **a,*b;
    a=init_matrix_d (3,3);
    b=malloc((3+1)*sizeof(double));

    a[1][1]=sxx;a[1][2]=syx;a[1][3]=szx;b[1]=-sx;
    a[2][1]=sxy;a[2][2]=syy;a[2][3]=szy;b[2]=-sy;
    a[3][1]=sxz;a[3][2]=syz;a[3][3]=szz;b[3]=-sz;
    
    double *xx;
    xx=malloc((3+1)*sizeof(double));

    agaus (a, b, 3, xx);

    double kk;

    kk=sqrt(xx[1]*xx[1]+xx[2]*xx[2]+xx[3]*xx[3]);
     *nx=xx[1]/kk;  *ny=xx[2]/kk; *nz=xx[3]/kk; *p=1/kk;


    free(xx);
    free_matrix_d(a,3);
    free(b);

}

static double rseed=55.0;

double random_gen ()

{
 rseed=rseed-65536*(int)(rseed/65536);
 rseed=rseed*2053+13849;
 rseed=rseed-65536*(int)(rseed/65536);
 return rseed/65536;  
   



}

#define DIM 3

void fit_line_simp (float *pos_pl, int ncomp, double *a, double *b)
{
    int i,j,k;
    double d[3],sl,t;
    for (i=0;i<DIM;i++)
     {a[i]=0;b[i]=0;}
   int  fl,numvec=0;
     for(i=0;i<ncomp;i++)
      for (j=0;j<ncomp;j++)
       if (i!=j)
        {
          sl=0;
          for(k=0;k<DIM;k++) {d[k]=pos_pl[i*DIM+k]-pos_pl[j*DIM+k]; sl+=d[k]*d[k];}
          if (sl>25.0)
            {
              if (numvec==0)
               {
                 sl=sqrt(sl);
                  for (k=0;k<DIM;k++)
                  { d[k]/=sl;b[k]+=d[k];}
                 numvec++;

               }
              else
               {
                  if (d[1]*b[1]+d[2]*b[2]+d[3]*b[3]>0)
                    {
                          sl=sqrt(sl);
                          for (k=0;k<DIM;k++)
                          { d[k]/=sl;b[k]+=d[k];}
                          numvec++;

                    }


               }
            }
        }


       sl=0;fl=0;
       for (i=0;i<DIM;i++)
         if (fabs(b[i])>sl)
          {
             sl=fabs(b[i]);
              fl=i;
          }

        if (sl+1.0==1.0) {printf("Singular Problem!\n");exit(1);}
          sl=b[fl];
         for (i=0;i<DIM;i++) {a[i]=0;b[i]/=sl;}

       double kk[3];
      for (i=0;i<DIM;i++) kk[i]=0;
      for (i=0;i<ncomp; i++)
       {
          for (j=0;j<DIM;j++)
           kk[j]+=pos_pl[i*DIM+j];
       }

      for (i=0;i<DIM;i++) kk[i]/=ncomp;

         t=kk[fl];

       for (i=0;i<DIM;i++)
        if (i!=fl) a[i]=kk[i]-b[i]*t;

//x=a[0]+b[0]*t
//y=a[1]+b[1]*t
//z=a[2]+b[2]*t



}

void line_slope_point (double *a, double *b, double *kk)
//solve line equation with known slope b[DIM] passing know points kk[DIM]
{
double sl,t;
int fl;
int i;

       sl=0;fl=0;
       for (i=0;i<DIM;i++)
         if (fabs(b[i])>sl)
          {
             sl=fabs(b[i]);
              fl=i;
          }

        if (sl+1.0==1.0) {printf("Singular Problem!\n");exit(1);}
          sl=b[fl];
         for (i=0;i<DIM;i++) {a[i]=0;b[i]/=sl;}



         t=kk[fl];

       for (i=0;i<DIM;i++)
        if (i!=fl) a[i]=kk[i]-b[i]*t;

//x=a[0]+b[0]*t
//y=a[1]+b[1]*t
//z=a[2]+b[2]*t



}


void dist2line(double *a, double *b, double *c, double *sl)
//distance from c[0..DIM-1] to line x[i]=a[i]+b[i]*t;

{
   
   double s=0;
   int i;
//   printf("%f \n",s);
//   printf("%f \n",s);
   for (i=0;i<DIM;i++) s=s+b[i]*b[i];
   s=sqrt(fabs(s));
//   printf("%f \n",s);
   double u[DIM];
   for (i=0;i<DIM;i++) u[i]=b[i]/s;

   double aa=0;
   for (i=0;i<DIM;i++)
    aa=aa+u[i]*(c[i]-a[i]);

   double ddd[DIM];
    for (i=0;i<DIM;i++)
     ddd[i]=c[i]-a[i]-aa*u[i];
   (*sl)=sqrt(fabs(ddd[0]*ddd[0]+ddd[1]*ddd[1]+ddd[2]*ddd[2]));

}


#undef DIM
