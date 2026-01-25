#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "geom_op.h"

#define DIM 3



float fsgn(float a)
{
  if (a>0) return 1.0;
  if (a<0) return -1.0;
  return 0.0;

}



float bond(float *ci,float *c1)

{
   float r14[DIM];
   
   VECSUB(r14,ci,c1);
   
   return VECLEN(r14);

}



float angle(float *ci,float *c1,float *c2)
{
   float r14[DIM],r21[DIM];
   float r14v,r21v,cs,ca;

    VECSUB(r14,ci,c1);
    r14v=VECLEN(r14);

    VECSUB(r21,c1,c2);
    r21v=VECLEN(r21);
      
    cs=VECDOT(r14,r21)/r14v/r21v;


    
      ca=atan(sqrt(1-cs*cs)/fabs(cs));
    
    if (cs>0) ca=PI-ca;
      return ca/PI*180;

}



float dihedral(float *ci,float *c1,float *c2,float *c3)

{
    float r[DIM],r14[DIM],r21[DIM],r23[DIM];
    float ra[DIM],rb[DIM],rc[DIM],rd[DIM];
    float cd,r21v,cs,rbv,rdv;

     
    VECSUB(r14, ci,c1);
  
     


   VECSUB(r21,c1,c2);
      

      r21v=VECLEN(r21);

     VECSUB(r23,c3,c2);
      
       VECOUT(r,r21,r23);

       

       cs=VECDOT(r21,r14)/r21v/r21v;
       VECELONG(ra,r21,cs);

      
     
      

       cs=VECDOT(r21,r23)/r21v/r21v;
       VECELONG(rc,r21,cs);

       VECSUB(rb,r14,ra);

       VECSUB(rd,r23,rc);
       
       rdv=VECLEN(rd);
       rbv=VECLEN(rb);

       

       cs=VECDOT(rb,rd)/rdv/rbv;



       cd=atan(sqrt(fabs(1-cs*cs))/fabs(cs));

       if (cs<0) cd=PI-cd;

       cd*=fsgn(VECDOT(r14,r));
       return cd/PI*180;

}



void calxyz(float *ci, float *c1, float *c2, float *c3, float b, float a, float d)
{
   
    float ca,cb,cd;
    float r[DIM],h[DIM];
    float r21[DIM],r23[DIM];
    float rb[DIM],rb21[DIM],rc[DIM];
    float rv,hv,r21v,rbv;


  
   
      ca=a;
      cb=b;
      cd=d;
      ca=ca/180*PI;
      cd=cd/180*PI;
   /*   case ci of
      1:begin
         atm[1].x:=0;
         atm[1].y:=0;
         atm[1].z:=0;

        end;*/

   /*   2:begin
        atm[2].x:=0;
        atm[2].y:=0;
        atm[2].z:=cb;

        end;*/

   /*   3:begin
        c1:=atm[3].zmatrix[1];
        c2:=atm[3].zmatrix[2];

        atm[3].y:=0;

        rz21:=atm[c1].z-atm[c2].z;
        rz21:=rz21/abs(rz21);

        atm[3].z:=atm[c1].z+cos(Pi-ca)*cb*rz21;
        atm[3].x:=sin(ca)*cb


        end;*/

    

        /*     c1:=atm[ci].zmatrix[1];
             c2:=atm[ci].zmatrix[2];
             c3:=atm[ci].zmatrix[3];*/
             
            VECSUB(r21, c1, c2);
            VECSUB(r23, c3, c2);
             
            VECOUT(r, r21, r23);

            VECOUT(h, r21, r);

             

             rv=VECLEN(r);
             hv=VECLEN(h);
             r21v=VECLEN(r21);


             VECELONG(r,r,1/rv);
             
             VECELONG(h,h,1/hv);

             VECELONG(r21,r21,1/r21v);


             rbv=cb*sin(ca);

             VECELONG(rb21, r21, cos(PI-ca)*cb); 
             
             VECELONG(rc, h, cos(PI-fabs(cd)) );
             VECELONG(rb, r, sin(fabs(cd))*fsgn(cd) );
             VECADD (rb, rb, rc);
             VECELONG (rb, rb, rbv);


             VECADD( ci, rb, rb21);
             VECADD( ci, ci, c1);

            

        

      


   







}



void  addH_typ0(int *list_A,float *cord, float *hyg)
{
  
  float *ip;
   ip=&cord[list_A[1]*DIM];
   VECCPY(hyg, ip);

}

void addH_typ1(int *list_A, float *cord, float *hyg)
{
  float *ci,*c1, *c2, *c3;
  float ci1[DIM],ci2[DIM],ci3[DIM];
  float c1v, c2v, c3v;
  c1=&cord[list_A[2]*DIM];
  c2=&cord[list_A[3]*DIM];
  c3=&cord[list_A[4]*DIM];
  ci=&cord[list_A[1]*DIM];

  VECSUB(ci1,c1,ci);
  VECSUB(ci2,c2,ci);
  VECSUB(ci3,c3,ci);
  c1v=VECLEN(ci1);
  c2v=VECLEN(ci2);
  c3v=VECLEN(ci3);

  VECELONG(ci1,ci1,1/c1v);
  VECELONG(ci2,ci2, 1/c2v);
  VECELONG(ci3,ci3, 1/c3v);

  VECADD(ci1,ci1,ci2);
  VECADD(ci1,ci1,ci3);
  VECSUB(hyg,ci,ci1);




}


void addH_typ2(int *list_A, float *cord, float *hyg)
{
float *c1,*c2,*c3;
  c1=&cord[list_A[1]*DIM];
  c2=&cord[list_A[2]*DIM];
  c3=&cord[list_A[3]*DIM];

  if (list_A[4]==1)
    calxyz(hyg, c1, c2, c3, 1.0, 109.5, 120.0);
   else 
    calxyz(hyg, c1, c2, c3, 1.0 , 109.5, -120.0);


}




void addH_typ4(int *list_A, float *cord, float *hyg)
{
float *c1,*c2,*c3;
  c1=&cord[list_A[1]*DIM];
  c2=&cord[list_A[2]*DIM];
  c3=&cord[list_A[3]*DIM];

  
    calxyz(hyg, c1, c2, c3, 1.0, 120.0, 180.0);


}





void addH_typ3(int *list_A, float *cord, float *hyg, float dih)
{
float *c1,*c2,*c3;
  c1=&cord[list_A[1]*DIM];
  c2=&cord[list_A[2]*DIM];
  c3=&cord[list_A[3]*DIM];

  
    calxyz(hyg, c1, c2, c3, 1.0, 109.5, dih);


}







double cal_NOE(int *list_A, int *list_B,float *cord)
//cal a NOE as min distance^(-1/6) between two NOW hydrogens from list_A and list_B


{

int typA, typB;
 typA=list_A[0]; typB=list_B[0];
 
 if ((typA==3) && (typB==3))
 {

   float hygA[3*DIM], hygB[3*DIM];
   addH_typ3(list_A, cord, &hygA[0*DIM] , 180.0);
   addH_typ3(list_A, cord, &hygA[1*DIM] , -60.0);
   addH_typ3(list_A, cord, &hygA[2*DIM] , 60.0);
   addH_typ3(list_B, cord, &hygB[0*DIM] , 180.0);
   addH_typ3(list_B, cord, &hygB[1*DIM] , -60.0);
   addH_typ3(list_B, cord, &hygB[2*DIM] , 60.0);

   float min=10000;
   float ci[DIM],civ,*c1,*c2;
   
   int i,j;
     for(i=0;i<3;i++)
      for (j=0;j<3;j++)
       {
         c1=&hygA[i*DIM];
         c2=&hygB[j*DIM];

         VECSUB(ci,c1,c2);
         civ=VECLEN(ci);

         if (civ<min) {min=civ; }
//       printf("%f \n",min);
       }

  return pow(min,-6);
  }


  if ((typA==3) && (typB!=3))
  {

      float hygA[3*DIM], hygB[DIM];
   addH_typ3(list_A, cord, &hygA[0*DIM] , 180.0);
   addH_typ3(list_A, cord, &hygA[1*DIM] , -60.0);
   addH_typ3(list_A, cord, &hygA[2*DIM] , 60.0);
    if (list_B[0]==0) addH_typ0(list_B, cord, hygB);
    if (list_B[0]==1) addH_typ1(list_B, cord, hygB);
    if (list_B[0]==2) addH_typ2(list_B, cord, hygB);
    if (list_B[0]==4) addH_typ4(list_B, cord, hygB);
   int i;
   float min=10000,ci[DIM],civ,*c1;
   for (i=0;i<3;i++)
    {
     c1=&hygA[i*DIM];
     VECSUB(ci,c1,hygB);
     civ=VECLEN(ci);
     if (civ<min) min=civ;
  //    printf("%f \n",min);
    }


   return pow(min,-6);
  }

  if ((typA!=3) && (typB==3))
  {

      float hygB[3*DIM], hygA[DIM];
   addH_typ3(list_B, cord, &hygB[0*DIM] , 180.0);
   addH_typ3(list_B, cord, &hygB[1*DIM] , -60.0);
   addH_typ3(list_B, cord, &hygB[2*DIM] , 60.0);
    if (list_A[0]==0) addH_typ0(list_A, cord, hygA);
    if (list_A[0]==1) addH_typ1(list_A, cord, hygA);
    if (list_A[0]==2) addH_typ2(list_A, cord, hygA);
    if (list_A[0]==4) addH_typ4(list_A, cord, hygA);
   int i;
   float min=10000,ci[DIM],civ,*c1;
   for (i=0;i<3;i++)
    {
     c1=&hygB[i*DIM];
     VECSUB(ci,c1,hygA);
     civ=VECLEN(ci);
     if (civ<min) min=civ;
    //       printf("%f \n",min);
    }

   

  return pow(min,-6);
  }

      float hygB[DIM], hygA[DIM];
    if (list_A[0]==0) addH_typ0(list_A, cord, hygA);
    if (list_A[0]==1) addH_typ1(list_A, cord, hygA);
    if (list_A[0]==2) addH_typ2(list_A, cord, hygA);
    if (list_A[0]==4) addH_typ4(list_A, cord, hygA);
    if (list_B[0]==0) addH_typ0(list_B, cord, hygB);
    if (list_B[0]==1) addH_typ1(list_B, cord, hygB);
    if (list_B[0]==2) addH_typ2(list_B, cord, hygB);
    if (list_B[0]==4) addH_typ4(list_B, cord, hygB);
   float ci[DIM],civ;
   VECSUB(ci,hygA,hygB);
   civ=VECLEN(ci);
   
  return pow(civ,-6);



}











#undef DIM

