//general header for open a data file and construct its parse according a PDB file
//argv1:data file  argv2 : num of res  argv3: PDB file


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "string_op.h"
#include "file_op.h"
#include "parse_cg.h"
#include "geom_op.h"


#define DIM 3
#define MAXBUF 5000

int main (int agrc, char * argv[])
{

int list_A[MAXNOE*MAX_LIST], list_B[MAXNOE*MAX_LIST];
double dis[MAXNOE];
int natm, nres, nNOE, nframe;
FILE *fp,*fdat;
char buffer[MAXBUF],buf[MAXBUF];
char **resnm;
natm =str2int(argv[1]);
int parseMod;
int i,j,k;
parseMod=str2int(argv[4]);
nres=str2int(argv[2]);
int *HD21,*HD22,*ND2,*HD1,*ND1,*OD1,*OD2,*OG1, *CD1,*CD2, *OE, *NE, *AUX,*OG,*CG,*CG1,*CG2,*H, *CA,  *C, *NH, *CB, *O;

resnm=malloc((nres+1)*sizeof(char*));
for (i=0;i<nres+1;i++) resnm[i]=malloc(11*sizeof(char));

CA=malloc ((nres+1)*sizeof(int));
C=malloc ((nres+1)*sizeof(int));
NH=malloc ((nres+1)*sizeof(int));
CB=malloc ((nres+1)*sizeof(int));
O=malloc ((nres+1)*sizeof(int));
H=malloc ((nres+1)*sizeof(int));
CG1=malloc ((nres+1)*sizeof(int));
CG2=malloc ((nres+1)*sizeof(int));
CG=malloc ((nres+1)*sizeof(int));
OG=malloc ((nres+1)*sizeof(int));
CD1=malloc ((nres+1)*sizeof(int));
CD2=malloc ((nres+1)*sizeof(int));
OE=malloc ((nres+1)*sizeof(int));
NE=malloc ((nres+1)*sizeof(int));
AUX=malloc ((nres+1)*sizeof(int));
OG1=malloc ((nres+1)*sizeof(int));
OD1=malloc ((nres+1)*sizeof(int));
OD2=malloc ((nres+1)*sizeof(int));
HD1=malloc ((nres+1)*sizeof(int));
ND1=malloc ((nres+1)*sizeof(int));
ND2=malloc ((nres+1)*sizeof(int));
HD21=malloc ((nres+1)*sizeof(int));
HD22=malloc ((nres+1)*sizeof(int));



struct atm_P **pres, **atm_idx;

pres=malloc((nres+1)*sizeof(struct atm_P *));
atm_idx=malloc(natm*sizeof(struct atm_P *));

int t_resid=0;
int last_resid=-1;

fp=fopen(argv[3],"r");
for (i=0;i<natm;i++)
{
get_line(fp,buffer);
trcspace(buffer);
parstr(buffer,1,buf);

while (strcmp(buf,"ATOM")!=0)
{
get_line(fp,buffer);
trcspace(buffer);
parstr(buffer,1,buf);



}

parstr(buffer,5, buf);
if ((buf[0]<'0')||(buf[0]>'9')) parstr(buffer,6,buf); /*in case of chain identifier*/
j=str2int(buf);
if (j!=last_resid)
{
  last_resid=j;
  t_resid++;
}
parstr(buffer,4,buf);
strcpy(resnm[t_resid],buf);

}


fclose(fp);

fp=fopen(argv[3],"r");
construct_parse(natm, nres, pres, atm_idx, fp,parseMod);
/*parseMod
  1  for  no halt for error
  0  for sensative mod
*/
 
fclose(fp);


 for (i=1;i<=nres;i++)
 {
  CA[i]= find_atm_noerr (pres[i], "CA" );
  C[i]= find_atm_noerr (pres[i], "C");
  O[i]= find_atm_noerr (pres[i], "O");
  NH[i]= find_atm_noerr (pres[i], "N");
  H[i]= find_atm_noerr (pres[i], "H");
  CB[i]= find_atm_noerr (pres[i], "CB");
  CG[i]= find_atm_noerr (pres[i], "CG");
  CG1[i]= find_atm_noerr (pres[i], "CG1");
  CG2[i]= find_atm_noerr (pres[i], "CG2");
  OG[i]= find_atm_noerr (pres[i], "OG");
  OE[i]= find_atm_noerr (pres[i], "OE1");
  NE[i]= find_atm_noerr (pres[i], "NE2");
  CD1[i]= find_atm_noerr (pres[i], "CD1");
  CD2[i]= find_atm_noerr (pres[i], "CD2");
  AUX[i]= find_atm_noerr (pres[i], "AUXL");
  OG1[i]= find_atm_noerr (pres[i], "OG1");
  OD1[i]= find_atm_noerr (pres[i], "OD1");
  OD2[i]= find_atm_noerr (pres[i], "OD2");
  HD1[i]= find_atm_noerr (pres[i], "HD1");
  ND1[i]= find_atm_noerr (pres[i], "ND1");
  ND2[i]= find_atm_noerr (pres[i], "ND2");
  HD21[i]= find_atm_noerr (pres[i], "HD21");
  HD22[i]= find_atm_noerr (pres[i], "HD22");


  printf("%s CA:%d C:%d O:%d N:%d H:%d CB:%d CG1:%d CG2:%d CG:%d OG:%d OG1:%d OE:%d NE:%d CD1:%d CD2:%d OD1:%d OD2:%d HD1:%d ND1:%d ND2:%d\n", resnm[i],CA[i], C[i], O[i], NH[i], H[i],CB[i],CG1[i],CG2[i],CG[i],OG[i],OG1[i],OE[i],NE[i],CD1[i],CD2[i],OD1[i],OD2[i],HD1[i],ND1[i],ND2[i]);
  }

printf ("[ exclusions ]\n");

  for (i=1;i<nres;i++)
  {


    if (i<=nres-4)
     { 
      
      if ((CA[i]>=0) && (CB[i+4]>=0))
       { 
         printf("%d    %d\n", CA[i]+1,  CB[i+4]+1);
         
       }
       
      if ((CA[i]>=0) && (CA[i+4]>=0)&&(strcmp(resnm[i+4],"GLY")==0))
       {
         printf("%d    %d   ;  0.671  0.50 GLY nm\n", CA[i]+1,  CA[i+4]+1);

       }


     }

    if (i<=nres-3)
     { 
      if ((strcmp(resnm[i+3],"PRO")!=0)&&(strcmp(resnm[i+3],"DPR")!=0))
      if ((O[i]>=0) && (NH[i+3]>=0))
       {
         printf("%d    %d  ;  Oi-Ni+3 \n", O[i]+1,  NH[i+3]+1);

       }



     }




    if (i<=nres-4)
     { 
      if ((strcmp(resnm[i+4],"PRO")!=0)&&(strcmp(resnm[i+4],"DPR")!=0))
      if ((O[i]>=0) && (NH[i+4]>=0))
       {
         printf("%d    %d  ;  Oi-Ni+4 \n", O[i]+1,  NH[i+4]+1);
  
       }

  

     } 






// correct for gamma particles

    if ((OG[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d\n", OG[i]+1,  H[i]+1);

    }
    if ((OG[i]>=0) && (O[i]>=0))
    {
      printf("%d    %d\n", OG[i]+1,  O[i]+1);

    }
    if ((OG[i]>=0) && (NH[i+1]>=0))
    {
      printf("%d    %d\n", OG[i]+1,  NH[i+1]+1);

    }


    if ((OG[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d\n", OG[i]+1,  H[i+1]+1);

    }



    if (i>1)
    {

     if ((OG[i]>=0) && (O[i-1]>=0))
     {
       printf("%d    %d\n", OG[i]+1,  O[i-1]+1);

     }



    }



// og1
      if ((OG1[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d\n", OG1[i]+1,  H[i]+1);

    }
    if ((OG1[i]>=0) && (O[i]>=0))
    {
      printf("%d    %d\n", OG1[i]+1,  O[i]+1);

    }
    if ((OG1[i]>=0) && (NH[i+1]>=0))
    {
      printf("%d    %d\n", OG1[i]+1,  NH[i+1]+1);

    }


    if ((OG1[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d\n", OG1[i]+1,  H[i+1]+1);

    }



    if (i>1)
    {

     if ((OG1[i]>=0) && (O[i-1]>=0))
     {
       printf("%d    %d\n", OG1[i]+1,  O[i-1]+1);

     }



    }



 if ((strncmp(resnm[i],"ASP",3)==0)||(strcmp(resnm[i],"ASN")==0))
 {

    if ((CG[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d ;  CG-Hi asx\n", CG[i]+1,  H[i]+1);
      
    }
    if ((CG[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d   ;  CG-Hi+1 asx\n", CG[i]+1,  H[i+1]+1);
      
    }
    
    if ((CG[i]>=0) && (NH[i+1]>=0))
    {
      printf("%d    %d    ;  CG-Ni+1 asx \n", CG[i]+1,  NH[i+1]+1);
       
    }

    if ((ND2[i]>=0) && (NH[i]>=0))
    {
      printf("%d    %d    ;  ND2-Ni asx \n", ND2[i]+1,  NH[i]+1);

    }

    if ((ND2[i]>=0) && (NH[i+1]>=0))
    {
      printf("%d    %d    ;  ND2-Ni+1 asx \n", ND2[i]+1,  NH[i+1]+1);

    }


  if (strcmp(resnm[i],"ASN")==0)
{
    if ((HD21[i]>=0) && (O[i]>=0))
    {
      printf("%d    %d    ;  HD21-Oi asx\n", HD21[i]+1,  O[i]+1);

    }

    if ((OD1[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d    ;  OD1-Hi asx\n", OD1[i]+1,  H[i]+1);

    }

    if ((OD1[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d    ;  OD1-Hi+1 asx\n", OD1[i]+1,  H[i+1]+1);

    }




}


  if (strncmp(resnm[i],"ASP",3)==0)
   {


    if ((OD1[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d ;  OD1-Hi asx\n", OD1[i]+1,  H[i]+1);

    }

    if ((OD1[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d   ;  OD1-Hi+1 asx\n", OD1[i]+1,  H[i+1]+1);

    }

    if ((OD2[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d ;  OD2-Hi asx\n", OD2[i]+1,  H[i]+1);

    }

    if ((OD2[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d   ;  OD2-Hi+1 asx\n", OD2[i]+1,  H[i+1]+1);

    }

    if ((CG[i]>=0) && (O[i]>=0))
    {
      printf("%d    %d   ;  CG-Oi asx\n", CG[i]+1,  O[i]+1);

    }





   }


 }

 else

 if ((strcmp(resnm[i],"PHE")==0)||(strcmp(resnm[i],"TYR")==0)||(strcmp(resnm[i],"TRP")==0)||(strncmp(resnm[i],"HIS",3)==0))
  {


    if ((CG[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d\n", CG[i]+1,  H[i]+1);

    }
    if ((CG[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d\n", CG[i]+1,  H[i+1]+1);

    }


    if ((CG[i]>=0) && (NH[i+1]>=0))
    {
      printf("%d    %d\n", CG[i]+1,  NH[i+1]+1);

    }
    if ((CD1[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d\n", CD1[i]+1,  H[i]+1);

    }

    if ((ND1[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d\n", ND1[i]+1,  H[i]+1);

    }


    if ((CD2[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d\n", CD2[i]+1,  H[i]+1);

    }



    if ((CD1[i]>=0) && (NH[i]>=0))
    {
      printf("%d    %d\n", CD1[i]+1,  NH[i]+1);

    }

    if ((ND1[i]>=0) && (NH[i]>=0))
    {
      printf("%d    %d\n", ND1[i]+1,  NH[i]+1);
     
    }


    if ((CD2[i]>=0) && (NH[i]>=0))
    {
      printf("%d    %d\n", CD2[i]+1,  NH[i]+1);

    }


    if (i>1)
    {

     if ((CD1[i]>=0) && (O[i-1]>=0))
     {
       printf("%d    %d\n", CD1[i]+1,  O[i-1]+1);

     }



    }


    if (i>1)
    {

     if ((CD2[i]>=0) && (O[i-1]>=0))
     {
       printf("%d    %d\n", CD2[i]+1,  O[i-1]+1);

     }



    }

    if ((CD1[i]>=0) && (NH[i+1]>=0))
    {
      printf("%d    %d\n", CD1[i]+1,  NH[i+1]+1);

    }

    if ((ND1[i]>=0) && (NH[i+1]>=0))
    {
      printf("%d    %d\n", ND1[i]+1,  NH[i+1]+1);

    }


    if ((CD2[i]>=0) && (NH[i+1]>=0))
    {
      printf("%d    %d\n", CD2[i]+1,  NH[i+1]+1);

    }


   }
 else
 {
    if ((CG[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d\n", CG[i]+1,  H[i]+1);

    }
    if ((CG[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d\n", CG[i]+1,  H[i+1]+1);

    }




 }


    if ((CG1[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d\n", CG1[i]+1,  H[i]+1);

    }
    if ((CG1[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d\n", CG1[i]+1,  H[i+1]+1);

    }
  




    if ((CG2[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d\n", CG2[i]+1,  H[i]+1);

    }
    if ((CG2[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d\n", CG2[i]+1,  H[i+1]+1);

    }


  











   if ((NH[i]>=0) && (C[i+1]>=0))
    {
      printf("%d    %d\n", NH[i]+1,  C[i+1]+1);     

    }

    if ((O[i]>=0) && (O[i+1]>=0))
    {
      printf("%d    %d\n", O[i]+1,  O[i+1]+1);

    }




    if ((NH[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d\n", NH[i]+1,  H[i+1]+1);
    
    }


   if ((NH[i]>=0) && (O[i+1]>=0))
    {
      printf("%d    %d\n", NH[i]+1,  O[i+1]+1); 

    }

   if ((CA[i]>=0) && (O[i+1]>=0))
    {
      printf("%d    %d\n", CA[i]+1,  O[i+1]+1); 

    }

   if ((C[i]>=0) && (O[i+1]>=0))
    {
      printf("%d    %d\n", C[i]+1,  O[i+1]+1); 

    }

   if ((O[i]>=0) && (C[i+1]>=0))
    {
      printf("%d    %d\n", O[i]+1,  C[i+1]+1); 

    }


   if ((CA[i]>=0) && (C[i+1]>=0))
    {
      printf("%d    %d\n", CA[i]+1,  C[i+1]+1);

    }



   if (i<nres-1)
   {

   if ((NH[i]>=0) && (C[i+2]>=0))
    {
      printf("%d    %d\n", NH[i]+1,  C[i+2]+1); 

    }


   if ((NH[i]>=0) && (O[i+2]>=0))
    {
      printf("%d    %d\n", NH[i]+1,  O[i+2]+1); 

    }

   if ((CA[i]>=0) && (O[i+2]>=0))
    {
      printf("%d    %d\n", CA[i]+1,  O[i+2]+1); 

    }

   if ((C[i]>=0) && (O[i+2]>=0))
    {
      printf("%d    %d\n", C[i]+1,  O[i+2]+1); 

    }

   if ((C[i]>=0) && (NH[i+2]>=0))
    {
      printf("%d    %d\n", C[i]+1,  NH[i+2]+1); 

    }


   if ((O[i]>=0) && (CA[i+2]>=0))
    {
      printf("%d    %d\n", O[i]+1,  CA[i+2]+1); 

    }

   if ((O[i]>=0) && (NH[i+2]>=0))
    {
      printf("%d    %d\n", O[i]+1,  NH[i+2]+1); 

    }


   if ((CA[i]>=0) && (CA[i+2]>=0))
    {
      printf("%d    %d\n", CA[i]+1,  CA[i+2]+1);

    }


   if ((CA[i]>=0) && (C[i+2]>=0))
    {
      printf("%d    %d\n", CA[i]+1,  C[i+2]+1);

    }


   if ((C[i]>=0) && (CA[i+2]>=0))
    {
      printf("%d    %d\n", C[i]+1,  CA[i+2]+1);

    }


   if ((C[i]>=0) && (C[i+2]>=0))
    {
      printf("%d    %d\n", C[i]+1,  C[i+2]+1);

    }



    if (i>1)
      {
       if ((C[i-1]>=0) && (O[i+2]>=0))
       {
         printf("%d    %d\n", C[i-1]+1,  O[i+2]+1);

       }
      }







   }





  }

  printf ("[ pairs ]\n");


  for (i=1;i<nres;i++)
  {




//cai-->cbi+4   0.671  0.47

    if (i<=nres-4)
     {

      if ((CA[i]>=0) && (CB[i+4]>=0))
       {
         printf("%d    %d  1  2.89314E-02 3.11858E-04  ;  0.671  0.47 nm\n", CA[i]+1,  CB[i+4]+1);

       }


      if ((CA[i]>=0) && (CA[i+4]>=0)&&(strcmp(resnm[i+4],"GLY")==0))
       { 
         printf("%d    %d  1  4.19375E-02  6.55273E-04  ;  0.671  0.50 GLY nm\n", CA[i]+1,  CA[i+4]+1);
       
       }




     }






//reinforce Oi-Ni+3 HB for turn

    if (i<=nres-3)
     {
      if ((strcmp(resnm[i+3],"PRO")!=0)&&(strcmp(resnm[i+3],"DPR")!=0))
      if ((O[i]>=0) && (NH[i+3]>=0))
       {
         printf("%d    %d  1  1.123685E-02    2.147396E-06 ;  Oi-Ni+3 \n", O[i]+1,  NH[i+3]+1);

       }



     }




    if (i<=nres-4)
     {
      if ((strcmp(resnm[i+4],"PRO")!=0)&&(strcmp(resnm[i+4],"DPR")!=0))
      if ((O[i]>=0) && (NH[i+4]>=0))
       {
         printf("%d    %d  1  1.123685E-02    2.147396E-06 ;  Oi-Ni+4 \n", O[i]+1,  NH[i+4]+1);

       }



     }






// correct for gamma particles

    if ((OG[i]>=0) && (H[i]>=0))
    { 
      printf("%d    %d  1  1.13380E-03  1.28550E-07 ; 2.5 0.22 nm OG-H \n", OG[i]+1,  H[i]+1);
    
    }
    if ((OG[i]>=0) && (O[i]>=0))
    { 
      printf("%d    %d  1 7.74841E-04   1.50095E-06 ; 2.5  0.27 nm OG-Oi  rep\n", OG[i]+1,  O[i]+1);
    
    }
    if ((OG[i]>=0) && (NH[i+1]>=0))
    { 
      printf("%d    %d  1  3.83970E-03  4.12285E-06 ;   OG-Ni+1 \n", OG[i]+1,  NH[i+1]+1);

    }


    if ((OG[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d  1  1.13380E-03  1.28550E-07  ;  Og-Hi+1\n", OG[i]+1,  H[i+1]+1);

    }
    
    
      
    if (i>1)
    {

     if ((OG[i]>=0) && (O[i-1]>=0))
     { 
       printf("%d    %d  1  3.87420E-03   1.50095E-06 ;  Og-Oi-1 \n", OG[i]+1,  O[i-1]+1);
      
     }
    


    }



// og1
       if ((OG1[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d  1  2.26760E-04  1.28550E-07  ; OG1-Hi rep\n", OG1[i]+1,  H[i]+1);

    }
    if ((OG1[i]>=0) && (O[i]>=0))
    {
      printf("%d    %d  1  7.74841E-04   1.50095E-06 ; OG1-Oi rep\n", OG1[i]+1,  O[i]+1);

    }
    if ((OG1[i]>=0) && (NH[i+1]>=0))
    {
      printf("%d    %d  1  3.83970E-03  4.12285E-06 ;  OG1-Ni+1 \n", OG1[i]+1,  NH[i+1]+1);

    }


    if ((OG1[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d  1  1.13380E-03  1.28550E-07  ;  Og1-Hi+1\n", OG1[i]+1,  H[i+1]+1);
      
    }
    


    if (i>1)
    {
    
     if ((OG1[i]>=0) && (O[i-1]>=0))
     {  
       printf("%d    %d  1  3.87420E-03   1.50095E-06 ;  Og1-Oi-1 \n", OG1[i]+1,  O[i-1]+1);
       
     }
     


    }



 if ((strncmp(resnm[i],"ASP",3)==0)||(strcmp(resnm[i],"ASN")==0))
 {

    if ((CG[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d  1  0.0     8.5E-07  ;  CG-Hi asx\n", CG[i]+1,  H[i]+1);
    
    }
    if ((CG[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d  1  0.0     8.5E-07  ;  CG-Hi+1 asx\n", CG[i]+1,  H[i+1]+1);

    }

    if ((CG[i]>=0) && (NH[i+1]>=0))
    {
      printf("%d    %d  1  1.91985E-03   2.06142E-06  ;  CG-Ni+1 asx \n", CG[i]+1,  NH[i+1]+1);
       
    }

    if ((ND2[i]>=0) && (NH[i]>=0))
    {
      printf("%d    %d  1  1.30345E-03   9.50217E-07  ;  ND2-Ni asx \n", ND2[i]+1,  NH[i]+1);

    }

    if ((ND2[i]>=0) && (NH[i+1]>=0))
    {
      printf("%d    %d  1  2.76211E-03   4.26692E-06  ;  ND2-Ni+1 asx \n", ND2[i]+1,  NH[i+1]+1);

    }



  if (strcmp(resnm[i],"ASN")==0)
{
    if ((HD21[i]>=0) && (O[i]>=0))
    {
      printf("%d    %d  1  1.275068E-03  2.139210E-08  ;  HD21-Oi asx\n", HD21[i]+1,  O[i]+1);

    }
    if ((OD1[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d   1  1.275068E-03  2.139210E-08  ;  OD1-Hi asx\n", OD1[i]+1,  H[i]+1);

    }


    if ((OD1[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d   1  1.275068E-03  2.139210E-08 ;  OD1-Hi+1 asx\n", OD1[i]+1,  H[i+1]+1);

    }



}



  if (strncmp(resnm[i],"ASP",3)==0)
   {


    if ((OD1[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d   1  1.275068E-03  2.139210E-08 ;  OD1-Hi asx\n", OD1[i]+1,  H[i]+1);

    }

    if ((OD1[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d   1  1.275068E-03  2.139210E-08 ;;  OD1-Hi+1 asx\n", OD1[i]+1,  H[i+1]+1);

    }

    if ((OD2[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d   1  1.275068E-03  2.139210E-08 ;;  OD2-Hi asx\n", OD2[i]+1,  H[i]+1);

    }


    if ((OD2[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d   1  1.275068E-03  2.139210E-08 ;;  OD2-Hi+1 asx\n", OD2[i]+1,  H[i+1]+1);

    }


    if ((CG[i]>=0) && (O[i]>=0))
    {
      printf("%d    %d   1  7.41398E-03   2.29029E-06 ;;  CG-Oi asx\n", CG[i]+1,  O[i]+1);

    }




   }





 }

 else

if ((strcmp(resnm[i],"PHE")==0)||(strcmp(resnm[i],"TYR")==0)||(strcmp(resnm[i],"TRP")==0)||(strncmp(resnm[i],"HIS",3)==0))
 {

    if ((CG[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d  1  6.53184E-04  4.76171E-07  ;  CG-Hi aro\n", CG[i]+1,  H[i]+1);

    }
    if ((CG[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d  1  6.53184E-04  4.76171E-07  ;  CG-Hi+1 aro\n", CG[i]+1,  H[i+1]+1);

    }


    if ((CG[i]>=0) && (NH[i+1]>=0))
    {
      printf("%d    %d  1  1.91985E-03   2.06142E-06  ;  CG-Ni+1 aro \n", CG[i]+1,  NH[i+1]+1);

    }
    if ((CD1[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d 1  6.53184E-04  4.76171E-07  ;  Hi-CD1 aro\n", CD1[i]+1,  H[i]+1);
    
    }
      
    if ((ND1[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d 1  6.53184E-04  4.76171E-07  ;  Hi-ND1 aro\n", ND1[i]+1,  H[i]+1);
    
    }
      

    if ((CD2[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d 1  6.53184E-04  4.76171E-07  ;  Hi-CD2 aro\n", CD2[i]+1,  H[i]+1);
    
    }


    if ((CD1[i]>=0) && (NH[i]>=0))
    {
      printf("%d    %d   1  1.30345E-03   9.50217E-07  ;  CD1-Ni  aro\n", CD1[i]+1,  NH[i]+1);

    }

    if ((ND1[i]>=0) && (NH[i]>=0))
    {
      printf("%d    %d   1  1.30345E-03   9.50217E-07  ;  ND1-Ni  aro\n", ND1[i]+1,  NH[i]+1);

    }



    if ((CD2[i]>=0) && (NH[i]>=0))
    {
      printf("%d    %d   1  1.30345E-03   9.50217E-07  ;  CD2-Ni  aro\n", CD2[i]+1,  NH[i]+1);

    }



    if (i>1)
    {

     if ((CD1[i]>=0) && (O[i-1]>=0))
     {
       printf("%d    %d   1     1.06354E-03  6.32621E-07  ;   CD1-Oi-1  aro\n", CD1[i]+1,  O[i-1]+1);

     }



    }

    if (i>1)
    {

     if ((CD2[i]>=0) && (O[i-1]>=0))
     {
       printf("%d    %d    1     1.06354E-03  6.32621E-07  ;   CD2-Oi-1  aro\n", CD2[i]+1,  O[i-1]+1);

     }



    }

    if ((CD1[i]>=0) && (NH[i+1]>=0))
    {
      printf("%d    %d   1  2.76211E-03  4.26692E-06  ;  CD1-Ni+1 aro\n", CD1[i]+1,  NH[i+1]+1);

    }

    if ((ND1[i]>=0) && (NH[i+1]>=0))
    {
      printf("%d    %d   1  2.76211E-03  4.26692E-06  ;  ND1-Ni+1 aro\n", ND1[i]+1,  NH[i+1]+1);

    }



    if ((CD2[i]>=0) && (NH[i+1]>=0))
    {
      printf("%d    %d   1  2.76211E-03  4.26692E-06  ;  CD2-Ni+1 aro\n", CD2[i]+1,  NH[i+1]+1);

    }

 }
else
 {
    if ((CG[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d  1    1.74675E-03   1.70645E-06      ; H-CG H 0.2 0.24 nm rep\n", CG[i]+1,  H[i]+1);

    }
    if ((CG[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d  1    1.74675E-03   1.70645E-06      ; H-CG H 0.2 0.24 nm rep\n", CG[i]+1,  H[i+1]+1);

    }


  
  

  }


    if ((CG1[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d  1  1.06354E-03  6.32621E-07      ; H-CG H 0.2 0.24 nm rep\n", CG1[i]+1,  H[i]+1);

    }
    if ((CG1[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d  1  1.06354E-03  6.32621E-07      ; H-CG H 0.2 0.24 nm rep\n", CG1[i]+1,  H[i+1]+1);

    }




    if ((CG2[i]>=0) && (H[i]>=0))
    {
      printf("%d    %d  1  1.06354E-03  6.32621E-07      ; H-CG H 0.2 0.24 nm rep\n", CG2[i]+1,  H[i]+1);

    }
    if ((CG2[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d  1  1.06354E-03  6.32621E-07      ; H-CG H 0.2 0.24 nm rep\n", CG2[i]+1,  H[i+1]+1);

    }





   if ((NH[i]>=0) && (C[i+1]>=0))
    {
      printf("%d    %d  1   3.37244E-03   5.68668E-06 ; 0.5 0.345 nm\n", NH[i]+1,  C[i+1]+1);

    }


    if ((O[i]>=0) && (O[i+1]>=0))
    {
      printf("%d    %d  1   4.13270E-03   5.33725E-06 ; Oi-Oi+1 \n", O[i]+1,  O[i+1]+1);

    }






    if ((NH[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d  1   4.42E-04   1.64E-08  ;   Ni---Hi+1 \n", NH[i]+1,  H[i+1]+1);

    }





   if ((NH[i]>=0) && (O[i+1]>=0))
    {
      printf("%d    %d  1   2.87870E-03  2.31737E-06 ; 0.894 0.305 nm\n", O[i+1]+1, NH[i]+1);

    }

   if ((CA[i]>=0) && (O[i+1]>=0))
    {
      printf("%d    %d  1   2.26146E-03   3.19637E-06  ; 0.6  0.335 nm\n", CA[i]+1,  O[i+1]+1);

    }

   if ((C[i]>=0) && (O[i+1]>=0))
    {
      printf("%d    %d  1    1.91985E-03   2.06142E-06  ; 0.6  0.32 nm\n", C[i]+1,  O[i+1]+1);

    }

   if ((O[i]>=0) && (C[i+1]>=0))
    {
      printf("%d    %d  1    9.29809E-04   3.60227E-07  ; Oi-Ci+1\n", O[i]+1,  C[i+1]+1);

    }


   if ((CA[i]>=0) && (C[i+1]>=0))
    {
      printf("%d    %d  1    2.49170E-03   6.92920E-06  ;  0.45  0.375 nm\n", CA[i]+1,  C[i+1]+1);

    }

   if ((CB[i]>=0) && (H[i+1]>=0))
    {
      printf("%d    %d  1    1.74675E-03   1.70645E-06      ; H-Cb H 0.2 0.24 nm rep\n", CB[i]+1,  H[i+1]+1);

    }



   if (i<nres-1)
   {

   if ((NH[i]>=0) && (C[i+2]>=0))
    {
      printf("%d    %d  1    3.37244E-03   5.68668E-06 ; 0.671 0.345 nm\n", NH[i]+1,  C[i+2]+1);

    }




   if ((NH[i]>=0) && (O[i+2]>=0))
    {
      printf("%d    %d  1    2.87870E-03  2.31737E-06 ; 0.894 0.305 nm\n", O[i+2]+1, NH[i]+1);

    }

   if ((CA[i]>=0) && (O[i+2]>=0))
    {
      printf("%d    %d  1    2.26146E-03   3.19637E-06   ; 0.6  0.335 nm\n", CA[i]+1,  O[i+2]+1);

    }

   if ((C[i]>=0) && (O[i+2]>=0))
    {
      printf("%d    %d  1    1.91985E-03   2.06142E-06   ; 0.6  0.32 nm\n", C[i]+1,  O[i+2]+1);

    }

   if ((C[i]>=0) && (NH[i+2]>=0))
    {
      printf("%d    %d  1     3.37244E-03   5.68668E-06 ; 0.671 0.345 nm\n", C[i]+1,  NH[i+2]+1);

    }



   if ((O[i]>=0) && (CA[i+2]>=0))
    {
      printf("%d    %d  1    2.26146E-03   3.19637E-06  ; 0.6  0.335 nm\n", CA[i+2]+1, O[i]+1);

    }

   if ((O[i]>=0) && (NH[i+2]>=0))
    {
      printf("%d    %d  1    1.10468E-03   3.41254E-07 ; Oi-Ni+2\n", O[i]+1,  NH[i+2]+1);

    }


   if ((CA[i]>=0) && (CA[i+2]>=0))
    {
      printf("%d    %d  1    2.81500E-03   9.90525E-06 ; 0.45  0.39 nm\n", CA[i]+1,  CA[i+2]+1);

    }


   if ((CA[i]>=0) && (C[i+2]>=0))
    {
      printf("%d    %d  1    2.49170E-03   6.92920E-06  ;  0.45  0.375 nm\n", CA[i]+1,  C[i+2]+1);

    }


   if ((C[i]>=0) && (CA[i+2]>=0))
    {
      printf("%d    %d  1    2.49170E-03   6.92920E-06  ;  0.45  0.375 nm\n", C[i]+1,  CA[i+2]+1);

    }


   if ((C[i]>=0) && (C[i+2]>=0))
    {
      printf("%d    %d  1    2.17678E-03   4.73838E-06  ;  0.45  0.36 nm\n", C[i]+1,  C[i+2]+1);

    }



    if (i>1)
      {
       if ((C[i-1]>=0) && (O[i+2]>=0))
       {
         printf("%d    %d  1     1.91985E-03   2.06142E-06   ; 0.6  0.32 nm\n", C[i-1]+1,  O[i+2]+1);

       }
      }






   }





  





  }

printf("[ impropers ]\n");
  for (i=1;i<nres;i++)
  {
      if ((C[i]>=0) && (O[i]>=0) && (NH[i+1]>=0) && (CA[i+1]>=0))
        printf("%d  %d  %d  %d  2  0.0 75.0\n",O[i]+1, C[i]+1, NH[i+1]+1, CA[i+1]+1);

  }


free_parse(nres, pres);

free(pres);
free(atm_idx);


for (i=0;i<nres+1;i++) free (resnm[i]);
free(resnm);
free(C);
free(CA);
free(NH);
free(CB);
free(O);
free(H);
free(CG1);
free(CG2);
free(CG);
free (OG);
free(OE);
free(NE);
free(AUX);
free(CD1);
free(CD2);
free (OG1);
free(OD1);
free(OD2);
free(HD1);
free(ND1);
free(ND2);
free(HD21);
free(HD22);

return 0;
}
