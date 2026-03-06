#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "string_op.h"


#define  MAXBUFFER 50000


int sgn(float snum)
{
  if (snum<0)  {return -1;}
  else if (snum>0) { return 1;}
    else {return 0;}

}




int str2int(char *numstr)

{

  int tempval;
  int vsign;
  int vi;

   if (strlen(numstr)==0) { tempval=0;vsign=0;}
    else
    {
      vsign=1;
      tempval=0;
       for (vi=0;vi<(int)strlen(numstr);vi++){
         if (numstr[vi]=='-') { vsign=-1;}
          else
           {
             if ((numstr[vi]>='0') && (numstr[vi]<='9'))
                {
                  tempval*=10;
                  tempval+=numstr[vi]-'0';
                }

           }
       }

     }
   tempval*=vsign;
  return tempval;

}


void int2str(int valu, char *tempstr)
{

  int tempval;
  int i=0;

  
   tempval=(int)fabs(valu);
   do
  {
   

   tempstr[i]=tempval-(int)(tempval/10)*10+'0';i++;

   tempval/=10;

   } while (tempval!=0);

   if (valu<0) {tempstr[i]='-';i++;}
   tempstr[i]=0;
   i--;
   int j;
   char p;
   for (j=0;j<=(int)(i/2);j++) 
   {
    p=tempstr[j];
    tempstr[j]=tempstr[i-j];
    tempstr[i-j]=p;

    }
 
  
}

#define FALSE 0
#define TRUE  1

float str2float(char * numstr)
{

  int fi,fexp;
  float fv,fv1;
  int fsign;
  int fnumbgn,fscfcnt,fint,fsiexp;

  fv=0;
  fnumbgn=FALSE;
  fscfcnt=FALSE;
  fint=TRUE;
  fv1=1;
  fexp=0;
  for(fi=0;fi<(int)strlen(numstr);fi++) 
   {
//     printf("%c %10.3f%10.3f %d %d %d\n",numstr[fi],fv ,fv1, fnumbgn, fint, fscfcnt);
     
     if (numstr[fi]!=' ') 
      {
      if (fnumbgn==FALSE) 
       {
          if (numstr[fi]=='-') 
          {
          fsign=-1;
          }
          else
          {
          fsign=1;
           if (numstr[fi]=='.') 
           {
           fint=FALSE;
           }
            else
               {
                 fv=(float)(numstr[fi]-'0');
               }

          }
         fnumbgn=TRUE;
       }
      else

      {
         if (fscfcnt==FALSE) 
         {
           if ((numstr[fi]=='E') || (numstr[fi]=='e')) 
           {
           fscfcnt=TRUE;
           fv1=1;
           fsiexp=TRUE;
           }
           else
           {
            if (fint) 
              {
                 if (numstr[fi]=='.') 
                 {
                 fint=FALSE;
                 }
                 else
                 {
                   fv=fv*10;
                   fv=fv+(float)(numstr[fi]-'0');
                 }
              }
              else
              {
                 fv1/=10;
                 fv+=fv1*(float)(numstr[fi]-'0');
              }
           }
         }
         else
         {
            if (numstr[fi]=='-') 
            {
            fsiexp=FALSE;
            }
            else
            {
              if (numstr[fi]!='+')  {
              fexp*=10;
              fexp+=numstr[fi]-'0';
                                             }
            }



         }



      }
     }

   }

   if (fexp>0) 
   {
     for( fi=1;fi<= fexp;fi++) 
      if (fsiexp) { fv*=10;} else {fv/=10;}

   }
  fv*=(float)fsign;
  return fv;

}


#undef FALSE
#undef TRUE


void float2str(float valu, char *sstr, int ndigit)
{
    int i;
    int sint,sdigit,snum;
    float sfrac;

  i=0;
 sint=(int)(fabs(valu));
 sfrac=fabs(valu)-(float)sint;
 do
 {
  sdigit=sint-((int)(sint / 10))*10;
  sstr[i]='0'+sdigit;
  i++;
  sint/=10;
 } while (sint!=0);


 if (valu<0) { sstr[i]='-';i++;}
i--;

   int j;
   char p;
   for (j=0;j<=(int)(i/2);j++)
   {
    p=sstr[j];
    sstr[j]=sstr[i-j];
    sstr[i-j]=p;

    }


i++;



 if (sfrac>1e-8) 
 {
  sstr[i]='.';
  i++;
  snum=1;
  do
  {
  sdigit=(int)((float)(sfrac*10));
  sstr[i]='0'+sdigit;
  i++;
  sfrac=sfrac*10-(float)sdigit;
  snum++;
  } while ((sfrac>1e-8) && (snum<=ndigit));


 }


}





void trcspace(char *tstr)
{
   char   tempstr[MAXBUFFER],temp[MAXBUFFER];
   int ti,tend,tspcnt;
   int i;


  if (strlen(tstr)==0) {
    tempstr[0]='|';
    tempstr[1]='n';
    tempstr[2]='o';
    tempstr[3]='n';
    tempstr[4]='e';
     i=5;
  
    }
  else
  {
    i=0;
    tend=(int)strlen(tstr);
    tspcnt=0;
    for( ti=0;ti< tend;ti++)  if (tstr[ti]==' ')  tspcnt++;
    if (tspcnt==tend) 
     {
    tempstr[0]='|';
    tempstr[1]='n';
    tempstr[2]='o';
    tempstr[3]='n';
    tempstr[4]='e';
     i=5;
      
     }
     else
     {
      int j;
      for (j=1;j<=tend;j++) 

       {temp[j]=tstr[j-1];

         if (temp[j]==9) temp[j]=' ';
        }
      temp[0]=' ';
      
      temp[tend+1]=0;     
     for( ti=0;ti<tend+1;ti++) 
          if (temp[ti]!=' ') 
           {
            if (temp[ti-1]==' ') { tempstr[i]='|';i++;}
            tempstr[i]=temp[ti];
            i++;

           }




     }


  }

  tempstr[i]='|';
  i++;
  tempstr[i]=0;

  for (ti=0;ti<=i;ti++) tstr[ti]=tempstr[ti]; 


}

#undef MAXBUFFER




int parcnt(char *pstr)
{
  int  ti,pj;

   pj=0;
   for(ti=0;ti<(int)strlen(pstr);ti++)  if (pstr[ti]=='|')  pj++;
   pj--;
   return pj;




}






void parstr(char *pstr, int pnum, char *tempstr)
{
   int pai,pj;
   


   if (pnum>parcnt(pstr)) 
   {tempstr[0]=0;exit;}
   else
   {
      int i=0;

       
      pj=0;
      for( pai=0;pai<(int)strlen(pstr);pai++) 
       if (pstr[pai]=='|') { pj++;}
         else
          {
           if (pj==pnum) { tempstr[i]=pstr[pai];i++;}

          }

      tempstr[i]=0;
   }


}


void cpstrfrag (char *a,  char *b, int be, int en)
{
 int slen=strlen(a);

 if (be<1)  be=1;
 if (en>slen) en=slen;
 if (en<be)  en=be;
 if (be>en)  be=en;

  be--;
  en--;

  int i=0;
  int j;


   for (j=be;j<=en;j++)
   {

     b[i]=a[j];
     i++;

    }


   b[i]=0;

   



}
