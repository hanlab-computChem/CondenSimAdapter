#include <stdlib.h>
#include <stdio.h>
#include "file_op.h"


#define DIM 3



void get_x_data (FILE *fp, float *x,int  natom)
{
int i,j,k;

  for(i=0;i<natom;i++){

    for (j=0;j<DIM;j++){
     READ(fp, k);
     x[i*DIM+j]=(float)k/1000;
     

      }


   }


}

void put_x_data (FILE *fp, float *x,int natom)
{


int i,j,k;


   for (i=0;i<natom;i++){
      for (j=0;j<DIM;j++){
        k=(int)((float)(x[i*DIM+j]*1000.0));
        WRITE(fp, k);

          }



    }




}


void get_head_data (FILE *fp, int *natom, int *nframe )
{
int k;
 READ(fp,k);
 *natom=k;
 READ(fp,k);
 *nframe=k;

}


void put_head_data (FILE *fp, int natom , int nframe)
{

WRITE(fp,natom);
WRITE(fp,nframe);

}

#define MAXCHAR_LINE  10000

void get_line (FILE *fp, char *buffer)
{


int i=0;
char gc;
  



  if (feof(fp)==0) READ(fp,gc);
     while ((gc!=10) && (feof(fp)==0)) 
      {
        if (gc!=13) { buffer[i]=gc;i++;}
        READ(fp,gc);
       if (i>MAXCHAR_LINE) {i=0;break;}       
      }

     buffer[i]=0;



}


void seek_data (FILE *fp, int setFrame, int natm)

{
if (fseek(fp, (2+DIM*natm*(setFrame-1))*sizeof(int), SEEK_SET)
       !=0)
    {printf("File has not enough frames\n");exit(1);}

}






#undef DIM
