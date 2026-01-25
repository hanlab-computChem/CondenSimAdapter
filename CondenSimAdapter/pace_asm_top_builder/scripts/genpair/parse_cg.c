#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "file_op.h"
#include "string_op.h"
#include "parse_cg.h"










#define MAXBUF 5000
#define DIM 3
#define MAXRES      120
#define MAXCHAR_RES 500






void init_parse_gmx (char *init)
{
int i=0;
//for cg
//strcpy(&init[i*MAXCHAR_RES],
//       "|ALA|N|1|CA|2|C|0|CB|0|");
//i++;


//for modified cg
strcpy(&init[i*MAXCHAR_RES],
       "|ALA|N|2|H|0|CA|2|C|1|CB|0|O|0|");
i++;



strcpy(&init[i*MAXCHAR_RES],
       "|TYR|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|2|CD1|1|CD2|1|CE1|1|CE2|0|CZ|1|OH|0|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|GLU|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|1|CD|2|OE1|0|OE2|0|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|GLUH|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|1|CD|2|OE1|0|OE2|0|");
i++;


strcpy(&init[i*MAXCHAR_RES],
       "|VAL|N|2|H|0|CA|2|C|1|CB|2|O|0|CG1|0|CG2|0|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|HISA|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|2|ND1|2|CD2|1|HD1|0|CE1|0|NE2|0|");
i++;


strcpy(&init[i*MAXCHAR_RES],
       "|HISD|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|2|ND1|2|CD2|1|HD1|0|CE1|0|NE2|0|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|HIS|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|2|ND1|2|CD2|1|HD1|0|CE1|0|NE2|0|");
i++;


strcpy(&init[i*MAXCHAR_RES],
       "|HISB|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|2|ND1|1|CD2|1|CE1|0|NE2|1|HE2|0|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|HISE|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|2|ND1|1|CD2|1|CE1|0|NE2|1|HE2|0|");
i++;


strcpy(&init[i*MAXCHAR_RES],
       "|HISH|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|2|ND1|2|CD2|1|HD1|0|CE1|0|NE2|1|HE2|0|");
i++;



strcpy(&init[i*MAXCHAR_RES],
       "|GLN|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|1|CD|2|OE1|0|NE2|2|HE21|0|HE22|0|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|LYS|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|1|CD|1|CE|1|NZ|0|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|LYSH|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|1|CD|1|CE|1|NZ|0|");
i++;



strcpy(&init[i*MAXCHAR_RES],
       "|ORN|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|1|CD|1|NE|0|");
i++;   


strcpy(&init[i*MAXCHAR_RES],
       "|DAB|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|1|ND|0|");
i++;   



strcpy(&init[i*MAXCHAR_RES],
       "|NLE|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|1|CD|1|CE|0|");
i++;


//strcpy(&init[i*MAXCHAR_RES],
//       "|LYS|N|1|CA|2|C|1|CB|1|O|0|CD|1|NZ|0|");
//i++;




strcpy(&init[i*MAXCHAR_RES],
       "|LEU|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|2|CD1|0|CD2|0|");
i++;

//strcpy(&init[i*MAXCHAR_RES],
//       "|LEU|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|0|");
//i++;




strcpy(&init[i*MAXCHAR_RES],
       "|PHE|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|2|CD1|1|CD2|1|CE1|1|CE2|0|CZ|0|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|ASP|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|2|OD1|0|OD2|0|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|ASPH|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|2|OD1|0|OD2|0|");
i++;



strcpy(&init[i*MAXCHAR_RES],
       "|GLY|N|2|H|0|CA|1|C|1|O|0|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|SER|N|2|H|0|CA|2|C|1|CB|1|O|0|OG|0|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|ASN|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|2|OD1|0|ND2|2|HD21|0|HD22|0|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|ILE|N|2|H|0|CA|2|C|1|CB|2|O|0|CG1|1|CG2|0|CD|0|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|MET|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|1|SD|1|CE|0|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|NH2|N|2|H1|0|H2|0|");
i++;

//for cg
//strcpy(&init[i*MAXCHAR_RES],
//       "|ACE|CA|1|C|0|");
//i++;

//for modified cg
strcpy(&init[i*MAXCHAR_RES],
       "|ACE|CA|1|C|1|O|0|");
i++;



strcpy(&init[i*MAXCHAR_RES],
       "|CU2|CU|0|");
i++;


strcpy(&init[i*MAXCHAR_RES],
       "|ZN2|ZN|0|");
i++;


strcpy(&init[i*MAXCHAR_RES],
       "|TRP|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|2|CD1|1|CD2|2|NE1|1|CE2|1|CE3|1|HE1|0|CZ2|0|CZ3|1|CH2|0|");
i++;


strcpy(&init[i*MAXCHAR_RES],
       "|PRO|N|1|CA|2|C|1|CB|1|O|0|CG|1|CD|0|");
// modified for ffNEW
//strcpy(&init[i*MAXCHAR_RES],
//       "|PRO|N|1|CA|2|C|1|CB|0|O|0|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|DPR|N|1|CA|2|C|1|CB|1|O|0|CG|1|CD|0|");
// modified for ffNEW
//strcpy(&init[i*MAXCHAR_RES],
//       "|PRO|N|1|CA|2|C|1|CB|0|O|0|");
i++;


strcpy(&init[i*MAXCHAR_RES],
       "|ARG|N|2|H|0|CA|2|C|1|CB|1|O|0|CG|1|CD|1|NE|1|CZ|2|NH1|0|NH2|0|");
i++;



strcpy(&init[i*MAXCHAR_RES],
       "|CYS|N|2|H|0|CA|2|C|1|CB|1|O|0|SG|1|HG|0|");
i++;




strcpy(&init[i*MAXCHAR_RES],
       "|CYSH|N|2|H|0|CA|2|C|1|CB|1|O|0|SG|1|HG|0|");
i++;


strcpy(&init[i*MAXCHAR_RES],
       "|CYS2|N|2|H|0|CA|2|C|1|CB|1|O|0|SG|0|");
i++;


strcpy(&init[i*MAXCHAR_RES],
       "|THR|N|2|H|0|CA|2|C|1|CB|2|O|0|OG1|1|CG2|0|HG2|0|");
i++;


//for cg and its modification
strcpy(&init[i*MAXCHAR_RES],
       "|NMA|N|2|H|0|CA|0|");
i++;






strcpy(&init[i*MAXCHAR_RES],
       "|none|");


}



void init_hdb_gmx (char *init)
{

int i=0;
strcpy(&init[i*MAXCHAR_RES],
       "|ALA|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|TYR|CA|-1|CB|-2|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|GLU|CA|-1|CB|-2|CG|-2|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|VAL|CA|-1|CB|-1|CG1|-3|CG2|-3|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|HIS|CA|-1|CB|-2|CD2|-4|CE1|-4|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|GLN|CA|-1|CB|-2|CG|-2|");
i++;

//strcpy(&init[i*MAXCHAR_RES],
//       "|LYS|CA|-1|CB|-2|CG|-2|CD|-2|CE|-2|NZ|-3|");
//i++;

strcpy(&init[i*MAXCHAR_RES],
       "|LYS|");
i++;



strcpy(&init[i*MAXCHAR_RES],
       "|LEU|CA|-1|CB|-2|CG|-1|CD1|-3|CD2|-3|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|PHE|CA|-1|CB|-2|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|ASP|CA|-1|CB|-2|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|GLY|CA|-2|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|SER|CA|-1|CB|-2|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|ASN|CA|-1|CB|-2|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|ILE|CA|-1|CB|-1|CG1|-2|CG2|-3|CD|-3|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|MET|CA|-1|CB|-2|CG|-2|CE|-3|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|NH2|");
i++;


strcpy(&init[i*MAXCHAR_RES],
       "|ACE|");
i++;



strcpy(&init[i*MAXCHAR_RES],
       "|CU2|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|ZN2|");
i++;

strcpy(&init[i*MAXCHAR_RES],
       "|TRP|CA|-1|CB|-2|");
i++;


strcpy(&init[i*MAXCHAR_RES],
       "|PRO|CA|-1|CB|-2|CG|-2|CD|-2|");
i++;




strcpy(&init[i*MAXCHAR_RES],
       "|ARG|CA|-1|CB|-2|CG|-2|CD|-2|");
i++;



strcpy(&init[i*MAXCHAR_RES],
       "|CYS|CA|-1|CB|-2|");
i++;


strcpy(&init[i*MAXCHAR_RES],
       "|THR|CA|-1|CB|-1|CG2|-3|");
i++;


strcpy(&init[i*MAXCHAR_RES],
       "|NMA|");
i++;


strcpy(&init[i*MAXCHAR_RES],
       "|none|");

}


void construct_parse (int natm, int nres,  struct atm_P **pres,struct atm_P **atm_idx, FILE *fp, int parseMod )
//pres, pres[i] --> atm_P struct of atm N of res i
//atm_idx[i] --> atm_P struct of atm i

{

//should allocate the memory for  pres before enter this function

char buffer[MAXBUF],bufstr[MAXBUF],bufchar;

int *atm_res;
 atm_res=malloc(natm*sizeof(int));

char *atm_str,*res_str;
 atm_str=malloc(natm*10*sizeof(char));
 res_str=malloc((nres+1)*10*sizeof(char));

int i=0,j;
while (i<natm)
{


do{

get_line(fp, &buffer[0]);
trcspace (&buffer[0]);
parstr(&buffer[0],1,&bufstr[0]);
} while (strcmp(bufstr, "ATOM")!=0);


parstr (buffer,3,bufstr);
strcpy (&atm_str[i*10], bufstr);


parstr (buffer,5,bufstr);

if ((bufstr[0]<'0') || (bufstr[0]>'9')) parstr(buffer, 6, bufstr);

atm_res[i]=str2int(bufstr);

parstr (buffer,4,bufstr);
strcpy (&res_str[atm_res[i]*10], bufstr);


/*for (j=0;j<DIM;j++){
parstr(buffer,j+s_phase,bufstr);
cord[i*DIM+j]=str2float(bufstr);

}*/






i++;
}


int head,reer,k,con,l;
struct atm_P *seq[MAXBUF],*ip,*ipx;
char parse_gmx[MAXRES][MAXCHAR_RES], hdb_gmx[MAXRES][MAXCHAR_RES];

init_parse_gmx (&parse_gmx[0][0]);
init_hdb_gmx(&hdb_gmx[0][0]);

for (i=1;i<=nres;i++)
{

j=0;
parstr ( &parse_gmx[j][0],1,bufstr);

 while (strcmp(&res_str[i*10], bufstr)!=0)
 {
 j++;
 parstr ( &parse_gmx[j][0],1,bufstr);
 if (strcmp(bufstr,"none")==0) { printf("%S%d is not defined in parse_gmx.\n",&res_str[i*10],i);exit(1);}

 }



 head=reer=2;
  k=0;
    parstr( &parse_gmx[j][0], head, bufstr);


     while ((strcmp(&atm_str[k*10],bufstr)!=0) || (atm_res[k]!=i))
     {
      k++;
      if (k>=natm) {printf("k overflow!\n");exit(1);}
      }
     
 
      ip=malloc(sizeof(struct atm_P));
  


      strcpy(ip->nm,bufstr);
      ip->idx=k;
      ip->isH=0;
     pres[i]=ip;
      ip->res=i;

      ip->up=NULL;
      ip->down=NULL;
      ip->left=NULL;
      atm_idx[k]=ip;
      seq[head]=ip;



 
  do
  {
     parstr(&parse_gmx[j][0],head+1,buffer);
     con=str2int(buffer);



      
     ipx=NULL;     
    while (con>0)
    {

     reer+=2;
     parstr(&parse_gmx[j][0],reer,buffer);
 
     k=0;


     while ((strcmp(&atm_str[k*10],buffer)!=0) || (atm_res[k]!=i) )
     {
       k++;
      if (parseMod==0) 
       if ((k>=natm) && (buffer[0]!='H')) {

          
          printf("In reisude %s%d, atom %s is missing. \n",&res_str[i*10],i,buffer);exit(1);



          }

        if (k>=natm) break;


      }

      if (k<natm)
      {
      ip=malloc(sizeof(struct atm_P));
      strcpy(ip->nm,buffer);
      
      ip->idx=k;
      ip->isH=0;
      ip->res=i;
      if (buffer[0]=='H') ip->isH=1;
      l=2;
      while (l<parcnt(&hdb_gmx[j][0]))
      {
       parstr(&hdb_gmx[j][0],l,bufstr);
       if (strcmp(buffer,bufstr)==0)
        {parstr(&hdb_gmx[j][0],l+1,bufstr);ip->isH=str2int(bufstr); }
       l+=2;
      }

      ip->down=NULL;
      ip->up=seq[head];
      seq[reer]=ip;
      atm_idx[k]=ip;
      ip->left=ipx;
      ipx=ip;
      }
      else
      {
      seq[reer]=NULL;

       }
      con--;

     }

    // printf("%d %d %p\n",head,reer,ipx);
   
    
    seq[head]->down=ipx;
    do
    {
    head+=2;
    if (head>reer) continue;

     }  while (seq[head]==NULL);

  }while (head<=reer);

}





free(atm_res);
free(atm_str);
free(res_str);


}




#undef DIM
#undef MAXRES
#undef MAXCHAR_RES


void free_parse (int nres, struct atm_P **pres)
{


struct atm_P *seq[MAXBUF], *ip;
int head,reer;
int i;

for (i=1;i<=nres;i++)
{

head=reer=1;
seq[1]=pres[i];
do
{

ip=seq[head];
 if (ip->down!=NULL)
  {
    ip=ip->down;
    reer++;
    seq[reer]=ip;
    while (ip->left!=NULL)
    {
      ip=ip->left;
      reer++;
      seq[reer]=ip;
    }

   }


head++;
} while (head<=reer);

}


}



void bIOSYM_format (char *a, int *res_id, char *atm_str)
//BIOSTM formatted NOE token
{
  int i,j=0,k=0;
  char buffer[MAXBUF];
  for (i=0;i<(int)strlen(a);i++)
  {
    if (k==1) {buffer[j]=a[i];j++; }
    if (k==2) {atm_str[j]=a[i];j++; }

     if (a[i]=='_') {k=1;}
     if ((k==1) && (a[i]==':')) { buffer[j]=0; *res_id=str2int(buffer); j=0;k=2;}


  }
   
  if (k<2) {printf("%s is not BIOSTM format\n",a);exit(1);}
  atm_str[j]=0;



}




void init_NOE_list (int *list_A, int *list_B, struct atm_P ** atm_idx,int natm, int *nNOE,float *up_bond, FILE *fp)

{
 int *list_AA[MAXNOE], *list_BB[MAXNOE];
 float dis[MAXNOE];
 char buffer[MAXBUF],bufstr[MAXBUF],bufstr1[MAXBUF],bufstr2[MAXBUF],bufstr3[MAXBUF],bufstr4[MAXBUF],bufstr5[MAXBUF];
   
  do{
      get_line(fp, buffer);
      trcspace(buffer);
      parstr(buffer,1,bufstr);
  }while (strcmp(bufstr,"#NOE_distance")!=0);   

int i,j,k,m,n,n1;
char l;
struct atm_P *ip;
 *nNOE=0;


  while (strcmp(bufstr,"!")!=0)
 {
   get_line(fp,buffer);
   trcspace (buffer);

   parstr(buffer,1,bufstr);

  if ((strcmp(bufstr,"!")!=0) && ( strcmp(bufstr,"none")))
    {


//for list_AA
   parstr(buffer,1,bufstr);
      bIOSYM_format ( bufstr, &i, bufstr1);


     if (strcmp(bufstr1,"HN")==0) strcpy(bufstr1,"H");
     j=-1;
      for (k=0;k<natm;k++)
      {
      // printf("%p \n",atm_idx[k]);
     
       if (atm_idx[k]!=NULL)
       {
         
         if ((atm_idx[k]->res==i) && (strcmp(atm_idx[k]->nm,bufstr1)==0)) 
        {
        j=k;
        ip=atm_idx[j]->up;
        if (ip->isH<0) j=-1;
        break;

         }}

      }    
  //   printf("ok\n");
  //  printf ("%d  %s \n ", j, bufstr1);
     if (j==-1)
        
      {   


           
          for(k=(int)strlen(bufstr1)-1;k>=0;k--)
              {
               if ((bufstr1[k]<'0') || (bufstr1[k]>'9')) 
                {
                  m=0;n1=0;
                  bufstr3[m]=bufstr1[k];
                  m++;
                  bufstr5[n1]=bufstr1[k];
                  n1++;
                  if (k<(int)strlen(bufstr1)-1) {bufstr3[m]=bufstr1[k+1];m++;}
                  if (k<(int)strlen(bufstr1)-2) {bufstr5[n1]=bufstr1[k+1];n1++;}

                  bufstr5[n1]=0;
                  bufstr3[m]=0;
                  break;

                } 
              }




             if (k<0) {printf("%s is not valid restaint\n",bufstr);exit(0);}
             l=bufstr1[k];
                 
           for(m=0;m<natm;m++)
            if (atm_idx[m]!=NULL)
              {
                  strcpy(bufstr2, atm_idx[m]->nm);


                for(k=(int)strlen(bufstr2)-1;k>=0;k--)
                    {
                     if ((bufstr2[k]<'0') || (bufstr2[k]>'9'))
                      {
                        n=0;
                        bufstr4[n]=bufstr2[k];
                        
                        n++;
                        if (k<(int)strlen(bufstr2)-1) {bufstr4[n]=bufstr2[k+1];n++;}
                        bufstr4[n]=0;
                       
                        
                        break;
                      } 
                    }




                  if (((strcmp(bufstr3,bufstr4)==0) || (strcmp(bufstr5, bufstr4)==0)) && 
                     ( atm_idx[m]->res==i) && (atm_idx[m]->isH!=1)) {j=m;break;}
                  
              }

          if (j==-1) {printf("%s is not valid restaint\n",bufstr);exit(0);}

          m=-atm_idx[j]->isH;
          if ((m==-1)||(m==0)) {printf("conflict H\n");exit(0);}
          
         list_AA[*nNOE]=malloc(MAX_LIST*sizeof(int));
          list_AA[*nNOE][0]=m;
          if (m==1)
          {
            list_AA[*nNOE][1]=j;
            if (atm_idx[j]->up==NULL) {printf("wrong connection for atom %d\n",j);exit(0);}
            list_AA[*nNOE][2]=atm_idx[j]->up->idx;
            

            if (atm_idx[j]->down==NULL) {printf("wrong connection child  for atom %d\n",j); exit(0);}
            ip=atm_idx[j]->down;
            for (k=3;k<=4;k++)
            {
              while ((ip->isH==1) && (ip->left!=NULL))
              {
                ip=ip->left;

               }

              if (ip->isH!=1) 
                 list_AA[*nNOE][k]=ip->idx;
                 else {printf("not enough hevay atoms for atom %d\n",j);exit(0); }
             if (ip->left!=NULL) ip=ip->left;

            }


          }
          else if ((m==2)||(m==4))
               {
                                   



                list_AA[*nNOE][1]=j;
                if (atm_idx[j]->up==NULL) {printf("wrong connection for atom %d\n",j);exit(0);}
                list_AA[*nNOE][2]=atm_idx[j]->up->idx;
            

                if (atm_idx[j]->down==NULL) {printf("wrong connection child  for atom %d\n",j); exit(0);}
                ip=atm_idx[j]->down;
                for (k=3;k<=3;k++)
                {
                  while ((ip->isH==1) && (ip->left!=NULL))
                  {
                    ip=ip->left;

                   }

                  if (ip->isH!=1) 
                     list_AA[*nNOE][k]=ip->idx;
                     else {printf("not enough hevay atoms for atom %d\n",j);exit(0); }
                  if (ip->left!=NULL) ip=ip->left;

                }


                   if (m==2) 
                       {
                          l=bufstr1[(int)strlen(bufstr1)-1];
                          if (l=='1') list_AA[*nNOE][4]=1; else list_AA[*nNOE][4]=2;

                       }







               }
               else 
                    {



                      list_AA[*nNOE][1]=j;

                      if (atm_idx[j]->up==NULL) {printf("wrong connection for atom %d\n",j);exit(0);}
                      list_AA[*nNOE][2]=atm_idx[j]->up->idx;
                      ip=atm_idx[j]->up;
                      if (ip->up==NULL) {printf("wrong up connection for atom %d\n",j);exit(0);} 
                      list_AA[*nNOE][3]=ip->up->idx;





















                    }
                    
                    
                    

      }
      else
      {
        list_AA[*nNOE]=malloc(MAX_LIST*sizeof(int));
        list_AA[*nNOE][0]=0;
        list_AA[*nNOE][1]=j;


      }
     



 //    for (j=0;j<5;j++) printf("%d  ",list_AA[*nNOE][j]);
  //   printf("        ");




//for list_BB
   parstr(buffer,2,bufstr);
      bIOSYM_format ( bufstr, &i, bufstr1);

     if (strcmp(bufstr1,"HN")==0) strcpy(bufstr1,"H");
     j=-1;
      for (k=0;k<natm;k++)
      {

      if (atm_idx[k]!=NULL)

      { if ((atm_idx[k]->res==i) && (strcmp(atm_idx[k]->nm,bufstr1)==0))
        {
        j=k;
        ip=atm_idx[j]->up;
        if (ip->isH<0) j=-1;
        break;

         }}

      }

   // printf ("%s %d %s ", bufstr, j, bufstr1);
  //  printf("%d\n",j); 
     if (j==-1)

      {



          for(k=(int)strlen(bufstr1)-1;k>=0;k--)
              {
               if ((bufstr1[k]<'0') || (bufstr1[k]>'9'))
                {
                  m=0;n1=0;
                  bufstr3[m]=bufstr1[k];
                  m++;
                  bufstr5[n1]=bufstr1[k];
                  n1++;
                  if (k<(int)strlen(bufstr1)-1) {bufstr3[m]=bufstr1[k+1];m++;}
                  if (k<(int)strlen(bufstr1)-2) {bufstr5[n1]=bufstr1[k+1];n1++;}

                  bufstr5[n1]=0;
                  bufstr3[m]=0;
                  break;

                }
              }




             if (k<0) {printf("%s is not valid restaint\n",bufstr);exit(0);}
             l=bufstr1[k];

           for(m=0;m<natm;m++)
            if (atm_idx[m]!=NULL)
              {
                  strcpy(bufstr2, atm_idx[m]->nm);


                for(k=(int)strlen(bufstr2)-1;k>=0;k--)
                    {
                     if ((bufstr2[k]<'0') || (bufstr2[k]>'9'))
                      {
                        n=0;
                        bufstr4[n]=bufstr2[k];

                        n++;
                        if (k<(int)strlen(bufstr2)-1) {bufstr4[n]=bufstr2[k+1];n++;}
                        bufstr4[n]=0;


                        break;
                      }
                    }




                  if (((strcmp(bufstr3,bufstr4)==0) || (strcmp(bufstr5, bufstr4)==0)) &&
                     ( atm_idx[m]->res==i) && (atm_idx[m]->isH!=1)) {j=m;break;}

              }

          if (j==-1) {printf("%s is not valid restaint\n",bufstr);exit(0);}

          m=-atm_idx[j]->isH;
          if ((m==-1)||(m==0)) {printf("conflict H\n");exit(0);}

         list_BB[*nNOE]=malloc(MAX_LIST*sizeof(int));
          list_BB[*nNOE][0]=m;
          if (m==1)
          {
            list_BB[*nNOE][1]=j;
            if (atm_idx[j]->up==NULL) {printf("wrong connection for atom %d\n",j);exit(0);}
            list_BB[*nNOE][2]=atm_idx[j]->up->idx;


            if (atm_idx[j]->down==NULL) {printf("wrong connection child  for atom %d\n",j); exit(0);}
            ip=atm_idx[j]->down;
            for (k=3;k<=4;k++)
            {
              while ((ip->isH==1) && (ip->left!=NULL))
              {
                ip=ip->left;

               }

              if (ip->isH!=1)
                 list_BB[*nNOE][k]=ip->idx;
                 else {printf("not enough hevay atoms for atom %d\n",j);exit(0); }
             if (ip->left!=NULL) ip=ip->left;

            }


          }
          else if ((m==2)||(m==4))
               {




                list_BB[*nNOE][1]=j;
                if (atm_idx[j]->up==NULL) {printf("wrong connection for atom %d\n",j);exit(0);}
                list_BB[*nNOE][2]=atm_idx[j]->up->idx;


                if (atm_idx[j]->down==NULL) {printf("wrong connection child  for atom %d\n",j); exit(0);}
                ip=atm_idx[j]->down;
                for (k=3;k<=3;k++)
                {
                  while ((ip->isH==1) && (ip->left!=NULL))
                  {
                    ip=ip->left;

                   }

                  if (ip->isH!=1)
                     list_BB[*nNOE][k]=ip->idx;
                     else {printf("not enough hevay atoms for atom %d\n",j);exit(0); }
                  if (ip->left!=NULL) ip=ip->left;

                }


                   if (m==2)
                       {
                          l=bufstr1[(int)strlen(bufstr1)-1];
                          if (l=='1') list_BB[*nNOE][4]=1; else list_BB[*nNOE][4]=2;

                       }







               }
               else
                    {



                      list_BB[*nNOE][1]=j;

                      if (atm_idx[j]->up==NULL) {printf("wrong connection for atom %d\n",j);exit(0);}
                      list_BB[*nNOE][2]=atm_idx[j]->up->idx;
                      ip=atm_idx[j]->up;
                      if (ip->up==NULL) {printf("wrong up connection for atom %d\n",j);exit(0);}
                      list_BB[*nNOE][3]=ip->up->idx;





















                    }




      }
      else
      {
        list_BB[*nNOE]=malloc(MAX_LIST*sizeof(int));
        list_BB[*nNOE][0]=0;
        list_BB[*nNOE][1]=j;


      }





//     for (j=0;j<5;j++) printf("%d  ",list_BB[*nNOE][j]);
  //   printf("\n");

     
     parstr(buffer,4,bufstr);
     dis[*nNOE]=str2float(bufstr);




     (*nNOE)++;



//   printf ("%d \n",*nNOE);



    }


  }


 /*  up-bond=malloc(*nNOE*sizeof(float));
   list_A=malloc(*nNOE*sizeof(int*));
   list_B=malloc(*nNOE*sizeof(int*));  */ 


   for (i=0;i<*nNOE;i++)
    {
    for (j=0;j<MAX_LIST;j++)
     { 
        list_A[i*MAX_LIST+j]=list_AA[i][j];
        list_B[i*MAX_LIST+j]=list_BB[i][j];
       
      }
         free(list_AA[i]);
         free(list_BB[i]);
         up_bond[i]=dis[i];

     }




}



/*void free_NOE_list(int nNOE, int **list_A, int **list_B, float *up-bond)
{
 int i;
  for (i=0;i<nNOE;i++)
   {
     free(list_A[i]);
     free(list_B[i]);
   
    }

   free(list_A);
   free(list_B);
   free(up-bond);


}*/


int find_atm_name (struct atm_P *en, char *cnm )
//input : pres[i]  out : atm idx with name of cnm
{ 
  
  
struct atm_P *seq[MAXBUF], *ip;
int head,reer;


    

head=reer=1;
seq[1]=en;

if (strcmp(en->nm,cnm)==0) return en->idx;

do
{ 
  
ip=seq[head];
 if (ip->down!=NULL)
  {
    ip=ip->down;
    reer++;
    if (strcmp(ip->nm, cnm)==0) return ip->idx;
    seq[reer]=ip;
    while (ip->left!=NULL)
    {
      
      ip=ip->left; 
      reer++;
      if (strcmp(ip->nm,cnm)==0) return ip->idx;
      seq[reer]=ip;
    }
 
   }

  
head++;
} while (head<=reer);
      
 
   printf("no match atom for %s in res %d\n",cnm, en->res);
   exit(0);

}





int find_atm_noerr (struct atm_P *en, char *cnm )
//input : pres[i]  out : atm idx with name of cnm
{


struct atm_P *seq[MAXBUF], *ip;
int head,reer;




head=reer=1;
seq[1]=en;

if (strcmp(en->nm,cnm)==0) return en->idx;

do
{

ip=seq[head];
 if (ip->down!=NULL)
  {
    ip=ip->down;
    reer++;
    if (strcmp(ip->nm, cnm)==0) return ip->idx;
    seq[reer]=ip;
    while (ip->left!=NULL)
    {

      ip=ip->left;
      reer++;
      if (strcmp(ip->nm,cnm)==0) return ip->idx;
      seq[reer]=ip;
    }

   }


head++;
} while (head<=reer);

    return -1;
//   printf("no match atom for %s in res %d\n",cnm, en->res);
//   exit(0);

}






















void bbdih_list (int *list_A, int *list_B, int resN,int nres, struct atm_P **pres)
{

if (resN<=1) {printf("First Res has no bb dih\n");exit(0);}
if (resN>=nres) {printf("Last Res has no bb dih\n");exit(0);}

 
 
   
   list_A[0]=find_atm_noerr ( pres[resN-1], "C");
   list_A[1]=find_atm_noerr ( pres[resN], "N");
   list_A[2]=find_atm_noerr (pres[resN], "CA");
   list_A[3]=find_atm_noerr (pres[resN], "C");

   list_B[0]=find_atm_noerr (pres[resN],"N");
   list_B[1]=find_atm_noerr (pres[resN],"CA");
   list_B[2]=find_atm_noerr (pres[resN],"C");
   list_B[3]=find_atm_noerr (pres[resN+1],"N");






}



#undef MAXBUF
