#ifndef FILE_OP
#define FILE_OP

#define READ(a,b)   fread(&(b),sizeof(b),1,(a))
#define WRITE(a,b)  fwrite(&(b),sizeof(b),1,(a))


void get_x_data (FILE *fp, float *x,int natom);
void put_x_data (FILE *fp, float *x,int natom);
void get_head_data (FILE *fp, int *natom, int *nframe );
void put_head_data (FILE *fp, int natom , int nframe);
void get_line (FILE *fp, char *buffer);
void seek_data (FILE *fp, int setFrame, int natm);
#endif
