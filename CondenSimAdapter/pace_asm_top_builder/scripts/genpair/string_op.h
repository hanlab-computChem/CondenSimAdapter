#ifndef STRING_OP
#define STRING_OP


int sgn(float snum);
int str2int(char *numstr);
void int2str(int valu, char *tempstr);
float str2float(char *numstr);
void float2str(float valu, char *sstr, int ndigit);
void trcspace(char *tstr);
int parcnt(char *pstr);
void parstr(char *pstr, int pnum, char *tempstr);
void cpstrfrag (char *a,  char *b, int be, int en);


#endif
