#ifndef PARSE_H_CG
#define PARSE_H_CG

#define MAXNOE 5000
#define MAX_LIST 50


 struct   atm_P {

   int idx,isH,res;
  //isH 1:is H 0:heavy atm without missing H
  //isH :-1 lack of 1H :-2 lack of 2H -3: lack of 3H
  //isH:-4  lack of aromatic H
   char nm[10];
   struct atm_P *left, * down, *up;


};







 



void init_hdb_gmx (char *init);
void init_parse_gmx (char *init);
void construct_parse (int natm, int nres,  struct atm_P **pres, struct atm_P **atm_idx, FILE *fp,int parseMod );
void free_parse (int nres, struct atm_P **pres);
void bIOSYM_format (char *a, int *res_id, char *atm_str);
void init_NOE_list (int *list_A, int *list_B, struct atm_P ** atm_idx,int natm, int *nNOE,float * up_bond, FILE *fp);
void bbdih_list (int *list_A, int *list_B, int resN,int nres, struct atm_P **pres);
int find_atm_name (struct atm_P *en, char *cnm );
int find_atm_noerr (struct atm_P *en, char *cnm );

#endif
