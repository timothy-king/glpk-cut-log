/* cutlog01.c (api extension routines) */


#include "glpapi.h"

int glp_get_it_cnt(glp_prob *P){
  if(P == NULL){
    return 0;
  }else{
    return P->it_cnt;
  }
}
