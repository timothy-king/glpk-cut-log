/* cutlog01.c (api extension routines) */


#include "glpapi.h"
#include "glpios.h"

int glp_get_it_cnt(glp_prob *P){
  if(P == NULL){
    return 0;
  }else{
    return P->it_cnt;
  }
}


/* int glp_ios_cut_klass(glp_tree *T){ */
/*   if(T == NULL){ */
/*     return 0; */
/*   }else{ */
/*     return T->cut_klass; */
/*   } */
/* } */
/* /\* determine the type of cut routine *\/ */


int glp_ios_get_cut(glp_tree *T, int i, int* ind, double* val, int* klass, int* type, double* rhs){
  xassert(T != NULL);

  IOSCUT* cut;
  int len;
  IOSAIJ* aij;
  glp_prob* prob;

  if (T->reason != GLP_ICUTADDED){
    xerror("glp_ios_get_cut: not called during cut added.\n");
  }
  cut = ios_find_row(T->local, i);
  len = 0;
  for(len = 0, aij = cut->ptr; aij != NULL; aij = aij->next)
  {
    len++;
    if(ind != NULL){ ind[len] = aij->j; }
    if(val != NULL){ val[len] = aij->val; }
  }
  if(klass != NULL){ *klass = cut->klass; }
  if(type != NULL){ *type = cut->type; }
  if(rhs != NULL){ *rhs = cut->rhs; }
  return len;
}


IOSAUX *ios_create_aux(){
  IOSAUX *aux;
  aux = xmalloc(sizeof(IOSAUX));
  aux->r = 0;
  aux->r_mult = 0.0;
  aux->c = 0;
  aux->c_mult = 0.0;
  return aux;
}

void ios_delete_aux(IOSAUX *aux){
  xassert(aux != NULL);
  xfree(aux);
  return;
}

static void cut_set_gmi_aux(IOSCUT *cut, int j){
  xassert(cut != NULL);
  if(cut->aux == NULL){
    cut->aux = ios_create_aux();
  }
  cut->aux->r = j;
  cut->aux->r_mult = +1.0;
  cut->aux->c = 0;
  cut->aux->c_mult = 0.0;
}

void ios_cut_set_gmi_aux(glp_tree *T, int ord, int j){
  IOSCUT *cut;
  cut = ios_find_row(T->local, ord);
  xassert(cut != NULL);
  cut_set_gmi_aux(cut, j);
}

void glp_ios_cut_get_aux(glp_tree *tree, int ord, int *r, double *rm, int *c, double *cm){
  IOSCUT* cut;
  if (tree->reason != GLP_ICUTADDED){
    xerror("glp_ios_cut_get_gmi_aux: not called during cut added.\n");
  }
  cut = ios_find_row(tree->local, ord);
  xassert(cut != NULL);
  if(cut->aux != NULL){
    if(r  != NULL){ *r  = cut->aux->r; }
    if(rm != NULL){ *rm = cut->aux->r_mult; }
    if(c  != NULL){ *c  = cut->aux->c; }
    if(cm != NULL){ *cm = cut->aux->c_mult; }
  }
  return;
}

