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
  if ( cut == NULL ) {
    xerror("glp_ios_get_cut: called with an invalid index.");
  }
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


IOSAUX *ios_create_aux(int n){
  IOSAUX *aux;
  aux = xmalloc(sizeof(IOSAUX));
  aux->nrows = n;
  aux->rows = xcalloc(1+n, sizeof(int));
  aux->mult = xcalloc(1+n, sizeof(double));
  return aux;
}

void ios_delete_aux(IOSAUX *aux){
  xassert(aux != NULL);
  xfree(aux->rows);
  xfree(aux->mult);
  xfree(aux);
  return;
}

static void cut_set_aux(IOSCUT *cut, int n,
                        const int rows[], const double coeffs[]){
  int i;
  xassert( cut != NULL );
  if( cut->aux != NULL ) {
    ios_delete_aux(cut-> aux);
  }

  cut->aux = ios_create_aux(n);
  xassert( cut->aux->nrows == n );
  for ( i = 1; i <= n; i++)
  {  cut->aux->rows[i] = rows[i];
     cut->aux->mult[i] = coeffs[i];
  }
}

void ios_cut_set_single_aux(glp_tree *T, int ord, int j){
  IOSCUT *cut;
  cut = ios_find_row(T->local, ord);
  xassert(cut != NULL);

  /* set up arrays */
  int ind[1+1];
  double coeffs[1+1];
  ind[1] = j;
  coeffs[1] = +1.0;

  /* call general procedure */
  cut_set_aux(cut, 1, ind, coeffs);
}

void ios_cut_set_aux(glp_tree *T, int ord, int n,
                     const int rows[], const double coeffs[]){
  IOSCUT *cut;
  cut = ios_find_row(T->local, ord);
  xassert(cut != NULL);
  cut_set_aux(cut, n, rows, coeffs);
}

int glp_ios_cut_get_aux_nrows(glp_tree *tree, int ord){
  IOSCUT *cut;
  IOSAUX *aux;
  if (tree->reason != GLP_ICUTADDED){
    xerror("glp_ios_cut_get_aux_nrows: not called during cut added.\n");
  }
  cut = ios_find_row(tree->local, ord);
  if ( cut == NULL ){
    xerror("glp_ios_cut_get_aux_nrows: not called on a valid cut.\n");
  }
  aux = cut->aux;
  return (aux == NULL) ? 0 : aux->nrows;
}

void glp_ios_cut_get_aux_rows(glp_tree *tree, int ord,
                              int rows[], double coeffs[]){
  IOSCUT *cut;
  IOSAUX *aux;
  int j, nrows;
  if (tree->reason != GLP_ICUTADDED){
    xerror("glp_ios_cut_get_aux_rows: not called during cut added.\n");
  }
  cut = ios_find_row(tree->local, ord);
  if ( cut == NULL ){
    xerror("glp_ios_cut_get_aux_rows: not called on a valid cut.\n");
  }
  aux = cut->aux;
  if( aux != NULL ){
    nrows = aux->nrows;
    for ( j = 1; j <= nrows; j++ )
    {  if ( rows != NULL ) { rows[j] = aux->rows[j]; }
       if ( coeffs != NULL ) { coeffs[j] = aux->mult[j]; }
    }
  }
  return;
}

