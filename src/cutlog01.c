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


IOSAUX *ios_create_aux(int n, const int rows[], const double coeffs[]){
  IOSAUX *aux;
  int i;
  aux = xmalloc(sizeof(IOSAUX));
  aux->nrows = n;
  aux->rows = xcalloc(1+n, sizeof(int));
  aux->mult = xcalloc(1+n, sizeof(double));
  aux->selected = -1;
  aux->mir_cset = NULL;

  for ( i = 1; i <= n; i++)
  {  aux->rows[i] = rows[i];
     aux->mult[i] = coeffs[i];
  }

  return aux;
}

void ios_delete_aux(IOSAUX *aux){
  xassert(aux != NULL);
  xfree(aux->rows);
  xfree(aux->mult);
  if( aux->mir_cset != NULL ){
    xfree(aux->mir_cset);
  }
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

  cut->aux = ios_create_aux(n, rows, coeffs);
  xassert( cut->aux->nrows == n );
}

static void cut_set_aux_mir(IOSAUX *aux, double delta,
                            int n, const char cset[]){
  int i;
  xassert( aux != NULL );
  if ( aux->mir_cset != NULL )
  {  xfree(aux->mir_cset);
  }

  aux->mir_cset = xcalloc(1+n, sizeof(char));
  for ( i = 1; i <= n; i++)
  {  aux->mir_cset[i] = cset[i];
  }

  aux->mir_delta = delta;
}

void ios_cut_set_aux_mir(glp_tree *T, int ord, double delta,
                         int n, const char cset[]){
  IOSCUT *cut;
  cut = ios_find_row(T->local, ord);
  xassert(cut != NULL);
  cut_set_aux_mir(cut->aux, delta, n, cset);
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


void glp_ios_cut_get_aux_mir(glp_tree *tree, int ord,
                             char *cset, double *delta){
  IOSCUT *cut;
  IOSAUX *aux;
  int j, n;
  if ( tree->reason != GLP_ICUTADDED ){
    xerror("glp_ios_cut_get_aux_mir: not called during cut added.\n");
  }
  cut = ios_find_row(tree->local, ord);
  if ( cut == NULL ){
    xerror("glp_ios_cut_get_aux_mir: not called on a cut.\n");
  }
  if ( cut->klass != GLP_RF_MIR ){
    xerror("glp_ios_cut_get_aux_mir: not called on a mir cut.\n");
  }
  aux = cut->aux;
  if ( delta != NULL ){
    (*delta) = aux->mir_delta;
  }
  n = tree->n;
  if( cset != NULL ){
    for ( j=1; j <= n; j++ ){
      cset[j] = (aux->mir_cset == NULL) ? 0 : aux->mir_cset[j];
    }
  }
}

void ios_cut_set_selected(IOSCUT *cut, int sel){
#ifdef CUT_DEBUG
  static int i = 0;
  ++i;
  printf("ios_cut_set_selected: %d %d %p\n", i, sel, cut);
#endif

  IOSAUX *aux;
  aux = cut->aux;
  if ( aux != NULL ){
    aux->selected = sel;
  }
}

int glp_ios_selected_cuts(glp_tree *tree, int ords[], int sel[]){
  int len, j, N, s;
  IOSPOOL* pool;
  IOSCUT* cut;
  IOSAUX* aux;
  if ( tree == NULL ){
    xerror("glp_ios_selected_cuts: not called with a valid tree.\n");
  }
  if ( tree->reason != GLP_ICUTSELECT ){
    xerror("glp_ios_selected_cuts: not called during cut selected.\n");
  }

  pool = tree->local;
  if ( pool == NULL ){
    xerror("glp_ios_selected_cuts: called on a malformed tree.\n");
  }

  for (len = 0, j = 1, cut = pool->head; cut != NULL; cut = cut->next, j++)
  {  aux = cut->aux;
#ifdef CUT_DEBUG
     printf("glp_ios_selected_cuts: %d %p\n", j, cut);
#endif
     if ( aux != NULL )
     { s = aux->selected;
       if ( s >= 0 )
       {  len++;
          if (ords != NULL) { ords[len] = j; }
          if (sel != NULL)  { sel[len] = s; }
       }
     }
  }
  return len;
}
