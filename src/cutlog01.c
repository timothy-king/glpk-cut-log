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
     xfree(aux->mir_subst);
     xfree(aux->mir_vlb_rows);
     xfree(aux->mir_vub_rows);
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

static void cut_set_aux_mir(IOSAUX *aux, double delta, int m, int n,
                            const char cset[], const char subst[],
                            const int vlbrs[], const int vubrs[]){
  int i;
  xassert( aux != NULL );
  if ( aux->mir_cset != NULL )
  {  xfree(aux->mir_cset);
     xfree(aux->mir_subst);
     xfree(aux->mir_vlb_rows);
     xfree(aux->mir_vub_rows);
  }

  aux->mir_cset     = xcalloc(1+n+m, sizeof(char));
  aux->mir_subst    = xcalloc(1+n+m, sizeof(char));
  aux->mir_vlb_rows = xcalloc(1+n+m, sizeof(int));
  aux->mir_vub_rows = xcalloc(1+n+m, sizeof(int));

  for ( i = 1; i <= n+m; i++)
  {  aux->mir_cset[i] = cset[i];
     aux->mir_subst[i] = subst[i];
     aux->mir_vlb_rows[i] = vlbrs[i];
     aux->mir_vub_rows[i] = vubrs[i];
  }

  aux->mir_delta = delta;
}

void ios_cut_set_aux_mir(glp_tree *T, int ord, double delta,
                         const char cset[], const char subst[],
                         const int vlbrs[], const int vubrs[]){
  int m, n;
  IOSCUT *cut;
  glp_prob *mip;
  mip = T->mip;
  m = mip->m;
  n = mip->n;
  cut = ios_find_row(T->local, ord);
  xassert(cut != NULL);
  cut_set_aux_mir(cut->aux, delta, m, n, cset, subst, vlbrs, vubrs);
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


void glp_ios_cut_get_mir_subst(glp_tree *tree, int ord, char subst[]);
/* gets mir cut substition information. */
void glp_ios_cut_get_mir_virtual_rows(glp_tree *tree, int ord,
                                      int vlb[], int vub[]);
/* gets mir cut virtual bounds rows. */

void glp_ios_cut_get_mir_cset(glp_tree *tree, int ord, char *cset){
  glp_prob *mip;
  IOSCUT *cut;
  IOSAUX *aux;
  int j, n, m;
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
  mip = tree->mip;
  m = mip->m;
  n = mip->n;

  if( cset != NULL ){
    for ( j=1; j <= n+m; j++ ){
      cset[j] = (aux->mir_cset == NULL) ? 0 : aux->mir_cset[j];
    }
  }
}
void glp_ios_cut_get_mir_subst(glp_tree *tree, int ord, char *subst){
  glp_prob *mip;
  IOSCUT *cut;
  IOSAUX *aux;
  int j, n, m;
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
  mip = tree->mip;
  m = mip->m;
  n = mip->n;

  if( subst != NULL ){
    for ( j=1; j <= n+m; j++ ){
      subst[j] = (aux->mir_subst == NULL) ? 0 : aux->mir_subst[j];
    }
  }
}
void glp_ios_cut_get_mir_virtual_rows(glp_tree *tree, int ord, int vlb_rows[], int vub_rows[]){
  glp_prob *mip;
  IOSCUT *cut;
  IOSAUX *aux;
  int j, n, m;
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
  mip = tree->mip;
  m = mip->m;
  n = mip->n;

  for ( j=1; j <= n+m; j++ ){
    vlb_rows[j] = (aux->mir_vlb_rows == NULL) ? 0 : aux->mir_vlb_rows[j];
    vub_rows[j] = (aux->mir_vub_rows == NULL) ? 0 : aux->mir_vub_rows[j];
  }
}
double glp_ios_cut_get_mir_delta(glp_tree *tree, int ord){
  glp_prob *mip;
  IOSCUT *cut;
  IOSAUX *aux;
  int j, n, m;
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
  return aux->mir_delta;
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

int glp_ios_branch_log(glp_tree *tree, double *br_val, int* parent, int* dn, int* up){
  IOSNPD *node;
  int br_result, br_var;
  int p, d, u;
  double v;
  glp_prob *mip;
  if ( tree == NULL ){
    xerror("glp_ios_branch_log: not called with a valid tree \n");
  }
  if ( tree->reason != GLP_LI_BRANCH ){
    xerror("glp_ios_branch_log: not called during GLP_LI_BRANCH \n");
  }
  mip = tree->mip;
  if ( mip == NULL ){
    xerror("glp_ios_branch_log: not called with a valid tree\n");
  }
  br_result = tree->br_result;
  br_var = tree->br_var;
  switch(br_result){
  case 0:
    p = tree->br_node;
    node = tree->slot[p].node;
    break;
  case 1:
  case 2:
    node = tree->curr;
    p = node->p;
    break;
  default:
    xerror("glp_ios_branch_log: br_result is not properly set \n");
  }
  if( node == NULL ){
    xerror("glp_ios_branch_log: not called with a valid tree \n");
  }
  switch(br_result){
  case 0:
    v = node->br_val;
    d = tree->dn_child;
    u = tree->up_child;
    break;
  case 1:
    v = mip->col[br_var]->prim;
    if(tree->br_to_up){
      d = -1;
      u = p;
    }else{
      d = p;
      u = -1;
    }
    break;
  case 2:
    v = mip->col[br_var]->prim;
    d = -1;
    u = -1;
    break;
  default:
    xerror("glp_ios_branch_log: not called with a valid tree \n");
  }

  if(br_val != NULL){ *br_val = v; }
  if(parent != NULL){ *parent = p; }
  if(dn != NULL){ *dn = d; }
  if(up != NULL){ *up = u; }

  return br_var;
}

int glp_ios_node_ord(glp_tree *tree, int p){
  IOSNPD *node;
  if ( tree == NULL ){
    xerror("glp_ios_node_ord: not called with a valid tree.\n");
  }
  if (!(1 <= p && p <= tree->nslots)){
    xerror("glp_ios_node_ord: not called with a valid p.\n");
  }
  node = tree->slot[p].node;
  return node->ord;
}

void ios_cb_rows_deleted(glp_tree *T, int n, const int* rows){
  if (T->parm->cb_func != NULL)
  {
    xassert(T->reason == 0);
    xassert(T->deleting_rows == NULL);
    xassert(T->num_deleting_rows == 0);
    T->num_deleting_rows = n;
    T->deleting_rows = rows;

    T->reason = GLP_LI_DELROW;
    T->parm->cb_func(T, T->parm->cb_info);
    T->reason = 0;
    T->num_deleting_rows = 0;
    T->deleting_rows = NULL;
  }
}

int glp_ios_rows_deleted(glp_tree *tree, int* rows){
  if ( tree == NULL ){
    xerror("glp_ios_rows_deleted: not called with a valid tree.\n");
  }
  if ( tree->reason != GLP_LI_DELROW ){
    xerror("glp_ios_rows_deleted: not called with a valid reason.\n");
  }

  int j;
  if(rows != NULL){
    for(j=1; j <= tree->num_deleting_rows; j++){
      rows[j] = tree->deleting_rows[j];
    }
  }
  return tree->num_deleting_rows;
}
