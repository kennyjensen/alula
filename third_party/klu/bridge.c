// SPDX-License-Identifier: GPL-2.0-or-later
// Project-owned ABI wrapper. The SuiteSparse sources are unmodified.
#include <stdlib.h>
#include "klu.h"
typedef struct { klu_symbolic *symbolic; klu_numeric *numeric; klu_common common; } Factor;
static int last_status;
int mses_klu_status(void) { return last_status; }
void mses_klu_free(Factor *f) {
  if (!f) return;
  if (f->numeric) klu_free_numeric(&f->numeric, &f->common);
  if (f->symbolic) klu_free_symbolic(&f->symbolic, &f->common);
  free(f);
}
Factor *mses_klu_factor_with_pivot(int n, int *ap, int *ai, double *ax, int ordering, double pivot_tolerance) {
  if (!(pivot_tolerance > 0 && pivot_tolerance <= 1)) { last_status=KLU_INVALID; return NULL; }
  Factor *f = calloc(1, sizeof(Factor));
  if (!f) { last_status=KLU_OUT_OF_MEMORY; return NULL; }
  klu_defaults(&f->common); f->common.ordering=ordering;
  f->common.tol=pivot_tolerance;
  f->symbolic=klu_analyze(n,ap,ai,&f->common);
  if(f->symbolic) f->numeric=klu_factor(ap,ai,ax,f->symbolic,&f->common);
  last_status=f->common.status;
  if(!f->numeric || last_status!=KLU_OK) { mses_klu_free(f); return NULL; }
  return f;
}
// Preserve the original ABI and upstream diagonal preference.
Factor *mses_klu_factor(int n, int *ap, int *ai, double *ax, int ordering) {
  return mses_klu_factor_with_pivot(n,ap,ai,ax,ordering,0.001);
}
int mses_klu_solve(Factor *f, double *rhs) {
  int ok=klu_solve(f->symbolic,f->numeric,f->symbolic->n,1,rhs,&f->common);
  last_status=f->common.status; return ok;
}
int mses_klu_nnz(Factor *f) { return f->numeric->lnz+f->numeric->unz+f->numeric->nzoff; }

// Explicit symbolic ordering; KLU owns original RHS and solution maps.
Factor *mses_klu_factor_given(int n, int *ap, int *ai, double *ax,
                            int *p, int *q, double pivot, int btf) {
  if (!(pivot > 0 && pivot <= 1) || (btf != 0 && btf != 1)) {
    last_status=KLU_INVALID; return NULL;
  }
  Factor *f=calloc(1,sizeof(Factor));
  if (!f) { last_status=KLU_OUT_OF_MEMORY; return NULL; }
  klu_defaults(&f->common); f->common.tol=pivot; f->common.btf=btf;
  f->symbolic=klu_analyze_given(n,ap,ai,p,q,&f->common);
  if (f->symbolic) f->numeric=klu_factor(ap,ai,ax,f->symbolic,&f->common);
  last_status=f->common.status;
  if (!f->numeric || last_status!=KLU_OK) { mses_klu_free(f); return NULL; }
  return f;
}
int mses_klu_btf(Factor *f) { return f->symbolic->do_btf; }
int mses_klu_blocks(Factor *f) { return f->symbolic->nblocks; }
