// SPDX-License-Identifier: GPL-2.0-or-later
#include "sparse.h"
#include "klu.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct AlulaSparse {
    int n, nnz;
    int *starts, *rows;
    double *values, *work, *scale;
    klu_common common;
    klu_symbolic *symbolic;
    klu_numeric *numeric;
    AlulaSparseStats stats;
};

typedef struct { int row, column, original; } Key;
static int order(const void *x, const void *y)
{
    const Key *a=x, *b=y;
    return a->column == b->column ? a->row-b->row : a->column-b->column;
}

AlulaSparse *alula_sparse_create(int n, const AlulaSparseEntry *entries, size_t count, int *slots)
{
    if (!entries || n < 1 || count < 1 || count > INT_MAX) return NULL;
    AlulaSparse *m=calloc(1,sizeof *m);
    Key *keys=calloc(count,sizeof *keys);
    if (!m || !keys) { free(m); free(keys); return NULL; }
    klu_defaults(&m->common);
    m->n=n;
    m->starts=calloc((size_t)n+1,sizeof *m->starts);
    m->rows=calloc(count,sizeof *m->rows);
    m->values=calloc(count,sizeof *m->values);
    m->work=calloc((size_t)n,sizeof *m->work);
    m->scale=calloc((size_t)n,sizeof *m->scale);
    if (!m->starts || !m->rows || !m->values || !m->work || !m->scale) goto fail;
    for (int i=0; i<(int)count; ++i) {
        int row=entries[i].row, col=entries[i].column;
        if (row < 0 || row >= n || col < 0 || col >= n) goto fail;
        keys[i]=(Key){row,col,i};
    }
    qsort(keys,count,sizeof *keys,order);
    int k=0;
    for (int i=0; i<(int)count; ++i) {
        if (!i || order(keys+i-1,keys+i)) {
            m->starts[keys[i].column+1]++;
            m->rows[k++]=keys[i].row;
        }
        if (slots) slots[keys[i].original]=k-1;
    }
    m->nnz=m->stats.nonzeros=k;
    for (int i=1; i<=n; ++i) m->starts[i]+=m->starts[i-1];
    m->symbolic=klu_analyze(n,m->starts,m->rows,&m->common);
    if (!m->symbolic) goto fail;
    free(keys);
    return m;
fail:
    free(keys);
    alula_sparse_free(m);
    return NULL;
}

double *alula_sparse_values(AlulaSparse *m) { return m->values; }
int alula_sparse_nonzeros(const AlulaSparse *m) { return m->nnz; }
AlulaSparseStats alula_sparse_stats(const AlulaSparse *m) { return m->stats; }

void alula_sparse_product(const AlulaSparse *m, const double *x, double *out)
{
    memset(out,0,(size_t)m->n*sizeof *out);
    for (int j=0; j<m->n; ++j)
        for (int k=m->starts[j]; k<m->starts[j+1]; ++k)
            out[m->rows[k]]+=m->values[k]*x[j];
}

/* Recompute the residual on the unfactored matrix with compensated sums.
 * Keep the correction RHS in work for iterative refinement. */
static double relative_residual(AlulaSparse *m, const double *rhs, const double *x)
{
    double bnorm=0, error=0;
    memset(m->scale,0,(size_t)m->n*sizeof *m->scale);
    memcpy(m->work,rhs,(size_t)m->n*sizeof *rhs);
    for (int j=0; j<m->n; ++j) {
        if (!isfinite(x[j])) return INFINITY;
        bnorm=fmax(bnorm,fabs(rhs[j]));
        for (int k=m->starts[j]; k<m->starts[j+1]; ++k) {
            int i=m->rows[k];
            double term=-m->values[k]*x[j]-m->scale[i], next=m->work[i]+term;
            m->scale[i]=(next-m->work[i])-term;
            m->work[i]=next;
        }
    }
    for (int i=0; i<m->n; ++i) {
        if (!isfinite(m->work[i])) return INFINITY;
        error=fmax(error,fabs(m->work[i]));
    }
    return bnorm > 0 ? error/bnorm : error;
}

int alula_sparse_solve(AlulaSparse *m, const double *rhs, double *x)
{
    if (!m || !rhs || !x || rhs == x) return 0;
    for (int k=0; k<m->nnz; ++k) if (!isfinite(m->values[k])) return 0;
    for (int i=0; i<m->n; ++i) if (!isfinite(rhs[i])) return 0;
    int reuse=m->numeric != NULL;
    for (int attempt=0; attempt<=reuse; ++attempt) {
        int refactor=reuse && attempt == 0;
        if (refactor) {
            m->stats.refactor_count++;
            if (!klu_refactor(m->starts,m->rows,m->values,m->symbolic,m->numeric,&m->common)) continue;
        } else {
            klu_free_numeric(&m->numeric,&m->common);
            m->stats.factor_count++;
            m->numeric=klu_factor(m->starts,m->rows,m->values,m->symbolic,&m->common);
            if (!m->numeric) return 0;
        }
        if (m->common.status != KLU_OK) continue;
        memcpy(x,rhs,(size_t)m->n*sizeof *x);
        if (!klu_solve(m->symbolic,m->numeric,m->n,1,x,&m->common)) continue;
        double e=relative_residual(m,rhs,x);
        for (int refine=0; refine<3 && isfinite(e) && e>1e-10; ++refine) {
            if (!klu_solve(m->symbolic,m->numeric,m->n,1,m->work,&m->common)) break;
            for (int i=0; i<m->n; ++i) x[i]+=m->work[i];
            e=relative_residual(m,rhs,x);
        }
        if (e <= 1e-10) {
            m->stats.max_relative_residual=fmax(m->stats.max_relative_residual,e);
            m->stats.factor_nonzeros=m->numeric->lnz+m->numeric->unz-m->n+m->numeric->nzoff;
            return 1;
        }
    }
    /* The next caller gets a fresh pivot search after an unusable solve. */
    klu_free_numeric(&m->numeric,&m->common);
    return 0;
}

void alula_sparse_free(AlulaSparse *m)
{
    if (!m) return;
    klu_free_numeric(&m->numeric,&m->common);
    klu_free_symbolic(&m->symbolic,&m->common);
    free(m->starts); free(m->rows); free(m->values); free(m->work); free(m->scale);
    free(m);
}
