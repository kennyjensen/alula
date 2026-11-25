// SPDX-License-Identifier: GPL-2.0-or-later
#include "sparse.h"
#include "klu.h"
#include "klu_version.h" /* GET_POINTER for the vendored real, int-indexed factors. */

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct AlulaSparse {
    int n, nnz;
    int *starts, *rows;
    double *values, *work, *scale, *corrections;
    double (*lanes)[ALULA_SPARSE_BATCH];
    long double *wide;
    klu_common common;
    klu_symbolic *symbolic, *matched_order;
    klu_numeric *numeric;
    AlulaSparseStats stats;
};

typedef struct { int row, column, original; } Key;
static int order(const void *x, const void *y)
{
    const Key *a=x, *b=y;
    return a->column == b->column ? a->row-b->row : a->column-b->column;
}

/* KLU's block triangular ordering can permute the equation rows even when
 * it finds just one block. For expensive factorizations, also estimate AMD
 * on the original row/column pairing. Keep the matched order as a numeric
 * fallback if the cheaper symbolic prediction leads to poor pivots. */
static void compare_single_block_order(AlulaSparse *m)
{
    if (m->symbolic->nblocks!=1 || m->symbolic->est_flops<1e7) return;
    m->common.btf=0;
    klu_symbolic *aligned=klu_analyze(m->n,m->starts,m->rows,&m->common);
    m->common.btf=1;
    if (aligned && aligned->est_flops<m->symbolic->est_flops && aligned->lnz<=m->symbolic->lnz) {
        m->matched_order=m->symbolic;
        m->symbolic=aligned;
    } else {
        klu_free_symbolic(&aligned,&m->common);
    }
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
    compare_single_block_order(m);
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

/* Compensation corrects addition error, but products rounded to double can
 * still swamp the Newton RHS when large terms cancel. Refine those solves
 * using wider products and sums of the original entries. */
static double wide_residual(AlulaSparse *m,const double *rhs,const double *x)
{
    double bnorm=0,error=0;
    for (int i=0;i<m->n;++i) {
        m->wide[i]=rhs[i];
        bnorm=fmax(bnorm,fabs(rhs[i]));
    }
    for (int j=0;j<m->n;++j) {
        if (!isfinite(x[j])) return INFINITY;
        for (int k=m->starts[j];k<m->starts[j+1];++k)
            m->wide[m->rows[k]]-=(long double)m->values[k]*x[j];
    }
    for (int i=0;i<m->n;++i) {
        m->work[i]=(double)m->wide[i];
        if (!isfinite(m->work[i])) return INFINITY;
        error=fmax(error,fabs(m->work[i]));
    }
    return bnorm>0?error/bnorm:error;
}

static const double linear_tolerance=1e-10;

#if defined(FP_FAST_FMA) && LDBL_MANT_DIG == 113
/* TwoSum and FMA recover the rounding errors of the addition and product.
 * On targets with hardware FMA and software binary128, use this cheaper
 * residual for correction steps. Verify accepted solutions in binary128. */
static double fused_residual(AlulaSparse *m,const double *rhs,const double *x)
{
    double bnorm=0,error=0;
    memset(m->scale,0,(size_t)m->n*sizeof *m->scale);
    memcpy(m->work,rhs,(size_t)m->n*sizeof *m->work);
    for (int j=0;j<m->n;++j) {
        if (!isfinite(x[j])) return INFINITY;
        bnorm=fmax(bnorm,fabs(rhs[j]));
        for (int k=m->starts[j];k<m->starts[j+1];++k) {
            int i=m->rows[k];
            double a=m->work[i], b=-m->values[k]*x[j];
            double sum=a+b, z=sum-a;
            double sum_error=(a-(sum-z))+(b-z);
            double product_error=fma(-m->values[k],x[j],-b);
            m->scale[i]+=sum_error+product_error;
            m->work[i]=sum;
        }
    }
    for (int i=0;i<m->n;++i) {
        m->work[i]+=m->scale[i];
        if (!isfinite(m->work[i])) return INFINITY;
        error=fmax(error,fabs(m->work[i]));
    }
    return bnorm>0?error/bnorm:error;
}
#endif

static double refined_residual(AlulaSparse *m,const double *rhs,const double *x)
{
#if defined(FP_FAST_FMA) && LDBL_MANT_DIG == 113
    double e=fused_residual(m,rhs,x);
    if (isfinite(e) && e>linear_tolerance) return e;
#endif
    return wide_residual(m,rhs,x);
}

static int prepare_wide(AlulaSparse *m)
{
    if (LDBL_MANT_DIG<=DBL_MANT_DIG) return 0;
    if (!m->wide) m->wide=malloc((size_t)m->n*sizeof *m->wide);
    return m->wide!=NULL;
}

static int correct(AlulaSparse *m,double *x)
{
    if (!klu_solve(m->symbolic,m->numeric,m->n,1,m->work,&m->common)) return 0;
    for (int i=0;i<m->n;++i) x[i]+=m->work[i];
    return 1;
}

static void record_solve(AlulaSparse *m,double e)
{
    m->stats.max_relative_residual=fmax(m->stats.max_relative_residual,e);
    m->stats.factor_nonzeros=m->numeric->lnz+m->numeric->unz-m->n+m->numeric->nzoff;
}

static double reuse_residual(AlulaSparse *m,const double *rhs,const double *x)
{
    double e=relative_residual(m,rhs,x);
    return e>linear_tolerance && e<1e-7 && prepare_wide(m)?refined_residual(m,rhs,x):e;
}

/* Use the previous factors as a preconditioner for the current matrix.
 * A small number of contracting corrections can be much cheaper than
 * refactorization. The residual always uses the newly assembled matrix. */
static int reuse_factors(AlulaSparse *m,const double *rhs,double *x)
{
    memcpy(x,rhs,(size_t)m->n*sizeof *x);
    if (!klu_solve(m->symbolic,m->numeric,m->n,1,x,&m->common)) return 0;
    double previous=INFINITY;
    for (int i=0;i<=6;++i) {
        double e=reuse_residual(m,rhs,x);
        if (e<=linear_tolerance) {
            record_solve(m,e);
            ++m->stats.factor_reuses;
            return 1;
        }
        if (i==6 || !isfinite(e) || e>.5*previous || !correct(m,x)) break;
        previous=e;
    }
    return 0;
}

static int factor_and_solve(AlulaSparse *m,const double *rhs,double *x)
{
    int reuse=m->numeric != NULL;
    if (reuse && reuse_factors(m,rhs,x)) return 1;
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
        int use_wide=e>linear_tolerance && prepare_wide(m);
        if (use_wide) e=refined_residual(m,rhs,x);
        for (int refine=0; refine<3 && isfinite(e) && e>linear_tolerance; ++refine) {
            if (!correct(m,x)) break;
            e=use_wide?refined_residual(m,rhs,x):relative_residual(m,rhs,x);
        }
#if defined(FP_FAST_FMA) && LDBL_MANT_DIG == 113
        if (use_wide && e>linear_tolerance) e=wide_residual(m,rhs,x);
#endif
        if (e<=linear_tolerance) {
            record_solve(m,e);
            return 1;
        }
        m->stats.failed_relative_residual=e;
    }
    /* The next caller gets a fresh pivot search after an unusable solve. */
    klu_free_numeric(&m->numeric,&m->common);
    return 0;
}

static int finite_system(AlulaSparse *m,const double *rhs,int count)
{
    m->stats.failed_relative_residual=0;
    for (int k=0;k<m->nnz;++k) if (!isfinite(m->values[k])) return 0;
    for (size_t i=0;i<(size_t)m->n*count;++i) if (!isfinite(rhs[i])) return 0;
    return 1;
}

int alula_sparse_solve(AlulaSparse *m,const double *rhs,double *x)
{
    if (!m || !rhs || !x || rhs==x || !finite_system(m,rhs,1)) return 0;
    if (factor_and_solve(m,rhs,x)) return 1;
    if (!m->matched_order) return 0;
    klu_free_symbolic(&m->symbolic,&m->common);
    m->symbolic=m->matched_order;
    m->matched_order=NULL;
    return factor_and_solve(m,rhs,x);
}

/* Each factor entry updates a contiguous row of independent right-hand sides.
 * Read KLU's current factors directly, including after a numerical refactor. */
static void substitute_column(double (*restrict lanes)[ALULA_SPARSE_BATCH],int column,
    Unit *block,int *offsets,int *lengths,const double *diagonal)
{
    double pivot[ALULA_SPARSE_BATCH];
    memcpy(pivot,lanes[column],sizeof pivot);
    if (diagonal) {
        for (int j=0;j<ALULA_SPARSE_BATCH;++j) pivot[j]/=diagonal[column];
        memcpy(lanes[column],pivot,sizeof pivot);
    }
    int any=0;
    for (int j=0;j<ALULA_SPARSE_BATCH;++j) any|=(pivot[j]!=0);
    if (!any) return;
    int length,*rows;
    double *values;
    GET_POINTER(block,offsets,lengths,rows,values,column,length);
    for (int k=0;k<length;++k) {
        double *row=lanes[rows[k]];
        for (int j=0;j<ALULA_SPARSE_BATCH;++j) row[j]-=values[k]*pivot[j];
    }
}

static int solve_factor_batch(AlulaSparse *m,int count,double *vectors)
{
    size_t n=(size_t)m->n;
    int single_block=count>4 && m->n>1 && m->symbolic->nblocks==1;
    if (single_block && !m->lanes) m->lanes=calloc(n,sizeof *m->lanes);
    if (!single_block || !m->lanes)
        return klu_solve(m->symbolic,m->numeric,m->n,count,vectors,&m->common);

    klu_numeric *f=m->numeric;
    for (int i=0;i<m->n;++i) {
        double scale=f->Rs?f->Rs[i]:1;
        for (int j=0;j<count;++j) m->lanes[i][j]=vectors[j*n+f->Pnum[i]]/scale;
        for (int j=count;j<ALULA_SPARSE_BATCH;++j) m->lanes[i][j]=0;
    }
    for (int i=0;i<m->n;++i)
        substitute_column(m->lanes,i,f->LUbx[0],f->Lip,f->Llen,NULL);
    for (int i=m->n-1;i>=0;--i)
        substitute_column(m->lanes,i,f->LUbx[0],f->Uip,f->Ulen,f->Udiag);
    for (int i=0;i<m->n;++i)
        for (int j=0;j<count;++j) vectors[j*n+m->symbolic->Q[i]]=m->lanes[i][j];
    return 1;
}

int alula_sparse_solve_batch(AlulaSparse *m,int count,const double *rhs,double *x)
{
    if (!m || count<1 || count>ALULA_SPARSE_BATCH || !rhs || !x || rhs==x) return 0;
    if (count==1) return alula_sparse_solve(m,rhs,x);
    size_t n=(size_t)m->n,bytes=n*sizeof *x;
    if (!m->numeric) {
        if (!alula_sparse_solve(m,rhs,x)) return 0;
        return alula_sparse_solve_batch(m,count-1,rhs+n,x+n);
    }
    if (!finite_system(m,rhs,count)) return 0;
    if (!m->corrections) m->corrections=malloc(ALULA_SPARSE_BATCH*bytes);
    int active[ALULA_SPARSE_BATCH],verified[ALULA_SPARSE_BATCH]={0};
    double previous[ALULA_SPARSE_BATCH];
    memcpy(x,rhs,(size_t)count*bytes);
    int pending=count;
    for (int j=0;j<count;++j) {active[j]=j;previous[j]=INFINITY;}
    for (int pass=0;pending && pass<=6;++pass) {
        double *vectors=pass?m->corrections:x;
        if (!solve_factor_batch(m,pending,vectors)) break;
        int next=0;
        for (int k=0;k<pending;++k) {
            int j=active[k];
            double *column=x+j*n;
            if (pass) for (size_t i=0;i<n;++i) column[i]+=vectors[k*n+i];
            double e=reuse_residual(m,rhs+j*n,column);
            if (e<=linear_tolerance) {
                verified[j]=1;
                record_solve(m,e);
                ++m->stats.factor_reuses;
            } else if (pass<6 && m->corrections && isfinite(e) && e<=.5*previous[j]) {
                /* Only unconverged columns need another triangular solve.
                 * Compact toward the front after consuming this correction. */
                previous[j]=e;
                active[next]=j;
                memcpy(m->corrections+next*n,m->work,bytes);
                ++next;
            }
        }
        pending=next;
    }
    for (int j=0;j<count;++j)
        if (!verified[j] && !alula_sparse_solve(m,rhs+j*n,x+j*n)) return 0;
    return 1;
}

void alula_sparse_free(AlulaSparse *m)
{
    if (!m) return;
    klu_free_numeric(&m->numeric,&m->common);
    klu_free_symbolic(&m->symbolic,&m->common);
    klu_free_symbolic(&m->matched_order,&m->common);
    free(m->starts); free(m->rows); free(m->values); free(m->work); free(m->scale);
    free(m->wide); free(m->corrections); free(m->lanes); free(m);
}
