// SPDX-License-Identifier: GPL-2.0-or-later
#include "dense.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef ALULA_LAPACK
/* Fortran LP64 LAPACK ABI. */
extern void dgetrf_(const int *, const int *, double *, const int *, int *, int *);
extern void dgetrs_(const char *, const int *, const int *, const double *,
    const int *, const int *, double *, const int *, int *);

int alula_dense_factor(double *a,size_t n,int *pivot)
{
    if (!a || !pivot || !n || n>INT_MAX) return 0;
    int size=(int)n, info;
    /* Present the physical equations as rows so LAPACK pivots on the same
     * columns as the scalar solver, including the trailing-edge condition. */
    for (size_t i=0;i<n;++i) for (size_t j=i+1;j<n;++j) {
        double value=a[i*n+j];
        a[i*n+j]=a[j*n+i];
        a[j*n+i]=value;
    }
    dgetrf_(&size, &size, a, &size, pivot, &info);
    return info==0;
}

int alula_dense_apply_rows(const double *lu,size_t n,const int *pivot,double *rhs,size_t columns)
{
    if (!lu || !pivot || !rhs || !n || n>INT_MAX || !columns || columns>INT_MAX || n>SIZE_MAX/sizeof(double)/columns) return 0;
    int size=(int)n, count=(int)columns, info;
    double *b=columns==1?rhs:malloc(n*columns*sizeof *b);
    if (!b) return 0;
    if (columns>1)
        for (size_t i=0;i<n;++i) for (size_t j=0;j<columns;++j) b[j*n+i]=rhs[i*columns+j];
    dgetrs_("N",&size,&count,lu,&size,pivot,b,&size,&info);
    if (columns>1) {
        for (size_t i=0;i<n;++i) for (size_t j=0;j<columns;++j) rhs[i*columns+j]=b[j*n+i];
        free(b);
    }
    return info==0;
}

int alula_dense_solve(double *a, double *b, size_t n, int columns)
{
    if (!a || !b || !n || n>INT_MAX || (columns!=1 && columns!=2)) return 0;
    int size=(int)n, info=-1;
    int *pivot=malloc(n*sizeof *pivot);
    if (!pivot) return 0;
    if (alula_dense_factor(a,n,pivot)) dgetrs_("N",&size,&columns,a,&size,pivot,b,&size,&info);
    free(pivot);
    return info == 0;
}
#else
/* The pivots in [first,last) are ready. Finish their upper rows in order,
 * then apply the same row operations to the rest of the matrix. Four
 * successive pivots share one load and store per trailing entry. */
static void update_trailing(double *a,size_t n,size_t first,size_t last)
{
    for (size_t i=first+1;i<n;++i) {
        double *row=a+i*n;
        size_t stop=i<last?i:last, k=first;
        for (;stop-k>=4;k+=4) {
            const double *p0=a+k*n, *p1=p0+n, *p2=p1+n, *p3=p2+n;
            double s0=row[k], s1=row[k+1], s2=row[k+2], s3=row[k+3];
            for (size_t j=last;j<n;++j) {
                double v=row[j]-s0*p0[j];
                v-=s1*p1[j];
                v-=s2*p2[j];
                row[j]=v-s3*p3[j];
            }
        }
        for (;k<stop;++k) {
            const double *pivot=a+k*n;
            double scale=row[k];
            for (size_t j=last;j<n;++j) row[j]-=scale*pivot[j];
        }
    }
}

static int factor(double *a,double *b,size_t n,int columns,int *pivots)
{
    size_t block = n < 96 ? n : 24;
    for (size_t first=0;first<n;first+=block) {
        size_t last=n-first<block?n:first+block;
        for (size_t k=first;k<last;++k) {
            size_t pivot = k;
            for (size_t i = k+1; i < n; ++i)
                if (fabs(a[i*n+k]) > fabs(a[pivot*n+k])) pivot = i;
            if (!(fabs(a[pivot*n+k]) > DBL_MIN)) return 0;
            if (pivots) pivots[k]=(int)pivot;
            if (pivot != k) {
                /* Carry this block's saved multipliers with the swapped rows. */
                for (size_t j=pivots?0:first;j<n;++j) {
                    double tmp = a[k*n+j];
                    a[k*n+j] = a[pivot*n+j];
                    a[pivot*n+j] = tmp;
                }
                for (int rhs = 0; rhs < columns; ++rhs) {
                    double *column = b+rhs*n;
                    double tmp = column[k];
                    column[k] = column[pivot];
                    column[pivot] = tmp;
                }
            }
            const double *row = a+k*n;
            for (size_t i = k+1; i < n; ++i) {
                double *target = a+i*n;
                double scale = target[k]/row[k];
                target[k]=scale;
                for (size_t j=k+1;j<last;++j) target[j]-=scale*row[j];
                if (columns) b[i] -= scale*b[k];
                if (columns==2) b[n+i] -= scale*b[n+k];
            }
        }
        if (last<n) update_trailing(a,n,first,last);
    }
    return 1;
}

int alula_dense_factor(double *a,size_t n,int *pivots)
{
    return a && pivots && n && n<=INT_MAX && factor(a,NULL,n,0,pivots);
}

/* Apply a triangular row to a batch of independent right-hand sides.
 * Accumulating four factors at a time cuts the target's load/store traffic
 * while preserving the order of subtraction within each right-hand side. */
static void subtract_rows(double *restrict target,const double *restrict rhs,
    const double *factor,size_t first,size_t end,size_t columns)
{
    size_t j=first;
    for (;end-j>=4;j+=4) {
        const double *r0=rhs+j*columns,*r1=r0+columns,*r2=r1+columns,*r3=r2+columns;
        double a=factor[j],b=factor[j+1],c=factor[j+2],d=factor[j+3];
        for (size_t k=0;k<columns;++k) {
            double v=target[k]-a*r0[k];
            v-=b*r1[k];
            v-=c*r2[k];
            target[k]=v-d*r3[k];
        }
    }
    for (;j<end;++j) {
        const double *row=rhs+j*columns;
        double scale=factor[j];
        for (size_t k=0;k<columns;++k) target[k]-=scale*row[k];
    }
}

int alula_dense_apply_rows(const double *lu,size_t n,const int *pivots,double *rhs,size_t columns)
{
    if (!lu || !pivots || !rhs || !n || n>INT_MAX || !columns || columns>INT_MAX || n>SIZE_MAX/sizeof(double)/columns) return 0;
    for (size_t i=0;i<n;++i) {
        size_t p=(size_t)pivots[i];
        if (p!=i) for (size_t c=0;c<columns;++c) {
            double tmp=rhs[i*columns+c];
            rhs[i*columns+c]=rhs[p*columns+c]; rhs[p*columns+c]=tmp;
        }
    }
    for (size_t i=0;i<n;++i) subtract_rows(rhs+i*columns,rhs,lu+i*n,0,i,columns);
    for (size_t i=n;i-- >0;) {
        double *target=rhs+i*columns;
        subtract_rows(target,rhs,lu+i*n,i+1,n,columns);
        for (size_t c=0;c<columns;++c) target[c]/=lu[i*n+i];
    }
    return 1;
}

int alula_dense_solve(double *a,double *b,size_t n,int columns)
{
    if (!a || !b || !n || n>INT_MAX || (columns!=1 && columns!=2) || !factor(a,b,n,columns,NULL)) return 0;
    for (size_t i = n; i-- > 0;) {
        for (size_t j = i+1; j < n; ++j) {
            b[i] -= a[i*n+j]*b[j];
            if (columns==2) b[n+i] -= a[i*n+j]*b[n+j];
        }
        b[i] /= a[i*n+i];
        if (columns==2) b[n+i] /= a[i*n+i];
    }
    return 1;
}
#endif
