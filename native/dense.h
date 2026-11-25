// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_DENSE_H
#define ALULA_DENSE_H
#include <stddef.h>
/* In-place pivoted LU. A is row-major. B contains one or two contiguous columns. */
int alula_dense_solve(double *a,double *b,size_t n,int columns);
/* Keep a pivoted LU for repeated solves. Factor layout is backend-specific.
 * Each RHS matrix here is row-major, including the one-column case. */
int alula_dense_factor(double *a,size_t n,int *pivots);
int alula_dense_apply_rows(const double *lu,size_t n,const int *pivots,double *rhs,size_t columns);
#endif
