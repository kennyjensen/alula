// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_SPARSE_H
#define ALULA_SPARSE_H

#include <stddef.h>

typedef struct AlulaSparse AlulaSparse;
typedef struct { int row, column; } AlulaSparseEntry;
typedef struct {
    int nonzeros, factor_nonzeros, factor_count, refactor_count;
    double max_relative_residual;
} AlulaSparseStats;

/* Canonicalizes duplicates. If supplied, slots receives one value index per input entry. */
AlulaSparse *alula_sparse_create(int n, const AlulaSparseEntry *entries, size_t count, int *slots);
double *alula_sparse_values(AlulaSparse *matrix);
int alula_sparse_nonzeros(const AlulaSparse *matrix);
void alula_sparse_product(const AlulaSparse *matrix, const double *x, double *out);
int alula_sparse_solve(AlulaSparse *matrix, const double *rhs, double *solution);
AlulaSparseStats alula_sparse_stats(const AlulaSparse *matrix);
void alula_sparse_free(AlulaSparse *matrix);

#endif
