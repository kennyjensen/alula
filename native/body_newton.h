// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_BODY_NEWTON_H
#define ALULA_BODY_NEWTON_H

#include "body.h"
#include "sparse.h"

typedef struct AlulaBodyNewton AlulaBodyNewton;

/* An owned output buffer passed to the chart callback. Layout is point x/y
 * pairs, density, mass, circulation/source/doublets, then wake widths. */
typedef struct {
    int density, mass, strength, wake, size;
} AlulaBodyPhysicalLayout;

typedef int (*AlulaBodyDecode)(void *context,const double *packed,double *physical,double *derivatives);

typedef struct {
    /* CSR derivative support from each physical scalar to packed columns.
     * Starts has physical_layout.size + 1 entries, columns is sorted and
     * unique per scalar. Decode fills the corresponding derivative values
     * when derivatives is non-NULL. Neither pointer is retained. */
    const int *starts,*columns;
    AlulaBodyDecode decode;
    void *context;
    /* Optional permutation from original residual rows to KLU rows.
     * Pairing rows with their corresponding unknown reduces fill in AMD. */
    const int *equation_order;
    /* Optional accepted-step recentering. Must preserve physical state and
     * the derivative support. On failure, leave the chart unchanged. */
    int (*rebase)(void *context,double *packed);
} AlulaBodyChart;

typedef enum {
    ALULA_BODY_SOLVED, ALULA_BODY_ITERATION_LIMIT, ALULA_BODY_BAD_STATE,
    ALULA_BODY_LINEAR_FAILURE, ALULA_BODY_LINE_SEARCH_FAILURE, ALULA_BODY_INVALID_FINAL_GRID
} AlulaBodySolveStatus;

typedef struct {
    AlulaBodySolveStatus status;
    int iterations, backtracks, colors;
    double norm, max_mach;
    AlulaSparseStats linear;
} AlulaBodySolveResult;

AlulaBodyPhysicalLayout alula_body_physical_layout(AlulaBodySize size);
/* The body must outlive the solver. The chart context remains caller-owned. */
AlulaBodyNewton *alula_body_newton_create(AlulaBody *body,const AlulaBodyChart *chart);
/* Product and residual use the original equation order. */
int alula_body_newton_product(AlulaBodyNewton *solver,const double *packed,
    const double *direction,double *product);
AlulaBodySolveResult alula_body_newton_solve(AlulaBodyNewton *solver,double *packed,
    int max_iterations,double tolerance);
void alula_body_newton_free(AlulaBodyNewton *solver);

#endif
