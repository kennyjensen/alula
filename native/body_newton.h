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
    /* Optional faster recentering with the same contract. Newton calls it
     * immediately after a successful decode of the same trial vector. */
    int (*rebase_decoded)(void *context,double *packed);
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
    int rejected_chart, rejected_flow, rejected_residual;
    AlulaBodyStatus last_rejected_flow;
    double residual_by_family[ALULA_BODY_FAMILIES];
} AlulaBodySolveResult;

AlulaBodyPhysicalLayout alula_body_physical_layout(AlulaBodySize size);
/* The body must outlive the solver. The chart context remains caller-owned. */
AlulaBodyNewton *alula_body_newton_create(AlulaBody *body,const AlulaBodyChart *chart);
/* Product and residual use the original equation order. */
int alula_body_newton_product(AlulaBodyNewton *solver,const double *packed,
    const double *direction,double *product);
/* Linearize at a converged Euler state, then predict the packed response
 * to changes in prescribed displacement on the same coordinate chart.
 * The Jacobian is shared by all response columns until the next Euler residual. */
int alula_body_newton_linearize(AlulaBodyNewton *s,const double *packed);
/* Also compute the inhomogeneous Euler correction -J^-1 R at a trial state. */
int alula_body_newton_correction(AlulaBodyNewton *s,const double *packed,double *correction);
int alula_body_newton_response(AlulaBodyNewton *s,const double *packed,double *prediction);
/* Collect forcing vectors for independent perturbations of displacement,
 * then solve them on the same linearized chart. Arrays are column-major. */
int alula_body_newton_displacement_rhs(AlulaBodyNewton *s,const double *packed,double *rhs);
int alula_body_newton_predict_batch(AlulaBodyNewton *s,const double *packed,int count,
    const double *rhs,double *predictions);
/* Prepare the body's work arrays for an exact decoded physical state.
 * The vector uses AlulaBodyPhysicalLayout. */
int alula_body_newton_prepare(AlulaBodyNewton *solver,const double *physical);
AlulaBodySolveResult alula_body_newton_solve(AlulaBodyNewton *solver,double *packed,
    int max_iterations,double tolerance);
typedef void (*AlulaEulerProgress)(void *context,int iteration,double residual);
/* Invoked on the evaluated seed and each accepted update, while its cell
 * pressures are current. The observer must not modify or reenter the solve. */
AlulaBodySolveResult alula_body_newton_solve_observed(AlulaBodyNewton *solver,double *packed,
    int max_iterations,double tolerance,AlulaEulerProgress progress,void *context);
void alula_body_newton_free(AlulaBodyNewton *solver);

#endif
