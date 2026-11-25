// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_AIRFOIL_COUPLED_H
#define ALULA_AIRFOIL_COUPLED_H

#include "airfoil_bl.h"

typedef enum {
    ALULA_COUPLED_OK, ALULA_COUPLED_INVALID, ALULA_COUPLED_ALLOCATION,
    ALULA_COUPLED_BL_START, ALULA_COUPLED_LINE_SEARCH, ALULA_COUPLED_ITERATION_LIMIT
} AlulaCoupledStatus;

typedef struct {
    double reynolds,ncrit,displacement_tolerance,euler_tolerance,relaxation;
    const double (*trips)[2];
    int max_iterations,euler_iterations;
} AlulaCoupledControls;

typedef struct {
    AlulaCoupledStatus status;
    int iterations,backtracks,euler_updates;
    double displacement_residual,bl_residual;
    AlulaBodySolveResult euler;
    AlulaAirfoilBL *profiles; /* Owned by the caller. */
} AlulaCoupledResult;

/* Converge the displacement feedback by alternating Euler and BL solves.
 * The caller supplies an already converged Euler grid. Displacement
 * tolerance is in reference chords. Failed attempts restore the last
 * accepted geometry, Euler state, and displacement. */
AlulaCoupledResult alula_airfoil_coupled_solve(AlulaAirfoilGrid *grid,const AlulaCoupledControls *controls);

#endif
