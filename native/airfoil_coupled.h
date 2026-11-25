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
    double displacement_scale; /* Optional thickness scaling in reference chords. */
    /* Accepted outer iterates. The grid and profiles are borrowed during the callback. */
    void (*progress)(AlulaAirfoilGrid *,const AlulaAirfoilBL *,int iteration,double residual,void *context);
    void *progress_context;
} AlulaCoupledControls;

typedef struct {
    AlulaCoupledStatus status;
    int iterations,backtracks,euler_updates,used_interaction;
    int used_local_interaction,euler_grid_repairs,euler_rejections,bl_rejections,merit_rejections;
    int coupled_newton_updates; /* Accepted fallback steps, including work rolled back on failure. */
    double displacement_residual,bl_residual,edge_residual;
    AlulaBodySolveResult euler;
    AlulaAirfoilBL *profiles; /* Owned by the caller. */
    int mach_stages,mach_retries;
} AlulaCoupledResult;

/* Converge the displacement feedback by alternating Euler and BL solves.
 * The caller supplies an already converged Euler grid. Displacement
 * tolerance is in reference chords, and gives the edge matching tolerance
 * in freestream speed units. A stalled interaction march can switch to
 * a Newton solve of the BL stations. Failed attempts restore the last accepted
 * geometry, Euler state, and displacement. */
AlulaCoupledResult alula_airfoil_coupled_solve(AlulaAirfoilGrid *grid,const AlulaCoupledControls *controls);

/* Solve the starting condition, then advance on the same grid to target_mach.
 * A failed continuation keeps the last accepted condition and returns a failure
 * status. The caller owns the returned profiles as in alula_airfoil_coupled_solve. */
AlulaCoupledResult alula_airfoil_coupled_solve_to_mach(AlulaAirfoilGrid *grid,
    const AlulaCoupledControls *controls,double target_mach);

#endif
