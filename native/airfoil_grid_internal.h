// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_AIRFOIL_GRID_INTERNAL_H
#define ALULA_AIRFOIL_GRID_INTERNAL_H

#include "airfoil_grid.h"
#include "body_chart.h"

typedef struct { AlulaPoint *points; size_t count; } AlulaMaterialBase;

struct AlulaAirfoilGrid {
    AlulaBody *body;
    AlulaBodyCoordinates *coordinates;
    AlulaBodyNewton *solver;
    AlulaBodyChart chart;
    AlulaBodyCut cuts[6];
    int tubes[7],contour_index[6],elements;
    double *packed,*physical;
    AlulaPoint *nodes;
    double base_width[6]; /* material gap, before BL displacement */
    AlulaMaterialBase bases[6];
    AlulaBodySolveResult last_solve;
    AlulaEulerProgress progress;
    void *progress_context;
    int converged;
    double *stations;
    double mach,lower_y,upper_y;
};

AlulaBodyStatus alula_airfoil_grid_admit_seed(AlulaAirfoilGrid *grid,double *residual);
/* Decode the current converged state and make its Euler stencils available
 * to postprocessing. */
int alula_airfoil_grid_prepare_flow(AlulaAirfoilGrid *grid,AlulaBodyState *state);
/* Integrate already evaluated cell pressures for an accepted iterate. */
int alula_airfoil_pressure_prepared(AlulaAirfoilGrid *grid,AlulaAirfoilPressureLoads *loads,
    AlulaAirfoilPressureVisitor visit,void *context);
/* The caller checkpoints the chart and packed state before attempting this
 * repair. A successful repair must also pass the Euler seed checks. */
int alula_airfoil_grid_repair_seed(AlulaAirfoilGrid *grid,double *residual);

typedef struct AlulaAirfoilGridSeed AlulaAirfoilGridSeed;
/* A seed borrows the contours until it is freed. Grids own their state and
 * can outlive the seed. Build grids sequentially with any mesh resolution.
 * keep_retry_mesh retains one potential mesh for multi-element retries. */
AlulaAirfoilGridSeed *alula_airfoil_grid_prepare(const AlulaContour *contours,size_t count,
    const AlulaAirfoilGridInput *input,int keep_retry_mesh,AlulaAirfoilGridStatus *status);
/* Interval counts, tube counts, surface spacing and crossline method may vary. */
AlulaAirfoilGrid *alula_airfoil_grid_create_seeded(AlulaAirfoilGridSeed *seed,
    const AlulaAirfoilGridInput *mesh,int wall_following,AlulaAirfoilGridStatus *status,
    AlulaAirfoilGridReport *report);
void alula_airfoil_grid_seed_free(AlulaAirfoilGridSeed *seed);
/* Try a bounded local repair of a stalled grid. Retain it only if a
 * fresh Newton solve meets the original residual and convexity checks. */
int alula_airfoil_grid_repair(AlulaAirfoilGrid *grid,int iterations,double tolerance,int prefer_tangent,
    AlulaBodySolveResult *result);
/* Advance an inviscid target on its own grid toward a coarse coupled
 * displacement. Return the verified fraction retained on the target. */
double alula_airfoil_grid_continue_displacement(AlulaAirfoilGrid *target,AlulaAirfoilGrid *source);

/* Packed reference displacement in grid length units, as in body_chart.h. */
struct AlulaAirfoilBL;
struct AlulaAirfoilBL *alula_airfoil_bl_create_interaction(AlulaAirfoilGrid *grid,double reynolds,double ncrit,
    const double (*trips)[2],const double *displacement,int local);

#endif
