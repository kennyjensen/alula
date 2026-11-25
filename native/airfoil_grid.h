// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_AIRFOIL_GRID_H
#define ALULA_AIRFOIL_GRID_H

#include "body_newton.h"

typedef struct AlulaAirfoilGrid AlulaAirfoilGrid;
typedef enum {
    ALULA_CROSSLINE_AUTO, ALULA_CROSSLINE_ALGEBRAIC, ALULA_CROSSLINE_POTENTIAL
} AlulaCrosslinePlacement;
typedef struct {
    double mach, alpha, chord;
    AlulaPoint moment;
    double inlet_x, outlet_x, lower_y, upper_y;
    int surface_intervals, inlet_intervals, outlet_intervals;
    const int *tubes; /* count+1 passages, from lower farfield to upper */
    AlulaCrosslinePlacement crosslines;
    AlulaBodyMode mode;
    int ismom; /* 1 through 4 for hybrid mode */
    double surface_uniformity; /* 0: cosine station spacing, 1: uniform */
} AlulaAirfoilGridInput;

typedef enum {
    ALULA_GRID_OK, ALULA_GRID_INVALID, ALULA_GRID_PANEL, ALULA_GRID_GUIDE,
    ALULA_GRID_ALLOCATION, ALULA_GRID_STATIONS, ALULA_GRID_FOLDED
} AlulaAirfoilGridStatus;

typedef struct {
    int nx, elements, invalid_cells;
    int first_group, first_interval, first_tube;
    double minimum_corner_sine;
    AlulaCrosslinePlacement crosslines;
    AlulaBodyStatus flow_seed; /* actual Euler residual check on accepted grids */
} AlulaAirfoilGridReport;

/* Build a panel-guided grid and a native moving-coordinate Euler system
 * directly from closed contours. Returned element indices are ordered by
 * inlet streamfunction. Finite bases use separated material wake banks. */
AlulaAirfoilGrid *alula_airfoil_grid_create(const AlulaContour *contours, size_t count,
    const AlulaAirfoilGridInput *input, AlulaAirfoilGridStatus *status, AlulaAirfoilGridReport *report);
AlulaBodySize alula_airfoil_grid_size(const AlulaAirfoilGrid *grid);
double alula_airfoil_grid_mach(const AlulaAirfoilGrid *grid);
const AlulaBodyCut *alula_airfoil_grid_cuts(const AlulaAirfoilGrid *grid);
const int *alula_airfoil_grid_tubes(const AlulaAirfoilGrid *grid);
/* Input contour index for a body in grid passage order, or -1. */
int alula_airfoil_grid_contour_index(const AlulaAirfoilGrid *grid,int body);
/* View the current physical state. Arrays are valid until the next solve,
 * state query, or free. */
int alula_airfoil_grid_state(AlulaAirfoilGrid *grid, AlulaBodyState *state);
/* Measure the current grid, including its invalid cell locations. This
 * populates geometry fields only. */
int alula_airfoil_grid_current_quality(AlulaAirfoilGrid *grid,AlulaAirfoilGridReport *report);
typedef struct { double cl,cm,pressure_drag,cx,cy; } AlulaAirfoilPressureLoads;
typedef enum {
    ALULA_PRESSURE_UPPER, ALULA_PRESSURE_LOWER, ALULA_PRESSURE_BASE
} AlulaAirfoilPressurePart;
typedef void (*AlulaAirfoilPressureVisitor)(void *context,int contour,
    AlulaAirfoilPressurePart part,AlulaPoint point,double cp);
/* Integrate the current converged Euler interface pressures on the
 * material surface. Retained finite-base segments interpolate TE pressure. */
int alula_airfoil_grid_pressure_loads(AlulaAirfoilGrid *grid,AlulaAirfoilPressureLoads *loads);
/* Optionally visit the same samples used for the load integral. Surface
 * samples run from stagnation to trailing edge, base samples from lower to
 * upper edge. Contour indices refer to the caller's input order. The visitor
 * must not reenter the grid API on the same grid. */
int alula_airfoil_grid_pressure_results(AlulaAirfoilGrid *grid,AlulaAirfoilPressureLoads *loads,
    AlulaAirfoilPressureVisitor visit,void *context);
AlulaBodySolveResult alula_airfoil_grid_solve(AlulaAirfoilGrid *grid, int iterations, double tolerance);
/* Seed a finer grid from a converged grid with the same geometry, operating
 * point and domain. Successful transfers have a convex grid and an
 * admissible gas state. Failure preserves the target's original packed
 * state and chart. */
int alula_airfoil_grid_transfer(AlulaAirfoilGrid *target,AlulaAirfoilGrid *source,double geometry_weight);
void alula_airfoil_grid_free(AlulaAirfoilGrid *grid);

#endif
