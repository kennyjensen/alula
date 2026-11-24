// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_PANEL_H
#define ALULA_PANEL_H

#include <stddef.h>

typedef struct { double x, y; } AlulaPoint;
#define ALULA_PANEL_MAX_WETTED 700
#define ALULA_PANEL_MAX_TOTAL 4096

typedef struct {
    const AlulaPoint *points;
    size_t count;
    int finite_base;
    size_t upper_index, lower_index;
} AlulaContour;
typedef struct AlulaPanelSystem AlulaPanelSystem;
typedef struct {
    double cl, cm, pressure_drag, circulation_lift, residual_bound;
} AlulaPanelResult;

/* Contours must be disjoint, counterclockwise, and closed. Sharp contours start and end at the TE.
 * For a finite base, set finite_base and give the two distinct TE vertex indices
 * in the original point array (excluding the repeated closing vertex).
 * An assembly containing a finite base uses the nodal streamfunction system.
 * The prepared system is immutable and does not borrow the contour arrays.
 * Angles are in degrees. cp, when supplied, has one entry per panel. */
AlulaPanelSystem *alula_panel_prepare(const AlulaContour *contours,
    size_t count, double chord, AlulaPoint moment);
int alula_panel_solve(const AlulaPanelSystem *system, double alpha,
    AlulaPanelResult *result, double *cp);
void alula_panel_free(AlulaPanelSystem *system);

#endif
