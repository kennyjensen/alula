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
    int tubes[7];
    double *packed,*physical;
    AlulaPoint *nodes;
    double base_width[6]; /* material gap, before BL displacement */
    AlulaMaterialBase bases[6];
    AlulaBodySolveResult last_solve;
    int converged;
};

#endif
