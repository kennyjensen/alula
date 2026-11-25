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
    int converged;
    double *stations;
    double mach,lower_y,upper_y;
};

AlulaBodyStatus alula_airfoil_grid_admit_seed(AlulaAirfoilGrid *grid,double *residual);

/* Packed reference displacement in grid length units, as in body_chart.h. */
struct AlulaAirfoilBL;
struct AlulaAirfoilBL *alula_airfoil_bl_create_interaction(AlulaAirfoilGrid *grid,double reynolds,double ncrit,
    const double (*trips)[2],const double *displacement,int local);

#endif
