// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_AIRFOIL_BL_INTERNAL_H
#define ALULA_AIRFOIL_BL_INTERNAL_H

#include "airfoil_coupled.h"
#include "body_internal.h"

struct AlulaAirfoilBL {
    AlulaAirfoilBLBranch branches[MAX_BODIES][3];
    AlulaBLStation *stations;
    AlulaBLRegime *regimes;
    AlulaBLEdge *edges;
    double trip[MAX_BODIES][2];
    AlulaBLContext conditions;
    int bodies,complete;
};

/* A snapshot owns its edge samples and can be reevaluated at different BL states. */
AlulaAirfoilBL *alula_airfoil_bl_sample(AlulaAirfoilGrid *g,double reynolds,double ncrit,
    const double (*trips)[2],int probe);
/* A probe samples wall cells into scratch storage, preserving the base
 * Euler stencils for the next linear response. */
int alula_airfoil_bl_resample(AlulaAirfoilBL *p,AlulaAirfoilGrid *g,const double (*trips)[2],int probe);

/* Attempt a coupled correction from a complete interaction march. On failure
 * keep the original accepted grid and profiles paired. */
void alula_airfoil_coupled_newton(AlulaAirfoilGrid *g,const AlulaCoupledControls *controls,
    AlulaCoupledResult *result,int coupled_flow);

#endif
