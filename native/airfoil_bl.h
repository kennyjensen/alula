// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_AIRFOIL_BL_H
#define ALULA_AIRFOIL_BL_H

#include "airfoil_grid.h"
#include "bl.h"

typedef struct AlulaAirfoilBL AlulaAirfoilBL;
typedef enum { ALULA_BL_UPPER, ALULA_BL_LOWER, ALULA_BL_MERGED_WAKE } AlulaBLBranch;

typedef struct {
    const AlulaBLStation *stations;
    const AlulaBLRegime *regimes;
    size_t count; /* Only result.completed entries are initialized. */
    AlulaBLMarchResult result;
} AlulaAirfoilBLBranch;

/* Initialize prescribed-edge profiles from a converged native Euler grid.
 * Trips are upper/lower material fractions in grid passage order. A fraction
 * of one or a NULL trips pointer selects automatic transition with terminal
 * shear initialization. The returned object may contain partial profiles. */
AlulaAirfoilBL *alula_airfoil_bl_create(AlulaAirfoilGrid *grid,double reynolds,double ncrit,
    const double (*trips)[2]);
const AlulaAirfoilBLBranch *alula_airfoil_bl_branch(const AlulaAirfoilBL *profiles,
    int body,AlulaBLBranch branch);
int alula_airfoil_bl_complete(const AlulaAirfoilBL *profiles);
/* Compressible Squire–Young wake drag estimate, summed over the wakes.
 * Uses the BL edge gas convention and momentum/displacement thickness. */
int alula_airfoil_bl_wake_drag(const AlulaAirfoilBL *profiles,double *coefficient);
void alula_airfoil_bl_free(AlulaAirfoilBL *profiles);

#endif
