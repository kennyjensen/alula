// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_PANEL_BL_H
#define ALULA_PANEL_BL_H
#include "bl.h"
#include "panel_displacement.h"

/* Shared-field initialization for simultaneous panel/BL coupling. Stations
 * follow the displacement operator's material order. The march can use an
 * inverse interval at separation. Its edge mismatch is measured explicitly. */
typedef struct {
    size_t bodies,count;
    double stagnation[6],trip_arc[6][2],join_gap[6],fluid_gap[6],tangent_cross[6];
    size_t upper_first[6],lower_first[6];
    AlulaBLStation *stations;
    AlulaBLRegime *regimes;
    double *speed; /* Signed BL edge speed */
    AlulaBLContext physics;
    double interval_residual,outer_residual;
    int complete;
    size_t failed_body,failed_station;
    AlulaBLStatus status;
} AlulaPanelBL;

typedef struct {
    int iterations,backtracks;
    int converged; /* BL equations and common edge matching on this fixed wake. */
    double residual;
    AlulaBLStatus status;
} AlulaPanelBLNewton;

typedef void (*AlulaPanelBLProgress)(void *context,const AlulaPanelBL *profile,int iteration);

/* Improve a complete initial profile on its displacement operator. On
 * failure the last accepted iterate is retained. The wake path is fixed. */
AlulaPanelBLNewton alula_panel_bl_solve(const AlulaPanelDisplacement *outer,AlulaPanelBL *profile,
    int max_iterations,double tolerance,AlulaPanelBLProgress progress,void *context);

AlulaPanelBL *alula_panel_bl_initialize(const AlulaPanelDisplacement *outer,double reynolds,double ncrit,
    const double (*trips)[2]);
void alula_panel_bl_free(AlulaPanelBL *profile);
#endif
