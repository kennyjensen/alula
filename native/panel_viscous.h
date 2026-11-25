// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_PANEL_VISCOUS_H
#define ALULA_PANEL_VISCOUS_H
#include "panel_bl.h"

typedef enum {
    ALULA_VISCOUS_OK, ALULA_VISCOUS_INPUT, ALULA_VISCOUS_GEOMETRY,
    ALULA_VISCOUS_SEED, ALULA_VISCOUS_NEWTON, ALULA_VISCOUS_WAKE_REVERSE,
    ALULA_VISCOUS_WAKE_LIMIT, ALULA_VISCOUS_OUTSIDE_MODEL
} AlulaViscousStatus;

typedef struct {
    size_t wake_count;
    double wake_length,alpha,reynolds,ncrit;
    const double (*trips)[2];
    int newton_limit,wake_limit;
    double equation_tolerance,wake_tolerance,wake_relaxation;
} AlulaPanelViscousOptions;

/* Input coordinates use the reference chord. wake_length is in local
 * element chords. Moments are about (0.25, 0). Snapshots borrow their arrays
 * during the callback. wake_residual is unavailable until the BL converges. */
typedef struct {
    const AlulaPanelDisplacement *outer;
    const AlulaPanelBL *profile;
    int iteration,wake_iteration;
    double alpha,target_alpha,equation_residual,wake_residual,cl,cm,pressure_drag,wake_drag,max_displacement;
} AlulaPanelSnapshot;
typedef void (*AlulaPanelProgress)(void *,const AlulaPanelSnapshot *);
typedef struct {
    AlulaViscousStatus status;
    AlulaBLStatus bl_status;
    AlulaPanelSnapshot final;
} AlulaPanelViscousResult;

AlulaPanelViscousOptions alula_panel_viscous_defaults(void);
AlulaPanelViscousResult alula_panel_viscous_solve(const AlulaContour *contours,size_t count,
    const AlulaPanelViscousOptions *options,AlulaPanelProgress progress,void *context);
/* Release the final snapshot after any solve, including unsuccessful runs. */
void alula_panel_viscous_free(AlulaPanelViscousResult *result);
#endif
