// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_PANEL_DISPLACEMENT_H
#define ALULA_PANEL_DISPLACEMENT_H

#include "panel.h"

typedef struct AlulaPanelDisplacement AlulaPanelDisplacement;
typedef struct {
    size_t bodies,stations,surface_stations,wake_count,surface_start[7];
    /* Upper TE to lower TE for each body, followed by the wakes in body order.
     * Arc length starts at zero independently for each surface and wake. */
    const AlulaPoint *points;
    const double *arc,*inviscid,*influence; /* row-major d(signed edge speed)/d(mass) */
} AlulaDisplacementData;

typedef enum { ALULA_WAKE_STRAIGHT, ALULA_WAKE_INVISCID } AlulaWakeStart;

/* Contours must be disjoint, closed and counterclockwise, with trailing
 * edges specified as in panel.h. Unit freestream, coordinates in reference
 * chords and alpha in degrees. wake_length is in local chords. Optional prescribed paths
 * contain bodies*wake_count points, including each physical wake origin.
 * With explicit paths, initialization and wake_length are unused. */
AlulaPanelDisplacement *alula_panel_displacement_create(const AlulaContour *contours,
    size_t bodies,double alpha,size_t wake_count,double wake_length,AlulaWakeStart initialization,
    const AlulaPoint *paths);
const AlulaDisplacementData *alula_panel_displacement_data(const AlulaPanelDisplacement *system);
/* Reuse the identical mass influence and fixed wake geometry at a new
 * incidence. Updates the borrowed inviscid array. Call between solves. */
int alula_panel_displacement_set_incidence(AlulaPanelDisplacement *system,double alpha);
/* Signed Ue*delta* on each surface, positive downstream on each wake.
 * Input and output buffers must be disjoint. */
int alula_panel_displacement_evaluate(const AlulaPanelDisplacement *system,const double *mass,
    size_t stations,double *speed);
/* Mean velocity at each wake segment midpoint (body-major). */
int alula_panel_displacement_wake_velocity(const AlulaPanelDisplacement *system,
    const double *mass,size_t stations,AlulaPoint *velocity);
void alula_panel_displacement_free(AlulaPanelDisplacement *system);

#endif
