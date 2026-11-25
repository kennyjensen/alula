// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_STREAMTUBE_INTERNAL_H
#define ALULA_STREAMTUBE_INTERNAL_H

#include "streamtube.h"

/* Initialize physical sections before applying the shared speed filter.
 * A previous stencil from the same residual evaluation may be supplied
 * when both cells use the same geometry policy and share the four points,
 * mass and gas inputs on their common interval. */
AlulaCellStatus alula_cell_initialize(AlulaCell *cell, const AlulaCellInput *input,
    const AlulaCell *previous);
AlulaCellStatus alula_cell_flux(AlulaCell *cell);
int alula_cell_geometry_only(AlulaCell *cell, const AlulaCellInput *input);

typedef struct {
    double coefficient, scale, slope, difference, ratio, delta, spacing;
    int second, kink;
} AlulaUpwindStencil;

AlulaCellStatus alula_upwind_prepare(const AlulaUpwindInput *input, AlulaUpwindValue *value,
    AlulaUpwindStencil *stencil);

/* The body Jacobian keeps the dissipation controls fixed across all colors. */
static inline double alula_upwind_difference_direction(const AlulaUpwindStencil *s,
    const AlulaUpwindInput *d)
{
    double dr=(d->spacing[1]-s->ratio*d->spacing[0])/s->spacing;
    return -(d->speed[2]-d->speed[1])+(s->second ?
        s->ratio*(d->speed[1]-d->speed[0])+s->delta*dr : 0);
}

static inline double alula_upwind_speed_direction(const AlulaUpwindStencil *s,
    const AlulaUpwindInput *d)
{
    if (s->coefficient==0 && s->slope==0) return d->speed[2];
    double dc=s->scale*(s->slope*(.5*d->mach_sq[0]+.5*d->mach_sq[1]));
    return d->speed[2]+(dc*s->difference+s->coefficient*alula_upwind_difference_direction(s,d));
}

#endif
