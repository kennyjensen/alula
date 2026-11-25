// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_TRACE_H
#define ALULA_TRACE_H

#include "panel.h"

typedef enum { ALULA_TRACE_POTENTIAL, ALULA_TRACE_ARC } AlulaTraceParameter;
typedef enum {
    ALULA_TRACE_END, ALULA_TRACE_PLANE, ALULA_TRACE_INVALID,
    ALULA_TRACE_DOMAIN, ALULA_TRACE_MIN_STEP, ALULA_TRACE_STEP_LIMIT,
    ALULA_TRACE_PLANE_NOT_REACHED, ALULA_TRACE_PLANE_UNRESOLVED
} AlulaTraceStatus;

typedef struct {
    AlulaPoint position;
    double parameter, potential;
} AlulaTracePoint;

typedef struct {
    AlulaTraceParameter parameter;
    double end, tolerance, potential_tolerance, max_step, min_step, max_spatial_step;
    int max_steps, stop_at_plane;
    AlulaPoint plane_normal;
    double plane_offset;
    void *context;
    int (*admissible_point)(void *context, AlulaPoint point);
    int (*admissible_segment)(void *context, AlulaPoint start, AlulaPoint end);
    void (*visit)(void *context, const AlulaTracePoint *point);
} AlulaTraceControls;

typedef struct {
    AlulaTraceStatus status;
    int accepted, rejected;
    AlulaTracePoint last;
} AlulaTraceResult;

/* Trace forward or backward in velocity potential or physical arc. With
 * potential parameterization, seed.parameter and seed.potential must agree.
 * max_spatial_step is required for potential traces. potential_tolerance
 * controls the integrated potential error in arc traces. Optional callbacks
 * guard the stages and accepted chords or record the seed and accepted path. */
AlulaTraceResult alula_trace_panel(const AlulaPanelField *field,
    AlulaTracePoint seed, const AlulaTraceControls *controls);

#endif
