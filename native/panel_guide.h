// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_PANEL_GUIDE_H
#define ALULA_PANEL_GUIDE_H

#include "contour.h"
#include "trace.h"

typedef struct AlulaPanelGuides AlulaPanelGuides;
typedef struct AlulaPanelStreamline AlulaPanelStreamline;
typedef struct {
    double inlet_x, outlet_x, inlet_lower_y, inlet_upper_y;
    double length_scale, relative_tolerance, seed_fraction;
    /* Optional per-element upper bound on the physical seed offset. */
    const double *seed_limit;
} AlulaPanelGuideControls;

typedef struct {
    const AlulaTracePoint *points;
    size_t count;
} AlulaGuidePath;

typedef struct {
    size_t element;
    const AlulaContourSpline *curve;
    double stagnation, strain, wall_level, inlet_level, seed_distance;
    AlulaPoint stagnation_point, trailing_edge, trailing_vector;
    /* Both paths start at the off-body seed and end on their named plane.
     * The potential origin of each path is local to its seed. */
    AlulaGuidePath incoming, wake;
} AlulaPanelBodyGuide;

typedef enum {
    ALULA_GUIDE_OK, ALULA_GUIDE_INVALID, ALULA_GUIDE_ALLOCATION,
    ALULA_GUIDE_STAGNATION, ALULA_GUIDE_SEED, ALULA_GUIDE_TRACE, ALULA_GUIDE_ORDER,
    ALULA_GUIDE_UNRESOLVED, ALULA_GUIDE_AMBIGUOUS
} AlulaPanelGuideStatus;

/* The field must outlive the guides. Body guides are ordered by captured
 * inlet streamfunction and retain the original element indices. The
 * capture levels have body_count+2 entries and include both farfields. */
AlulaPanelGuides *alula_panel_guides_create(const AlulaPanelField *field,
    const AlulaPanelGuideControls *controls, AlulaPanelGuideStatus *status);
const AlulaPanelBodyGuide *alula_panel_guides_bodies(const AlulaPanelGuides *guides, size_t *count);
const double *alula_panel_guides_levels(const AlulaPanelGuides *guides);
AlulaGuidePath alula_panel_guides_farfield(const AlulaPanelGuides *guides, int upper);
/* Find one resolved x crossing and reintegrate from its bracket's first
 * point. Reject multiple crossings and intervals with unresolved turns. */
AlulaPanelGuideStatus alula_panel_guides_sample_x(const AlulaPanelGuides *guides,
    AlulaGuidePath path, double x, AlulaTracePoint *point);
/* Trace an interior streamfunction level within a passage from inlet to
 * outlet. The returned line and its points belong to the caller. */
AlulaPanelStreamline *alula_panel_guides_trace_level(const AlulaPanelGuides *guides,
    double level, AlulaPanelGuideStatus *status);
AlulaGuidePath alula_panel_streamline_path(const AlulaPanelStreamline *line);
void alula_panel_streamline_free(AlulaPanelStreamline *line);
void alula_panel_guides_free(AlulaPanelGuides *guides);

#endif
