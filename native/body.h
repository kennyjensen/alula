// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_BODY_H
#define ALULA_BODY_H

#include "streamtube.h"

typedef enum {
    ALULA_BODY_MOMENTUM, ALULA_BODY_ISENTROPIC, ALULA_BODY_HYBRID
} AlulaBodyMode;

typedef enum {
    ALULA_BODY_STREAMWISE, ALULA_BODY_INLET_DENSITY, ALULA_BODY_INTERNAL_PRESSURE,
    ALULA_BODY_FARFIELD_PRESSURE, ALULA_BODY_CUT_PRESSURE, ALULA_BODY_END_TANGENCY,
    ALULA_BODY_WAKE_GAP, ALULA_BODY_TRAILING_KUTTA, ALULA_BODY_LEADING_KUTTA,
    ALULA_BODY_FARFIELD_MATCH, ALULA_BODY_FAMILIES
} AlulaBodyFamily;

typedef struct {
    int leading, trailing, rounded;
    /* ISMOM3 inherited region; use through == 0 for the default region. */
    int through, lower_tube, upper_tube;
} AlulaBodyCut;

typedef struct {
    int nx, bodies;
    const int *tubes;          /* bodies+1 passages, lower to upper */
    const AlulaBodyCut *cuts;   /* bodies entries */
    double gamma, mach, alpha, length_scale, correction;
    AlulaPoint center;
    AlulaBodyMode mode;
    int upwind, ismom, positive_simple, displaced, independent_wakes, bank_outlet;
    double mucon, mcrit, epsilon;
} AlulaBodyInput;

typedef struct {
    /* Group major, then station, then bank. Each passage has nx+1 rows of
     * tubes[group]+1 banks. Density is section major, then global tube. */
    const AlulaPoint *nodes;
    const double *density, *mass; /* nx*nt and nt physical values */
    /* Circulation uses positive lift sign. Other strengths are dimensional
     * in the same normalized coordinates as the nodes. */
    double circulation, source, doublet_x, doublet_y;
    /* Body major, then downstream station TE+1 through nx, when enabled. */
    const double *wake_width;
} AlulaBodyState;

typedef struct {
    int rows, nodes, density, tubes, wake_widths;
} AlulaBodySize;

typedef struct {
    double residual[ALULA_BODY_FAMILIES], max_mach;
} AlulaBodyDiagnostics;

typedef enum {
    ALULA_BODY_OK, ALULA_BODY_INVALID, ALULA_BODY_CELL, ALULA_BODY_BOUNDARY,
    ALULA_BODY_GAP, ALULA_BODY_NONFINITE
} AlulaBodyStatus;

typedef struct AlulaBody AlulaBody;

/* Prescribed configuration is copied. Evaluations reuse workspace, so use a
 * separate object for each concurrent evaluation. */
AlulaBody *alula_body_create(const AlulaBodyInput *input);
AlulaBodySize alula_body_size(const AlulaBody *body);
AlulaBodyStatus alula_body_residual(AlulaBody *body, const AlulaBodyState *state,
    double *residual, AlulaBodyDiagnostics *diagnostics);
/* Apply the analytic derivative in physical coordinates. The direction
 * uses the same layout as the state; NULL direction arrays are zero.
 * Residual and product are separate buffers of size rows. */
AlulaBodyStatus alula_body_linearize(AlulaBody *body, const AlulaBodyState *state,
    const AlulaBodyState *direction, double *residual, double *product,
    AlulaBodyDiagnostics *diagnostics);
void alula_body_free(AlulaBody *body);

#endif
