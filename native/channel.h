// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_CHANNEL_H
#define ALULA_CHANNEL_H

#include "sparse.h"

typedef enum {
    ALULA_CHANNEL_MOMENTUM,
    ALULA_CHANNEL_ISENTROPIC,
    ALULA_CHANNEL_UPWIND,
    ALULA_CHANNEL_HYBRID
} AlulaChannelMode;

typedef struct {
    int nx, nt;
    const double *x, *lower, *upper; /* nx+1 entries */
    const double *mass, *h0, *rho_total; /* nt entries */
    const double *inlet_slope, *outlet_slope; /* nt-1 entries, NULL means zero */
    const double *outlet_pressure; /* nt entries for reservoir modes */
    double gamma, rho_ref, p_ref, correction, mucon, mcrit, epsilon;
    AlulaChannelMode mode;
} AlulaChannelInput;

typedef struct AlulaChannel AlulaChannel;
typedef enum {
    ALULA_CHANNEL_CONVERGED, ALULA_CHANNEL_LIMIT, ALULA_CHANNEL_INVALID,
    ALULA_CHANNEL_LINEAR_FAILURE, ALULA_CHANNEL_LINE_SEARCH_FAILURE
} AlulaChannelStatus;

typedef struct {
    AlulaChannelStatus status;
    int iterations, backtracks;
    double residual, max_mach;
    AlulaSparseStats linear;
} AlulaChannelResult;

/* State ordering matches the JS verification systems: log densities,
 * interior node offsets, and (in reservoir modes) log tube masses.
 * Prescribed inputs are copied. The channel owns mutable workspace, so
 * concurrent solves must use separate channels. Linear stats are cumulative
 * for a channel. Only accepted iterations are copied into the output state. */
AlulaChannel *alula_channel_create(const AlulaChannelInput *input);
int alula_channel_size(const AlulaChannel *channel);
int alula_channel_residual(AlulaChannel *channel, const double *state, double *residual);
int alula_channel_jacobian_product(AlulaChannel *channel, const double *state,
    const double *direction, double *product);
AlulaChannelResult alula_channel_solve(AlulaChannel *channel, double *state,
    int max_iterations, double tolerance);
void alula_channel_free(AlulaChannel *channel);

#endif
