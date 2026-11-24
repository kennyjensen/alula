// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_STREAMTUBE_H
#define ALULA_STREAMTUBE_H

#include "panel.h"

typedef enum {
    ALULA_CELL_OK,
    ALULA_CELL_INPUT,
    ALULA_CELL_GEOMETRY,
    ALULA_CELL_ENTHALPY,
    ALULA_CELL_PRESSURE,
    ALULA_CELL_TRANSPORT,
    ALULA_CELL_NONFINITE
} AlulaCellStatus;

typedef struct {
    AlulaPoint lower[3], upper[3];
    double density[2], mass, h0, gamma, correction;
    int positive_simple, biased;
    double transport[2];
} AlulaCellInput;

typedef struct {
    AlulaPoint section[2], direction[2], side[2], streamwise, transverse;
    double area, normal_area[2], length[2], curvature;
} AlulaCellGeometry;

typedef struct {
    double density, speed, pressure, enthalpy, mach_sq;
} AlulaSection;

typedef struct {
    AlulaCellGeometry geometry;
    AlulaSection state[2];
    double transport[2], artificial_h[2];
    double correction, pressure[2], momentum, entropy, artificial_entropy, isentropic;
} AlulaCellValue;

/* Reusable state for applying multiple columns of the local Jacobian. */
typedef struct {
    AlulaCellInput input;
    AlulaCellValue value;
    AlulaPoint gap[3], segment[2];
} AlulaCell;

typedef struct {
    AlulaPoint lower[3], upper[3];
    double density[2], mass, h0, correction, transport[2];
} AlulaCellDirection;

/* On success the output owns all data. Outputs on failure must be discarded. */
AlulaCellStatus alula_cell_prepare(AlulaCell *cell, const AlulaCellInput *input);
AlulaCellStatus alula_cell_apply(const AlulaCell *cell, const AlulaCellDirection *direction,
    AlulaCellValue *derivative);

typedef struct {
    double speed[3], mach_sq[2], spacing[2], mucon, mcrit, gamma;
} AlulaUpwindInput;

typedef struct {
    double speed, coefficient, correction;
} AlulaUpwindValue;

/* An optional direction differentiates all controls and stencil inputs. */
AlulaCellStatus alula_streamtube_upwind(const AlulaUpwindInput *input,
    AlulaUpwindValue *value, const AlulaUpwindInput *direction, AlulaUpwindValue *derivative);

typedef struct {
    double fraction, residual, loss[2], compression;
} AlulaBlendValue;

AlulaCellStatus alula_streamtube_blend(const AlulaCellValue *cell, double epsilon,
    AlulaBlendValue *value, const AlulaCellValue *direction, double epsilon_direction,
    AlulaBlendValue *derivative);

#endif
