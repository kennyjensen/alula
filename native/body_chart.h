// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_BODY_COORDINATES_H
#define ALULA_BODY_COORDINATES_H

#include "body_newton.h"

typedef struct AlulaBodyCoordinates AlulaBodyCoordinates;

typedef struct {
    const AlulaPoint *surface; /* upper TE -> leading edge -> lower TE */
    int count;
    double stagnation;
    const double *lower_fraction,*upper_fraction; /* LE through TE, inclusive */
    const double *lower_thickness,*upper_thickness; /* required if displaced */
    const double *wake_width; /* nx - TE entries if displaced */
    int material_wake;
} AlulaBodySurface;

typedef struct {
    const AlulaPoint *raw_nodes; /* prepared mesh before displacement */
    const double *mass_weights; /* one per tube, grouped lower to upper */
    const double *capture_levels; /* bodies+2 increasing streamfunction levels */
    const AlulaBodySurface *surfaces;
    int primary_body, walls_only, body_station_normals, te_center_wake;
} AlulaBodyCoordinateInput;

/* The body must outlive coordinates and any solver created from the chart.
 * Inputs are copied. Initial packed zero means the supplied raw grid, contour
 * stagnation, capture levels, unit density, and zero multipole strengths. */
AlulaBodyCoordinates *alula_body_coordinates_create(AlulaBody *body,
    const AlulaBodyCoordinateInput *input);
AlulaBodyChart alula_body_coordinates_chart(AlulaBodyCoordinates *coordinates);
int alula_body_coordinates_stagnation_column(const AlulaBodyCoordinates *coordinates,int body);
void alula_body_coordinates_free(AlulaBodyCoordinates *coordinates);

#endif
