// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_BODY_COORDINATES_H
#define ALULA_BODY_COORDINATES_H

#include "body_newton.h"
#include "contour.h"

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
void alula_body_coordinates_seed_circulation(const AlulaBodyCoordinates *coordinates,
    double *packed,double circulation);
/* The caller checkpoints packed variables and origins before installing a
 * candidate seed. The material wall locations keep their contour constraint. */
int alula_body_coordinates_seed_state(AlulaBodyCoordinates *coordinates,double *packed,
    const double *stagnation,const double *passage_mass,const AlulaBodyState *source);
/* Decode only the cells touching either passage bank for BL edge probes.
 * nodes and gas need size.nodes points and size.density + size.tubes doubles.
 * The returned state is valid only at those banks. */
int alula_body_coordinates_edge_state(AlulaBodyCoordinates *coordinates,const double *packed,
    AlulaPoint *nodes,double *gas,AlulaBodyState *state);
int alula_body_coordinates_seed_nodes(AlulaBodyCoordinates *coordinates,double *packed,const AlulaPoint *nodes);
typedef struct {
    const AlulaContourSpline *curve;
    const double *fractions;
    double stagnation,initial_stagnation;
    const double *arc_lengths; /* Cumulative spline length for displacement-enabled bodies. */
} AlulaBodySurfaceView;
int alula_body_coordinates_surface(const AlulaBodyCoordinates *coordinates,const double *packed,
    int element,int upper,AlulaBodySurfaceView *surface);
/* Each body contributes lower then upper LE-to-TE thickness (inclusive),
 * followed by wake station widths from TE+1. Units are grid lengths.
 * Updating these values preserves the sparse chart graph. */
size_t alula_body_coordinates_displacement_count(const AlulaBodyCoordinates *coordinates);
int alula_body_coordinates_get_displacement(const AlulaBodyCoordinates *coordinates,double *values);
int alula_body_coordinates_set_displacement(AlulaBodyCoordinates *coordinates,const double *values);
/* Save and restore the chart origin and normals independently of the packed
 * state and prescribed displacements. Storage is checkpoint_count points. */
size_t alula_body_coordinates_checkpoint_count(const AlulaBodyCoordinates *coordinates);
int alula_body_coordinates_save(const AlulaBodyCoordinates *coordinates,AlulaPoint *checkpoint);
int alula_body_coordinates_restore(AlulaBodyCoordinates *coordinates,const AlulaPoint *checkpoint);
void alula_body_coordinates_free(AlulaBodyCoordinates *coordinates);

#endif
