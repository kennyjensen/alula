// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_PANEL_INTERNAL_H
#define ALULA_PANEL_INTERNAL_H

#include "panel.h"
#include <math.h>

#define ALULA_PI 3.1415926535897932384626433832795

typedef struct {
    AlulaPoint a, center, tangent, normal;
    double length;
    size_t node;
} Panel;

typedef struct {
    Panel geometry;
    AlulaPoint end, cut_origin, cut_direction;
    AlulaPoint tangential;
    size_t output, body, last_node;
    double source, vortex, log_length;
} Segment;

typedef struct { double cx, cy, cm; } Loads;
typedef struct { double xx, xy, yy; } Pressure;

struct AlulaPanelSystem {
    size_t count;
    double residual_bound;
    double gx, gy;
    Loads loads[3];
    Pressure *pressure;
    Segment *segments;
    double *responses;
    AlulaPoint boundary_psi[6];
    size_t wetted, response_size, bodies;
    int dirichlet;
};

struct AlulaPanelField {
    const AlulaPanelSystem *system;
    AlulaPoint freestream;
    struct PanelGroup *groups;
    AlulaPoint strengths[];
};

static inline double dot(AlulaPoint a, AlulaPoint b)
{
    return a.x*b.x + a.y*b.y;
}

static inline AlulaPoint subtract(AlulaPoint a, AlulaPoint b)
{
    return (AlulaPoint){a.x-b.x, a.y-b.y};
}

static inline double cross(AlulaPoint a, AlulaPoint b)
{
    return a.x*b.y - a.y*b.x;
}

static inline Pressure pressure_between(double ax, double ay, double bx, double by)
{
    return (Pressure){1-.5*(ax*ax+bx*bx), -(ax*ay+bx*by), 1-.5*(ay*ay+by*by)};
}

typedef struct {
    size_t start, end, offset, wetted;
    AlulaPoint origin, direction, wake_direction;
    int finite;
} PanelRange;

typedef struct {
    AlulaPoint point, tangent;
    size_t body;
    int interior;
} PanelRow;

int alula_panel_dimensions(const AlulaContour *contours,size_t count,size_t *total,size_t *wet);
int alula_panel_boundary_setup(const AlulaContour *contours,size_t count,size_t total,
    Segment *segments,PanelRange *ranges,PanelRow *rows,size_t wetted);
void alula_panel_boundary_assemble(double *a,double *bx,double *by,size_t count,
    const Segment *segments,size_t wet,size_t total,const PanelRow *rows,const PanelRange *ranges,double chord);

AlulaPanelSystem *alula_panel_allocate(size_t count);
void alula_panel_integrate(AlulaPanelSystem *system, size_t index, const Panel *panel,
    Pressure p, Pressure gradient, double chord, AlulaPoint moment);
void alula_panel_vortex(AlulaPoint point, const Panel *p, int self, AlulaPoint q[2]);
void alula_panel_vortex_psi(AlulaPoint point, const Segment *s, double q[2]);
double alula_panel_base_psi(AlulaPoint point, const Segment *s, double source, double vortex);
AlulaPanelSystem *alula_panel_prepare_streamfunction(const AlulaContour *contours,
    size_t count, size_t total, size_t wet, double chord, AlulaPoint moment);

#endif
