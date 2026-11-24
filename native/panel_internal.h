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

typedef struct { double cx, cy, cm; } Loads;
typedef struct { double xx, xy, yy; } Pressure;

struct AlulaPanelSystem {
    size_t count;
    double residual_bound;
    double gx, gy;
    Loads loads[3];
    Pressure *pressure;
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

AlulaPanelSystem *alula_panel_allocate(size_t count);
void alula_panel_integrate(AlulaPanelSystem *system, size_t index, const Panel *panel,
    Pressure p, Pressure gradient, double chord, AlulaPoint moment);
int alula_panel_factor(double *a, double *b, size_t n);
void alula_panel_vortex(AlulaPoint point, const Panel *p, int self, AlulaPoint q[2]);
AlulaPanelSystem *alula_panel_prepare_streamfunction(const AlulaContour *contours,
    size_t count, size_t total, double chord, AlulaPoint moment);

#endif
