// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_CONTOUR_H
#define ALULA_CONTOUR_H

#include "panel.h"

typedef struct {
    int count;
    double *arc;
    AlulaPoint *points, *second;
} AlulaContourSpline;

/* Natural cubic through an ordered upper-TE to lower-TE surface. */
int alula_contour_spline_create(AlulaContourSpline *curve, const AlulaPoint *points, int count);
void alula_contour_spline_free(AlulaContourSpline *curve);
void alula_contour_spline_value(const AlulaContourSpline *curve, double parameter,
    AlulaPoint *point, AlulaPoint *tangent, AlulaPoint *curvature);
/* Eight-point Gaussian integration over each cubic spline interval.
 * cumulative has curve->count entries and may be reused across queries. */
void alula_contour_arc_build(const AlulaContourSpline *curve,double *cumulative);
double alula_contour_arc(const AlulaContourSpline *curve,const double *cumulative,double parameter);

#endif
