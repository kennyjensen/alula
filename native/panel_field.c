// SPDX-License-Identifier: GPL-2.0-or-later
#include "panel_internal.h"

#include <stdlib.h>

AlulaPanelField *alula_panel_field_create(const AlulaPanelSystem *system, double alpha)
{
    if (!system || !isfinite(alpha) || fabs(alpha) > 90) return NULL;
    AlulaPanelField *field = malloc(sizeof *field + system->count*sizeof *field->strengths);
    if (!field) return NULL;
    double angle = alpha*ALULA_PI/180;
    field->system = system;
    field->freestream = (AlulaPoint){cos(angle), sin(angle)};
    const double *x = system->responses, *y = x+system->response_size;
    for (size_t i = 0; i < system->count; ++i) {
        const Segment *s = system->segments+i;
        size_t a = s->geometry.node, b = i < system->wetted ? a+1 : s->last_node;
        double ga = field->freestream.x*x[a]+field->freestream.y*y[a];
        double gb = field->freestream.x*x[b]+field->freestream.y*y[b];
        field->strengths[i] = i < system->wetted ? (AlulaPoint){ga, gb} :
            (AlulaPoint){s->source*(ga-gb), s->vortex*(ga-gb)};
    }
    return field;
}

void alula_panel_field_free(AlulaPanelField *field)
{
    free(field);
}

static void add_psi(double *psi, double *correction, double value)
{
    double term = value-*correction;
    double next = *psi+term;
    *correction = (next-*psi)-term;
    *psi = next;
}

static int endpoint(AlulaPoint point, const Segment *s)
{
    return (point.x == s->geometry.a.x && point.y == s->geometry.a.y) ||
        (point.x == s->end.x && point.y == s->end.y);
}

int alula_panel_field_sample(const AlulaPanelField *field, AlulaPoint point,
    AlulaPoint *velocity, double *streamfunction)
{
    if (!field || (!velocity && !streamfunction) || !isfinite(point.x) || !isfinite(point.y)) return 0;
    const AlulaPanelSystem *system = field->system;
    AlulaPoint q = field->freestream;
    double psi = cross(q, point), correction = 0;
    for (size_t i = 0; i < system->wetted; ++i) {
        const Segment *s = system->segments+i;
        AlulaPoint g = field->strengths[i];
        if (velocity) {
            if (endpoint(point, s)) return 0;
            AlulaPoint v[2];
            alula_panel_vortex(point, &s->geometry, 0, v);
            q.x += g.x*v[0].x+g.y*v[1].x;
            q.y += g.x*v[0].y+g.y*v[1].y;
        }
        if (streamfunction) {
            double basis[2];
            alula_panel_vortex_psi(point, s, basis);
            add_psi(&psi, &correction, g.x*basis[0]);
            add_psi(&psi, &correction, g.y*basis[1]);
        }
    }
    for (size_t i = system->wetted; i < system->count; ++i) {
        const Segment *s = system->segments+i;
        AlulaPoint g = field->strengths[i];
        if (velocity) {
            if (endpoint(point, s)) return 0;
            AlulaPoint v[2];
            alula_panel_vortex(point, &s->geometry, 0, v);
            double u = v[0].x+v[1].x, w = v[0].y+v[1].y;
            q.x += g.x*w+g.y*u;
            q.y += -g.x*u+g.y*w;
        }
        if (streamfunction) add_psi(&psi, &correction, alula_panel_base_psi(point, s, g.x, g.y));
    }
    if ((velocity && (!isfinite(q.x) || !isfinite(q.y))) ||
        (streamfunction && !isfinite(psi))) return 0;
    if (velocity) *velocity = q;
    if (streamfunction) *streamfunction = psi;
    return 1;
}

double alula_panel_field_branch_increment(const AlulaPanelField *field,
    AlulaPoint a, AlulaPoint b)
{
    if (!field || !isfinite(a.x) || !isfinite(a.y) || !isfinite(b.x) || !isfinite(b.y)) return NAN;
    const AlulaPanelSystem *system = field->system;
    double increment = 0;
    for (size_t i = system->wetted; i < system->count; ++i) {
        const Segment *s = system->segments+i;
        double source = field->strengths[i].x;
        if (source == 0) continue;
        AlulaPoint ap = subtract(a, s->cut_origin), bp = subtract(b, s->cut_origin);
        double sa = cross(s->cut_direction, ap), sb = cross(s->cut_direction, bp);
        if ((sa < 0) == (sb < 0)) continue;
        double t = sa/(sa-sb);
        AlulaPoint intersection = {ap.x+t*(bp.x-ap.x), ap.y+t*(bp.y-ap.y)};
        if (dot(s->cut_direction, intersection) > 0)
            increment += (sa < 0 ? -source : source)*s->geometry.length;
    }
    return increment;
}
