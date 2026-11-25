// SPDX-License-Identifier: GPL-2.0-or-later
#include "panel_internal.h"
#include "dense.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void alula_panel_vortex(AlulaPoint point, const Panel *p, int self, AlulaPoint q[2])
{
    AlulaPoint d = {point.x - p->a.x, point.y - p->a.y};
    double x = dot(d, p->tangent);
    double y = -dot(d, p->normal);
    double sx = 0, sy = -.5;
    if (!self) {
        double r1 = x*x + y*y;
        double r2 = (x-p->length)*(x-p->length) + y*y;
        sx = log(r1/r2)/(4*ALULA_PI);
        sy = atan2(y*p->length, x*(x-p->length) + y*y)/(2*ALULA_PI);
    }
    double end_x = (x*sx + y*sy)/p->length - 1/(2*ALULA_PI);
    double end_y = (x*sy - y*sx)/p->length;
    q[0] = (AlulaPoint){
        -(sy-end_y)*p->tangent.x - (sx-end_x)*p->tangent.y,
        -(sy-end_y)*p->tangent.y + (sx-end_x)*p->tangent.x
    };
    q[1] = (AlulaPoint){
        -end_y*p->tangent.x - end_x*p->tangent.y,
        -end_y*p->tangent.y + end_x*p->tangent.x
    };
}

AlulaPanelSystem *alula_panel_allocate(size_t count)
{
    AlulaPanelSystem *system = calloc(1, sizeof *system);
    if (!system) return NULL;
    system->count = count;
    system->pressure = calloc(count, sizeof *system->pressure);
    if (!system->pressure) {
        free(system);
        return NULL;
    }
    return system;
}

void alula_panel_integrate(AlulaPanelSystem *system, size_t index, const Panel *panel,
    Pressure p, Pressure gradient, double chord, AlulaPoint moment)
{
    system->pressure[index] = p;
    double coefficients[] = {p.xx, p.xy, p.yy};
    double slopes[] = {gradient.xx, gradient.xy, gradient.yy};
    double fx = -panel->normal.x*panel->length/chord;
    double fy = -panel->normal.y*panel->length/chord;
    double cm = -cross(subtract(panel->center, moment), (AlulaPoint){fx, fy})/chord;
    double slope_moment = -pow(panel->length/chord, 2)/12;
    for (size_t j = 0; j < 3; ++j) {
        system->loads[j].cx += coefficients[j]*fx;
        system->loads[j].cy += coefficients[j]*fy;
        system->loads[j].cm += coefficients[j]*cm + slopes[j]*slope_moment;
    }
}

void alula_panel_free(AlulaPanelSystem *system)
{
    if (!system) return;
    free(system->pressure);
    free(system->segments);
    free(system->responses);
    free(system);
}

int alula_panel_dimensions(const AlulaContour *contours,size_t count,size_t *total,size_t *wet)
{
    if (!contours || !count || count>6) return 0;
    size_t n = 0, wetted = 0;
    for (size_t e = 0; e < count; ++e) {
        const AlulaContour *c = contours+e;
        if (!c->points || c->count < 4 || c->count > ALULA_PANEL_MAX_TOTAL+1) return 0;
        size_t panels = c->count-1;
        n += panels;
        if (c->finite_base) {
            if (c->upper_index >= panels || c->lower_index >= panels ||
                c->upper_index == c->lower_index) return 0;
            size_t surface = (c->lower_index+panels-c->upper_index)%panels;
            if (surface < 3 || c->count < 9) return 0;
            wetted += surface;
        } else {
            wetted += panels;
        }
    }
    if (n > ALULA_PANEL_MAX_TOTAL || wetted > ALULA_PANEL_MAX_WETTED) return 0;
    *total=n;
    *wet=wetted;
    return 1;
}

AlulaPanelSystem *alula_panel_prepare(const AlulaContour *contours,
    size_t count, double chord, AlulaPoint moment)
{
    size_t n,wetted;
    if (!(chord>0) || !isfinite(chord) || !isfinite(moment.x) || !isfinite(moment.y) ||
        !alula_panel_dimensions(contours,count,&n,&wetted)) return NULL;
    if (wetted<n) return alula_panel_prepare_streamfunction(contours, count, n, wetted, chord, moment);
    size_t size = n+count;
    AlulaPanelSystem *system = alula_panel_allocate(n);
    if (!system) return NULL;
    Panel *panels = calloc(n, sizeof *panels);
    double *a = calloc(size*size, sizeof *a);
    double *factors = malloc(size*size*sizeof *factors);
    double *tangential = calloc(n*size, sizeof *tangential);
    double *b = calloc(2*size, sizeof *b);
    system->segments = calloc(n, sizeof *system->segments);
    if (!a || !b || !factors || !tangential || !panels || !system->segments) goto fail;

    size_t first = 0;
    for (size_t e = 0; e < count; ++e) {
        const AlulaPoint *points = contours[e].points;
        size_t last = contours[e].count-1;
        if (points[0].x != points[last].x || points[0].y != points[last].y) goto fail;
        double area = 0;
        for (size_t k = 0; k < last; ++k) {
            AlulaPoint p = points[k], q = points[k+1];
            double dx = q.x-p.x, dy = q.y-p.y, length = hypot(dx, dy);
            if (!(length > 0) || !isfinite(length)) goto fail;
            area += (p.x-points[0].x)*(q.y-points[0].y)
                - (p.y-points[0].y)*(q.x-points[0].x);
            panels[first+k] = (Panel){
                p, {p.x+.5*dx, p.y+.5*dy}, {dx/length, dy/length},
                {dy/length, -dx/length}, length, first+e+k
            };
            system->segments[first+k] = (Segment){
                .geometry=panels[first+k], .end=q, .body=e, .log_length=log(length)
            };
        }
        if (!(area > 0) || !isfinite(area)) goto fail;
        a[(n+e)*size+first+e] = 1;
        a[(n+e)*size+first+e+last] = 1;
        first += last;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(n > 128)
#endif
    for (size_t i = 0; i < n; ++i) {
        const Panel *p = panels+i;
        b[i] = -p->normal.x;
        b[size+i] = -p->normal.y;
        for (size_t j = 0; j < n; ++j) {
            const Panel *source = panels+j;
            AlulaPoint q[2];
            alula_panel_vortex(p->center, source, i == j, q);
            for (size_t k = 0; k < 2; ++k) {
                a[i*size+source->node+k] += dot(q[k], p->normal);
                tangential[i*size+source->node+k] += dot(q[k], p->tangent);
            }
        }
    }
    memcpy(factors, a, size*size*sizeof *a);
    if (!alula_dense_solve(factors, b, size, 2)) goto fail;
    for (size_t i = 0; i < 2*size; ++i)
        if (!isfinite(b[i])) goto fail;

    double rx = 0, ry = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(max:rx,ry) if(n > 128)
#endif
    for (size_t i = 0; i < size; ++i) {
        double nx = i < n ? panels[i].normal.x : 0;
        double ny = i < n ? panels[i].normal.y : 0;
        for (size_t j = 0; j < size; ++j) {
            nx += a[i*size+j]*b[j];
            ny += a[i*size+j]*b[size+j];
        }
        rx = isfinite(nx) ? fmax(rx, fabs(nx)) : INFINITY;
        ry = isfinite(ny) ? fmax(ry, fabs(ny)) : INFINITY;
        if (i < n) {
            double qx = panels[i].tangent.x, qy = panels[i].tangent.y;
            for (size_t j = 0; j < size; ++j) {
                qx += tangential[i*size+j]*b[j];
                qy += tangential[i*size+j]*b[size+j];
            }
            system->pressure[i] = pressure_between(qx, qy, qx, qy);
            system->segments[i].tangential = (AlulaPoint){qx, qy};
        }
    }
    system->residual_bound = hypot(rx, ry);
    if (!(system->residual_bound < 1e-9)) goto fail;
    // Cp = 1 - (qx*u + qy*v)^2. Integrate its three coefficients once
    // using u*u + v*v = 1, then evaluate forces in constant time per angle.
    for (size_t i = 0; i < n; ++i) {
        const Panel *p = panels+i;
        alula_panel_integrate(system, i, p, system->pressure[i], (Pressure){0}, chord, moment);
        size_t k = p->node;
        system->gx += .5*(b[k]+b[k+1])*p->length/chord;
        system->gy += .5*(b[size+k]+b[size+k+1])*p->length/chord;
    }
    free(panels);
    free(a);
    free(factors);
    free(tangential);
    system->wetted = n;
    system->bodies = count;
    system->response_size = size;
    system->responses = b;
    return system;
fail:
    free(panels);
    free(a);
    free(factors);
    free(tangential);
    free(b);
    alula_panel_free(system);
    return NULL;
}

int alula_panel_solve(const AlulaPanelSystem *s, double alpha,
    AlulaPanelResult *result, double *cp)
{
    if (!s || !result || !isfinite(alpha) || fabs(alpha) > 90) return 0;
    double angle = alpha*ALULA_PI/180;
    double u = cos(angle), v = sin(angle), cx = 0, cy = 0, cm = 0;
    double weights[] = {u*u, u*v, v*v};
    for (size_t i = 0; i < 3; ++i) {
        cx += weights[i]*s->loads[i].cx;
        cy += weights[i]*s->loads[i].cy;
        cm += weights[i]*s->loads[i].cm;
    }
    if (cp) for (size_t i = 0; i < s->count; ++i) {
        Pressure p = s->pressure[i];
        cp[i] = weights[0]*p.xx + weights[1]*p.xy + weights[2]*p.yy;
    }
    *result = (AlulaPanelResult){-cx*v+cy*u, cm, cx*u+cy*v,
        -2*(u*s->gx+v*s->gy), s->residual_bound};
    return isfinite(result->cl) && isfinite(cm) && isfinite(result->pressure_drag)
        && isfinite(result->circulation_lift);
}
