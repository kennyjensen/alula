// SPDX-License-Identifier: GPL-2.0-or-later
#include "trace.h"

#include <float.h>
#include <math.h>

typedef struct {
    AlulaPoint position;
    double potential;
} State;

static int finite_point(AlulaPoint p)
{
    return isfinite(p.x) && isfinite(p.y);
}

static AlulaPoint advance(AlulaPoint p, AlulaPoint direction, double step)
{
    return (AlulaPoint){p.x+step*direction.x, p.y+step*direction.y};
}

static double distance(AlulaPoint a, AlulaPoint b)
{
    return hypot(a.x-b.x, a.y-b.y);
}

static int slope(const AlulaPanelField *field, const AlulaTraceControls *c,
    AlulaPoint p, State *result)
{
    AlulaPoint v;
    if (!finite_point(p) || (c->admissible_point && !c->admissible_point(c->context, p)) ||
        !alula_panel_field_sample(field, p, &v, NULL)) return 0;
    double speed_squared = v.x*v.x+v.y*v.y;
    if (!(speed_squared > 0) || !isfinite(speed_squared)) return 0;
    double speed = c->parameter == ALULA_TRACE_ARC ? sqrt(speed_squared) : 1;
    double scale = c->parameter == ALULA_TRACE_ARC ? speed : speed_squared;
    *result = (State){{v.x/scale, v.y/scale}, speed};
    return finite_point(result->position);
}

static int rk4(const AlulaPanelField *field, const AlulaTraceControls *c,
    State p, State a, double step, State *out)
{
    State b, d, e;
    if (!slope(field, c, advance(p.position, a.position, .5*step), &b) ||
        !slope(field, c, advance(p.position, b.position, .5*step), &d) ||
        !slope(field, c, advance(p.position, d.position, step), &e)) return 0;
    AlulaPoint sum = {a.position.x+2*b.position.x+2*d.position.x+e.position.x,
        a.position.y+2*b.position.y+2*d.position.y+e.position.y};
    *out = (State){advance(p.position, sum, step/6),
        p.potential + step*(a.potential+2*b.potential+2*d.potential+e.potential)/6};
    return 1;
}

static int trial(const AlulaPanelField *field, const AlulaTraceControls *c,
    State p, State a, double step, State *out, double *error)
{
    State whole, mid, mid_slope, half;
    if (!rk4(field, c, p, a, step, &whole) ||
        !rk4(field, c, p, a, .5*step, &mid) ||
        !slope(field, c, mid.position, &mid_slope) ||
        !rk4(field, c, mid, mid_slope, .5*step, &half)) return 0;
    State correction = {{(half.position.x-whole.position.x)/15,
        (half.position.y-whole.position.y)/15}, (half.potential-whole.potential)/15};
    *out = (State){advance(half.position, correction.position, 1), half.potential+correction.potential};
    if (!finite_point(out->position) || !isfinite(out->potential) ||
        (c->admissible_point && !c->admissible_point(c->context, out->position)) ||
        (c->admissible_segment && !c->admissible_segment(c->context, p.position, out->position))) return 0;
    *error = hypot(correction.position.x, correction.position.y)/c->tolerance;
    if (c->parameter == ALULA_TRACE_ARC)
        *error = fmax(*error, fabs(correction.potential)/c->potential_tolerance);
    return 1;
}

static double plane_distance(AlulaPoint p, const AlulaTraceControls *c, double length)
{
    return (p.x*c->plane_normal.x+p.y*c->plane_normal.y-c->plane_offset)/length;
}

static AlulaTraceStatus event_step(const AlulaPanelField *field, const AlulaTraceControls *c,
    State p, State a, double plane_length, double side, double *step, State *next, double error)
{
    double lower = 0, upper = 1, full = *step;
    double plane = plane_distance(next->position, c, plane_length);
    double low_value = side*plane_distance(p.position,c,plane_length), high_value = side*plane;
    int last_side = 0;
    for (int i = 0; i < 60 && fabs(plane) > c->tolerance; ++i) {
        /* Interpolate the signed distance, halving a repeatedly retained
         * endpoint's weight to keep the bracket contracting on curved paths. */
        double fraction = lower+(upper-lower)*low_value/(low_value-high_value);
        if (!(fraction>lower && fraction<upper)) fraction = .5*(lower+upper);
        *step = full*fraction;
        if (!trial(field, c, p, a, *step, next, &error)) return ALULA_TRACE_DOMAIN;
        plane = plane_distance(next->position, c, plane_length);
        if (side*plane>0) {
            lower=fraction;
            low_value=side*plane;
            if (last_side>0) high_value*=.5;
            last_side=1;
        } else {
            upper=fraction;
            high_value=side*plane;
            if (last_side<0) low_value*=.5;
            last_side=-1;
        }
    }
    if (!(fabs(plane) <= c->tolerance) || error > 1) return ALULA_TRACE_PLANE_UNRESOLVED;
    return ALULA_TRACE_PLANE;
}

AlulaTraceResult alula_trace_panel(const AlulaPanelField *field,
    AlulaTracePoint seed, const AlulaTraceControls *c)
{
    AlulaTraceResult r = {.status=ALULA_TRACE_INVALID, .last=seed};
    if (!field || !c || !finite_point(seed.position) || !isfinite(seed.parameter) ||
        !isfinite(seed.potential) || !isfinite(c->end) || seed.parameter == c->end ||
        !(c->tolerance > 0) || !isfinite(c->tolerance) || !(c->min_step > 0) ||
        !(c->max_step >= c->min_step) || !isfinite(c->max_step) || c->max_steps < 1 ||
        (c->parameter != ALULA_TRACE_ARC && c->parameter != ALULA_TRACE_POTENTIAL)) return r;
    int by_arc = c->parameter == ALULA_TRACE_ARC;
    if (by_arc ? (!(c->potential_tolerance > 0) || !isfinite(c->potential_tolerance)) :
        (!(c->max_spatial_step > 0) || !isfinite(c->max_spatial_step) || seed.parameter != seed.potential)) return r;
    double plane_length = 1, side = 0;
    if (c->stop_at_plane) {
        plane_length = hypot(c->plane_normal.x, c->plane_normal.y);
        if (!(plane_length > 0) || !isfinite(plane_length) || !isfinite(c->plane_offset)) return r;
        side = plane_distance(seed.position, c, plane_length);
        if (!isfinite(side)) return r;
    }
    if (c->admissible_point && !c->admissible_point(c->context, seed.position)) return r;
    if (c->visit) c->visit(c->context, &r.last);
    if (c->stop_at_plane && fabs(side) <= c->tolerance) {
        r.status = ALULA_TRACE_PLANE;
        return r;
    }
    side = copysign(1, side);
    double sign = copysign(1, c->end-seed.parameter);
    double small = 4*DBL_EPSILON*(by_arc ? fmax(DBL_MIN, fmax(fabs(seed.parameter), fabs(c->end))) :
        fmax(1, fabs(c->end)));
    double step_size = c->max_step;
    r.status = ALULA_TRACE_STEP_LIMIT;
    AlulaTraceStatus last_rejection = ALULA_TRACE_MIN_STEP;
    for (int i = 0; i < c->max_steps; ++i) {
        double remaining = fabs(c->end-r.last.parameter);
        if (remaining <= small) {
            r.status = c->stop_at_plane ? ALULA_TRACE_PLANE_NOT_REACHED : ALULA_TRACE_END;
            break;
        }
        State p = {r.last.position, r.last.potential}, a, next;
        if (!slope(field, c, p.position, &a)) { r.status = ALULA_TRACE_DOMAIN; break; }
        double magnitude = fmin(step_size, remaining);
        if (!by_arc) magnitude = fmin(magnitude, c->max_spatial_step/hypot(a.position.x, a.position.y));
        if (magnitude < c->min_step) { r.status = last_rejection; break; }
        double step = sign*magnitude, error;
        if (!trial(field, c, p, a, step, &next, &error)) {
            step_size = .5*magnitude;
            last_rejection = ALULA_TRACE_DOMAIN;
            ++r.rejected;
            continue;
        }
        if (error > 1 || (!by_arc && distance(p.position, next.position) > c->max_spatial_step)) {
            step_size = magnitude*fmax(.1, fmin(.5, .8*pow(error, -.2)));
            last_rejection = ALULA_TRACE_MIN_STEP;
            ++r.rejected;
            continue;
        }
        int at_plane = c->stop_at_plane && side*plane_distance(next.position, c, plane_length) <= 0;
        if (at_plane) {
            r.status = event_step(field, c, p, a, plane_length, side, &step, &next, error);
            if (r.status != ALULA_TRACE_PLANE) break;
        }
        r.last.parameter += step;
        r.last.position = next.position;
        r.last.potential = by_arc ? next.potential : r.last.parameter;
        ++r.accepted;
        if (c->visit) c->visit(c->context, &r.last);
        if (at_plane) break;
        if (fabs(c->end-r.last.parameter) <= small) {
            r.status = c->stop_at_plane ? ALULA_TRACE_PLANE_NOT_REACHED : ALULA_TRACE_END;
            break;
        }
        step_size = fmin(c->max_step, magnitude*fmax(.5, fmin(2, .9*pow(fmax(error, 1e-30), -.2))));
    }
    return r;
}
