// SPDX-License-Identifier: GPL-2.0-or-later
#include "panel_internal.h"
#include "dense.h"

#include <float.h>
#include <stdlib.h>
#include <string.h>

static AlulaPoint scale(AlulaPoint p, double x) { return (AlulaPoint){x*p.x, x*p.y}; }
static AlulaPoint add(AlulaPoint a, AlulaPoint b) { return (AlulaPoint){a.x+b.x, a.y+b.y}; }
static double norm(AlulaPoint a) { return hypot(a.x, a.y); }
static int same(AlulaPoint a, AlulaPoint b) { return a.x == b.x && a.y == b.y; }

static int inside(AlulaPoint p, const AlulaContour *c)
{
    int result = 0;
    for (size_t i = 1; i < c->count; ++i) {
        AlulaPoint a = c->points[i-1], b = c->points[i];
        if ((a.y > p.y) != (b.y > p.y) &&
            p.x < a.x + (b.x-a.x)*(p.y-a.y)/(b.y-a.y)) result = !result;
    }
    return result;
}

static int side(AlulaPoint a, AlulaPoint b, AlulaPoint p)
{
    AlulaPoint d = subtract(b, a), r = subtract(p, a);
    double v = cross(d, r);
    double coord = fmax(fmax(fabs(a.x), fabs(a.y)), fmax(fabs(b.x), fabs(b.y)));
    coord = fmax(coord, fmax(fabs(p.x), fabs(p.y)));
    double tol = 32*DBL_EPSILON*(fabs(d.x*r.y)+fabs(d.y*r.x)+coord*norm(d));
    return v > tol ? 1 : v < -tol ? -1 : 0;
}

static int chart_origin(const AlulaContour *c, const Segment *base, size_t nb,
    AlulaPoint center, double gap, AlulaPoint direction, AlulaPoint *origin)
{
    double tol = 32*DBL_EPSILON*fmax(gap, fmax(fabs(center.x), fabs(center.y))), distance = 0;
    for (size_t k = 0; k < nb; ++k) {
        const Panel *p = &base[k].geometry;
        double departure = dot(subtract(center, p->center), p->normal);
        if (departure <= tol) continue;
        double projection = dot(direction, p->normal);
        if (!(projection > 0)) return 0;
        distance = fmax(distance, (departure+tol)/projection);
    }
    *origin = subtract(center, scale(direction, distance));
    if (distance > 0 && !inside(*origin, c)) return 0;
    for (size_t k = 0; k < nb; ++k) {
        const Panel *p = &base[k].geometry;
        if (dot(subtract(*origin, p->center), p->normal) > tol) return 0;
        for (size_t i = 1; i < c->count; ++i) {
            AlulaPoint a = c->points[i-1], b = c->points[i];
            if (side(*origin, p->center, a)*side(*origin, p->center, b) == -1 &&
                side(a, b, *origin)*side(a, b, p->center) == -1) return 0;
        }
        if (distance > 0 && !inside(scale(add(*origin, p->center), .5), c)) return 0;
    }
    return 1;
}

static int blocked(const AlulaContour *contours, size_t count, size_t body,
    AlulaPoint origin, AlulaPoint direction)
{
    for (size_t e = 0; e < count; ++e) if (e != body) {
        const AlulaContour *c = contours+e;
        /* Include the repeated endpoint's degenerate edge as in the reference. */
        for (size_t i = 0; i < c->count; ++i) {
            AlulaPoint a = c->points[i], b = c->points[(i+1)%c->count];
            AlulaPoint edge = subtract(b, a), offset = subtract(a, origin);
            double denominator = cross(direction, edge);
            double tol = 64*DBL_EPSILON*fmax(norm(offset), norm(edge));
            if (fabs(denominator) <= tol) {
                if (fabs(cross(offset, direction)) <= tol &&
                    fmax(dot(offset, direction), dot(subtract(b, origin), direction)) >= 0) return 1;
            } else {
                double ray = cross(offset, edge)/denominator;
                double segment = cross(offset, direction)/denominator;
                if (ray >= -tol && segment >= -64*DBL_EPSILON && segment <= 1+64*DBL_EPSILON) return 1;
            }
        }
    }
    return 0;
}

static int order_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static int chart(const AlulaContour *contours, size_t count, size_t body,
    const Segment *base, size_t nb, AlulaPoint center, double gap, PanelRange *r)
{
    const AlulaContour *c = contours+body;
    AlulaPoint preferred = r->direction;
    if (!chart_origin(c, base, nb, center, gap, preferred, &r->origin)) return 0;
    if (!blocked(contours, count, body, r->origin, preferred)) return 1;

    double angles[ALULA_PANEL_MAX_TOTAL+8];
    size_t n = 0;
    angles[n++] = -ALULA_PI/2;
    angles[n++] = ALULA_PI/2;
    for (size_t e = 0; e < count; ++e) if (e != body)
        for (size_t i = 0; i < contours[e].count; ++i) {
            AlulaPoint d = subtract(contours[e].points[i], center);
            double angle = atan2(cross(preferred, d), dot(preferred, d));
            if (fabs(angle) < ALULA_PI/2) angles[n++] = angle;
        }
    qsort(angles, n, sizeof *angles, order_double);
    double angle = atan2(preferred.y, preferred.x), best = INFINITY;
    for (size_t i = 1; i < n;) {
        double lo = angles[i-1], hi = lo;
        for (; i < n; ++i) {
            if (angles[i] == angles[i-1]) continue;
            double a = angle+.5*(angles[i]+angles[i-1]);
            AlulaPoint d = {cos(a), sin(a)}, origin;
            if (!chart_origin(c, base, nb, center, gap, d, &origin) ||
                blocked(contours, count, body, origin, d)) break;
            hi = angles[i];
        }
        if (hi > lo && fabs(lo+hi) < best) {
            double a = angle+.5*(lo+hi);
            AlulaPoint d = {cos(a), sin(a)}, origin;
            if (chart_origin(c, base, nb, center, gap, d, &origin) &&
                !blocked(contours, count, body, origin, d)) {
                best = fabs(lo+hi);
                r->origin = origin;
                r->direction = d;
            }
        }
        ++i;
    }
    return isfinite(best);
}

/* Solve the parametric spline for the endpoint derivatives used by TECALC.
 * The sharp-TE probe uses the natural cubic; finite bases use SPLIND's
 * zero-third-derivative endpoints on the open wetted contour. */
static void endpoint_derivatives(const Segment *segments, size_t n, int finite,
    AlulaPoint *left, AlulaPoint *right)
{
    double diag[ALULA_PANEL_MAX_WETTED+1], upper[ALULA_PANEL_MAX_WETTED+1];
    AlulaPoint rhs[ALULA_PANEL_MAX_WETTED+1];
    diag[0] = 1;
    upper[0] = finite ? 1 : 0;
    rhs[0] = finite ? scale(segments[0].geometry.tangent, 2) : (AlulaPoint){0};
    for (size_t i = 1; i < n; ++i) {
        const Panel *a = &segments[i-1].geometry, *b = &segments[i].geometry;
        double lower = finite ? b->length : a->length, factor = lower/diag[i-1];
        diag[i] = 2*(a->length+b->length)-factor*upper[i-1];
        upper[i] = finite ? a->length : b->length;
        AlulaPoint value = finite ?
            scale(add(scale(b->tangent, a->length), scale(a->tangent, b->length)), 3) :
            scale(subtract(b->tangent, a->tangent), 6);
        rhs[i] = subtract(value, scale(rhs[i-1], factor));
    }
    if (finite) {
        double factor = 1/diag[n-1];
        diag[n] = 1-factor*upper[n-1];
        rhs[n] = subtract(scale(segments[n-1].geometry.tangent, 2), scale(rhs[n-1], factor));
    } else {
        diag[n] = 1;
        rhs[n] = (AlulaPoint){0};
    }
    rhs[n] = scale(rhs[n], 1/diag[n]);
    for (size_t i = n; i-- > 0;)
        rhs[i] = scale(subtract(rhs[i], scale(rhs[i+1], upper[i])), 1/diag[i]);
    if (finite) {
        *left = rhs[0];
        *right = rhs[n];
    } else {
        const Panel *a = &segments[0].geometry, *b = &segments[n-1].geometry;
        *left = subtract(a->tangent, scale(rhs[1], a->length/6));
        *right = add(b->tangent, scale(rhs[n-1], b->length/6));
    }
}

static double primitive(double x, double y)
{
    double r2 = x*x+y*y, ay = fabs(y);
    return (r2 == 0 ? 0 : x*log(r2))-2*x+(ay == 0 ? 0 : 2*ay*atan(x/ay));
}

static double xlogx(double x) { return x == 0 ? 0 : x*log(x); }

/* Powers include the half-panel scaling. The far-field series terminates
 * once the remaining geometric tail is below double precision. */
static void log_series(double x, double y, double r2, double *zeroth,
    double *odd, double *angle)
{
    double ir = .5*x/r2, ii = -.5*y/r2;
    double re = ir, im = ii;
    for (int k = 1; k <= 48; ++k) {
        if (k%2) {
            if (odd) *odd += .5*re/(k*(k+2));
        } else {
            *zeroth -= re/(k*(k+1));
            if (angle) *angle -= im/(k*(k+1));
        }
        double next = re*ir-im*ii;
        im = re*ii+im*ir;
        re = next;
        if (fabs(re)+fabs(im) < DBL_EPSILON/8) break;
    }
}

void alula_panel_vortex_psi(AlulaPoint point, const Segment *s, double q[2])
{
    const Panel *p = &s->geometry;
    AlulaPoint d = subtract(point, p->a);
    double x = dot(d, p->tangent)/p->length, y = -dot(d, p->normal)/p->length;
    double zx = x-.5, r2 = zx*zx+y*y, length_log = s->log_length;
    double coefficient = -p->length/(2*ALULA_PI);
    if (r2 >= 1) {
        double zeroth = length_log+.5*log(r2), odd = 0;
        log_series(zx, y, r2, &zeroth, &odd, NULL);
        q[0] = coefficient*(.5*zeroth+odd);
        q[1] = coefficient*(.5*zeroth-odd);
    } else {
        double zeroth = primitive(x, y)-primitive(x-1, y);
        double a2 = x*x+y*y, b2 = (x-1)*(x-1)+y*y;
        double first = x*zeroth-.5*(xlogx(a2)-a2-xlogx(b2)+b2);
        q[0] = .5*coefficient*(length_log+zeroth-first);
        q[1] = .5*coefficient*(length_log+first);
    }
}

double alula_panel_base_psi(AlulaPoint point, const Segment *s, double source, double vortex)
{
    const Panel *p = &s->geometry;
    AlulaPoint d = subtract(point, p->a), common = subtract(point, s->cut_origin);
    int at_start = same(point, p->a), at_end = same(point, s->end);
    double cut_cross = cross(s->cut_direction, common), cut_dot = dot(s->cut_direction, common);
    double coordinate_max = fmax(fmax(fabs(point.x), fabs(point.y)),
        fmax(fabs(s->cut_origin.x), fabs(s->cut_origin.y)));
    double tol = 16*DBL_EPSILON*fmax(p->length, coordinate_max);
    if (!at_start && !at_end && source != 0 && cut_dot >= 0 && fabs(cut_cross) <= tol)
        return NAN;
    double x = at_start ? 0 : at_end ? 1 : dot(d, p->tangent)/p->length;
    double y = at_start || at_end ? 0 : -dot(d, p->normal)/p->length;
    double zx = x-.5, r2 = zx*zx+y*y, midpoint_angle = atan2(y, zx);
    double real = s->log_length, imaginary;
    if (r2 >= 1) {
        real += .5*log(r2);
        imaginary = midpoint_angle;
        log_series(zx, y, r2, &real, NULL, &imaginary);
    } else {
        real += .5*(primitive(x, y)-primitive(x-1, y));
        imaginary = x*atan2(y, x)-(x-1)*atan2(y, x-1);
        if (y != 0) imaginary += .5*y*log((x*x+y*y)/((x-1)*(x-1)+y*y));
    }
    AlulaPoint relative = subtract(point, p->center);
    double angle = atan2(cut_cross, cut_dot);
    if (angle < 0) angle += 2*ALULA_PI;
    angle += atan2(cross(common, relative), dot(common, relative));
    return p->length*(source*(imaginary+angle-midpoint_angle)-vortex*real)/(2*ALULA_PI);
}

int alula_panel_boundary_setup(const AlulaContour *contours, size_t count, size_t total,
    Segment *segments, PanelRange *ranges, PanelRow *rows, size_t wetted)
{
    size_t wet = 0, base = total, offset = 0;
    for (size_t e = 0; e < count; ++e) {
        const AlulaContour *c = contours+e;
        size_t n = c->count-1, upper = c->finite_base ? c->upper_index : 0;
        size_t nw = c->finite_base ? (c->lower_index+n-upper)%n : n;
        size_t nb = n-nw;
        if (!same(c->points[0], c->points[n])) return 0;
        PanelRange *r = ranges+e;
        *r = (PanelRange){.start=wet+e, .end=wet+e+nw, .offset=wet, .wetted=nw, .finite=c->finite_base};
        base -= nb;
        double area = 0;
        for (size_t k = 0; k < n; ++k) {
            size_t idx = (upper+k)%n;
            AlulaPoint a = c->points[idx], b = c->points[idx+1], d = subtract(b, a);
            double length = norm(d);
            if (!(length > 0) || !isfinite(length)) return 0;
            area += cross(subtract(a, c->points[0]), subtract(b, c->points[0]));
            AlulaPoint t = scale(d, 1/length);
            Segment *s = segments+(k < nw ? wet+k : base+k-nw);
            *s = (Segment){
                .geometry={a, scale(add(a, b), .5), t, {t.y, -t.x}, length, r->start+k},
                .end=b, .output=offset+idx, .body=e, .log_length=log(length)
            };
            if (k < nw) rows[r->start+k] = (PanelRow){.point=a, .body=e};
        }
        if (!(area > 0) || !isfinite(area)) return 0;
        AlulaPoint left, right;
        if (wetted<total) endpoint_derivatives(segments+wet,nw,r->finite,&left,&right);
        else {
            left=segments[wet].geometry.tangent;
            right=segments[wet+nw-1].geometry.tangent;
        }
        r->wake_direction=subtract(segments[wet+nw-1].geometry.tangent,segments[wet].geometry.tangent);
        AlulaPoint top = segments[wet].geometry.a, bottom = segments[wet+nw-1].end;
        if (r->finite) {
            AlulaPoint t = scale(subtract(right, left), .5);
            double length = norm(t);
            if (!(length > 0) || !isfinite(length)) return 0;
            r->wake_direction = r->direction = scale(t, 1/length);
            AlulaPoint gap = subtract(top, bottom);
            AlulaPoint bisector = subtract(segments[wet+nw-1].geometry.tangent, segments[wet].geometry.tangent);
            if (!(cross(bisector, gap) > 0)) return 0;
            if (!chart(contours, count, e, segments+base, nb, scale(add(top, bottom), .5), norm(gap), r))
                return 0;
            for (size_t k = 0; k < nb; ++k) {
                Segment *s = segments+base+k;
                s->source = -.5*cross(t, s->geometry.tangent);
                s->vortex = -.5*dot(t, s->geometry.tangent);
                s->cut_origin = r->origin;
                s->cut_direction = r->direction;
                s->geometry.node = r->start;
                s->last_node = r->end;
            }
            rows[r->end] = (PanelRow){.point=bottom, .body=e};
        } else {
            AlulaPoint t = subtract(scale(right, 1/norm(right)), scale(left, 1/norm(left)));
            double length = norm(t);
            if (!(length > 0) || !isfinite(length)) return 0;
            t = scale(t, 1/length);
            double distance = .1*fmin(segments[wet].geometry.length, segments[wet+nw-1].geometry.length);
            AlulaPoint probe = subtract(top, scale(t, distance));
            if (!inside(probe, c)) return 0;
            rows[r->end] = (PanelRow){probe, t, e, 1};
        }
        wet += nw;
        offset += n;
    }
    return wet==wetted;
}

static void assemble_row(double *a, double *bx, double *by, const PanelRow *row, size_t strength,
    const Segment *segments, size_t wet, size_t total, const PanelRange *ranges, double chord)
{
    double inv_chord = 1/chord;
    if (row->interior) {
        *bx = -row->tangent.x;
        *by = -row->tangent.y;
    } else {
        *bx = -row->point.y*inv_chord;
        *by = row->point.x*inv_chord;
        a[strength+row->body] = -1;
    }
    for (size_t j = 0; j < wet; ++j) {
        const Segment *s = segments+j;
        double q[2];
        if (row->interior) {
            AlulaPoint velocity[2];
            alula_panel_vortex(row->point, &s->geometry, 0, velocity);
            for (size_t k = 0; k < 2; ++k) q[k] = dot(velocity[k], row->tangent);
        } else {
            alula_panel_vortex_psi(row->point, s, q);
            q[0] *= inv_chord;
            q[1] *= inv_chord;
        }
        a[s->geometry.node] += q[0];
        a[s->geometry.node+1] += q[1];
    }
    for (size_t j = wet; j < total; ++j) {
        const Segment *s = segments+j;
        const PanelRange *r = ranges+s->body;
        double q;
        if (row->interior) {
            AlulaPoint velocity[2];
            alula_panel_vortex(row->point, &s->geometry, 0, velocity);
            AlulaPoint v = add(velocity[0], velocity[1]), source = {v.y, -v.x};
            q = s->source*dot(source, row->tangent)+s->vortex*dot(v, row->tangent);
        } else {
            q = alula_panel_base_psi(row->point, s, s->source, s->vortex)*inv_chord;
        }
        a[r->start] += q;
        a[r->end] -= q;
    }
}

void alula_panel_boundary_assemble(double *a,double *bx,double *by,size_t count,
    const Segment *segments,size_t wet,size_t total,const PanelRow *rows,const PanelRange *ranges,double chord)
{
    size_t strength=wet+count, size=strength+count;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(wet > 128)
#endif
    for (size_t i = 0; i < strength; ++i)
        assemble_row(a+i*size, bx+i, by+i, rows+i, strength, segments, wet, total, ranges, chord);
    for (size_t e = 0; e < count; ++e) {
        a[(strength+e)*size+ranges[e].start] = 1;
        a[(strength+e)*size+ranges[e].end] = 1;
    }
}

AlulaPanelSystem *alula_panel_prepare_streamfunction(const AlulaContour *contours,
    size_t count, size_t total, size_t wet, double chord, AlulaPoint moment)
{
    AlulaPanelSystem *system = alula_panel_allocate(total);
    Segment *segments = calloc(total, sizeof *segments);
    PanelRange ranges[6];
    PanelRow rows[ALULA_PANEL_MAX_WETTED+6];
    double *a = NULL, *factors = NULL, *b = NULL;
    int success = 0;
    if (!system || !segments || !alula_panel_boundary_setup(contours,count,total,segments,ranges,rows,wet))
        goto done;
    size_t strength = wet+count, size = strength+count;
    a = calloc(size*size, sizeof *a);
    factors = malloc(size*size*sizeof *factors);
    b = calloc(2*size, sizeof *b);
    if (!a || !factors || !b) goto done;
    alula_panel_boundary_assemble(a,b,b+size,count,segments,wet,total,rows,ranges,chord);
    memcpy(factors, a, size*size*sizeof *a);
    if (!alula_dense_solve(factors, b, size, 2)) goto done;
    double rx = 0, ry = 0;
    for (size_t i = 0; i < size; ++i) {
        double nx = 0, ny = 0;
        if (i < strength) {
            const PanelRow *r = rows+i;
            nx = r->interior ? r->tangent.x : r->point.y/chord;
            ny = r->interior ? r->tangent.y : -r->point.x/chord;
        }
        for (size_t j = 0; j < size; ++j) {
            nx += a[i*size+j]*b[j];
            ny += a[i*size+j]*b[size+j];
        }
        if (!isfinite(nx) || !isfinite(ny)) goto done;
        rx = fmax(rx, fabs(nx));
        ry = fmax(ry, fabs(ny));
    }
    system->residual_bound = hypot(rx, ry);
    if (!(system->residual_bound < 1e-9)) goto done;
    for (size_t i = 0; i < total; ++i) {
        Segment *s = segments+i;
        const PanelRange *r = ranges+s->body;
        int base = i >= wet;
        size_t a0 = base ? r->start : s->geometry.node;
        size_t a1 = base ? r->end : a0+1;
        double ax = b[a0], ay = b[size+a0], bx = b[a1], by = b[size+a1];
        if (!base) s->tangential = (AlulaPoint){.5*(ax+bx), .5*(ay+by)};
        Pressure p = pressure_between(ax, ay, bx, by);
        Pressure gradient = base ? (Pressure){0} :
            (Pressure){ax*ax-bx*bx, 2*(ax*ay-bx*by), ay*ay-by*by};
        alula_panel_integrate(system, s->output, &s->geometry, p, gradient, chord, moment);
        double gx = base ? s->vortex*(ax-bx) : .5*(ax+bx);
        double gy = base ? s->vortex*(ay-by) : .5*(ay+by);
        system->gx += gx*s->geometry.length/chord;
        system->gy += gy*s->geometry.length/chord;
    }
    system->wetted = wet;
    system->bodies = count;
    system->dirichlet = 1;
    for (size_t e = 0; e < count; ++e)
        system->boundary_psi[e] = (AlulaPoint){chord*b[strength+e], chord*b[size+strength+e]};
    system->response_size = size;
    system->segments = segments;
    system->responses = b;
    segments = NULL;
    b = NULL;
    success = 1;
done:
    free(a);
    free(factors);
    free(b);
    free(segments);
    if (success) return system;
    alula_panel_free(system);
    return NULL;
}
