// SPDX-License-Identifier: GPL-2.0-or-later
#include "panel_internal.h"

#include <stdlib.h>

/* The complex velocity of a straight vortex sheet is proportional to
 * the integral of gamma(s)/(z-s). Expand the denominator only for
 * distant targets. Eight-point Gauss integration exactly integrates
 * each of the 15 moments of the linear sheet in exact arithmetic. */
enum { PANELS_PER_GROUP=16, MOMENT_COUNT=15 };
struct PanelGroup {
    AlulaPoint center, moments[MOMENT_COUNT];
    double radius, far_distance_sq;
};

static AlulaPoint product(AlulaPoint a,AlulaPoint b)
{
    return (AlulaPoint){a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x};
}

static void prepare_group(const AlulaPanelField *field,struct PanelGroup *group,size_t first,size_t end)
{
    const Segment *segments=field->system->segments;
    AlulaPoint lo=segments[first].geometry.a,hi=lo;
    for (size_t i=first;i<end;++i) {
        AlulaPoint points[]={segments[i].geometry.a,segments[i].end};
        for (int j=0;j<2;++j) {
            lo.x=fmin(lo.x,points[j].x); hi.x=fmax(hi.x,points[j].x);
            lo.y=fmin(lo.y,points[j].y); hi.y=fmax(hi.y,points[j].y);
        }
    }
    group->center=(AlulaPoint){.5*lo.x+.5*hi.x,.5*lo.y+.5*hi.y};
    double rx=fmax(hi.x-group->center.x,group->center.x-lo.x);
    double ry=fmax(hi.y-group->center.y,group->center.y-lo.y);
    group->radius=hypot(rx,ry);
    /* The omitted geometric-series tail is bounded by (1/16)^15/(1-1/16)
     * times the absolute integrated sheet strength divided by 2*pi*distance. */
    double cutoff=16*group->radius;
    group->far_distance_sq=cutoff*cutoff;
    if (!(group->far_distance_sq>0) || !isfinite(group->far_distance_sq)) {
        group->far_distance_sq=INFINITY;
        return;
    }
    static const double x[]={.0198550717512319,.1016667612931866,.2372337950418355,.4082826787521751,
        .5917173212478249,.7627662049581645,.8983332387068134,.9801449282487681};
    static const double weight[]={.0506142681451881,.1111905172266872,.1568533229389436,.1813418916891810,
        .1813418916891810,.1568533229389436,.1111905172266872,.0506142681451881};
    for (size_t i=first;i<end;++i) {
        AlulaPoint a=subtract(segments[i].geometry.a,group->center);
        AlulaPoint step=subtract(segments[i].end,segments[i].geometry.a);
        AlulaPoint gamma=field->strengths[i];
        for (int j=0;j<8;++j) {
            AlulaPoint pos={(a.x+x[j]*step.x)/group->radius,(a.y+x[j]*step.y)/group->radius};
            double g=(1-x[j])*gamma.x+x[j]*gamma.y;
            AlulaPoint term={g*segments[i].geometry.length*weight[j]/(2*ALULA_PI),0};
            for (int k=0;k<MOMENT_COUNT;++k) {
                group->moments[k].x+=term.x;
                group->moments[k].y+=term.y;
                term=product(term,pos);
            }
        }
    }
    for (int k=0;k<MOMENT_COUNT;++k)
        if (!isfinite(group->moments[k].x) || !isfinite(group->moments[k].y))
            group->far_distance_sq=INFINITY;
}

static int add_distant_velocity(const struct PanelGroup *group,AlulaPoint point,AlulaPoint *q)
{
    AlulaPoint d=subtract(point,group->center);
    double r2=dot(d,d);
    if (!(r2>group->far_distance_sq) || !isfinite(r2)) return 0;
    AlulaPoint inv={d.x/r2,-d.y/r2};
    AlulaPoint ratio={group->radius*inv.x,group->radius*inv.y};
    AlulaPoint sum=group->moments[MOMENT_COUNT-1];
    for (int k=MOMENT_COUNT-2;k>=0;--k) {
        sum=product(sum,ratio);
        sum.x+=group->moments[k].x;
        sum.y+=group->moments[k].y;
    }
    sum=product(sum,inv);
    q->x+=sum.y;
    q->y+=sum.x;
    return 1;
}

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
    size_t groups=(system->wetted+PANELS_PER_GROUP-1)/PANELS_PER_GROUP;
    field->groups=calloc(groups,sizeof *field->groups);
    if (field->groups) for (size_t k=0;k<groups;++k) {
        size_t first=k*PANELS_PER_GROUP,end=first+PANELS_PER_GROUP;
        prepare_group(field,field->groups+k,first,end<system->wetted?end:system->wetted);
    }
    return field;
}

void alula_panel_field_free(AlulaPanelField *field)
{
    if (field) free(field->groups);
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
    for (size_t first=0;first<system->wetted;first+=PANELS_PER_GROUP) {
        size_t end=first+PANELS_PER_GROUP;
        if (end>system->wetted) end=system->wetted;
        int distant=velocity && field->groups &&
            add_distant_velocity(field->groups+first/PANELS_PER_GROUP,point,&q);
        if (distant && !streamfunction) continue;
        for (size_t i=first;i<end;++i) {
            const Segment *s = system->segments+i;
            AlulaPoint g = field->strengths[i];
            if (velocity && !distant) {
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
