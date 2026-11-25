// SPDX-License-Identifier: GPL-2.0-or-later
#include "panel_guide.h"
#include "panel_internal.h"

#include <float.h>
#include <stdlib.h>

typedef struct { AlulaPoint lo, hi; } Box;
typedef struct AlulaPanelStreamline {
    AlulaTracePoint *points;
    size_t count, capacity;
    int failed;
} Path;

struct AlulaPanelGuides {
    const AlulaPanelField *field;
    double tolerance, psi_tolerance, scale, outlet_x;
    Box bounds[6];
    Box *edges;
    AlulaContourSpline curves[6];
    AlulaPanelBodyGuide bodies[6];
    double levels[8];
    Path incoming[6], wake[6], farfield[2];
};

typedef struct { const AlulaPanelGuides *guides; Path *path; } TraceContext;

static Box box(AlulaPoint a, AlulaPoint b)
{
    return (Box){{fmin(a.x,b.x), fmin(a.y,b.y)}, {fmax(a.x,b.x), fmax(a.y,b.y)}};
}

static int overlap(Box a, Box b, double tol)
{
    return a.lo.x <= b.hi.x+tol && b.lo.x <= a.hi.x+tol &&
        a.lo.y <= b.hi.y+tol && b.lo.y <= a.hi.y+tol;
}

static double length(AlulaPoint v) { return hypot(v.x,v.y); }

static AlulaPoint advance(AlulaPoint p, AlulaPoint t, double step)
{
    return (AlulaPoint){p.x+step*t.x, p.y+step*t.y};
}

static int side(AlulaPoint a, AlulaPoint b, AlulaPoint point, double tol)
{
    AlulaPoint d = subtract(b,a), r = subtract(point,a);
    double crossing = cross(d,r);
    double error = tol*(fabs(d.x)+fabs(d.y)) +
        8*DBL_EPSILON*(fabs(d.x*r.y)+fabs(d.y*r.x));
    return crossing>error ? 1 : crossing<-error ? -1 : 0;
}

static int exterior(const AlulaPanelGuides *g, AlulaPoint p)
{
    const AlulaPanelSystem *s = g->field->system;
    int active = 0, inside = 0;
    double tol = 1e-12*g->scale;
    Box point = {p,p};
    for (size_t e=0; e<s->bodies; ++e) if (overlap(g->bounds[e],point,tol)) active |= 1<<e;
    if (!active) return 1;
    for (size_t i=0; i<s->count; ++i) {
        const Segment *edge = s->segments+i;
        int bit = 1<<edge->body;
        if (!(active & bit)) continue;
        AlulaPoint a = edge->geometry.a, b = edge->end;
        if (overlap(g->edges[i],point,tol) && side(a,b,p,tol) == 0) return 0;
        if ((a.y>p.y) != (b.y>p.y) &&
            p.x < a.x+(b.x-a.x)*(p.y-a.y)/(b.y-a.y)) inside ^= bit;
    }
    return !inside;
}

static int clear_chord(const AlulaPanelGuides *g, AlulaPoint a, AlulaPoint b)
{
    const AlulaPanelSystem *s = g->field->system;
    Box span = box(a,b);
    double tol = 1e-12*g->scale;
    int active = 0;
    for (size_t e=0; e<s->bodies; ++e) if (overlap(g->bounds[e],span,tol)) active |= 1<<e;
    if (!active) return 1;
    for (size_t i=0; i<s->count; ++i) {
        const Segment *edge = s->segments+i;
        if (!(active & (1<<edge->body)) || !overlap(g->edges[i],span,tol)) continue;
        AlulaPoint c = edge->geometry.a, d = edge->end;
        if (side(a,b,c,tol)*side(a,b,d,tol)<=0 && side(c,d,a,tol)*side(c,d,b,tol)<=0) return 0;
    }
    return 1;
}

static int trace_point(void *context, AlulaPoint p)
{
    const TraceContext *c = context;
    return (!c->path || !c->path->failed) && exterior(c->guides,p);
}

static int trace_chord(void *context, AlulaPoint a, AlulaPoint b)
{
    const TraceContext *c = context;
    return clear_chord(c->guides,a,b);
}

static void record(void *context, const AlulaTracePoint *p)
{
    Path *path = ((TraceContext *)context)->path;
    if (path->count == path->capacity) {
        size_t capacity = path->capacity ? 2*path->capacity : 64;
        AlulaTracePoint *points = realloc(path->points,capacity*sizeof *points);
        if (!points) { path->failed=1; return; }
        path->points=points;
        path->capacity=capacity;
    }
    path->points[path->count++]=*p;
}

static AlulaPanelGuideStatus trace(const AlulaPanelGuides *g, AlulaPoint p,
    double plane, double potential_end, Path *path)
{
    TraceContext context = {g,path};
    AlulaTraceControls controls = {
        .parameter=ALULA_TRACE_POTENTIAL, .end=potential_end, .tolerance=g->tolerance,
        .max_step=.05*g->scale, .min_step=1e-12*g->scale, .max_spatial_step=.1*g->scale,
        .max_steps=20000, .stop_at_plane=1, .plane_normal={1,0}, .plane_offset=plane,
        .context=&context, .admissible_point=trace_point, .admissible_segment=trace_chord,
        .visit=record
    };
    AlulaTraceResult r = alula_trace_panel(g->field,(AlulaTracePoint){.position=p},&controls);
    if (path->failed) return ALULA_GUIDE_ALLOCATION;
    return r.status == ALULA_TRACE_PLANE && path->count>1 ? ALULA_GUIDE_OK : ALULA_GUIDE_TRACE;
}

static int project_level(const AlulaPanelGuides *g, AlulaPoint *p, double level, AlulaPoint direction)
{
    double psi;
    if (!exterior(g,*p) || !alula_panel_field_sample(g->field,*p,NULL,&psi)) return 0;
    for (int i=0; i<40; ++i) {
        double error=psi-level;
        if (fabs(error)<=g->psi_tolerance) return 1;
        AlulaPoint q;
        if (!alula_panel_field_sample(g->field,*p,&q,NULL)) return 0;
        double derivative=cross(q,direction);
        if (derivative==0) return 0;
        double step=-error/derivative;
        if (!isfinite(step)) return 0;
        step=copysign(fmin(fabs(step),.1*g->scale),step);
        int accepted=0;
        for (int j=0; j<=24; ++j,step*=.5) {
            AlulaPoint next=advance(*p,direction,step);
            double trial;
            if (!exterior(g,next) || !clear_chord(g,*p,next) ||
                alula_panel_field_branch_increment(g->field,*p,next)!=0 ||
                !alula_panel_field_sample(g->field,next,NULL,&trial) ||
                fabs(trial-level)>=fabs(error)) continue;
            *p=next;
            psi=trial;
            accepted=1;
            break;
        }
        if (!accepted) return 0;
    }
    return fabs(psi-level)<=g->psi_tolerance;
}

static int stagnation(const AlulaPanelField *field, size_t first, size_t end,
    const AlulaContourSpline *curve, double *parameter, double *slope)
{
    const Segment *segments=field->system->segments;
    double last=dot(field->freestream,segments[first].tangential);
    int crossings=0;
    if (!(last<0)) return 0;
    for (size_t i=first+1; i<end; ++i) {
        double next=dot(field->freestream,segments[i].tangential);
        if (last>=0 && next<0) return 0;
        if (last<0 && next>=0) {
            size_t k=i-first;
            double a=.5*(curve->arc[k-1]+curve->arc[k]);
            double b=.5*(curve->arc[k]+curve->arc[k+1]);
            *slope=(next-last)/(b-a);
            *parameter=a-last/ *slope;
            ++crossings;
        }
        last=next;
    }
    return crossings==1 && last>0 && isfinite(*slope);
}

static AlulaPanelGuideStatus prepare_body(AlulaPanelGuides *g, const AlulaPanelGuideControls *c,
    size_t e,size_t first,size_t end)
{
    const AlulaPanelSystem *system=g->field->system;
    AlulaPoint points[ALULA_PANEL_MAX_WETTED+1];
    for (size_t i=first; i<end; ++i) points[i-first]=system->segments[i].geometry.a;
    points[end-first]=system->segments[end-1].end;
    AlulaContourSpline *curve=g->curves+e;
    if (!alula_contour_spline_create(curve,points,(int)(end-first+1))) return ALULA_GUIDE_ALLOCATION;
    AlulaPanelBodyGuide *body=g->bodies+e;
    *body=(AlulaPanelBodyGuide){.element=e,.curve=curve};
    if (!stagnation(g->field,first,end,curve,&body->stagnation,&body->strain)) return ALULA_GUIDE_STAGNATION;
    AlulaPoint tangent;
    alula_contour_spline_value(curve,body->stagnation,&body->stagnation_point,&tangent,NULL);
    double speed=length(tangent);
    if (!(speed>0) || !isfinite(speed)) return ALULA_GUIDE_STAGNATION;
    body->strain/=speed;
    AlulaPoint unit={tangent.x/speed,tangent.y/speed}, normal={unit.y,-unit.x};
    double chord=0;
    for (int k=1; k<curve->count; ++k) chord=fmax(chord,length(subtract(points[k],points[0])));
    double seed=c->seed_fraction*chord;
    if (c->seed_limit) seed=fmin(seed,c->seed_limit[e]);
    body->seed_distance=seed;
    AlulaPoint upper=points[0], lower=points[end-first];
    body->trailing_edge=(AlulaPoint){.5*(upper.x+lower.x),.5*(upper.y+lower.y)};
    body->trailing_vector=subtract(upper,lower);
    if (!(body->stagnation_point.x>c->inlet_x && body->trailing_edge.x<c->outlet_x))
        return ALULA_GUIDE_INVALID;
    if (system->dirichlet) body->wall_level=dot(system->boundary_psi[e],g->field->freestream);
    else if (!alula_panel_field_sample(g->field,body->stagnation_point,NULL,&body->wall_level))
        return ALULA_GUIDE_SEED;
    AlulaPoint incoming=advance(body->stagnation_point,normal,seed);
    if (!project_level(g,&incoming,body->wall_level,unit)) return ALULA_GUIDE_SEED;
    double budget=4*(c->outlet_x-c->inlet_x);
    AlulaPanelGuideStatus status=trace(g,incoming,c->inlet_x,-budget,g->incoming+e);
    if (status!=ALULA_GUIDE_OK) return status;
    Path *path=g->incoming+e;
    body->incoming=(AlulaGuidePath){path->points,path->count};
    body->inlet_level=body->wall_level;
    for (size_t i=1; i<path->count; ++i)
        body->inlet_level+=alula_panel_field_branch_increment(g->field,path->points[i-1].position,path->points[i].position);
    AlulaPoint wake_seed;
    if (upper.x!=lower.x || upper.y!=lower.y) {
        AlulaPoint a,b;
        alula_contour_spline_value(curve,0,NULL,&a,NULL);
        alula_contour_spline_value(curve,curve->arc[curve->count-1],NULL,&b,NULL);
        AlulaPoint direction=subtract(b,a);
        double norm=length(direction);
        if (!(norm>0) || !isfinite(norm)) return ALULA_GUIDE_SEED;
        wake_seed=advance(body->trailing_edge,direction,seed/norm);
    } else {
        wake_seed=advance(body->trailing_edge,(AlulaPoint){1,0},seed);
        if (!project_level(g,&wake_seed,body->wall_level,(AlulaPoint){0,1})) return ALULA_GUIDE_SEED;
    }
    status=trace(g,wake_seed,c->outlet_x,budget,g->wake+e);
    if (status!=ALULA_GUIDE_OK) return status;
    body->wake=(AlulaGuidePath){g->wake[e].points,g->wake[e].count};
    return ALULA_GUIDE_OK;
}

static int order(const void *left,const void *right)
{
    double a=((const AlulaPanelBodyGuide *)left)->inlet_level;
    double b=((const AlulaPanelBodyGuide *)right)->inlet_level;
    return (a>b)-(a<b);
}

AlulaPanelGuides *alula_panel_guides_create(const AlulaPanelField *field,
    const AlulaPanelGuideControls *c,AlulaPanelGuideStatus *status)
{
    AlulaPanelGuideStatus result=ALULA_GUIDE_INVALID;
    AlulaPanelGuides *g=NULL;
    if (!field || !c || !isfinite(c->inlet_x) || !isfinite(c->outlet_x) ||
        !isfinite(c->inlet_lower_y) || !isfinite(c->inlet_upper_y) ||
        !(c->inlet_x<c->outlet_x && c->inlet_lower_y<c->inlet_upper_y) ||
        !(c->length_scale>0 && c->relative_tolerance>0 && c->seed_fraction>0) ||
        !isfinite(c->length_scale+c->relative_tolerance+c->seed_fraction)) goto fail;
    const AlulaPanelSystem *s=field->system;
    if (c->seed_limit) for (size_t e=0; e<s->bodies; ++e)
        if (!(c->seed_limit[e]>0) || !isfinite(c->seed_limit[e])) goto fail;
    g=calloc(1,sizeof *g);
    result=ALULA_GUIDE_ALLOCATION;
    if (!g) goto fail;
    g->edges=malloc(s->count*sizeof *g->edges);
    if (!g->edges) goto fail;
    g->field=field;
    g->outlet_x=c->outlet_x;
    g->scale=c->length_scale;
    g->tolerance=c->relative_tolerance*g->scale;
    g->psi_tolerance=fmin(1e-11*g->scale,g->tolerance);
    for (size_t e=0; e<s->bodies; ++e)
        g->bounds[e]=(Box){{INFINITY,INFINITY},{-INFINITY,-INFINITY}};
    for (size_t i=0; i<s->count; ++i) {
        const Segment *seg=s->segments+i;
        Box b=box(seg->geometry.a,seg->end);
        g->edges[i]=b;
        Box *bound=g->bounds+seg->body;
        bound->lo=(AlulaPoint){fmin(bound->lo.x,b.lo.x),fmin(bound->lo.y,b.lo.y)};
        bound->hi=(AlulaPoint){fmax(bound->hi.x,b.hi.x),fmax(bound->hi.y,b.hi.y)};
    }
    size_t first=0;
    for (size_t e=0; e<s->bodies; ++e) {
        size_t end=first;
        while (end<s->wetted && s->segments[end].body==e) ++end;
        if (end-first<3) { result=ALULA_GUIDE_INVALID; goto fail; }
        result=prepare_body(g,c,e,first,end);
        if (result!=ALULA_GUIDE_OK) goto fail;
        first=end;
    }
    qsort(g->bodies,s->bodies,sizeof *g->bodies,order);
    double budget=4*(c->outlet_x-c->inlet_x);
    for (int upper=0; upper<2; ++upper) {
        AlulaPoint p={c->inlet_x,upper?c->inlet_upper_y:c->inlet_lower_y};
        if (!alula_panel_field_sample(field,p,NULL,g->levels+(upper?s->bodies+1:0))) {
            result=ALULA_GUIDE_SEED; goto fail;
        }
        result=trace(g,p,c->outlet_x,budget,g->farfield+upper);
        if (result!=ALULA_GUIDE_OK) goto fail;
    }
    for (size_t e=0; e<s->bodies; ++e) g->levels[e+1]=g->bodies[e].inlet_level;
    for (size_t e=0; e<=s->bodies; ++e)
        if (!(g->levels[e+1]>g->levels[e])) { result=ALULA_GUIDE_ORDER; goto fail; }
    if (status) *status=ALULA_GUIDE_OK;
    return g;
fail:
    alula_panel_guides_free(g);
    if (status) *status=result;
    return NULL;
}

const AlulaPanelBodyGuide *alula_panel_guides_bodies(const AlulaPanelGuides *g,size_t *count)
{
    if (count) *count=g ? g->field->system->bodies : 0;
    return g ? g->bodies : NULL;
}

const double *alula_panel_guides_levels(const AlulaPanelGuides *g)
{
    return g ? g->levels : NULL;
}

AlulaGuidePath alula_panel_guides_farfield(const AlulaPanelGuides *g,int upper)
{
    if (!g || (upper!=0 && upper!=1)) return (AlulaGuidePath){0};
    const Path *p=g->farfield+upper;
    return (AlulaGuidePath){p->points,p->count};
}

static int directed(const AlulaPanelGuides *g, const AlulaTracePoint *a, const AlulaTracePoint *b)
{
    double orientation=(a->potential<b->potential?1:-1)*(a->position.x<b->position.x?1:-1);
    AlulaPoint qa,qb;
    return a->position.x!=b->position.x && a->potential!=b->potential &&
        alula_panel_field_sample(g->field,a->position,&qa,NULL) &&
        alula_panel_field_sample(g->field,b->position,&qb,NULL) &&
        orientation*qa.x>0 && orientation*qb.x>0;
}

AlulaPanelGuideStatus alula_panel_guides_sample_x(const AlulaPanelGuides *g,
    AlulaGuidePath path,double x,AlulaTracePoint *out)
{
    if (!g || !path.points || path.count<2 || !out || !isfinite(x)) return ALULA_GUIDE_INVALID;
    size_t lo=0,hi=0;
    int matches=0;
    const AlulaTracePoint *p=path.points;
    double direction=p[path.count-1].potential-p[0].potential;
    if (direction==0 || !isfinite(direction)) return ALULA_GUIDE_INVALID;
    for (size_t i=0;i<path.count;++i) {
        if (!isfinite(p[i].position.x) || !isfinite(p[i].position.y) || !isfinite(p[i].potential) ||
            (i && ((direction>0 && p[i].potential<=p[i-1].potential) ||
                (direction<0 && p[i].potential>=p[i-1].potential)))) return ALULA_GUIDE_INVALID;
        double b=p[i].position.x, a=i?p[i-1].position.x:b;
        if (x==b) {
            lo=hi=i;
            ++matches;
        } else if ((a<x && x<b) || (b<x && x<a)) {
            lo=i-1; hi=i;
            ++matches;
        }
        if (matches>1) return ALULA_GUIDE_AMBIGUOUS;
    }
    if (!matches) {
        /* A successfully integrated domain-plane endpoint carries a finite
         * event tolerance; preserve that endpoint when x matches its plane. */
        if (fabs(p[0].position.x-x)<=g->tolerance) {lo=hi=0;++matches;}
        size_t last=path.count-1;
        if (fabs(p[last].position.x-x)<=g->tolerance) {lo=hi=last;++matches;}
        if (!matches) return ALULA_GUIDE_UNRESOLVED;
        if (matches>1) return ALULA_GUIDE_AMBIGUOUS;
    }
    if (lo==hi) {
        if ((lo && !directed(g,p+lo-1,p+lo)) ||
            (hi+1<path.count && !directed(g,p+hi,p+hi+1))) return ALULA_GUIDE_UNRESOLVED;
        *out=p[lo];
        return ALULA_GUIDE_OK;
    }
    if (!directed(g,p+lo,p+hi)) return ALULA_GUIDE_UNRESOLVED;
    TraceContext context={g,NULL};
    AlulaTraceControls controls={
        .parameter=ALULA_TRACE_POTENTIAL,.end=p[hi].potential,.tolerance=g->tolerance,
        .max_step=.05*g->scale,.min_step=1e-12*g->scale,.max_spatial_step=.1*g->scale,
        .max_steps=20000,.stop_at_plane=1,.plane_normal={1,0},.plane_offset=x,
        .context=&context,.admissible_point=trace_point,.admissible_segment=trace_chord
    };
    AlulaTraceResult r=alula_trace_panel(g->field,p[lo],&controls);
    if (r.status!=ALULA_TRACE_PLANE) {
        if (r.status!=ALULA_TRACE_PLANE_NOT_REACHED || fabs(r.last.position.x-x)>g->tolerance)
            return ALULA_GUIDE_UNRESOLVED;
    }
    *out=r.last;
    return ALULA_GUIDE_OK;
}

static int inlet_seed(const AlulaPanelGuides *g,double level,AlulaPoint *p)
{
    size_t top=g->field->system->bodies+1;
    AlulaPoint low=g->farfield[0].points[0].position, high=g->farfield[1].points[0].position;
    double lower=low.y, upper=high.y;
    p->x=low.x;
    p->y=lower+(upper-lower)*(level-g->levels[0])/(g->levels[top]-g->levels[0]);
    AlulaPoint q;
    if (!alula_panel_field_sample(g->field,low,&q,NULL) || !(q.x>0) ||
        !alula_panel_field_sample(g->field,high,&q,NULL) || !(q.x>0)) return 0;
    for (int step=0;step<60;++step) {
        double psi;
        if (!alula_panel_field_sample(g->field,*p,&q,&psi) || !(q.x>0)) return 0;
        double error=psi-level;
        if (fabs(error)<=g->psi_tolerance) return 1;
        if (error<0) lower=p->y;
        else upper=p->y;
        double next=p->y-error/q.x;
        p->y=next>lower && next<upper ? next : lower+.5*(upper-lower);
    }
    return 0;
}

AlulaPanelStreamline *alula_panel_guides_trace_level(const AlulaPanelGuides *g,
    double level,AlulaPanelGuideStatus *status)
{
    AlulaPanelGuideStatus code=ALULA_GUIDE_INVALID;
    Path *line=NULL;
    if (!g || !(level>g->levels[0] && level<g->levels[g->field->system->bodies+1])) goto done;
    for (size_t i=0;i<g->field->system->bodies;++i) if (level==g->levels[i+1]) goto done;
    AlulaPoint seed;
    code=ALULA_GUIDE_SEED;
    if (!inlet_seed(g,level,&seed)) goto done;
    line=calloc(1,sizeof *line);
    code=ALULA_GUIDE_ALLOCATION;
    if (!line) goto done;
    code=trace(g,seed,g->outlet_x,4*(g->outlet_x-seed.x),line);
    if (code!=ALULA_GUIDE_OK) {
        alula_panel_streamline_free(line);
        line=NULL;
    }
done:
    if (status) *status=code;
    return line;
}

AlulaGuidePath alula_panel_streamline_path(const AlulaPanelStreamline *line)
{
    return line ? (AlulaGuidePath){line->points,line->count} : (AlulaGuidePath){0};
}

void alula_panel_streamline_free(AlulaPanelStreamline *line)
{
    if (!line) return;
    free(line->points);
    free(line);
}

void alula_panel_guides_free(AlulaPanelGuides *g)
{
    if (!g) return;
    for (int e=0; e<6; ++e) {
        alula_contour_spline_free(g->curves+e);
        free(g->incoming[e].points);
        free(g->wake[e].points);
    }
    free(g->farfield[0].points);
    free(g->farfield[1].points);
    free(g->edges);
    free(g);
}
