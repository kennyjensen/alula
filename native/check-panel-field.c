// SPDX-License-Identifier: GPL-2.0-or-later
#include "panel_internal.h"
#include "panel_guide.h"
#include "trace.h"

#include <stdio.h>

typedef struct { double x_limit; int visited; } TraceCheck;

static int chord(void *context, AlulaPoint start, AlulaPoint end)
{
    TraceCheck *test = context;
    return start.x < test->x_limit && end.x < test->x_limit;
}

static void visit(void *context, const AlulaTracePoint *point)
{
    TraceCheck *test = context;
    if (isfinite(point->parameter)) ++test->visited;
}

int main(void)
{
    AlulaPoint points[ALULA_PANEL_MAX_TOTAL+6];
    AlulaContour contours[6];
    size_t elements, used = 0;
    double alpha;
    int test_guides;
    if (scanf("%zu%lf%d", &elements, &alpha, &test_guides) != 3 || elements < 1 || elements > 6) return 2;
    for (size_t i = 0; i < elements; ++i) {
        size_t count, upper, lower;
        int finite;
        if (scanf("%zu%d%zu%zu", &count, &finite, &upper, &lower) != 4 ||
            count < 4 || count > ALULA_PANEL_MAX_TOTAL+6-used) return 2;
        for (size_t j = 0; j < count; ++j)
            if (scanf("%lf%lf", &points[used+j].x, &points[used+j].y) != 2) return 2;
        contours[i] = (AlulaContour){points+used, count, finite, upper, lower};
        used += count;
    }
    AlulaPanelSystem *system = alula_panel_prepare(contours, elements, 1, (AlulaPoint){.25, 0});
    AlulaPanelField *field = alula_panel_field_create(system, alpha);
    int status = 1;
    if (field && test_guides) {
        AlulaPanelGuideStatus reason;
        AlulaPanelGuideControls controls = {-3, 4, -3, 3, 1, 2e-9, .002, NULL};
        AlulaPanelGuides *guides = alula_panel_guides_create(field,&controls,&reason);
        printf("%d\n", (int)reason);
        if (guides) {
            size_t count;
            const AlulaPanelBodyGuide *bodies = alula_panel_guides_bodies(guides,&count);
            const double *levels = alula_panel_guides_levels(guides);
            for (size_t i=0; i<count+2; ++i) printf("%.17g%c", levels[i], i==count+1?'\n':' ');
            for (size_t i=0; i<count; ++i) {
                const AlulaPanelBodyGuide *b = bodies+i;
                const AlulaTracePoint *in = b->incoming.points, *out = b->wake.points;
                AlulaTracePoint up_mid,wake_mid;
                if (alula_panel_guides_sample_x(guides,b->incoming,-.6,&up_mid)!=ALULA_GUIDE_OK ||
                    alula_panel_guides_sample_x(guides,b->wake,1.6,&wake_mid)!=ALULA_GUIDE_OK ||
                    fabs(up_mid.position.x+.6)>2e-9 || fabs(wake_mid.position.x-1.6)>2e-9) {
                    alula_panel_guides_free(guides);
                    goto done;
                }
                printf("%zu %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                    b->element,b->stagnation,b->strain,b->wall_level,in[0].position.x,in[0].position.y,
                    in[b->incoming.count-1].position.y,out[0].position.x,out[0].position.y,
                    out[b->wake.count-1].position.x,out[b->wake.count-1].position.y,
                    up_mid.position.y,wake_mid.position.y);
            }
            AlulaTracePoint folded[] = {{{-3,-3},0,0},{{-2,-3},1,1},{{-3,-3},2,2}}, point;
            if (alula_panel_guides_sample_x(guides,(AlulaGuidePath){folded,3},-2.5,&point)!=ALULA_GUIDE_AMBIGUOUS) {
                alula_panel_guides_free(guides);
                goto done;
            }
            for (int upper=0; upper<2; ++upper) {
                AlulaGuidePath path = alula_panel_guides_farfield(guides,upper);
                if (path.count<2 || fabs(path.points[path.count-1].position.x-controls.outlet_x)>2e-9) {
                    alula_panel_guides_free(guides);
                    goto done;
                }
            }
            for (size_t i=0; i<=count; ++i) {
                double level=.5*(levels[i]+levels[i+1]);
                AlulaPanelStreamline *line=alula_panel_guides_trace_level(guides,level,&reason);
                AlulaGuidePath path=alula_panel_streamline_path(line);
                if (!line || path.count<2 || fabs(path.points[path.count-1].position.x-controls.outlet_x)>2e-9) {
                    fprintf(stderr,"Passage %zu trace status %d\n",i,(int)reason);
                    alula_panel_streamline_free(line);
                    alula_panel_guides_free(guides);
                    goto done;
                }
                printf("%.17g %.17g %.17g\n",level,path.points[0].position.y,path.points[path.count-1].position.y);
                alula_panel_streamline_free(line);
            }
            alula_panel_guides_free(guides);
        }
        status=0;
        goto done;
    }
    if (field) {
        for (size_t i = system->wetted; i < system->count; ++i) {
            const Segment *base = system->segments+i;
            if (base->source == 0) continue;
            AlulaPoint p = {base->cut_origin.x + .6*base->cut_direction.x,
                base->cut_origin.y + .6*base->cut_direction.y};
            double psi;
            if (alula_panel_field_sample(field, p, NULL, &psi)) goto done;
            break;
        }
        size_t probes;
        if (scanf("%zu", &probes) != 1) goto done;
        for (size_t i = 0; i < probes; ++i) {
            AlulaPoint a, b;
            if (scanf("%lf%lf%lf%lf", &a.x, &a.y, &b.x, &b.y) != 4) goto done;
            AlulaPoint velocity = {0};
            double psi = 0;
            int has_velocity = alula_panel_field_sample(field, b, &velocity, NULL);
            int has_psi = alula_panel_field_sample(field, b, NULL, &psi);
            if (has_velocity && has_psi) {
                AlulaPoint both;
                double both_psi;
                if (!alula_panel_field_sample(field, b, &both, &both_psi) ||
                    both.x != velocity.x || both.y != velocity.y || both_psi != psi)
                    goto done;
            }
            printf("%d %d %.17g %.17g %.17g %.17g\n", has_velocity, has_psi,
                velocity.x, velocity.y, psi, alula_panel_field_branch_increment(field, a, b));
        }
        AlulaTracePoint seed;
        int by_arc;
        AlulaTraceControls controls = {
            .tolerance=1e-9, .potential_tolerance=1e-9, .max_step=.05,
            .min_step=1e-12, .max_spatial_step=.08, .max_steps=20000,
            .plane_normal={1, 0}, .admissible_segment=chord, .visit=visit
        };
        TraceCheck trace;
        controls.context = &trace;
        while (scanf("%d%lf%lf%lf%lf%lf%d%lf%lf", &by_arc, &seed.position.x,
            &seed.position.y, &seed.parameter, &seed.potential, &controls.end,
            &controls.stop_at_plane, &controls.plane_offset, &trace.x_limit) == 9) {
            controls.parameter = by_arc ? ALULA_TRACE_ARC : ALULA_TRACE_POTENTIAL;
            trace.visited = 0;
            AlulaTraceResult r = alula_trace_panel(field, seed, &controls);
            if (trace.visited != r.accepted+1) goto done;
            printf("%d %d %d %.17g %.17g %.17g %.17g\n", (int)r.status, r.accepted, r.rejected,
                r.last.position.x, r.last.position.y, r.last.parameter, r.last.potential);
        }
        status = 0;
    }
done:
    alula_panel_field_free(field);
    alula_panel_free(system);
    return status;
}
