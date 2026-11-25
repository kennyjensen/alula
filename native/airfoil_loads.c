// SPDX-License-Identifier: GPL-2.0-or-later
#include "airfoil_grid_internal.h"
#include "body_internal.h"

#include <stdlib.h>

static double wall_cp(const AlulaBody *b,int element,int station,int upper)
{
    int passage=element+upper,tube=upper?0:b->tubes[passage]-1;
    double pressure=cell(b,passage,station,tube)->value.pressure[!upper];
    return 2*(pressure-b->p_inf);
}

static void integrate(const AlulaBody *b,AlulaAirfoilPressureLoads *loads,
    AlulaPoint a,AlulaPoint z,double cp_a,double cp_z)
{
    double chord=b->in.length_scale,average=.5*(cp_a+cp_z);
    AlulaPoint d=scale(sub(z,a),1/chord),r=scale(sub(mean(a,z),b->in.center),1/chord);
    loads->cx-=average*d.y;
    loads->cy+=average*d.x;
    loads->cm-=average*dot(d,r)+(cp_z-cp_a)*dot(d,d)/12;
}

static int integrate_surface(AlulaAirfoilGrid *g,int e,int upper,double leading_cp,
    AlulaAirfoilPressureLoads *loads)
{
    const AlulaBody *b=g->body;
    const AlulaBodyCut *cut=b->cuts+e;
    AlulaBodySurfaceView view;
    if (!alula_body_coordinates_surface(g->coordinates,g->packed,e,upper,&view)) return 0;
    const AlulaContourSpline *curve=view.curve;
    double arc=curve->arc[curve->count-1];
    AlulaPoint previous;
    alula_contour_spline_value(curve,view.stagnation,&previous,NULL,NULL);
    double previous_cp=leading_cp;
    for (int i=cut->leading+1;i<=cut->trailing;++i) {
        double f=view.fractions[i-cut->leading];
        double location=upper?(1-f)*view.stagnation:view.stagnation+f*(arc-view.stagnation);
        AlulaPoint next;
        alula_contour_spline_value(curve,location,&next,NULL,NULL);
        double cp=wall_cp(b,e,i,upper);
        if (upper) integrate(b,loads,next,previous,cp,previous_cp);
        else integrate(b,loads,previous,next,previous_cp,cp);
        previous=next;
        previous_cp=cp;
    }
    return 1;
}

static void integrate_base(const AlulaBody *b,const AlulaMaterialBase *base,int e,
    AlulaAirfoilPressureLoads *loads)
{
    int te=b->cuts[e].trailing;
    double cp_lower=wall_cp(b,e,te,0),cp_upper=wall_cp(b,e,te,1);
    double length=0;
    for (size_t i=1;i<base->count;++i) length+=norm(sub(base->points[i],base->points[i-1]));
    double distance=0,previous_cp=cp_lower;
    for (size_t i=1;i<base->count;++i) {
        distance+=norm(sub(base->points[i],base->points[i-1]));
        double cp=i+1==base->count?cp_upper:cp_lower+distance/length*(cp_upper-cp_lower);
        integrate(b,loads,base->points[i-1],base->points[i],previous_cp,cp);
        previous_cp=cp;
    }
}

int alula_airfoil_grid_pressure_loads(AlulaAirfoilGrid *g,AlulaAirfoilPressureLoads *out)
{
    if (!g || !g->converged || !out) return 0;
    AlulaBody *b=g->body;
    AlulaBodyState state;
    double *residual=malloc((size_t)b->size.rows*sizeof *residual);
    if (!residual) return 0;
    int valid=alula_airfoil_grid_state(g,&state) &&
        alula_body_residual(b,&state,residual,NULL)==ALULA_BODY_OK;
    free(residual);
    if (!valid) return 0;
    AlulaAirfoilPressureLoads loads={0};
    for (int e=0;e<b->in.bodies;++e) {
        int le=b->cuts[e].leading;
        double leading_cp=.5*(wall_cp(b,e,le,1)+wall_cp(b,e,le,0));
        if (!integrate_surface(g,e,1,leading_cp,&loads) || !integrate_surface(g,e,0,leading_cp,&loads))
            return 0;
        if (g->bases[e].count) integrate_base(b,g->bases+e,e,&loads);
    }
    loads.cl=loads.cy*b->freestream.x-loads.cx*b->freestream.y;
    loads.pressure_drag=loads.cx*b->freestream.x+loads.cy*b->freestream.y;
    if (!isfinite(loads.cl) || !isfinite(loads.pressure_drag) || !isfinite(loads.cm)) return 0;
    *out=loads;
    return 1;
}
