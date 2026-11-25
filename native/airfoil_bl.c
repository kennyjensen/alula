// SPDX-License-Identifier: GPL-2.0-or-later
#include "airfoil_bl.h"
#include "airfoil_grid_internal.h"
#include "body_internal.h"

#include <stdlib.h>

struct AlulaAirfoilBL {
    AlulaAirfoilBLBranch branches[MAX_BODIES][3];
    AlulaBLStation *stations;
    AlulaBLRegime *regimes;
    AlulaBLContext conditions;
    int bodies,complete;
};

static double edge_speed(const AlulaBody *b,int g,int i,int j)
{
    int station=i<b->in.nx?i:b->in.nx-1;
    const AlulaCellValue *v=&cell(b,g,station,j)->value;
    double q=.5*(v->state[0].speed+v->state[1].speed);
    double mach=.5*(v->state[0].mach_sq+v->state[1].mach_sq);
    return mach<1?q+.2*q*mach*(mach-1)*v->geometry.curvature:q;
}

static double parameter(const AlulaContourSpline *curve,double stagnation,double fraction,int upper)
{
    return upper?(1-fraction)*stagnation:stagnation+fraction*(curve->arc[curve->count-1]-stagnation);
}

static double wake_width(const AlulaBodySurfaceView *surface,double width,double distance)
{
    if (width<=0 || distance>=2.5*width) return 0;
    const AlulaContourSpline *c=surface->curve;
    AlulaPoint a,b;
    alula_contour_spline_value(c,0,NULL,&a,NULL);
    alula_contour_spline_value(c,c->arc[c->count-1],NULL,&b,NULL);
    double tangent=cross(a,b)/(norm(a)*norm(b)),limit=1.2;
    double slope=fabs(tangent)>=limit/hypot(1,limit)?copysign(limit,tangent):tangent/sqrt(1-tangent*tangent);
    double z=1-distance/(2.5*width),shape=3+2.5*slope;
    return width*(shape+(1-shape)*z)*z*z;
}

static void march_body(AlulaAirfoilBL *p,AlulaAirfoilGrid *g,const AlulaBLContext *ctx,
    const AlulaBodyState *state,int e,const double *trips,AlulaBLEdge *edges,size_t *offset,
    const double *displacement,double *interaction,int local)
{
    const AlulaBody *b=g->body;
    AlulaBodyCut cut=b->cuts[e];
    size_t ns=(size_t)(cut.trailing-cut.leading),nw=(size_t)(b->in.nx-cut.trailing+1);
    for (int side=0;side<3;++side) {
        size_t count=side==ALULA_BL_MERGED_WAKE?nw:ns;
        p->branches[e][side]=(AlulaAirfoilBLBranch){p->stations+*offset,p->regimes+*offset,count,
            {.status=ALULA_BL_INVALID}};
        *offset+=count;
    }
    AlulaBodySurfaceView s;
    if (!alula_body_coordinates_surface(g->coordinates,g->packed,e,1,&s)) return;
    double arc[ALULA_PANEL_MAX_WETTED+1];
    alula_contour_arc_build(s.curve,arc);
    double stagnation=alula_contour_arc(s.curve,arc,s.stagnation),scale=b->in.length_scale;
    for (int side=0;side<2;++side) {
        int upper=side==ALULA_BL_UPPER;
        if (!alula_body_coordinates_surface(g->coordinates,g->packed,e,upper,&s)) return;
        int passage=e+upper,index=upper?0:b->tubes[e]-1;
        double sign=upper?-1:1;
        for (size_t i=0;i<ns;++i) {
            double distance=parameter(s.curve,s.stagnation,s.fractions[i+1],upper);
            edges[i]=(AlulaBLEdge){sign*(alula_contour_arc(s.curve,arc,distance)-stagnation)/scale,
                edge_speed(b,passage,cut.leading+1+(int)i,index),0};
        }
        double fraction=trips?trips[side]:1;
        double trip_parameter=parameter(s.curve,s.initial_stagnation,fraction,upper);
        double trip=fraction<1?sign*(alula_contour_arc(s.curve,arc,trip_parameter)-stagnation)/scale:edges[ns-1].s;
        AlulaAirfoilBLBranch *surface=&p->branches[e][side];
        size_t start=(size_t)(surface->stations-p->stations);
        if (displacement) for (size_t i=0;i<ns;++i)
            interaction[i]=displacement[(upper?ns+1:0)+1+i]/scale;
        surface->result=alula_bl_march_surface_controlled(ctx,edges,ns,trip,1,
            displacement?interaction:NULL,local,p->stations+start,p->regimes+start);
    }
    const AlulaAirfoilBLBranch *upper=&p->branches[e][ALULA_BL_UPPER];
    const AlulaAirfoilBLBranch *lower=&p->branches[e][ALULA_BL_LOWER];
    if (upper->result.completed!=ns || lower->result.completed!=ns) return;
    AlulaAirfoilBLBranch *wake=&p->branches[e][ALULA_BL_MERGED_WAKE];
    double distance=lower->stations[ns-1].s,within=0;
    double width=g->base_width[e];
    AlulaPoint last={0};
    for (size_t k=0;k<nw;++k) {
        int i=cut.trailing+(int)k;
        AlulaPoint center=mean(point(b,state,e,i,b->tubes[e]),point(b,state,e+1,i,0));
        if (k) within+=norm(sub(center,last));
        last=center;
        edges[k]=(AlulaBLEdge){distance+within/scale,
            .5*(edge_speed(b,e,i,b->tubes[e]-1)+edge_speed(b,e+1,i,0)),wake_width(&s,width,within)/scale};
    }
    size_t start=(size_t)(wake->stations-p->stations);
    if (local && displacement) for (size_t j=1;j<nw;++j)
        interaction[j-1]=displacement[2*(ns+1)+j-1]/scale;
    wake->result=alula_bl_march_wake_interaction(ctx,upper->stations+ns-1,lower->stations+ns-1,
        edges,nw,local && displacement?interaction:NULL,p->stations+start);
    for (size_t k=0;k<wake->result.completed;++k) p->regimes[start+k]=ALULA_BL_WAKE;
}

AlulaAirfoilBL *alula_airfoil_bl_create_interaction(AlulaAirfoilGrid *g,double reynolds,double ncrit,
    const double (*trips)[2],const double *displacement,int local)
{
    AlulaBLContext ctx;
    if (!g || !g->converged || !alula_bl_context(&ctx,reynolds,g->body->in.mach,g->body->in.gamma,ncrit))
        return NULL;
    AlulaBody *b=g->body;
    size_t count=0;
    for (int i=0;i<b->in.bodies;++i) {
        if (trips) for (int s=0;s<2;++s) if (!(trips[i][s]>0 && trips[i][s]<=1)) return NULL;
        count+=(size_t)(2*(b->cuts[i].trailing-b->cuts[i].leading)+b->in.nx-b->cuts[i].trailing+1);
    }
    AlulaAirfoilBL *p=calloc(1,sizeof *p);
    AlulaBLEdge *edges=malloc((size_t)b->in.nx*sizeof *edges);
    double *interaction=displacement?malloc((size_t)b->in.nx*sizeof *interaction):NULL;
    double *residual=malloc((size_t)b->size.rows*sizeof *residual);
    AlulaBodyState state;
    if (!p || !edges || !residual || (displacement && !interaction) || !alula_airfoil_grid_state(g,&state) ||
        alula_body_residual(b,&state,residual,NULL)!=ALULA_BODY_OK) goto fail;
    p->stations=malloc(count*sizeof *p->stations);
    p->regimes=malloc(count*sizeof *p->regimes);
    if (!p->stations || !p->regimes) goto fail;
    p->bodies=b->in.bodies;
    p->conditions=ctx;
    p->complete=1;
    size_t start=0;
    for (int i=0;i<b->in.bodies;++i) {
        march_body(p,g,&ctx,&state,i,trips?trips[i]:NULL,edges,&start,displacement,interaction,local);
        if (displacement) {
            AlulaBodyCut cut=b->cuts[i];
            displacement+=2*(cut.trailing-cut.leading+1)+b->in.nx-cut.trailing;
        }
        for (int side=0;side<3;++side) if (p->branches[i][side].result.status!=ALULA_BL_OK) p->complete=0;
    }
    free(edges);free(residual);free(interaction);
    return p;
fail:
    free(edges);free(residual);free(interaction);alula_airfoil_bl_free(p);
    return NULL;
}

AlulaAirfoilBL *alula_airfoil_bl_create(AlulaAirfoilGrid *g,double reynolds,double ncrit,
    const double (*trips)[2])
{
    return alula_airfoil_bl_create_interaction(g,reynolds,ncrit,trips,NULL,0);
}

const AlulaAirfoilBLBranch *alula_airfoil_bl_branch(const AlulaAirfoilBL *p,int body,AlulaBLBranch side)
{
    return p && body>=0 && body<p->bodies && side>=ALULA_BL_UPPER && side<=ALULA_BL_MERGED_WAKE?
        &p->branches[body][side]:NULL;
}

int alula_airfoil_bl_complete(const AlulaAirfoilBL *p) {return p && p->complete;}

int alula_airfoil_bl_wake_drag(const AlulaAirfoilBL *p,double *coefficient)
{
    if (!p || !p->complete || !coefficient) return 0;
    const AlulaBLContext *c=&p->conditions;
    double sum=0;
    for (int i=0;i<p->bodies;++i) {
        const AlulaAirfoilBLBranch *wake=&p->branches[i][ALULA_BL_MERGED_WAKE];
        const AlulaBLStation *s=wake->stations+wake->count-1;
        AlulaBLProperties properties;
        if (alula_bl_station(c,s,ALULA_BL_WAKE,&properties)!=ALULA_BL_OK) return 0;
        double h_inf=1+c->gm1*c->mach*c->mach;
        double h=(s->delta_star-s->wake_gap)/s->theta;
        double drag=2*properties.rho*s->theta*pow(s->ue,2+.5*(h+h_inf));
        if (!positive(drag)) return 0;
        sum+=drag;
    }
    *coefficient=sum;
    return isfinite(sum);
}

void alula_airfoil_bl_free(AlulaAirfoilBL *p)
{
    if (!p) return;
    free(p->stations);free(p->regimes);free(p);
}
