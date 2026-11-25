// SPDX-License-Identifier: GPL-2.0-or-later
#include "airfoil_bl_internal.h"
#include "airfoil_grid_internal.h"

#include <stdlib.h>

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
    return alula_bl_dead_air(width,cross(a,b)/(norm(a)*norm(b)),distance);
}

static int sample_body(AlulaAirfoilBL *p,AlulaAirfoilGrid *g,const AlulaBodyState *state,
    int e,const double *trips,int probe)
{
    const AlulaBody *b=g->body;
    AlulaBodyCut cut=b->cuts[e];
    size_t ns=p->branches[e][0].count,nw=p->branches[e][ALULA_BL_MERGED_WAKE].count;
    AlulaEdgeCursor cursor[2];
    cursor[0].station=cursor[1].station=0;
    AlulaBodySurfaceView s;
    if (!alula_body_coordinates_surface(g->coordinates,g->packed,e,1,&s)) return 0;
    const double *arc=s.arc_lengths;
    if (!arc) return 0;
    double stagnation=alula_contour_arc(s.curve,arc,s.stagnation),scale=b->in.length_scale;
    for (int side=0;side<2;++side) {
        int upper=side==ALULA_BL_UPPER;
        if (!alula_body_coordinates_surface(g->coordinates,g->packed,e,upper,&s)) return 0;
        int passage=e+upper,index=upper?0:b->tubes[e]-1;
        double sign=upper?-1:1;
        AlulaBLEdge *edges=p->edges+(p->branches[e][side].stations-p->stations);
        for (size_t i=0;i<ns;++i) {
            double distance=parameter(s.curve,s.stagnation,s.fractions[i+1],upper);
            edges[i]=(AlulaBLEdge){sign*(alula_contour_arc(s.curve,arc,distance)-stagnation)/scale,
                alula_body_edge_speed(b,probe?state:NULL,passage,cut.leading+1+(int)i,index,&cursor[side]),0};
        }
        double fraction=trips?trips[side]:1;
        double trip_parameter=parameter(s.curve,s.initial_stagnation,fraction,upper);
        p->trip[e][side]=fraction<1?sign*(alula_contour_arc(s.curve,arc,trip_parameter)-stagnation)/scale:edges[ns-1].s;
    }
    const AlulaAirfoilBLBranch *lower=&p->branches[e][ALULA_BL_LOWER];
    AlulaBLEdge *edges=p->edges+(p->branches[e][ALULA_BL_MERGED_WAKE].stations-p->stations);
    double distance=p->edges[(lower->stations-p->stations)+ns-1].s,within=0;
    double width=g->base_width[e];
    AlulaPoint last={0};
    for (size_t k=0;k<nw;++k) {
        int i=cut.trailing+(int)k;
        AlulaPoint center=mean(point(b,state,e,i,b->tubes[e]),point(b,state,e+1,i,0));
        if (k) within+=norm(sub(center,last));
        last=center;
        edges[k]=(AlulaBLEdge){distance+within/scale,
            .5*(alula_body_edge_speed(b,probe?state:NULL,e,i,b->tubes[e]-1,&cursor[ALULA_BL_LOWER])+
                alula_body_edge_speed(b,probe?state:NULL,e+1,i,0,&cursor[ALULA_BL_UPPER])),wake_width(&s,width,within)/scale};
    }
    return 1;
}

int alula_airfoil_bl_resample(AlulaAirfoilBL *p,AlulaAirfoilGrid *g,const double (*trips)[2],int probe)
{
    AlulaBodyState state;
    if (probe) {
        if (!alula_body_coordinates_edge_state(g->coordinates,g->packed,g->nodes,g->physical,&state)) return 0;
    } else if (!alula_airfoil_grid_prepare_flow(g,&state)) return 0;
    for (int e=0;e<p->bodies;++e)
        if (!sample_body(p,g,&state,e,trips?trips[e]:NULL,probe)) return 0;
    return 1;
}

AlulaAirfoilBL *alula_airfoil_bl_sample(AlulaAirfoilGrid *g,double reynolds,double ncrit,
    const double (*trips)[2],int probe)
{
    AlulaBLContext ctx;
    if (!g || (!probe && !g->converged) || !alula_bl_context(&ctx,reynolds,g->body->in.mach,g->body->in.gamma,ncrit))
        return NULL;
    const AlulaBody *b=g->body;
    size_t count=0;
    for (int i=0;i<b->in.bodies;++i) {
        if (trips) for (int s=0;s<2;++s) if (!(trips[i][s]>0 && trips[i][s]<=1)) return NULL;
        count+=(size_t)(2*(b->cuts[i].trailing-b->cuts[i].leading)+b->in.nx-b->cuts[i].trailing+1);
    }
    AlulaAirfoilBL *p=calloc(1,sizeof *p);
    if (!p) return NULL;
    p->stations=malloc(count*sizeof *p->stations);
    p->regimes=malloc(count*sizeof *p->regimes);
    p->edges=malloc(count*sizeof *p->edges);
    if (!p->stations || !p->regimes || !p->edges) goto fail;
    p->bodies=b->in.bodies;
    p->conditions=ctx;
    size_t start=0;
    for (int e=0;e<p->bodies;++e) for (int side=0;side<3;++side) {
        AlulaBodyCut cut=b->cuts[e];
        size_t n=(size_t)(side==ALULA_BL_MERGED_WAKE?b->in.nx-cut.trailing+1:cut.trailing-cut.leading);
        p->branches[e][side]=(AlulaAirfoilBLBranch){p->stations+start,p->regimes+start,n,
            {.status=ALULA_BL_INVALID}};
        start+=n;
    }
    if (alula_airfoil_bl_resample(p,g,trips,probe)) return p;
fail:
    alula_airfoil_bl_free(p);
    return NULL;
}

AlulaAirfoilBL *alula_airfoil_bl_create_interaction(AlulaAirfoilGrid *g,double reynolds,double ncrit,
    const double (*trips)[2],const double *displacement,int local)
{
    AlulaAirfoilBL *p=alula_airfoil_bl_sample(g,reynolds,ncrit,trips,0);
    if (!p) return NULL;
    double *interaction=displacement?malloc((size_t)g->body->in.nx*sizeof *interaction):NULL;
    if (displacement && !interaction) {alula_airfoil_bl_free(p);return NULL;}
    p->complete=1;
    for (int e=0;e<p->bodies;++e) {
        AlulaAirfoilBLBranch *branches=p->branches[e];
        size_t ns=branches[0].count;
        for (int side=0;side<2;++side) {
            size_t start=(size_t)(branches[side].stations-p->stations);
            if (displacement) for (size_t j=0;j<ns;++j)
                interaction[j]=displacement[(side==ALULA_BL_UPPER?ns+1:0)+1+j]/g->body->in.length_scale;
            branches[side].result=alula_bl_march_surface_controlled(&p->conditions,p->edges+start,ns,
                p->trip[e][side],1,interaction,local,p->stations+start,p->regimes+start);
        }
        AlulaAirfoilBLBranch *wake=branches+ALULA_BL_MERGED_WAKE;
        size_t start=(size_t)(wake->stations-p->stations);
        if (branches[0].result.completed==ns && branches[1].result.completed==ns) {
            if (local && displacement) for (size_t j=1;j<wake->count;++j)
                interaction[j-1]=displacement[2*(ns+1)+j-1]/g->body->in.length_scale;
            wake->result=alula_bl_march_wake_interaction(&p->conditions,branches[0].stations+ns-1,
                branches[1].stations+ns-1,p->edges+start,wake->count,local?interaction:NULL,p->stations+start);
            for (size_t j=0;j<wake->result.completed;++j) p->regimes[start+j]=ALULA_BL_WAKE;
        }
        if (displacement) displacement+=2*(ns+1)+wake->count-1;
        for (int side=0;side<3;++side) if (branches[side].result.status!=ALULA_BL_OK) p->complete=0;
    }
    free(interaction);
    return p;
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
    free(p->edges);free(p->stations);free(p->regimes);free(p);
}
