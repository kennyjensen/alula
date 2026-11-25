// SPDX-License-Identifier: GPL-2.0-or-later
// Material trip geometry follows LEFIND, SEGSPL and XIFSET from XFOIL,
// Copyright (C) 2000 Mark Drela, via the Vibefoil JS port by Kenny Jensen.
#include "panel_bl.h"
#include "bl_internal.h"
#include "contour.h"
#include "panel_internal.h"
#include <stdlib.h>
#include <string.h>

static int trip_coordinates(const AlulaContourSpline *c,const double trips[2],double arc[2])
{
    const int end=c->count-1;
    const AlulaPoint *p=c->points;
    AlulaPoint te={.5*(p[0].x+p[end].x),.5*(p[0].y+p[end].y)};
    double total=c->arc[end];
    arc[0]=0;arc[1]=total;
    if (trips[0]==1 && trips[1]==1) return 1;
    int nose=1;
    for (int j=2;j<end;++j) {
        if (dot(subtract(p[j],te),subtract(p[j],te))>dot(subtract(p[nose],te),subtract(p[nose],te))) nose=j;
    }
    double sle=c->arc[nose];
    AlulaPoint le,t,bend;
    for (int k=0;k<50;++k) {
        alula_contour_spline_value(c,sle,&le,&t,&bend);
        AlulaPoint d=subtract(le,te);
        double ds=-dot(d,t)/(dot(t,t)+dot(d,bend));
        if (!isfinite(ds)) return 0;
        double bound=.02*hypot(d.x,d.y);
        ds=fmax(-bound,fmin(bound,ds));
        sle+=ds;
        if (fabs(ds)<1e-12*total) break;
        if (!(sle>0 && sle<total) || k==49) return 0;
    }
    alula_contour_spline_value(c,sle,&le,NULL,NULL);
    AlulaPoint chord=subtract(te,le);
    /* The chord projection is linear in the two SEGSPL cubics. Bracket
     * each surface root so an overhanging nose cannot exchange trips. */
    for (int side=0;side<2;++side) if (trips[side]!=1) {
        double lo=side?sle:0,hi=side?total:sle,target=trips[side]*dot(chord,chord);
        for (int step=0;step<55;++step) {
            double mid=.5*(lo+hi);
            AlulaPoint point;
            alula_contour_spline_value(c,mid,&point,NULL,NULL);
            if ((dot(subtract(point,le),chord)<target)==side) lo=mid;
            else hi=mid;
        }
        arc[side]=.5*(lo+hi);
    }
    return 1;
}

static AlulaBLStatus seed_body(AlulaPanelBL *p,const AlulaDisplacementData *d,size_t e,
    const double trips[2],AlulaBLEdge *edge,AlulaBLStation *march,AlulaBLRegime *regimes,double *inverse)
{
    size_t start=d->surface_start[e],end=d->surface_start[e+1]-1;
    size_t first_wake=d->surface_stations+e*d->wake_count;
    const double *q=d->inviscid,*s=d->arc;
    size_t left=start;
    while (left<end && q[left+1]<0) ++left;
    if (left==start || left+1>=end) return ALULA_BL_DOMAIN;
    size_t right=left+1+(q[left+1]==0);
    if (right>=end || !(q[left]<0 && q[right]>0)) return ALULA_BL_DOMAIN;
    p->upper_first[e]=left;p->lower_first[e]=right;
    double da=s[right]-s[left],du=q[right]-q[left];
    double upper_distance=-q[left]/du*da,lower_distance=q[right]/du*da;
    double stag=right==left+1?s[left]+upper_distance:s[left+1];
    p->stagnation[e]=stag;
    AlulaContourSpline curve;
    if (!alula_contour_spline_parabolic(&curve,d->points+start,(int)(end-start+1))) return ALULA_BL_INVALID;
    int ok=trip_coordinates(&curve,trips,p->trip_arc[e]);
    AlulaPoint top,bottom;
    alula_contour_spline_value(&curve,0,NULL,&top,NULL);
    alula_contour_spline_value(&curve,curve.arc[curve.count-1],NULL,&bottom,NULL);
    alula_contour_spline_free(&curve);
    if (!ok) return ALULA_BL_DOMAIN;
    for (int side=0;side<2;++side) {
        size_t count=side?end-right+1:left-start+1;
        double sign=side?1:-1;
        for (size_t k=0;k<count;++k) {
            size_t id=side?right+k:left-k;
            double distance=sign*(s[id]-stag);
            if (!k && right==left+1) distance=side?lower_distance:upper_distance;
            edge[k]=(AlulaBLEdge){.s=distance,.ue=sign*q[id]};
        }
        double trip=sign*(p->trip_arc[e][side]-stag);
        if (trip<0 || trip>edge[count-1].s) trip=edge[count-1].s;
        AlulaBLMarchResult r=alula_bl_march_surface_controlled(&p->physics,edge,count,trip,1,
            inverse,1,march,regimes);
        for (size_t k=0;k<r.completed;++k) {
            size_t id=side?right+k:left-k;
            p->stations[id]=march[k];p->regimes[id]=regimes[k];p->speed[id]=sign*march[k].ue;
        }
        if (r.status!=ALULA_BL_OK) {
            p->failed_station=side?right+r.completed:left-r.completed;
            return r.status;
        }
        p->interval_residual=fmax(p->interval_residual,r.maximum_residual);
    }
    if (right>left+1) {
        p->stations[left+1]=(AlulaBLStation){.theta=.5*(p->stations[left].theta+p->stations[right].theta),
            .delta_star=.5*(p->stations[left].delta_star+p->stations[right].delta_star)};
        p->regimes[left+1]=ALULA_BL_SIMILARITY;
    }
    AlulaPoint gap=subtract(d->points[start],d->points[end]),bisector=subtract(bottom,top);
    double width=.5*cross(bisector,gap);
    p->join_gap[e]=width;
    double tangent=cross(top,bottom)/(hypot(top.x,top.y)*hypot(bottom.x,bottom.y));
    double xmin=d->points[start].x,xmax=xmin;
    for (size_t i=start+1;i<=end;++i) { xmin=fmin(xmin,d->points[i].x);xmax=fmax(xmax,d->points[i].x); }
    double fluid_gap=hypot(gap.x,gap.y)<1e-4*fmax(1e-12,xmax-xmin)?0:width;
    p->fluid_gap[e]=fluid_gap;p->tangent_cross[e]=tangent;
    for (size_t j=0;j<d->wake_count;++j) {
        size_t id=first_wake+j;
        edge[j]=(AlulaBLEdge){s[end]-stag+s[id],q[id],alula_bl_dead_air(fluid_gap,tangent,s[id])};
    }
    AlulaBLMarchResult r=alula_bl_march_wake_join(&p->physics,p->stations+start,p->stations+end,
        edge,d->wake_count,width,inverse,p->stations+first_wake);
    for (size_t j=0;j<r.completed;++j) {
        size_t id=first_wake+j;
        p->regimes[id]=ALULA_BL_WAKE;p->speed[id]=p->stations[id].ue;
    }
    if (r.status!=ALULA_BL_OK) p->failed_station=first_wake+r.completed;
    p->interval_residual=fmax(p->interval_residual,r.maximum_residual);
    return r.status;
}

AlulaPanelBL *alula_panel_bl_initialize(const AlulaPanelDisplacement *outer,double reynolds,double ncrit,
    const double (*trips)[2])
{
    const AlulaDisplacementData *d=alula_panel_displacement_data(outer);
    AlulaBLContext physics;
    if (!d || !alula_bl_context(&physics,reynolds,0,1.4,ncrit)) return NULL;
    for (size_t e=0;trips && e<d->bodies;++e) for (int side=0;side<2;++side)
        if (!(trips[e][side]>0 && trips[e][side]<=1)) return NULL;
    size_t n=d->stations;
    AlulaPanelBL *p=calloc(1,sizeof *p);
    if (!p) return NULL;
    p->bodies=d->bodies;p->count=n;p->physics=physics;p->outer_residual=INFINITY;
    p->stations=calloc(n,sizeof *p->stations);p->regimes=calloc(n,sizeof *p->regimes);
    p->speed=calloc(n,sizeof *p->speed);
    AlulaBLEdge *edge=malloc(n*sizeof *edge);
    AlulaBLStation *march=malloc(n*sizeof *march);
    AlulaBLRegime *regime=malloc(n*sizeof *regime);
    double *zero=calloc(n,sizeof *zero),*mass=malloc(n*sizeof *mass),*speed=malloc(n*sizeof *speed);
    if (!p->stations || !p->regimes || !p->speed || !edge || !march || !regime || !zero || !mass || !speed) {
        alula_panel_bl_free(p);p=NULL;
        goto done;
    }
    for (size_t e=0;e<d->bodies;++e) {
        const double natural[]={1,1};
        p->status=seed_body(p,d,e,trips?trips[e]:natural,edge,march,regime,zero);
        if (p->status!=ALULA_BL_OK) { p->failed_body=e;goto done; }
    }
    for (size_t i=0;i<n;++i) mass[i]=p->speed[i]*p->stations[i].delta_star;
    if (!alula_panel_displacement_evaluate(outer,mass,n,speed)) { p->status=ALULA_BL_DOMAIN;goto done; }
    p->outer_residual=0;
    for (size_t i=0;i<n;++i) p->outer_residual=fmax(p->outer_residual,fabs(speed[i]-p->speed[i]));
    p->complete=1;
done:
    free(edge);free(march);free(regime);free(zero);free(mass);free(speed);
    return p;
}

void alula_panel_bl_free(AlulaPanelBL *p)
{
    if (!p) return;
    free(p->stations);free(p->regimes);free(p->speed);free(p);
}
