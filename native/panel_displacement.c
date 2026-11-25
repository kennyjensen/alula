// SPDX-License-Identifier: GPL-2.0-or-later
#include "panel_displacement.h"
#include "panel_sources.h"
#include "dense.h"

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* One Hermite wake source depends on at most three wake masses and four
 * trailing-edge surface masses at a sharp TE, or four interior wake masses. */
typedef struct { size_t count,id[7]; double value[7]; } Weights;
typedef struct { Panel p; Weights polynomial[3]; } WakeSegment;

struct AlulaPanelDisplacement {
    AlulaDisplacementData data;
    PanelRange bodies[6];
    size_t wet,solid;
    AlulaPoint freestream,*points,*basis;
    double *arc,*q0,*influence;
    Segment *panels;
    WakeSegment *wake;
};

static AlulaPoint advance(AlulaPoint p,AlulaPoint d,double step)
{
    return (AlulaPoint){p.x+step*d.x,p.y+step*d.y};
}

static double length(AlulaPoint p) { return hypot(p.x,p.y); }

static int unit(AlulaPoint *p)
{
    double l=length(*p);
    if (!(l>0) || !isfinite(l)) return 0;
    p->x/=l; p->y/=l;
    return 1;
}

static void weight(Weights *w,size_t id,double value)
{
    size_t i=0;
    while (i<w->count && w->id[i]!=id) ++i;
    if (i==w->count) { w->id[i]=id; ++w->count; }
    w->value[i]+=value;
}

static void combine(Weights *out,const Weights *w,double scale)
{
    for (size_t i=0;i<w->count;++i) weight(out,w->id[i],scale*w->value[i]);
}

static void accumulate(double *row,const Weights *w,double scale)
{
    for (size_t i=0;i<w->count;++i) row[w->id[i]]+=scale*w->value[i];
}

/* Projection of a prescribed piecewise-linear vortex sheet, including the
 * strength jumps that drive the finite base's source and vortex panels. */
static AlulaPoint vortex_velocity(const AlulaPanelDisplacement *o,const double *gamma,
    AlulaPoint point,const AlulaPoint *tangent,double *influence)
{
    AlulaPoint v=o->freestream;
    for (size_t i=0;i<o->solid;++i) {
        const Segment *s=o->panels+i;
        AlulaPoint b[2];
        alula_panel_vortex(point,&s->geometry,0,b);
        size_t node=s->geometry.node;
        if (i<o->wet) {
            for (size_t k=0;k<2;++k) {
                v=advance(v,b[k],gamma[node+k]);
                if (influence) influence[node+k]+=dot(b[k],*tangent);
            }
        } else {
            AlulaPoint sum={b[0].x+b[1].x,b[0].y+b[1].y};
            AlulaPoint base={s->source*sum.y+s->vortex*sum.x,-s->source*sum.x+s->vortex*sum.y};
            v=advance(v,base,gamma[node]-gamma[s->last_node]);
            if (influence) {
                double c=dot(base,*tangent);
                influence[node]+=c; influence[s->last_node]-=c;
            }
        }
    }
    return v;
}

static double wake_extent(double ratio,size_t intervals,double step)
{
    return ratio==1?intervals*step:step*expm1(intervals*log(ratio))/(ratio-1);
}

static int setup_wake(AlulaPanelDisplacement *o,size_t e,double extent,AlulaWakeStart start,
    const AlulaPoint *path)
{
    const PanelRange *body=o->bodies+e;
    size_t nw=o->data.wake_count, first=o->data.surface_stations+e*nw;
    AlulaPoint *points=o->points+first;
    WakeSegment *segments=o->wake+e*(nw-1);
    const Panel *top=&o->panels[body->offset].geometry;
    const Segment *bottom=o->panels+body->offset+body->wetted-1;
    AlulaPoint center={.5*(top->a.x+bottom->end.x),.5*(top->a.y+bottom->end.y)};
    double chord=0;
    for (size_t j=body->start;j<=body->end;++j)
        chord=fmax(chord,length(subtract(o->points[j],center)));
    AlulaPoint direction=body->wake_direction;
    if (!unit(&direction)) return 0;
    points[0]=advance(center,direction,body->finite?1e-4*chord:0);
    if (path) {
        if (!(length(subtract(path[0],points[0]))<=1e-10)) return 0;
        memcpy(points,path,nw*sizeof *path);
    } else {
        double step=.5*(top->length+bottom->geometry.length),lo=1,hi=4;
        extent*=chord;
        if (!isfinite(extent) || step*(nw-1)>extent) return 0;
        while (wake_extent(hi,nw-1,step)<extent) hi*=2;
        for (int k=0;k<70;++k) {
            double mid=.5*(lo+hi);
            if (wake_extent(mid,nw-1,step)>extent) hi=mid;
            else lo=mid;
        }
        double ratio=.5*(lo+hi),distance=0;
        for (size_t j=1;j<nw;++j) {
            distance+=step;
            if (start==ALULA_WAKE_INVISCID) {
                AlulaPoint prev=points[j-1],t=j==1?direction:vortex_velocity(o,o->q0,prev,NULL,NULL);
                if (!unit(&t)) return 0;
                t=vortex_velocity(o,o->q0,advance(prev,t,.5*step),NULL,NULL);
                if (!unit(&t)) return 0;
                points[j]=advance(prev,t,step);
            } else points[j]=advance(points[0],direction,distance);
            step*=ratio;
        }
    }
    for (size_t j=0;j+1<nw;++j) {
        AlulaPoint a=points[j],b=points[j+1],t=subtract(b,a);
        double h=length(t);
        if (!unit(&t)) return 0;
        segments[j].p=(Panel){.a=a,.center={.5*(a.x+b.x),.5*(a.y+b.y)},
            .length=h,.tangent=t,.normal={t.y,-t.x},.node=first+j};
        o->arc[first+j+1]=o->arc[first+j]+h;
    }
    Weights left={0};
    if (body->finite) {
        double h=segments[0].p.length;
        weight(&left,first,-1/h);weight(&left,first+1,1/h);
    } else {
        weight(&left,body->start,-1/top->length);weight(&left,body->start+1,1/top->length);
        weight(&left,body->end-1,-1/bottom->geometry.length);weight(&left,body->end,1/bottom->geometry.length);
    }
    for (size_t j=0;j+1<nw;++j) {
        Weights right={0},*coeff=segments[j].polynomial;
        double h=segments[j].p.length;
        if (j+2<nw) {
            double hr=segments[j+1].p.length,l=hr/(h*(h+hr)),r=h/(hr*(h+hr));
            weight(&right,first+j,-l); weight(&right,first+j+1,l-r); weight(&right,first+j+2,r);
        }
        coeff[0]=left;
        combine(coeff+1,&left,-4); combine(coeff+1,&right,-2);
        weight(coeff+1,first+j,-6/h);weight(coeff+1,first+j+1,6/h);
        combine(coeff+2,&left,3); combine(coeff+2,&right,3);
        weight(coeff+2,first+j,6/h);weight(coeff+2,first+j+1,-6/h);
        left=right;
    }
    return 1;
}

static int side(AlulaPoint a,AlulaPoint b,AlulaPoint p)
{
    AlulaPoint edge=subtract(b,a);
    double c=cross(edge,subtract(p,a)),tol=1e-11*length(edge);
    return (c>tol)-(c< -tol);
}

static int touch(AlulaPoint a,AlulaPoint b,AlulaPoint c,AlulaPoint d)
{
    const double eps=1e-11;
    if (fmax(a.x,b.x)+eps<fmin(c.x,d.x) || fmax(c.x,d.x)+eps<fmin(a.x,b.x) ||
        fmax(a.y,b.y)+eps<fmin(c.y,d.y) || fmax(c.y,d.y)+eps<fmin(a.y,b.y)) return 0;
    return side(a,b,c)*side(a,b,d)<=0 && side(c,d,a)*side(c,d,b)<=0;
}

static int valid_wakes(const AlulaPanelDisplacement *o)
{
    size_t nw=o->data.wake_count, ns=o->data.surface_stations;
    for (size_t e=0;e<o->data.bodies;++e) for (size_t j=0;j+1<nw;++j) {
        const AlulaPoint *ab=o->points+ns+e*nw+j;
        const Panel *p=&o->wake[e*(nw-1)+j].p;
        int inside[6]={0};
        for (size_t i=0;i<o->solid;++i) {
            const Segment *s=o->panels+i;
            AlulaPoint a=s->geometry.a,b=s->end,c=p->center;
            if ((a.y>c.y)!=(b.y>c.y) && c.x<a.x+(b.x-a.x)*(c.y-a.y)/(b.y-a.y))
                inside[s->body]^=1;
            const PanelRange *body=o->bodies+e;
            if (!j && !body->finite && (i==body->offset || i==body->offset+body->wetted-1)) continue;
            if (touch(ab[0],ab[1],a,b)) return 0;
        }
        for (size_t f=0;f<o->data.bodies;++f) if (inside[f]) return 0;
        for (size_t f=0;f<=e;++f) {
            size_t stop=f<e?nw-1:(j>1?j-1:0);
            for (size_t k=0;k<stop;++k) {
                const AlulaPoint *cd=o->points+ns+f*nw+k;
                if (touch(ab[0],ab[1],cd[0],cd[1])) return 0;
            }
        }
    }
    return 1;
}

/* Multiply the local dm/ds stencils into the RHS as each influence is
 * evaluated. This avoids allocating the dense zero-filled source map. */
static void source_row(const AlulaPanelDisplacement *o,AlulaPoint target,
    const AlulaPoint *tangent,const double *heading,const double *chart,double sign,double *row)
{
    size_t nw=o->data.bodies*(o->data.wake_count-1);
    for (size_t i=0;i<o->wet+nw;++i) {
        int wake=i>=o->wet, degree=wake?2:0;
        const WakeSegment *w=wake?o->wake+i-o->wet:NULL;
        const Panel *p=wake?&w->p:&o->panels[i].geometry;
        double coeff[3];
        if (tangent) {
            AlulaPoint v[3];
            alula_source_velocity(p,target,degree,wake && sign>0,0,v);
            for (int k=0;k<=degree;++k) coeff[k]=dot(v[k],*tangent);
        } else alula_source_psi(p,target,heading[i],chart[i],degree,coeff);
        if (wake) {
            for (int k=0;k<=degree;++k) accumulate(row,w->polynomial+k,sign*coeff[k]);
        } else {
            double v=sign*coeff[0]/p->length;
            row[p->node]-=v;row[p->node+1]+=v;
        }
    }
}

static int build_influence(AlulaPanelDisplacement *o,const PanelRow *rows,const double *lu,
    const int *pivots,double *rhs)
{
    const size_t ns=o->data.surface_stations,n=o->data.stations,nw=o->data.wake_count;
    size_t sources=o->wet+o->data.bodies*(nw-1);
    /* Each source has one physical heading and a fixed argument chart per
     * target body. Reuse them for every boundary point on that body. */
    double *angles=malloc((1+o->data.bodies)*sources*sizeof *angles);
    if (!angles) return 0;
    double *charts=angles+sources;
    for (size_t i=0;i<sources;++i) {
        const Panel *p=i<o->wet?&o->panels[i].geometry:&o->wake[i-o->wet].p;
        angles[i]=atan2(p->tangent.y,p->tangent.x);
    }
    for (size_t e=0;e<o->data.bodies;++e) {
        const PanelRange *b=o->bodies+e;
        AlulaPoint centroid={0};
        for (size_t i=b->start;i<=b->end;++i) centroid=advance(centroid,o->points[i],1);
        centroid.x/=b->wetted+1; centroid.y/=b->wetted+1;
        for (size_t i=0;i<sources;++i) {
            const Panel *p=i<o->wet?&o->panels[i].geometry:&o->wake[i-o->wet].p;
            AlulaPoint d=subtract(centroid,p->a);
            charts[e*sources+i]=i<o->wet && e==o->panels[i].body?ALULA_PI/2:
                atan2(cross(p->tangent,d),dot(p->tangent,d));
        }
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(ns>128)
#endif
    for (size_t i=0;i<ns;++i)
        source_row(o,rows[i].point,rows[i].interior?&rows[i].tangent:NULL,
            angles,charts+rows[i].body*sources,-1,rhs+i*n);
    free(angles);
    if (!alula_dense_apply_rows(lu,ns+o->data.bodies,pivots,rhs,n)) return 0;
    memcpy(o->influence,rhs,ns*n*sizeof *rhs);
    double *vortex=calloc(ns,sizeof *vortex);
    if (!vortex) return 0;
    int ok=1;
    for (size_t i=ns;i<n && ok;++i) {
        size_t e=(i-ns)/nw,j=(i-ns)%nw;
        double *row=o->influence+i*n;
        if (!j) {
            size_t te=o->bodies[e].end;
            o->q0[i]=o->q0[te];o->basis[i]=o->basis[te];
            memcpy(row,rhs+te*n,n*sizeof *row);
            continue;
        }
        WakeSegment *segments=o->wake+e*(nw-1);
        AlulaPoint tangent=advance(segments[j-1].p.tangent,segments[j==nw-1?j-1:j].p.tangent,1);
        if (length(tangent)<.1 || !unit(&tangent)) { ok=0; break; }
        memset(vortex,0,ns*sizeof *vortex);
        AlulaPoint p=o->points[i];
        vortex_velocity(o,o->q0,p,&tangent,vortex);
        o->basis[i]=tangent;
        source_row(o,p,&tangent,NULL,NULL,1,row);
        for (size_t k=0;k<ns;++k) {
            o->basis[i]=advance(o->basis[i],o->basis[k],vortex[k]);
            const double *gamma=rhs+k*n;
            for (size_t col=0;col<n;++col) row[col]+=vortex[k]*gamma[col];
        }
        o->q0[i]=dot(o->basis[i],o->freestream);
    }
    free(vortex);
    for (size_t i=0;i<n && ok;++i) if (!isfinite(o->q0[i])) ok=0;
    for (size_t i=0;i<n*n && ok;++i) if (!isfinite(o->influence[i])) ok=0;
    return ok;
}

AlulaPanelDisplacement *alula_panel_displacement_create(const AlulaContour *contours,
    size_t bodies,double alpha,size_t nw,double extent,AlulaWakeStart start,const AlulaPoint *paths)
{
    size_t solid,wet;
    if (!isfinite(alpha) || (!paths && (!(extent>0) || !isfinite(extent) ||
        (start!=ALULA_WAKE_STRAIGHT && start!=ALULA_WAKE_INVISCID))) || nw<4 ||
        !alula_panel_dimensions(contours,bodies,&solid,&wet)) return NULL;
    size_t ns=wet+bodies, size=ns+bodies;
    if (nw>((size_t)INT_MAX-ns)/bodies) return NULL;
    size_t n=ns+bodies*nw;
    if (n>SIZE_MAX/sizeof(double)/n) return NULL;
    AlulaPanelDisplacement *o=calloc(1,sizeof *o);
    if (!o) return NULL;
    o->wet=wet;o->solid=solid;
    o->freestream=(AlulaPoint){cos(alpha*(ALULA_PI/180)),sin(alpha*(ALULA_PI/180))};
    o->points=malloc(n*sizeof *o->points);o->arc=calloc(n,sizeof *o->arc);
    o->q0=malloc(n*sizeof *o->q0);o->influence=calloc(n*n,sizeof *o->influence);
    o->basis=malloc(n*sizeof *o->basis);
    o->panels=calloc(solid,sizeof *o->panels);o->wake=calloc(bodies*(nw-1),sizeof *o->wake);
    o->data=(AlulaDisplacementData){.stations=n,.bodies=bodies,.surface_stations=ns,.wake_count=nw,
        .points=o->points,.arc=o->arc,.inviscid=o->q0,.influence=o->influence};
    PanelRow *rows=calloc(ns,sizeof *rows);
    double *lu=calloc(size*size,sizeof *lu),*rhs=calloc(size*n,sizeof *rhs);
    int *pivots=malloc(size*sizeof *pivots);
    int ok=0;
    if (!rows || !pivots || !lu || !rhs || !o->points || !o->arc || !o->q0 || !o->influence ||
        !o->basis || !o->panels || !o->wake || !alula_panel_boundary_setup(contours,bodies,solid,o->panels,o->bodies,rows,wet)) goto done;
    alula_panel_boundary_assemble(lu,rhs,rhs+size,bodies,o->panels,wet,solid,rows,o->bodies,1);
    if (!alula_dense_factor(lu,size,pivots) || !alula_dense_apply_rows(lu,size,pivots,rhs,1) ||
        !alula_dense_apply_rows(lu,size,pivots,rhs+size,1)) goto done;
    for (size_t i=0;i<ns;++i) {
        o->basis[i]=(AlulaPoint){rhs[i],rhs[size+i]};
        o->q0[i]=dot(o->basis[i],o->freestream);
    }
    memset(rhs,0,2*size*sizeof *rhs);
    for (size_t e=0;e<bodies;++e) {
        const PanelRange *b=o->bodies+e;
        o->data.surface_start[e]=b->start;
        for (size_t j=0;j<b->wetted;++j) {
            const Segment *s=o->panels+b->offset+j;
            o->points[b->start+j]=s->geometry.a;
            o->arc[b->start+j+1]=o->arc[b->start+j]+s->geometry.length;
        }
        o->points[b->end]=o->panels[b->offset+b->wetted-1].end;
        if (!setup_wake(o,e,extent,start,paths?paths+e*nw:NULL)) goto done;
    }
    o->data.surface_start[bodies]=ns;
    ok=valid_wakes(o) && build_influence(o,rows,lu,pivots,rhs);
done:
    free(rows);free(lu);free(rhs);free(pivots);
    if (ok) return o;
    alula_panel_displacement_free(o);
    return NULL;
}

void alula_panel_displacement_free(AlulaPanelDisplacement *o)
{
    if (!o) return;
    free(o->points);free(o->arc);free(o->q0);free(o->influence);
    free(o->panels);free(o->wake);free(o->basis);free(o);
}

const AlulaDisplacementData *alula_panel_displacement_data(const AlulaPanelDisplacement *o)
{
    return o?&o->data:NULL;
}

int alula_panel_displacement_set_incidence(AlulaPanelDisplacement *o,double alpha)
{
    if (!o || !isfinite(alpha)) return 0;
    o->freestream=(AlulaPoint){cos(alpha*(ALULA_PI/180)),sin(alpha*(ALULA_PI/180))};
    for (size_t i=0;i<o->data.stations;++i) o->q0[i]=dot(o->basis[i],o->freestream);
    return 1;
}

int alula_panel_displacement_evaluate(const AlulaPanelDisplacement *o,const double *mass,size_t n,double *speed)
{
    if (!o || !mass || !speed || mass==speed || n!=o->data.stations) return 0;
    for (size_t i=0;i<n;++i) if (!isfinite(mass[i])) return 0;
    for (size_t i=0;i<n;++i) {
        const double *row=o->influence+i*n;
        double v=o->q0[i];
        for (size_t j=0;j<n;++j) v+=row[j]*mass[j];
        if (!isfinite(v)) return 0;
        speed[i]=v;
    }
    return 1;
}

int alula_panel_displacement_wake_velocity(const AlulaPanelDisplacement *o,
    const double *mass,size_t n,AlulaPoint *velocity)
{
    if (!o || !mass || !velocity || n!=o->data.stations) return 0;
    size_t nw=n-o->data.surface_stations-o->data.bodies;
    double *speed=malloc(n*sizeof *speed),*coeff=calloc(3*nw,sizeof *coeff);
    int ok=speed && coeff && alula_panel_displacement_evaluate(o,mass,n,speed);
    if (!ok) goto done;
    for (size_t j=0;j<nw;++j) for (int k=0;k<3;++k) {
        const Weights *w=o->wake[j].polynomial+k;
        for (size_t i=0;i<w->count;++i) coeff[3*j+k]+=w->value[i]*mass[w->id[i]];
    }
    for (size_t j=0;j<nw;++j) {
        AlulaPoint p=o->wake[j].p.center,v=vortex_velocity(o,speed,p,NULL,NULL);
        for (size_t i=0;i<o->wet+nw;++i) {
            int wake=i>=o->wet;
            size_t k=wake?i-o->wet:0;
            const Panel *source=wake?&o->wake[k].p:&o->panels[i].geometry;
            AlulaPoint basis[3];
            alula_source_velocity(source,p,wake?2:0,1,wake && k==j,basis);
            if (wake) {
                for (int l=0;l<3;++l) v=advance(v,basis[l],coeff[3*k+l]);
            } else v=advance(v,basis[0],(mass[source->node+1]-mass[source->node])/source->length);
        }
        if (!isfinite(v.x) || !isfinite(v.y)) { ok=0;break; }
        velocity[j]=v;
    }
done:
    free(speed);free(coeff);
    return ok;
}
