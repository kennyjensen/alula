// SPDX-License-Identifier: GPL-2.0-or-later
#include "grid_relax.h"
#include "panel_internal.h"

#include <float.h>
#include <stdlib.h>
#include <string.h>

#define CONSTRAINTS 12
#define POLYGON_VERTICES (CONSTRAINTS+4)

typedef struct {
    int id[3];
    double scale,target;
} Corner;

typedef struct {
    int count,corners[CONSTRAINTS];
    double radius;
} Vertex;

typedef struct { AlulaPoint gradient; double bound; } HalfPlane;

static double length(AlulaPoint p) { return hypot(p.x,p.y); }

static double determinant(AlulaPoint a,AlulaPoint b,AlulaPoint c)
{
    return cross(subtract(b,a),subtract(c,b));
}

static int clip(const AlulaPoint *p,int count,HalfPlane constraint,double slack,AlulaPoint *out)
{
    int next=0;
    for (int i=0;i<count;++i) {
        AlulaPoint a=p[i],b=p[i+1==count?0:i+1];
        double ca=dot(a,constraint.gradient)-constraint.bound+slack;
        double cb=dot(b,constraint.gradient)-constraint.bound+slack;
        if (ca>=0) out[next++]=a;
        if ((ca>=0)!=(cb>=0) && ca!=0 && cb!=0) {
            double s=ca/(ca-cb);
            out[next++]=(AlulaPoint){a.x+s*(b.x-a.x),a.y+s*(b.y-a.y)};
        }
    }
    return next;
}

static int feasible(const HalfPlane *h,int count,double slack,AlulaPoint *polygon)
{
    AlulaPoint p[2][POLYGON_VERTICES]={{{-1,-1},{1,-1},{1,1},{-1,1}}};
    int n=4;
    for (int i=0;i<count && n;++i) n=clip(p[i%2],n,h[i],slack,p[(i+1)%2]);
    if (n) memcpy(polygon,p[count%2],(size_t)n*sizeof *polygon);
    return n;
}

static AlulaPoint closest(const AlulaPoint *polygon,int n)
{
    AlulaPoint best=polygon[0];
    double squared=dot(best,best);
    for (int i=0;i<n;++i) {
        AlulaPoint a=polygon[i],v=subtract(polygon[i+1==n?0:i+1],a);
        double span=dot(v,v);
        double t=span>0?fmin(1,fmax(0,-dot(a,v)/span)):0;
        AlulaPoint p={a.x+t*v.x,a.y+t*v.y};
        double d=dot(p,p);
        if (d<squared) {squared=d;best=p;}
    }
    return best;
}

static AlulaPoint move_vertex(AlulaPoint *nodes,const Corner *corners,const Vertex *v,
    int k,AlulaPoint origin)
{
    AlulaPoint p=nodes[k],d=subtract(p,origin),current={d.x/v->radius,d.y/v->radius};
    HalfPlane planes[CONSTRAINTS]={0};
    double needed=0;
    for (int i=0;i<v->count;++i) {
        const Corner *r=corners+v->corners[i];
        AlulaPoint a=nodes[r->id[0]],b=nodes[r->id[1]],c=nodes[r->id[2]];
        AlulaPoint edge=r->id[0]==k?subtract(c,b):r->id[1]==k?subtract(a,c):subtract(b,a);
        double factor=v->radius/r->scale;
        AlulaPoint gradient={-edge.y*factor,edge.x*factor};
        double defect=r->target-determinant(a,b,c)/r->scale;
        planes[i]=(HalfPlane){gradient,defect+dot(gradient,current)};
        needed=fmax(needed,defect);
    }
    AlulaPoint polygon[POLYGON_VERTICES];
    double slack=0;
    int count=feasible(planes,v->count,slack,polygon);
    if (!count) {
        double lower=0,upper=needed+64*DBL_EPSILON*fmax(1,needed);
        count=feasible(planes,v->count,upper,polygon);
        if (!count) return p;
        for (int i=0;i<29;++i) {
            double mid=.5*(lower+upper);
            int n=feasible(planes,v->count,mid,polygon);
            if (n) {upper=mid;count=n;}
            else lower=mid;
        }
        slack=upper;
    }
    int contains_origin=1;
    for (int i=0;i<v->count;++i) if (planes[i].bound>slack) {contains_origin=0;break;}
    AlulaPoint local=contains_origin?(AlulaPoint){0}:closest(polygon,count);
    return (AlulaPoint){origin.x+v->radius*local.x,origin.y+v->radius*local.y};
}

int alula_repair_passage(AlulaPoint *nodes,int nx,int tubes)
{
    int width=tubes+1,n=(nx+1)*width,number=4*nx*tubes;
    Corner *corners=malloc((size_t)number*sizeof *corners);
    Vertex *vertices=calloc((size_t)n,sizeof *vertices);
    AlulaPoint *initial=malloc((size_t)n*sizeof *initial);
    unsigned char *active=malloc((size_t)n);
    int valid=0;
    if (!corners || !vertices || !initial || !active) goto done;
    memcpy(initial,nodes,(size_t)n*sizeof *nodes);
    int next=0;
    for (int i=0;i<nx;++i) for (int j=0;j<tubes;++j) {
        int k=i*width+j, q[]={k,k+width,k+width+1,k+1};
        for (int m=0;m<4;++m) {
            int a=q[m],b=q[(m+1)%4],c=q[(m+2)%4];
            double scale=length(subtract(nodes[b],nodes[a]))*length(subtract(nodes[c],nodes[b]));
            if (!(scale>0) || !isfinite(scale)) goto done;
            double v=determinant(nodes[a],nodes[b],nodes[c])/scale;
            corners[next]=(Corner){{a,b,c},scale,v>1e-12?fmin(.02,.5*v):.02};
            int ids[]={a,b,c};
            for (int t=0;t<3;++t) vertices[ids[t]].corners[vertices[ids[t]].count++]=next;
            ++next;
        }
    }
    for (int i=1;i<nx;++i) for (int j=1;j<tubes;++j) {
        int k=i*width+j;
        double d=INFINITY;
        int near[]={k-width,k+width,k-1,k+1};
        for (int l=0;l<4;++l) d=fmin(d,length(subtract(initial[k],initial[near[l]])));
        vertices[k].radius=6*d;
    }
    for (int sweep=0;sweep<=240;++sweep) {
        memset(active,0,(size_t)n);
        valid=1;
        int defective=0;
        for (int k=0;k<number;++k) {
            Corner r=corners[k];
            double value=determinant(nodes[r.id[0]],nodes[r.id[1]],nodes[r.id[2]])/r.scale;
            if (!(value>1e-12)) valid=0;
            if (value<r.target-1e-8) {
                defective=1;
                for (int t=0;t<3;++t) active[r.id[t]]=1;
            }
        }
        if ((valid && !defective) || sweep==240) break;
        double max_move=0;
        for (int i=0;i<n;++i) {
            int k=sweep%2?n-1-i:i;
            if (!active[k] || vertices[k].radius<=0) continue;
            AlulaPoint next_point=move_vertex(nodes,corners,vertices+k,k,initial[k]);
            max_move=fmax(max_move,length(subtract(next_point,nodes[k]))/vertices[k].radius);
            nodes[k]=next_point;
        }
        if (max_move<1e-11) break;
    }
done:
    free(active);free(initial);free(vertices);free(corners);
    return valid;
}
