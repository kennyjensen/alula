// SPDX-License-Identifier: GPL-2.0-or-later
#include "airfoil_grid_internal.h"
#include "body_internal.h"

#include <stdlib.h>
#include <string.h>

typedef struct { int index; double fraction; } Sample;

static double lerp(double a,double b,double t) { return a+t*(b-a); }

static AlulaPoint mix(AlulaPoint a,AlulaPoint b,double t)
{
    return (AlulaPoint){lerp(a.x,b.x,t),lerp(a.y,b.y,t)};
}

static Sample bracket(const double *x,int n,double value)
{
    if (n==1 || value<=x[0]) return (Sample){0,0};
    if (value>=x[n-1]) return (Sample){n-2,1};
    int lo=0,hi=n-1;
    while (hi-lo>1) {
        int mid=(lo+hi)/2;
        if (value<x[mid]) hi=mid;
        else lo=mid;
    }
    return (Sample){lo,(value-x[lo])/(x[hi]-x[lo])};
}

static int equal_points(const AlulaPoint *a,const AlulaPoint *b,size_t count)
{
    for (size_t i=0;i<count;++i)
        if (a[i].x!=b[i].x || a[i].y!=b[i].y) return 0;
    return 1;
}

static int compatible(const AlulaAirfoilGrid *f,const AlulaAirfoilGrid *p)
{
    const AlulaBody *a=f->body,*b=p->body;
    if (a->in.bodies!=b->in.bodies || a->in.mach!=b->in.mach || a->in.gamma!=b->in.gamma ||
        a->in.alpha!=b->in.alpha || a->in.length_scale!=b->in.length_scale ||
        a->in.center.x!=b->in.center.x || a->in.center.y!=b->in.center.y ||
        a->in.mode!=b->in.mode || a->in.ismom!=b->in.ismom ||
        f->stations[0]!=p->stations[0] || f->stations[a->in.nx]!=p->stations[b->in.nx] ||
        f->lower_y!=p->lower_y || f->upper_y!=p->upper_y) return 0;
    for (int g=0;g<a->in.bodies;++g) {
        AlulaBodySurfaceView af,bf;
        if (f->bases[g].count!=p->bases[g].count ||
            !equal_points(f->bases[g].points,p->bases[g].points,f->bases[g].count) ||
            !alula_body_coordinates_surface(f->coordinates,f->packed,g,0,&af) ||
            !alula_body_coordinates_surface(p->coordinates,p->packed,g,0,&bf) ||
            af.curve->count!=bf.curve->count ||
            !equal_points(af.curve->points,bf.curve->points,(size_t)af.curve->count)) return 0;
    }
    return 1;
}

/* Match material fractions on the wall and physical x along the two cuts.
 * Each bank can have its own station map in a passage between two bodies. */
static void stream_labels(const AlulaAirfoilGrid *g,int body,int upper,double *labels)
{
    int nx=g->body->in.nx;
    const double *x=g->stations;
    AlulaBodyCut cut=g->cuts[body];
    AlulaBodySurfaceView s;
    alula_body_coordinates_surface(g->coordinates,g->packed,body,upper,&s);
    for (int i=0;i<=nx;++i) {
        labels[i]=i<cut.leading?(x[i]-x[cut.leading])/(x[cut.leading]-x[0]):
            i<=cut.trailing?s.fractions[i-cut.leading]:
            1+(x[i]-x[cut.trailing])/(x[nx]-x[cut.trailing]);
    }
}

static void fractions(const double *mass,int count,double *x)
{
    x[0]=0;
    for (int i=0;i<count;++i) x[i+1]=x[i]+mass[i];
    double total=x[count];
    for (int i=1;i<=count;++i) x[i]/=total;
}

static void center(const double *edges,int count,double *x)
{
    for (int i=0;i<count;++i) x[i]=.5*(edges[i]+edges[i+1]);
}

static Sample combine(Sample a,Sample b,double eta)
{
    double low=a.index+a.fraction,high=b.index+b.fraction;
    if (eta==0 || low==high) return a;
    if (eta==1) return b;
    double i=lerp(low,high,eta);
    if (i==low) return a;
    if (i==high) return b;
    int k=(int)i;
    return (Sample){k,i-k};
}

static AlulaPoint interpolate_node(const AlulaBody *b,const AlulaBodyState *src,int g,Sample i,Sample j)
{
    int k=i.index,l=j.index;
    return mix(mix(point(b,src,g,k,l),point(b,src,g,k+1,l),i.fraction),
        mix(point(b,src,g,k,l+1),point(b,src,g,k+1,l+1),i.fraction),j.fraction);
}

static double interpolate_density(const AlulaBody *b,const AlulaBodyState *src,int g,Sample i,Sample j)
{
    int k=i.index,l=j.index,m=b->tubes[g]>1?l+1:l;
    double a=lerp(log(src->density[index_section(b,g,k,l)]),log(src->density[index_section(b,g,k+1,l)]),i.fraction);
    double z=lerp(log(src->density[index_section(b,g,k,m)]),log(src->density[index_section(b,g,k+1,m)]),i.fraction);
    return lerp(a,z,j.fraction);
}

static AlulaPoint wall_correction(const AlulaBody *f,const AlulaBody *p,const AlulaBodyState *target,
    const AlulaBodyState *source,int g,int i,int upper,Sample section)
{
    int body=upper?g:g-1;
    if (body<0 || body>=f->in.bodies || i<f->cuts[body].leading || i>f->cuts[body].trailing)
        return (AlulaPoint){0};
    int j=upper?f->tubes[g]:0,k=upper?p->tubes[g]:0;
    return sub(point(f,target,g,i,j),mix(point(p,source,g,section.index,k),
        point(p,source,g,section.index+1,k),section.fraction));
}

static void transfer_passage(AlulaAirfoilGrid *f,const AlulaAirfoilGrid *p,const AlulaBodyState *target,
    const AlulaBodyState *source,int g,double weight,AlulaPoint *nodes,double *work)
{
    const AlulaBody *fb=f->body,*pb=p->body;
    int pn=pb->in.nx,fn=fb->in.nx,pt=pb->tubes[g],ft=fb->tubes[g];
    double *p0=work,*p1=p0+pn+1,*f0=p1+pn+1,*f1=f0+fn+1,*pc0=f1+fn+1,*pc1=pc0+pn;
    double *py=pc1+pn,*fy=py+pt+1,*pcy=fy+ft+1;
    int lower=g?g-1:g,upper=g<fb->in.bodies?g:g-1;
    stream_labels(p,lower,g!=0,p0);stream_labels(f,lower,g!=0,f0);
    stream_labels(p,upper,g==fb->in.bodies,p1);stream_labels(f,upper,g==fb->in.bodies,f1);
    fractions(source->mass+pb->column[g],pt,py);
    fractions(target->mass+fb->column[g],ft,fy);
    center(p0,pn,pc0);center(p1,pn,pc1);center(py,pt,pcy);
    for (int i=0;i<=fn;++i) {
        Sample lo=bracket(p0,pn+1,f0[i]),hi=bracket(p1,pn+1,f1[i]);
        AlulaPoint dl=wall_correction(fb,pb,target,source,g,i,0,lo);
        AlulaPoint du=wall_correction(fb,pb,target,source,g,i,1,hi);
        for (int j=0;j<=ft;++j) {
            Sample y=bracket(py,pt+1,fy[j]);
            Sample x=combine(lo,hi,fy[j]);
            AlulaPoint q=interpolate_node(pb,source,g,x,y);
            q=add(q,add(scale(dl,1-fy[j]),scale(du,fy[j])));
            int k=fb->node[g]+i*(ft+1)+j;
            nodes[k]=mix(target->nodes[k],q,weight);
        }
    }
    for (int i=0;i<fn;++i) {
        Sample lo=bracket(pc0,pn,.5*(f0[i]+f0[i+1])),hi=bracket(pc1,pn,.5*(f1[i]+f1[i+1]));
        for (int j=0;j<ft;++j) {
            double eta=.5*(fy[j]+fy[j+1]);
            f->packed[index_section(fb,g,i,j)]=interpolate_density(pb,source,g,
                combine(lo,hi,eta),bracket(pcy,pt,eta));
        }
    }
}

int alula_airfoil_grid_transfer(AlulaAirfoilGrid *target,AlulaAirfoilGrid *source,double weight)
{
    if (!target || !source || target==source || target->converged || !source->converged ||
        !(weight>0 && weight<=1) || target->body->size.density<=source->body->size.density ||
        !compatible(target,source)) return 0;
    const AlulaBody *f=target->body,*p=source->body;
    int rows=f->size.rows;
    AlulaPoint *nodes=malloc((size_t)f->size.nodes*sizeof *nodes);
    AlulaPoint *chart=malloc(alula_body_coordinates_checkpoint_count(target->coordinates)*sizeof *chart);
    size_t n=(size_t)4*p->in.nx+2*f->in.nx+2*p->size.tubes+f->size.tubes+6;
    double *work=malloc((2*(size_t)rows+n)*sizeof *work);
    if (!nodes || !chart || !work) {free(nodes);free(chart);free(work);return 0;}
    double *packed=work,*residual=packed+rows,*interpolation=residual+rows;
    memcpy(packed,target->packed,(size_t)rows*sizeof *packed);
    alula_body_coordinates_save(target->coordinates,chart);

    double stagnation[MAX_BODIES],mass[MAX_BODIES+1]={0};
    AlulaBodyState ps,fs;
    int ok=0;
    if (!alula_airfoil_grid_state(source,&ps)) goto done;
    for (int g=0;g<=p->in.bodies;++g) {
        for (int j=p->column[g];j<p->column[g+1];++j) mass[g]+=ps.mass[j];
        if (g<p->in.bodies) {
            AlulaBodySurfaceView s;
            if (!alula_body_coordinates_surface(source->coordinates,source->packed,g,0,&s)) goto done;
            stagnation[g]=s.stagnation;
        }
    }
    if (!alula_body_coordinates_seed_state(target->coordinates,target->packed,stagnation,mass,&ps) ||
        !alula_airfoil_grid_state(target,&fs)) goto done;
    for (int g=0;g<=f->in.bodies;++g)
        transfer_passage(target,source,&fs,&ps,g,weight,nodes,interpolation);
    if (alula_body_coordinates_seed_nodes(target->coordinates,target->packed,nodes)) {
        AlulaAirfoilGridReport quality;
        if (alula_airfoil_grid_current_quality(target,&quality)) {
            if (quality.invalid_cells) {
                ok=alula_airfoil_grid_repair_seed(target,residual) &&
                    alula_airfoil_grid_current_quality(target,&quality) && !quality.invalid_cells;
            } else ok=alula_airfoil_grid_admit_seed(target,residual)==ALULA_BODY_OK;
        }
    }
done:
    if (!ok) {
        alula_body_coordinates_restore(target->coordinates,chart);
        memcpy(target->packed,packed,(size_t)rows*sizeof *packed);
    }
    free(nodes);free(chart);free(work);
    return ok;
}

static double sample_values(const double *x,const double *y,int n,double target)
{
    if (n==1) return y[0];
    Sample q=bracket(x,n,target);
    return lerp(y[q.index],y[q.index+1],q.fraction);
}

static int remap_displacement(AlulaAirfoilGrid *f,AlulaAirfoilGrid *p,const double *from,double *to)
{
    for (int g=0;g<f->elements;++g) {
        AlulaBodyCut coarse=p->cuts[g],fine=f->cuts[g];
        int nc=coarse.trailing-coarse.leading+1,nf=fine.trailing-fine.leading+1;
        for (int upper=0;upper<2;++upper) {
            AlulaBodySurfaceView ps,fs;
            if (!alula_body_coordinates_surface(p->coordinates,p->packed,g,upper,&ps) ||
                !alula_body_coordinates_surface(f->coordinates,f->packed,g,upper,&fs)) return 0;
            for (int i=0;i<nf;++i) to[i]=sample_values(ps.fractions,from,nc,fs.fractions[i]);
            from+=nc;
            to+=nf;
        }
        int nw=p->body->in.nx-coarse.trailing;
        for (int i=fine.trailing+1;i<=f->body->in.nx;++i)
            *to++=sample_values(p->stations+coarse.trailing+1,from,nw,f->stations[i]);
        from+=nw;
    }
    return 1;
}

double alula_airfoil_grid_continue_displacement(AlulaAirfoilGrid *f,AlulaAirfoilGrid *p)
{
    if (!f || !p || f==p || !f->converged || !p->converged || !compatible(f,p)) return 0;
    size_t n=alula_body_coordinates_displacement_count(f->coordinates);
    size_t source_count=alula_body_coordinates_displacement_count(p->coordinates);
    size_t packed_bytes=(size_t)f->body->size.rows*sizeof *f->packed;
    size_t chart_count=alula_body_coordinates_checkpoint_count(f->coordinates);
    double *work=malloc((3*n+source_count)*sizeof *work+packed_bytes);
    AlulaPoint *chart=malloc(chart_count*sizeof *chart);
    double fraction=0;
    if (!work || !chart) goto done;
    double *initial=work,*goal=initial+n,*trial=goal+n,*source=trial+n,*packed=source+source_count;
    if (!alula_body_coordinates_get_displacement(f->coordinates,initial) ||
        !alula_body_coordinates_get_displacement(p->coordinates,source) ||
        !remap_displacement(f,p,source,goal)) goto done;
    double step=.25;
    for (int it=0;it<100 && fraction<1 && step>=1e-5;++it) {
        double next=fmin(1,fraction+step);
        for (size_t i=0;i<n;++i) trial[i]=lerp(initial[i],goal[i],next);
        alula_body_coordinates_save(f->coordinates,chart);
        memcpy(packed,f->packed,packed_bytes);
        AlulaBodySolveResult previous=f->last_solve;
        if (!alula_body_coordinates_set_displacement(f->coordinates,trial)) break;
        AlulaBodySolveResult result=alula_airfoil_grid_solve(f,40,1e-9);
        if (result.status!=ALULA_BODY_SOLVED) alula_airfoil_grid_repair(f,40,1e-9,0,&result);
        if (result.status==ALULA_BODY_SOLVED) {
            fraction=next;
            step=fmin(.5,1.35*step);
        } else {
            alula_body_coordinates_restore(f->coordinates,chart);
            memcpy(f->packed,packed,packed_bytes);
            for (size_t i=0;i<n;++i) trial[i]=lerp(initial[i],goal[i],fraction);
            alula_body_coordinates_set_displacement(f->coordinates,trial);
            f->converged=1;
            f->last_solve=previous;
            step*=.5;
        }
    }
done:
    free(chart);
    free(work);
    return fraction;
}
