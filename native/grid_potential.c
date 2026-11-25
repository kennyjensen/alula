// SPDX-License-Identifier: GPL-2.0-or-later
#include "grid_potential.h"
#include "panel_internal.h"

#include <stdlib.h>
#include <string.h>

static AlulaPoint mix(AlulaPoint a,AlulaPoint b,double w)
{
    return (AlulaPoint){a.x+w*(b.x-a.x),a.y+w*(b.y-a.y)};
}

static double norm(AlulaPoint p) { return hypot(p.x,p.y); }

static int sample_potential(AlulaGuidePath path,double phi,AlulaPoint *point)
{
    if (path.count<2 || !isfinite(phi) || phi<path.points[0].potential ||
        phi>path.points[path.count-1].potential) return 0;
    size_t lo=0,hi=path.count-1;
    while (hi-lo>1) {
        size_t m=(hi+lo)/2;
        if (path.points[m].potential<phi) lo=m;
        else hi=m;
    }
    double a=path.points[lo].potential,b=path.points[hi].potential;
    if (!(b>a)) return 0;
    *point=mix(path.points[lo].position,path.points[hi].position,(phi-a)/(b-a));
    return isfinite(point->x) && isfinite(point->y);
}

static double surface_phi(const AlulaPanelField *field,const AlulaPanelBodyGuide *profile,
    double a,double b,double chord)
{
    /* Integrate the tangential speed just outside the material spline. */
    double p=.5*(a+b), h=.5*fabs(b-a);
    double total=0;
    for (int i=-1;i<=1;i+=2) {
        AlulaPoint q,t,v;
        alula_contour_spline_value(profile->curve,p+i*.5773502691896258*h,&q,&t,NULL);
        double mag=norm(t);
        if (!(mag>0)) return NAN;
        q.x+=1e-5*chord*t.y/mag;
        q.y-=1e-5*chord*t.x/mag;
        if (!alula_panel_field_sample(field,q,&v,NULL)) return NAN;
        total+=fabs(dot(t,v));
    }
    return h*total;
}

int alula_grid_guide_potential(const AlulaPanelField *field,
    const AlulaPanelBodyGuide *profiles,const AlulaAirfoilGrid *grid,const double *x,int nx,int count,double *phi)
{
    int rows=nx+1;
    for (int g=0;g<count;++g) {
        const AlulaPanelBodyGuide *b=profiles+g;
        AlulaBodyCut cut=grid->cuts[g];
        double *p=phi+(g+1)*rows;
        double inlet=b->incoming.points[b->incoming.count-1].potential;
        double to_stagnation=.5*b->strain*b->seed_distance*b->seed_distance;
        for (int i=0;i<cut.leading;++i) {
            double potential=p[i];
            if (x[i]>b->incoming.points[0].position.x) {
                double span=b->stagnation_point.x-b->incoming.points[0].position.x;
                double f=(b->stagnation_point.x-x[i])/span;
                potential=to_stagnation*(1-f*f);
            }
            if (!isfinite(potential)) return 0;
            p[i]=potential-inlet;
        }
        p[cut.leading]=to_stagnation-inlet;
        AlulaPoint q;
        if (!alula_panel_field_sample(field,b->wake.points[0].position,&q,NULL)) return 0;
        double wake_phi=norm(q)*norm(subtract(b->wake.points[0].position,b->trailing_edge));
        for (int i=cut.trailing+1;i<=nx;++i) {
            if (x[i]<b->wake.points[0].position.x) {
                p[i]=wake_phi*(x[i]-b->trailing_edge.x)/(b->wake.points[0].position.x-b->trailing_edge.x);
            } else {
                if (!isfinite(p[i])) return 0;
                p[i]=wake_phi+p[i];
            }
        }
    }
    return 1;
}

static int potential_boundaries(const AlulaPanelField *field,const AlulaPanelBodyGuide *profiles,
    const AlulaAirfoilGrid *grid,const double *guide_phi,int nx,int count,
    const double *fractions,double chord,double *phi)
{
    int rows=nx+1;
    memcpy(phi,guide_phi,(size_t)rows*sizeof *phi);
    memcpy(phi+(2*count+1)*rows,guide_phi+(count+1)*rows,(size_t)rows*sizeof *phi);
    for (int g=0;g<count;++g) {
        const AlulaPanelBodyGuide *b=profiles+g;
        AlulaBodyCut cut=grid->cuts[g];
        double *l=phi+(2*g+1)*rows,*u=l+rows;
        const double *p=guide_phi+(g+1)*rows;
        size_t incoming=(size_t)(cut.leading+1)*sizeof *phi;
        memcpy(l,p,incoming);
        memcpy(u,p,incoming);
        double end=b->curve->arc[b->curve->count-1],s=b->stagnation;
        const double *f=fractions+(size_t)g*4*rows;
        for (int i=cut.leading+1;i<=cut.trailing;++i) {
            double a=f[i-cut.leading-1],z=f[i-cut.leading];
            l[i]=l[i-1]+surface_phi(field,b,s+a*(end-s),s+z*(end-s),chord);
            u[i]=u[i-1]+surface_phi(field,b,s*(1-a),s*(1-z),chord);
        }
        for (int i=cut.trailing+1;i<=nx;++i) {
            l[i]=l[cut.trailing]+p[i];
            u[i]=u[cut.trailing]+p[i];
        }
    }
    for (int b=0;b<2*count+2;++b)
        for (int i=1;i<rows;++i) if (!(phi[b*rows+i]>phi[b*rows+i-1])) return 0;
    return 1;
}

typedef struct { int passage,tube; double fraction; } InteriorLine;

static int place_line(AlulaGuidePath path,double eta,const double *lower,const double *upper,
    AlulaPoint *nodes,int nx,int stride)
{
    double range=path.points[path.count-1].potential;
    /* Match relative potential along each bank while preserving the
     * traced streamline's inlet and outlet on the physical planes. */
    double extent=(1-eta)*(lower[nx]-lower[0])+eta*(upper[nx]-upper[0]);
    for (int i=0;i<=nx;++i) {
        double here=i==0?0:i==nx?range:range*((1-eta)*(lower[i]-lower[0])+eta*(upper[i]-upper[0]))/extent;
        if (!sample_potential(path,here,nodes+(size_t)i*stride)) return 0;
    }
    return 1;
}

int alula_fill_potential_grid(const AlulaPanelField *field,const AlulaPanelGuides *guides,
    const AlulaPanelBodyGuide *profiles,AlulaAirfoilGrid *grid,const AlulaPoint *boundaries,
    const double *guide_phi,const double *fractions,const double *weights,int nx,int count,double chord,
    AlulaPanelStreamline **const *lines)
{
    int rows=nx+1,ok=0;
    double *phi=malloc((size_t)(2*count+2)*rows*sizeof *phi);
    int capacity=alula_airfoil_grid_size(grid).tubes-count-1;
    InteriorLine *interior=capacity?malloc((size_t)capacity*sizeof *interior):NULL;
    if (!phi || (capacity && !interior) ||
        !potential_boundaries(field,profiles,grid,guide_phi,nx,count,fractions,chord,phi))
        goto done;
    const double *levels=alula_panel_guides_levels(guides);
    size_t offset[7]={0};
    int used=0,uncached=0;
    for (int g=0;g<=count;++g) {
        int tubes=grid->tubes[g];
        const AlulaPoint *a=boundaries+2*g*rows,*b=a+rows;
        double total=0,mass=0;
        for (int j=0;j<tubes;++j) total+=weights[j];
        for (int i=0;i<rows;++i) {
            grid->nodes[offset[g]+i*(tubes+1)]=a[i];
            grid->nodes[offset[g]+i*(tubes+1)+tubes]=b[i];
        }
        for (int j=1;j<tubes;++j) {
            mass+=weights[j-1];
            interior[used++]=(InteriorLine){g,j,mass/total};
            uncached+=lines[g][j]==NULL;
        }
        if (g<count) offset[g+1]=offset[g]+(size_t)rows*(tubes+1);
        weights+=tubes;
    }
    int bad=0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(|:bad) if(uncached>=4)
#endif
    for (int k=0;k<used;++k) {
        if (bad) continue;
        InteriorLine spec=interior[k];
        int g=spec.passage,j=spec.tube;
        double eta=spec.fraction;
        const double *lower=phi+2*g*rows;
        if (!lines[g][j])
            lines[g][j]=alula_panel_guides_trace_level(guides,levels[g]+eta*(levels[g+1]-levels[g]),NULL);
        if (!lines[g][j] || !place_line(alula_panel_streamline_path(lines[g][j]),eta,lower,lower+rows,
            grid->nodes+offset[g]+j,nx,grid->tubes[g]+1)) bad=1;
    }
    ok=!bad;
done:
    free(interior);
    free(phi);
    return ok;
}
