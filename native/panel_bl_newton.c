// SPDX-License-Identifier: GPL-2.0-or-later
#include "panel_bl.h"
#include "dense.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* A patch contributes three local integral equations. The fourth equation
 * at every station couples signed edge speed to all displacement masses. */
typedef enum {SURFACE, WAKE, JOIN, STAGNATION} Kind;
typedef struct {
    size_t id,body,up[2];
    int parents,side;
    Kind kind;
    AlulaBLRegime regime;
} Patch;
typedef struct {size_t left,right,start,end,wake;} Body;
typedef struct {
    const AlulaDisplacementData *d;
    const AlulaPanelBL *p;
    AlulaBLStation similarity; /* Unit stagnation velocity gradient. */
    double scale;
} System;
typedef struct {
    size_t columns[14],count;
    double value[14][3];
} Derivative;

static double origin(const System *sys,const Body *b,const double *x)
{
    const double *s=sys->d->arc;
    if (b->right>b->left+1) return s[b->left+1];
    double ul=-x[4*b->left+3],ur=x[4*b->right+3];
    return s[b->left]+(s[b->right]-s[b->left])*ul/(ul+ur);
}

static AlulaBLStation station(const System *sys,const Body *b,size_t i,const double *x)
{
    const double *arc=sys->d->arc;
    double s;
    if (i>=b->wake) s=arc[b->end]-origin(sys,b,x)+arc[i];
    else if (b->right==b->left+1 && (i==b->left || i==b->right))
        s=fabs(x[4*i+3])*(arc[b->right]-arc[b->left])/(x[4*b->right+3]-x[4*b->left+3]);
    else s=fabs(arc[i]-origin(sys,b,x));
    return (AlulaBLStation){s,x[4*i+1]*sys->scale,x[4*i+2]*sys->scale,
        fabs(x[4*i+3]),x[4*i],sys->p->stations[i].wake_gap};
}

static double trip_distance(const System *sys,const Body *b,const Patch *p,const double *x)
{
    double sign=p->side?1:-1,stag=origin(sys,b,x);
    double end=sign*(sys->d->arc[p->side?b->end:b->start]-stag);
    double trip=sign*(sys->p->trip_arc[p->body][p->side]-stag);
    return trip>=0?fmin(trip,end):end;
}

/* A frozen active set is used for derivatives. Rebuild it on each nonlinear
 * trial, including the stagnation bracket, material trips and natural e^N. */
static int classify(const System *sys,double *x,const AlulaBLRegime *old,
    Body *bodies,Patch *patch,AlulaBLRegime *active,size_t *order)
{
    const AlulaDisplacementData *d=sys->d;
    size_t next=0;
    for (size_t e=0;e<d->bodies;++e) {
        Body b={.start=d->surface_start[e],.end=d->surface_start[e+1]-1,
            .wake=d->surface_stations+e*d->wake_count};
        b.left=b.start;
        while (b.left<b.end && x[4*(b.left+1)+3]<-1e-13) ++b.left;
        b.right=b.left+1;
        if (b.right<b.end && fabs(x[4*b.right+3])<=1e-13) ++b.right;
        if (b.left==b.start || b.right>=b.end || x[4*b.left+3]>=0 || x[4*b.right+3]<=0) return 0;
        bodies[e]=b;
        for (int side=0;side<2;++side) {
            size_t count=side?b.end-b.right+1:b.left-b.start+1;
            int turbulent=0;
            AlulaBLStation prev={0};
            for (size_t k=0;k<count;++k) {
                size_t id=side?b.right+k:b.left-k;
                Patch p={.id=id,.body=e,.side=side,.kind=SURFACE,.parents=k!=0};
                if (k) p.up[0]=side?id-1:id+1;
                AlulaBLStation st=station(sys,&b,id,x);
                double trip=trip_distance(sys,&b,&p,x);
                if (!(trip>0 && (side?1:-1)*x[4*id+3]>0 && st.theta>0 && st.delta_star>st.theta)) return 0;
                if (turbulent) p.regime=ALULA_BL_TURBULENT;
                else if (k==0) {
                    p.regime=trip<=st.s?ALULA_BL_TRANSITION:ALULA_BL_SIMILARITY;
                    if (p.regime==ALULA_BL_SIMILARITY) x[4*id]=0;
                } else {
                    AlulaBLTransition tr;
                    if (alula_bl_check_transition(&sys->p->physics,&prev,&st,trip,&tr)!=ALULA_BL_OK) return 0;
                    p.regime=tr.active?ALULA_BL_TRANSITION:ALULA_BL_LAMINAR;
                    if (!tr.active) x[4*id]=tr.amplification;
                }
                if (p.regime>=ALULA_BL_TURBULENT) {
                    turbulent=1;
                    if (old[id]<ALULA_BL_TURBULENT) {
                        AlulaBLProperties props;
                        st.aux=.03;
                        if (alula_bl_station(&sys->p->physics,&st,ALULA_BL_TURBULENT,&props)!=ALULA_BL_OK) return 0;
                        x[4*id]=props.transition_shear;
                    }
                    if (!(x[4*id]>0)) return 0;
                }
                st.aux=x[4*id];prev=st;
                patch[id]=p;active[id]=p.regime;order[next++]=id;
            }
        }
        if (b.right!=b.left+1) {
            size_t id=b.left+1;
            patch[id]=(Patch){.kind=STAGNATION,.id=id,.body=e,.parents=2,.up={b.left,b.right}};
            active[id]=ALULA_BL_SIMILARITY;order[next++]=id;
        }
        for (size_t j=0;j<d->wake_count;++j) {
            size_t id=b.wake+j;
            if (!(x[4*id+3]>0 && x[4*id]>0 && x[4*id+1]>0 &&
                (x[4*id+2]-x[4*id+1])*sys->scale>sys->p->stations[id].wake_gap)) return 0;
            patch[id]=(Patch){.id=id,.body=e,.kind=j?WAKE:JOIN,.regime=ALULA_BL_WAKE,
                .parents=j?1:2,.up={j?id-1:b.start,b.end}};
            active[id]=ALULA_BL_WAKE;order[next++]=id;
        }
    }
    return next==d->stations;
}

static int local_residual(const System *sys,const Body *b,const Patch *p,const double *x,double f[3])
{
    size_t id=p->id;
    if (p->kind==STAGNATION) {
        f[0]=x[4*id];
        for (int k=1;k<=2;++k) f[k]=x[4*id+k]-.5*(x[4*b->left+k]+x[4*b->right+k]);
        return 1;
    }
    AlulaBLStation st=station(sys,b,id,x),previous={0};
    AlulaBLBlock block;
    AlulaBLStatus code;
    if (p->parents) previous=station(sys,b,p->up[0],x);
    if (p->kind==JOIN) {
        AlulaBLStation lower=station(sys,b,b->end,x);
        code=alula_bl_trailing_edge(&sys->p->physics,&previous,&lower,&st,sys->p->join_gap[p->body],&block);
    } else {
        double trip=p->kind==SURFACE?trip_distance(sys,b,p,x):DBL_MAX;
        if (!p->parents) {
            double slope=st.ue/st.s;
            if (!(slope>0 && isfinite(slope))) return 0;
            if (p->regime==ALULA_BL_TRANSITION) {
                previous=sys->similarity;
                previous.s=.5*trip;previous.ue=slope*previous.s;
                previous.theta/=sqrt(slope);previous.delta_star/=sqrt(slope);
            } else { st.s=1;st.ue=slope; }
        }
        code=alula_bl_interval(&sys->p->physics,&previous,&st,p->regime,trip,1,&block);
    }
    if (code!=ALULA_BL_OK) return 0;
    f[0]=block.residual[0]*(p->regime>=ALULA_BL_TURBULENT?20:1);
    double scale=p->kind==JOIN?1/sys->scale:1;
    f[1]=block.residual[1]*scale; f[2]=block.residual[2]*scale;
    return 1;
}

typedef struct {double bl,edge,rms;} Norm;
static double maximum(Norm norm) {return fmax(norm.bl,norm.edge);}

static int evaluate(const System *sys,const Body *body,const Patch *patch,const double *x,
    double *mass,double *f,Norm *norm)
{
    size_t n=sys->d->stations;
    *norm=(Norm){0};
    for (size_t i=0;i<n;++i) mass[i]=x[4*i+3]*x[4*i+2]*sys->scale;
    for (size_t i=0;i<n;++i) {
        const double *row=sys->d->influence+i*n;
        double q=sys->d->inviscid[i],*r=f+4*i;
        for (size_t j=0;j<n;++j) q+=row[j]*mass[j];
        r[3]=x[4*i+3]-q;
        if (!local_residual(sys,body+patch[i].body,patch+i,x,r)) return 0;
        for (int k=0;k<4;++k) norm->rms+=r[k]*r[k];
        norm->bl=fmax(norm->bl,fmax(fabs(r[0]),fmax(fabs(r[1]),fabs(r[2]))));
        norm->edge=fmax(norm->edge,fabs(r[3]));
    }
    norm->rms=sqrt(norm->rms/(4*n));
    return isfinite(norm->rms);
}

static void add_column(Derivative *d,size_t col)
{
    for (size_t k=0;k<d->count;++k) if (d->columns[k]==col) return;
    d->columns[d->count++]=col;
}

static int differentiate(const System *sys,const Body *body,const Patch *patch,double *x,
    const double *f,Derivative *derivatives)
{
    for (size_t i=0;i<sys->d->stations;++i) {
        const Patch *p=patch+i;
        const Body *b=body+p->body;
        Derivative *d=derivatives+i;
        d->count=0;
        for (int k=0;k<4;++k) add_column(d,4*i+k);
        for (int j=0;j<p->parents;++j) for (int k=0;k<4;++k) add_column(d,4*p->up[j]+k);
        add_column(d,4*b->left+3);add_column(d,4*b->right+3);
        for (size_t j=0;j<d->count;++j) {
            size_t col=d->columns[j];
            double v=x[col],h=4e-6*fmax(fabs(v),col%4==0?.01:.001),plus[3],minus[3];
            /* Keep finite differences inside the current stagnation bracket. */
            if (col%4==3 && (col/4==b->left || col/4==b->right)) h=fmin(h,.01*fabs(v));
            x[col]=v+h;
            int hi=local_residual(sys,b,p,x,plus);
            x[col]=v-h;
            int lo=local_residual(sys,b,p,x,minus);
            x[col]=v;
            if ((!hi && !lo) || !(h>0)) return 0;
            for (int k=0;k<3;++k)
                d->value[j][k]=hi&&lo?(plus[k]-minus[k])/(2*h):
                    hi?(plus[k]-f[4*i+k])/h:(f[4*i+k]-minus[k])/h;
        }
    }
    return 1;
}

/* A local mass response depends only on this body's upstream surface
 * stations and wake. Trim exact zero columns to keep small flaps and wakes
 * out of the large dense product. */
typedef struct {size_t begin,end;} Span;
static Span support(const double *row,size_t first,size_t last)
{
    while (first<last && row[first]==0) ++first;
    while (last>first && row[last-1]==0) --last;
    return (Span){first,last};
}

/* Retain one unknown per station. In a separated interval a fixed edge
 * speed can make the thickness block singular. Pivot among the four local
 * variables there. Both bracketing speeds must remain free since every
 * surface row can depend on them through the stagnation location. */
static double minor(const Derivative *d,int skip)
{
    int pick[3],count=0;
    for (int k=0;k<4;++k) if (k!=skip) pick[count++]=k;
    const double *a=d->value[pick[0]],*b=d->value[pick[1]],*c=d->value[pick[2]];
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]);
}

static int reduced_step(const System *sys,const Body *body,const Patch *patch,const size_t *order,
    const Derivative *der,const double *x,const double *f,double *dx)
{
    size_t n=sys->d->stations,w=n+1;
    if (n>SIZE_MAX/sizeof(double)/4/w) return 0;
    double *response=calloc(4*n*w,sizeof *response),*work=malloc(3*w*sizeof *work);
    double *mass=malloc(n*w*sizeof *mass),*matrix=malloc(n*n*sizeof *matrix),*rhs=malloc(n*sizeof *rhs);
    Span (*spans)[2]=malloc(n*sizeof *spans);
    int success=0;
    if (!spans || !response || !work || !mass || !matrix || !rhs) goto done;
    for (size_t e=0;e<sys->d->bodies;++e) {
        size_t left=body[e].left,right=body[e].right;
        response[(4*left+3)*w+left]=1;
        response[(4*right+3)*w+right]=1;
    }
    for (size_t l=0;l<n;++l) {
        size_t i=order[l];
        const Derivative *d=der+i;
        const Body *b=body+patch[i].body;
        int keep=3,dependent[3],count=0;
        if (i!=b->left && i!=b->right) {
            double best=0;
            for (int k=3;k>=0;--k) {
                double scale=fmax(fabs(x[4*i+k]),k==0?.05:.1);
                double quality=fabs(minor(d,k))/scale;
                if (quality>best) {keep=k;best=quality;}
            }
        }
        response[(4*i+keep)*w+i]=1;
        for (int k=0;k<4;++k) if (k!=keep) dependent[count++]=k;
        double a[9];
        for (int k=0;k<3;++k) {
            for (int j=0;j<3;++j) a[3*k+j]=d->value[dependent[j]][k];
            double *row=work+k*w;
            memset(row,0,w*sizeof *row);
            row[n]=-f[4*i+k];
            for (size_t j=0;j<d->count;++j) {
                if (j<4 && (int)j!=keep) continue;
                double coeff=d->value[j][k];
                if (coeff==0) continue;
                const double *source=response+d->columns[j]*w;
                for (size_t z=0;z<w;++z) row[z]-=coeff*source[z];
            }
        }
        int pivots[3];
        if (!alula_dense_factor(a,3,pivots) || !alula_dense_apply_rows(a,3,pivots,work,w)) goto done;
        for (int k=0;k<3;++k) memcpy(response+(4*i+dependent[k])*w,work+k*w,w*sizeof *work);
    }
    for (size_t j=0;j<n;++j) {
        const double *delta=response+(4*j+2)*w,*q=delta+w;
        for (size_t k=0;k<w;++k)
            mass[j*w+k]=sys->scale*(x[4*j+3]*delta[k]+x[4*j+2]*q[k]);
        const Body *b=body+patch[j].body;
        spans[j][0]=support(mass+j*w,b->start,b->end+1);
        spans[j][1]=support(mass+j*w,b->wake,b->wake+sys->d->wake_count);
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(n>=192)
#endif
    for (size_t i=0;i<n;++i) {
        double *row=matrix+i*n;
        const double *q=response+(4*i+3)*w;
        memcpy(row,q,n*sizeof *row);
        rhs[i]=-f[4*i+3]-q[n];
        for (size_t j=0;j<n;++j) {
            const double coeff=sys->d->influence[i*n+j],*m=mass+j*w;
            rhs[i]+=coeff*m[n];
            for (int part=0;part<2;++part) {
                Span span=spans[j][part];
                for (size_t k=span.begin;k<span.end;++k) row[k]-=coeff*m[k];
            }
        }
    }
    if (!alula_dense_solve(matrix,rhs,n,1)) goto done;
    for (size_t i=0;i<4*n;++i) {
        const double *row=response+i*w;
        double sum=row[n];
        for (size_t j=0;j<n;++j) sum+=row[j]*rhs[j];
        dx[i]=sum;
    }
    success=1;
done:
    free(response);free(work);free(mass);free(matrix);free(rhs);free(spans);
    return success;
}

static int full_step(const System *sys,const Derivative *der,const double *x,const double *f,double *dx)
{
    size_t n=sys->d->stations,m=4*n;
    if (m>SIZE_MAX/sizeof(double)/m) return 0;
    double *matrix=calloc(m*m,sizeof *matrix);
    if (!matrix) return 0;
    for (size_t i=0;i<n;++i) {
        const Derivative *d=der+i;
        for (int k=0;k<3;++k) for (size_t j=0;j<d->count;++j)
            matrix[(4*i+k)*m+d->columns[j]]=d->value[j][k];
        double *row=matrix+(4*i+3)*m;
        for (size_t j=0;j<n;++j) {
            double factor=sys->d->influence[i*n+j]*sys->scale;
            row[4*j+2]=-factor*x[4*j+3];
            row[4*j+3]=(i==j)-factor*x[4*j+2];
        }
    }
    for (size_t k=0;k<m;++k) dx[k]=-f[k];
    int ok=alula_dense_solve(matrix,dx,m,1);
    free(matrix);
    return ok;
}

/* Check the proposed direction against the original equations, including
 * cross-element interactions, before giving it to the nonlinear search. */
static int linear_balance(const System *sys,const Derivative *der,const double *x,const double *f,const double *dx)
{
    size_t n=sys->d->stations;
    for (size_t i=0;i<n;++i) {
        const Derivative *d=der+i;
        for (int k=0;k<4;++k) {
            double sum=f[4*i+k],scale=fabs(sum);
            if (k==3) {
                for (size_t j=0;j<n;++j) {
                    double factor=sys->d->influence[i*n+j]*sys->scale;
                    double a=-factor*x[4*j+3]*dx[4*j+2];
                    double b=((i==j)-factor*x[4*j+2])*dx[4*j+3];
                    sum+=a+b;scale+=fabs(a)+fabs(b);
                }
            } else for (size_t j=0;j<d->count;++j) {
                double term=d->value[j][k]*dx[d->columns[j]];
                sum+=term;scale+=fabs(term);
            }
            if (!isfinite(sum) || !isfinite(scale) || fabs(sum)>1e-9*scale+DBL_MIN) return 0;
        }
    }
    return 1;
}

static void publish(const System *sys,AlulaPanelBL *p,const Body *body,const Patch *patch,
    const AlulaBLRegime *active,const double *x,Norm norm)
{
    for (size_t e=0;e<sys->d->bodies;++e) {
        p->stagnation[e]=origin(sys,body+e,x);p->upper_first[e]=body[e].left;p->lower_first[e]=body[e].right;
    }
    for (size_t i=0;i<sys->d->stations;++i) {
        p->stations[i]=station(sys,body+patch[i].body,i,x);
        p->regimes[i]=active[i];p->speed[i]=x[4*i+3];
    }
    p->interval_residual=norm.bl;p->outer_residual=norm.edge;
}

AlulaPanelBLNewton alula_panel_bl_solve(const AlulaPanelDisplacement *outer,AlulaPanelBL *p,
    int max_iterations,double tolerance,AlulaPanelBLProgress progress,void *context)
{
    AlulaPanelBLNewton result={.status=ALULA_BL_INVALID,.residual=INFINITY};
    const AlulaDisplacementData *d=alula_panel_displacement_data(outer);
    if (!p || !p->complete || !d || p->count!=d->stations || p->bodies!=d->bodies ||
        max_iterations<0 || !(tolerance>0 && isfinite(tolerance))) return result;
    System sys={.d=d,.p=p,.scale=1/sqrt(p->physics.reynolds)};
    AlulaBLRegime regime;
    AlulaBLEdge unit={.s=1,.ue=1};
    if (alula_bl_march_surface(&p->physics,&unit,1,DBL_MAX,1,&sys.similarity,&regime).status!=ALULA_BL_OK)
        return result;
    size_t n=d->stations;
    double *x=malloc(4*n*sizeof *x),*trial=malloc(4*n*sizeof *trial),*dx=malloc(4*n*sizeof *dx);
    double *f=malloc(4*n*sizeof *f),*next=malloc(4*n*sizeof *next),*mass=malloc(n*sizeof *mass);
    AlulaBLRegime *active=malloc(n*sizeof *active),*changed=malloc(n*sizeof *changed);
    Patch *patch=malloc(n*sizeof *patch),*pending=malloc(n*sizeof *pending);
    Derivative *der=malloc(n*sizeof *der);
    size_t *order=malloc(n*sizeof *order),*next_order=malloc(n*sizeof *next_order);
    Body body[6],candidate[6];
    if (!x || !trial || !f || !next || !dx || !mass || !active || !changed || !patch || !pending ||
        !order || !next_order || !der) goto done;
    for (size_t i=0;i<n;++i) {
        x[4*i]=p->stations[i].aux;x[4*i+1]=p->stations[i].theta/sys.scale;
        x[4*i+2]=p->stations[i].delta_star/sys.scale;x[4*i+3]=p->speed[i];
    }
    Norm norm;
    result.status=ALULA_BL_DOMAIN;
    if (!classify(&sys,x,p->regimes,body,patch,active,order) || !evaluate(&sys,body,patch,x,mass,f,&norm)) goto done;
    for (;;) {
        result.residual=maximum(norm);
        if (progress) {
            publish(&sys,p,body,patch,active,x,norm);
            progress(context,p,result.iterations);
        }
        if (result.residual<=tolerance) {result.converged=1;result.status=ALULA_BL_OK;break;}
        if (result.iterations==max_iterations) {result.status=ALULA_BL_ITERATION_LIMIT;break;}
        if (!differentiate(&sys,body,patch,x,f,der)) break;
        if (!reduced_step(&sys,body,patch,order,der,x,f,dx) || !linear_balance(&sys,der,x,f,dx)) {
            if (!full_step(&sys,der,x,f,dx) || !linear_balance(&sys,der,x,f,dx)) {result.status=ALULA_BL_LINEAR;break;}
        }
        double limit=1;
        for (size_t i=0;i<n;++i) {
            for (int k=1;k<=2;++k) if (dx[4*i+k]<0) limit=fmin(limit,-.5*x[4*i+k]/dx[4*i+k]);
            if (fabs(dx[4*i+3])>.2) limit=fmin(limit,.2/fabs(dx[4*i+3]));
        }
        int accepted=0;
        for (int back=0;back<=20;++back) {
            double step=ldexp(limit,-back);
            for (size_t k=0;k<4*n;++k) trial[k]=x[k]+step*dx[k];
            for (size_t i=0;i<n;++i) {
                double gap=p->stations[i].wake_gap/sys.scale,excess=x[4*i+2]-x[4*i+1]-gap;
                trial[4*i+2]=fmax(trial[4*i+2],trial[4*i+1]+gap+.2*excess);
                if (active[i]>=ALULA_BL_TURBULENT) trial[4*i]=fmax(trial[4*i],.2*x[4*i]);
            }
            Norm updated;
            if (classify(&sys,trial,active,candidate,pending,changed,next_order) &&
                evaluate(&sys,candidate,pending,trial,mass,next,&updated) &&
                (updated.rms<=(1-1e-4*step)*norm.rms ||
                 (step<=limit*.5 && memcmp(active,changed,n*sizeof *active)!=0))) {
                memcpy(x,trial,4*n*sizeof *x);memcpy(f,next,4*n*sizeof *f);
                memcpy(body,candidate,d->bodies*sizeof *body);memcpy(patch,pending,n*sizeof *patch);
                memcpy(active,changed,n*sizeof *active);memcpy(order,next_order,n*sizeof *order);
                norm=updated;++result.iterations;accepted=1;break;
            }
            ++result.backtracks;
        }
        if (!accepted) {result.status=ALULA_BL_LINE_SEARCH;break;}
    }
    if (!progress) publish(&sys,p,body,patch,active,x,norm);
done:
    free(x);free(trial);free(dx);free(f);free(next);free(mass);free(active);free(changed);
    free(patch);free(pending);free(order);free(next_order);free(der);
    return result;
}
