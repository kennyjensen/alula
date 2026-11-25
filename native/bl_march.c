// SPDX-License-Identifier: GPL-2.0-or-later
#include "bl.h"

#include <float.h>
#include <math.h>

typedef struct {
    const AlulaBLContext *physics;
    const AlulaBLStation *previous;
    AlulaBLEdge edge;
    AlulaBLRegime regime;
    double trip,exponent,mach_sq,scale;
} March;

static int positive(double x) { return x>0 && isfinite(x); }

static double hk_to_h(double hk,double mach_sq) { return hk*(1+.113*mach_sq)+.29*mach_sq; }

static double edge_mach_sq(const AlulaBLContext *c,double speed)
{
    double energy=c->enthalpy_inverse*speed*speed;
    return energy/(c->gm1*(1-.5*energy));
}

static int initialize(March *m,const AlulaBLContext *c,const AlulaBLStation *previous,
    AlulaBLEdge edge,AlulaBLRegime regime,double trip,double exponent,double theta)
{
    double mach_sq=edge_mach_sq(c,edge.ue);
    if (!(mach_sq>=0) || !isfinite(mach_sq)) return 0;
    *m=(March){c,previous,edge,regime,trip,exponent,mach_sq,1/fmax(5*theta,1e-12)};
    return 1;
}

static int variables(const March *m) { return m->regime>=ALULA_BL_TURBULENT?3:2; }

static AlulaBLStation unpack(const March *m,const double x[3])
{
    double theta=exp(x[0]),hk=1+exp(x[1]),aux=variables(m)==3?exp(x[2]):0;
    return (AlulaBLStation){m->edge.s,theta,m->edge.wake_gap+theta*hk_to_h(hk,m->mach_sq),
        m->edge.ue,aux,m->edge.wake_gap};
}

static AlulaBLStatus residual(const March *m,const double x[3],double r[3],AlulaBLBlock *v)
{
    AlulaBLStation s=unpack(m,x);
    AlulaBLStatus code=alula_bl_interval(m->physics,m->previous,&s,m->regime,m->trip,m->exponent,v);
    if (code==ALULA_BL_OK) {
        r[0]=v->residual[1]; r[1]=v->residual[2];
        r[2]=variables(m)==3?m->scale*v->residual[0]:0;
    }
    return code;
}

static double max_norm(const double r[3],int n)
{
    double v=0;
    for (int i=0;i<n;++i) v=fmax(v,fabs(r[i]));
    return v;
}

static int linear(double a[3][3],const double b[3],int n,double x[3])
{
    for (int i=0;i<n;++i) x[i]=-b[i];
    for (int i=0;i<n;++i) {
        int largest=i;
        for (int j=i+1;j<n;++j) if (fabs(a[j][i])>fabs(a[largest][i])) largest=j;
        if (!(fabs(a[largest][i])>DBL_MIN) || !isfinite(a[largest][i])) return 0;
        if (largest!=i) {
            for (int k=i;k<n;++k) {
                double temp=a[i][k];a[i][k]=a[largest][k];a[largest][k]=temp;
            }
            double temp=x[i];x[i]=x[largest];x[largest]=temp;
        }
        for (int j=i+1;j<n;++j) {
            double f=a[j][i]/a[i][i];
            for (int k=i+1;k<n;++k) a[j][k]-=f*a[i][k];
            x[j]-=f*x[i];
        }
    }
    for (int i=n-1;i>=0;--i) {
        for (int j=i+1;j<n;++j) x[i]-=a[i][j]*x[j];
        x[i]/=a[i][i];
        if (!isfinite(x[i])) return 0;
    }
    return 1;
}

static int jacobian(const March *m,const double x[3],const double r[3],int n,double a[3][3])
{
    for (int j=0;j<n;++j) {
        double h=6e-6*fmax(1,fabs(x[j]));
        double hi[3]={x[0],x[1],x[2]},lo[3]={x[0],x[1],x[2]};
        double rp[3],rm[3];
        hi[j]+=h;lo[j]-=h;
        AlulaBLBlock scratch;
        int plus=residual(m,hi,rp,&scratch)==ALULA_BL_OK;
        int minus=residual(m,lo,rm,&scratch)==ALULA_BL_OK;
        if (!plus && !minus) return 0;
        for (int i=0;i<n;++i)
            a[i][j]=plus&&minus?(rp[i]-rm[i])/(hi[j]-lo[j]):
                plus?(rp[i]-r[i])/(hi[j]-x[j]):(r[i]-rm[i])/(x[j]-lo[j]);
    }
    return 1;
}

static AlulaBLStatus solve(const March *m,double theta,double hk,double aux,
    AlulaBLStation *s,AlulaBLBlock *v,AlulaBLMarchResult *report)
{
    double x[3]={log(theta),log(hk-1),log(fmax(aux,1e-8))},r[3];
    int n=variables(m);
    AlulaBLStatus status=residual(m,x,r,v);
    if (status!=ALULA_BL_OK) return status;
    for (int i=0;i<=40;++i) {
        double norm=max_norm(r,n);
        if (norm<=1e-10) {
            *s=unpack(m,x);
            return ALULA_BL_OK;
        }
        if (i==40) break;
        double matrix[3][3],direction[3]={0};
        if (!jacobian(m,x,r,n,matrix) || !linear(matrix,r,n,direction)) return ALULA_BL_LINEAR;
        double length=max_norm(direction,n);
        double step=length>1?1/length:1;
        int accepted=0;
        for (int trial=0;trial<=20;++trial,step*=.5) {
            double next[3]={x[0]+step*direction[0],x[1]+step*direction[1],x[2]+step*direction[2]};
            double test[3];
            AlulaBLBlock candidate;
            if (residual(m,next,test,&candidate)==ALULA_BL_OK &&
                max_norm(test,n)<(1-1e-4*step)*norm) {
                for (int k=0;k<n;++k) {x[k]=next[k];r[k]=test[k];}
                *v=candidate;
                ++report->iterations;
                accepted=1;
                break;
            }
            ++report->backtracks;
        }
        if (!accepted) break;
    }
    return ALULA_BL_ITERATION_LIMIT;
}

static AlulaBLStatus similarity(const AlulaBLContext *c,AlulaBLEdge edge,double exponent,
    AlulaBLStation *station,AlulaBLBlock *v,AlulaBLMarchResult *report)
{
    double seed=(exponent>.5?.292:.664)*sqrt(edge.s/(c->reynolds*edge.ue));
    if (!positive(seed)) return ALULA_BL_DOMAIN;
    March m;
    if (!initialize(&m,c,NULL,edge,ALULA_BL_SIMILARITY,DBL_MAX,exponent,seed)) return ALULA_BL_DOMAIN;
    return solve(&m,seed,exponent>.5?2.22:2.59,0,station,v,report);
}

static int profile_valid(const AlulaBLEdge *edges,size_t count,int wake)
{
    if (!edges || count==0) return 0;
    for (size_t i=0;i<count;++i)
        if (!positive(edges[i].s) || !positive(edges[i].ue) || !(edges[i].wake_gap>=0) ||
            !isfinite(edges[i].wake_gap) || (!wake && edges[i].wake_gap!=0) ||
            (i && !(edges[i].s>edges[i-1].s))) return 0;
    return 1;
}

static void publish(AlulaBLMarchResult *r,const AlulaBLBlock *v)
{
    r->maximum_residual=fmax(r->maximum_residual,max_norm(v->residual,3));
    ++r->completed;
}

AlulaBLMarchResult alula_bl_march_surface(const AlulaBLContext *c,
    const AlulaBLEdge *edges,size_t count,double trip,double exponent,
    AlulaBLStation *profile,AlulaBLRegime *regimes)
{
    AlulaBLMarchResult r={.status=ALULA_BL_INVALID};
    if (!c || !profile || !regimes || !positive(trip) || !isfinite(exponent) ||
        !profile_valid(edges,count,0)) return r;
    AlulaBLStation virtual_station;
    AlulaBLProperties prior={0};
    const AlulaBLStation *previous=NULL;
    if (trip<=edges[0].s) {
        double distance=.5*trip, speed=distance*edges[0].ue/edges[0].s;
        AlulaBLBlock v;
        r.status=similarity(c,(AlulaBLEdge){distance,speed,0},1,&virtual_station,&v,&r);
        if (r.status!=ALULA_BL_OK) return r;
        prior=v.downstream;
        previous=&virtual_station;
    }
    for (size_t i=0;i<count;++i) {
        AlulaBLBlock value;
        AlulaBLStation current;
        AlulaBLRegime regime=ALULA_BL_SIMILARITY;
        if (!previous) {
            r.status=similarity(c,edges[i],exponent,&current,&value,&r);
            if (r.status!=ALULA_BL_OK) return r;
        } else {
            regime=r.transition.active?ALULA_BL_TURBULENT:ALULA_BL_LAMINAR;
            int forced_here=!r.transition.active && trip>previous->s && trip<=edges[i].s;
            double theta=previous->theta,hk=prior.raw_hk,aux=previous->aux;
            March m;
            if (!initialize(&m,c,previous,edges[i],regime,trip,exponent,theta)) {
                r.status=ALULA_BL_DOMAIN;return r;
            }
            r.status=solve(&m,theta,hk,aux,&current,&value,&r);
            AlulaBLTransition selection={0};
            if (regime==ALULA_BL_LAMINAR && r.status==ALULA_BL_OK) {
                theta=current.theta;hk=value.downstream.raw_hk;
                r.status=alula_bl_check_transition(c,previous,&current,trip,&selection);
                if (r.status!=ALULA_BL_OK) return r;
                if (!selection.active) {
                    current.aux=selection.amplification;
                    r.status=alula_bl_interval(c,previous,&current,regime,trip,exponent,&value);
                }
            }
            if (regime==ALULA_BL_LAMINAR && (selection.active || forced_here)) {
                m.regime=regime=ALULA_BL_TRANSITION;
                r.status=solve(&m,theta,hk,.03,&current,&value,&r);
                if (r.status==ALULA_BL_OK) r.transition=value.transition;
            }
            if (r.status!=ALULA_BL_OK) return r;
        }
        profile[i]=current;
        regimes[i]=regime;
        previous=profile+i;
        prior=value.downstream;
        publish(&r,&value);
    }
    return r;
}

AlulaBLMarchResult alula_bl_march_wake(const AlulaBLContext *c,
    const AlulaBLStation *upper,const AlulaBLStation *lower,
    const AlulaBLEdge *edges,size_t count,AlulaBLStation *profile)
{
    AlulaBLMarchResult r={.status=ALULA_BL_INVALID};
    if (!c || !upper || !lower || !profile || !positive(upper->aux) || !positive(lower->aux) ||
        !profile_valid(edges,count,1)) return r;
    double theta=upper->theta+lower->theta;
    if (!positive(theta)) return r;
    AlulaBLEdge first=edges[0];
    AlulaBLStation start={first.s,theta,upper->delta_star+lower->delta_star+first.wake_gap,first.ue,
        (upper->aux*upper->theta+lower->aux*lower->theta)/theta,first.wake_gap};
    AlulaBLBlock v;
    r.status=alula_bl_trailing_edge(c,upper,lower,&start,first.wake_gap,&v);
    if (r.status!=ALULA_BL_OK) return r;
    profile[0]=start;
    publish(&r,&v);
    for (size_t i=1;i<count;++i) {
        const AlulaBLStation *previous=profile+i-1;
        double hk=v.downstream.raw_hk;
        March m;
        if (!initialize(&m,c,previous,edges[i],ALULA_BL_WAKE,DBL_MAX,1,previous->theta)) {
            r.status=ALULA_BL_DOMAIN;return r;
        }
        AlulaBLStation next;
        r.status=solve(&m,previous->theta,hk,previous->aux,&next,&v,&r);
        if (r.status!=ALULA_BL_OK) return r;
        profile[i]=next;
        publish(&r,&v);
    }
    return r;
}
