// SPDX-License-Identifier: GPL-2.0-or-later
// C port of Alula's integral BL kernels in src/viscous/integral.js and
// src/viscous/xfoil/xblsys.js (Vibefoil JS port by Kenny Jensen).
// Derived from XFOIL xblsys.f, Copyright (C) 2000 Mark Drela.
#include "bl.h"

#include <float.h>
#include <math.h>

typedef struct {
    AlulaBLStation input;
    AlulaBLProperties props;
    double d,hc,us,thickness;
} Station;

static double square(double x) { return x*x; }
static double cube(double x) { return x*x*x; }
static double mean(double a,double b) { return .5*(a+b); }
static double weighted(double a,double b,double f) { return (1-f)*a+f*b; }
static int positive(double v) { return v>0 && isfinite(v); }

int alula_bl_context(AlulaBLContext *ctx,double reynolds,double mach,double gamma,double ncrit)
{
    if (!ctx || !positive(reynolds) || !(mach>=0 && mach<1) || !(gamma>1) ||
        !isfinite(gamma) || !(ncrit>=0) || !isfinite(ncrit)) return 0;
    double gm1=gamma-1,mach_sq=mach*mach,t=1+.5*gm1*mach_sq;
    double hi=gm1*mach_sq/t,tr=1-.5*hi;
    double rho=pow(t,1/gm1),re=reynolds*(tr*sqrt(tr))*1.35/(tr+.35);
    if (!positive(rho) || !positive(re) || !isfinite(hi)) return 0;
    *ctx=(AlulaBLContext){reynolds,mach,gamma,ncrit,gm1,rho,hi,re};
    return 1;
}

static double laminar_shape(double hk)
{
    double t=hk-4.35;
    return hk<4.35 ? (.0111*square(t)-.0278*cube(t))/(hk+1)+1.528-.0002*square(t*hk) :
        .015*square(t)/hk+1.528;
}

static double turbulent_shape(double hk,double rt,double m)
{
    double ho=rt>400 ? 3+400/rt : 4,rz=fmax(rt,200);
    double hs=1.5+4/rz;
    if (hk<ho) {
        double hr=(ho-hk)/(ho-1);
        hs+=(2-hs)*square(hr)*1.5/(hk+.5);
    } else {
        double gr=log(rz),hd=hk-ho;
        hs+=square(hd)*(.007*gr/square(hd+4/gr)+.015/hk);
    }
    return (hs+.028*m)/(1+.014*m);
}

static double laminar_friction(double hk,double rt)
{
    double v=hk<5.5?.0727*cube(5.5-hk)/(hk+1):.015*square(1-1/(hk-4.5));
    return (v-.07)/rt;
}

static double turbulent_friction(double hk,double rt,double m)
{
    double fc=sqrt(1+.2*m),gr=fmax(3,log(rt/fc));
    double c=.3*exp(fmax(-20,-1.33*hk))*pow(gr/2.3026,-1.74-.31*hk);
    return (c+1.1e-4*(tanh(4-hk/.875)-1))/fc;
}

static double laminar_dissipation(double hk,double rt)
{
    double h=hk-4;
    double v=h<0?.00205*pow(-h,5.5):-.0016*square(h)/(1+.02*square(h));
    return (v+.207)/rt;
}

static AlulaBLStatus kinematics(const AlulaBLContext *c,const AlulaBLStation *in,Station *out)
{
    if (!in || !positive(in->s) || !positive(in->theta) || !positive(in->ue) ||
        !isfinite(in->aux) || !(in->wake_gap>=0) || !isfinite(in->wake_gap) ||
        !isfinite(in->delta_star) || !(in->delta_star-in->wake_gap>in->theta))
        return ALULA_BL_INVALID;
    double v2=square(in->ue)*c->enthalpy_inverse,temperature=1-.5*v2;
    if (!positive(temperature)) return ALULA_BL_DOMAIN;
    double mach_sq=v2/(c->gm1*temperature);
    double rho=c->rho_total*pow(1+.5*c->gm1*mach_sq,-1/c->gm1);
    double viscosity=(temperature*sqrt(temperature))*1.35/(temperature+.35)/c->re_total;
    double d=in->delta_star-in->wake_gap,h=d/in->theta,hk=(h-.29*mach_sq)/(1+.113*mach_sq);
    double rt=rho*in->ue*in->theta/viscosity;
    if (!positive(rho) || !positive(viscosity) || !(hk>1) || !isfinite(hk) ||
        !isfinite(mach_sq) || !positive(rt)) return ALULA_BL_DOMAIN;
    *out=(Station){.input=*in,.d=d,.props={.rho=rho,.viscosity=viscosity,
        .mach_sq=mach_sq,.h=h,.raw_hk=hk,.re_theta=rt}};
    return ALULA_BL_OK;
}

static int closure(Station *s,AlulaBLRegime regime)
{
    AlulaBLProperties *p=&s->props;
    int wake=regime==ALULA_BL_WAKE,turbulent=regime==ALULA_BL_TURBULENT;
    double hk=fmax(p->raw_hk,wake?1.00005:1.05),h=p->h,rt=p->re_theta,m=p->mach_sq;
    double hs=wake||turbulent?turbulent_shape(hk,rt,m):laminar_shape(hk);
    double us=.5*hs*(1-(hk-1)/(.75*h));
    if (wake) us=fmin(us,.99995);
    else if (us>.95) us=.98;
    double hc=m*(.064/(hk-.8)+.251);
    double hkc=turbulent?fmax(.01,hk-1-18/rt):hk-1;
    double equilibrium=sqrt((.5/(6.7*6.7*.75))*hs*(hk-1)*square(hkc)/((1-us)*h*square(hk)));
    double cf=0,dissipation=0;
    if (!wake) {
        double lc=laminar_friction(hk,rt);
        if (turbulent) {
            double tc=turbulent_friction(hk,rt,m);
            cf=fmax(lc,tc);
            dissipation=tc*us/hs*(.5+.5*tanh((hk-1)*log(rt)/2.1));
        } else {
            cf=lc;
            dissipation=laminar_dissipation(hk,rt);
        }
    }
    if (wake || turbulent) {
        double outer=.995-us;
        dissipation+=square(s->input.aux)*outer*2/hs + .15*square(outer)/rt*2/hs;
        double reference=wake ? 2.20*square(1-1/hk)/hk/(laminar_shape(hk)*rt) :
            laminar_dissipation(hk,rt);
        dissipation=fmax(dissipation,reference);
        if (wake) dissipation*=2;
    }
    p->hk=hk;
    p->hstar=hs;
    p->cf=cf;
    p->dissipation=dissipation;
    p->equilibrium_shear=equilibrium;
    p->transition_shear=1.8*exp(-3.3/(hk-1))*equilibrium;
    s->us=us;
    s->hc=hc;
    s->thickness=fmin((3.15+1.72/(hk-1))*s->input.theta+s->d,12*s->input.theta);
    return positive(hs) && positive(equilibrium) && isfinite(dissipation) && isfinite(cf) &&
        isfinite(hc) && positive(s->thickness);
}

AlulaBLStatus alula_bl_station(const AlulaBLContext *c,const AlulaBLStation *in,
    AlulaBLRegime regime,AlulaBLProperties *out)
{
    if (!c || !out || regime<ALULA_BL_SIMILARITY || regime>ALULA_BL_TRANSITION)
        return ALULA_BL_INVALID;
    Station s;
    AlulaBLStatus status=kinematics(c,in,&s);
    if (status!=ALULA_BL_OK) return status;
    if (!closure(&s,regime==ALULA_BL_TRANSITION?ALULA_BL_TURBULENT:regime)) return ALULA_BL_DOMAIN;
    *out=s.props;
    return ALULA_BL_OK;
}

/* XFOIL's default amplification correlation, DAMPL (IDAMPV=0). */
static double damping(double hk,double theta,double rt)
{
    double h=1/(hk-1),log_re=log10(rt);
    double critical=2.492*pow(h,.43)+.7*(tanh(14*h-9.24)+1);
    if (log_re<critical-.08) return 0;
    double r=fmin(1,(log_re-critical+.08)/.16);
    double slope=.028*(hk-1)-.0345*exp(-square(3.87*h-2.52));
    double factor=-.05+2.7*h-5.5*square(h)+3*cube(h);
    return (factor*slope/theta)*square(r)*(3-2*r);
}

static double base_rate(const Station *a,const Station *b,double hk_a,double hk_b)
{
    double ga=damping(hk_a,a->input.theta,a->props.re_theta);
    double gb=damping(hk_b,b->input.theta,b->props.re_theta);
    return sqrt(.5*(square(ga)+square(gb)));
}

static double axset(const AlulaBLContext *c,const Station *a,const Station *b,double hk_a,double hk_b,
    double amp_a,double amp_b)
{
    double arg=fmin(20,fmax(0,20*(c->ncrit-mean(amp_a,amp_b))));
    return base_rate(a,b,hk_a,hk_b)+.002*exp(-arg)/(a->input.theta+b->input.theta);
}

static double mid_friction(AlulaBLRegime regime,const AlulaBLProperties *a,const AlulaBLProperties *b)
{
    if (regime==ALULA_BL_WAKE) return 0;
    double hk=mean(a->hk,b->hk),rt=mean(a->re_theta,b->re_theta),m=mean(a->mach_sq,b->mach_sq);
    double laminar=laminar_friction(hk,rt);
    return regime==ALULA_BL_TURBULENT?fmax(laminar,turbulent_friction(hk,rt,m)):laminar;
}

static int differences(const AlulaBLContext *c,const Station *a,const Station *b,
    AlulaBLRegime regime,double exponent,double r[3])
{
    const AlulaBLProperties *p=&a->props,*q=&b->props;
    const AlulaBLStation *l=&a->input,*u=&b->input;
    int similarity=regime==ALULA_BL_SIMILARITY,wake=regime==ALULA_BL_WAKE;
    double dx=u->s-l->s;
    if (!similarity && dx<0) return 0;
    double xlog=similarity?1:log(u->s/l->s);
    double ulog=similarity?exponent:log(u->ue/l->ue);
    double tlog=similarity?.5*(1-exponent):log(u->theta/l->theta);
    double hlog=similarity?0:log(q->hstar/p->hstar);
    double hlogk=log((q->hk-1)/(p->hk-1));
    double upwind=1-.5*exp(-fmin(15,square(hlogk))*(wake?1:5)/square(q->hk));
    if (similarity) r[0]=u->aux;
    else if (regime==ALULA_BL_LAMINAR)
        r[0]=u->aux-l->aux-axset(c,a,b,p->hk,q->hk,l->aux,u->aux)*dx;
    else {
        if (!positive(l->aux) || !positive(u->aux)) return 0;
        double shift=wake?.9:1,sa=weighted(l->aux,u->aux,upwind);
        double cq=weighted(p->equilibrium_shear,q->equilibrium_shear,upwind);
        double cf=weighted(p->cf,q->cf,upwind),hk=weighted(p->hk,q->hk,upwind);
        double rt=mean(p->re_theta,q->re_theta),shape=wake?hk-1:fmax(.01,hk-1-18/rt);
        double hr=shape/(6.7*shift*hk);
        double uq=(.5*cf-square(hr))/(.75*mean(a->d,b->d));
        double scc=5.6*1.333/(1+mean(a->us,b->us)),de=mean(a->thickness,b->thickness);
        r[0]=scc*(cq-sa*shift)*dx-2*de*log(u->aux/l->aux)+2*de*(uq*dx-ulog);
    }
    double ha=mean(p->h,q->h),ma=mean(p->mach_sq,q->mach_sq),hsa=mean(p->hstar,q->hstar);
    double hwa=mean(l->wake_gap/l->theta,u->wake_gap/u->theta);
    double xot1=l->s/l->theta,xot2=u->s/u->theta;
    double cfx=.5*mid_friction(regime,p,q)*mean(l->s,u->s)/mean(l->theta,u->theta)
        +.25*(p->cf*xot1+q->cf*xot2);
    r[1]=tlog+(ha+2-ma+hwa)*ulog-.5*xlog*cfx;
    double dix=weighted(p->dissipation*xot1,q->dissipation*xot2,upwind);
    double cfx2=weighted(p->cf*xot1,q->cf*xot2,upwind);
    r[2]=hlog+(2*mean(a->hc,b->hc)/hsa+1-ha-hwa)*ulog+xlog*(.5*cfx2-dix);
    return isfinite(r[0]) && isfinite(r[1]) && isfinite(r[2]);
}

static AlulaBLStation interpolate(const AlulaBLStation *a,const AlulaBLStation *b,double f)
{
    return (AlulaBLStation){weighted(a->s,b->s,f),weighted(a->theta,b->theta,f),
        weighted(a->delta_star,b->delta_star,f),weighted(a->ue,b->ue,f),0,
        weighted(a->wake_gap,b->wake_gap,f)};
}

static double critical_residual(const AlulaBLContext *c,const Station *a,
    const Station *b,double fraction)
{
    Station at;
    const Station *end=b;
    if (fraction!=1) {
        AlulaBLStation station=interpolate(&a->input,&b->input,fraction);
        if (kinematics(c,&station,&at)!=ALULA_BL_OK) return NAN;
        end=&at;
    }
    double ax=axset(c,a,end,a->props.hk,end->props.raw_hk,a->input.aux,c->ncrit);
    return c->ncrit-a->input.aux-fraction*(b->input.s-a->input.s)*ax;
}

/* The solved amplification is the value extrapolated to the downstream
 * station by TRCHEK. If a trip shortens the growth calculation, use the
 * thermodynamic state at the trip and retain the original full interval
 * in the amplification equation. */
static double amplification(const AlulaBLContext *c,const Station *a,const Station *b,
    double dx,int clip_to_critical,double hk_b)
{
    double base=a->input.aux+dx*base_rate(a,b,a->props.hk,hk_b);
    double add=dx*.002/(a->input.theta+b->input.theta);
    double low=base,high=base+add;
    if (!isfinite(base) || !isfinite(high)) return NAN;
    double value=base;
    for (int i=0;i<60;++i) {
        double u=clip_to_critical?fmin(value,c->ncrit):value;
        double argument=20*(c->ncrit-mean(a->input.aux,u));
        double extra=add*exp(-fmin(20,fmax(0,argument)));
        double error=value-base-extra;
        if (fabs(error)<=2e-14*fmax(1,fabs(value))) return value;
        if (error>0) high=value;
        else low=value;
        double derivative=argument>0 && argument<20 && (!clip_to_critical || value<c->ncrit)?
            1-10*extra:1;
        double next=value-error/derivative;
        value=next>low && next<high ? next : mean(low,high);
    }
    return NAN;
}

static AlulaBLStatus transition(const AlulaBLContext *c,const Station *a,
    const Station *b,double trip,AlulaBLTransition *out)
{
    const AlulaBLStation *l=&a->input,*u=&b->input;
    if (!(u->s>l->s) || !isfinite(trip) || l->wake_gap!=0 || u->wake_gap!=0)
        return ALULA_BL_INVALID;
    double dx=u->s-l->s;
    double tolerance=1e-12*fmax(1,fabs(u->s));
    if (fabs(trip-u->s)<=tolerance) trip=u->s;
    int has_trip=trip>l->s && trip<=u->s;
    double limit=has_trip?(trip-l->s)/dx:1;
    double free_fraction=INFINITY,amp;
    double endpoint=critical_residual(c,a,b,limit);
    if (!isfinite(endpoint)) return ALULA_BL_DOMAIN;
    if (l->aux>=c->ncrit) free_fraction=0;
    else if (endpoint<=0) {
        double low=0,high=limit,rl=c->ncrit-l->aux,rh=endpoint;
        int last=0;
        for (int i=0;i<64;++i) {
            double f=low+(high-low)*rl/(rl-rh);
            if (!(f>low && f<high)) f=mean(low,high);
            double r=critical_residual(c,a,b,f);
            if (!isfinite(r)) return ALULA_BL_DOMAIN;
            if (f==low || f==high || r==0) {low=high=f;break;}
            /* Dampen a repeatedly retained endpoint. Each accepted
             * interpolation still bounds the same amplification root. */
            if (r>0) {
                low=f;rl=r;
                if (last>0) rh*=.5;
                last=1;
            } else {
                high=f;rh=r;
                if (last<0) rl*=.5;
                last=-1;
            }
        }
        free_fraction=mean(low,high);
        if (high-low>16*DBL_EPSILON*fmax(.001,free_fraction)) return ALULA_BL_TRANSITION_ROOT;
    }
    if (isfinite(free_fraction)) {
        amp=free_fraction>0?l->aux+(c->ncrit-l->aux)/free_fraction:
            amplification(c,a,a,dx,1,a->props.raw_hk);
        *out=(AlulaBLTransition){1,0,weighted(l->s,u->s,free_fraction),amp};
    } else if (has_trip) {
        AlulaBLStation s=interpolate(l,u,limit);
        Station at;
        if (kinematics(c,&s,&at)!=ALULA_BL_OK) return ALULA_BL_DOMAIN;
        amp=amplification(c,a,&at,dx,1,at.props.raw_hk);
        *out=(AlulaBLTransition){1,1,trip,amp};
    } else {
        amp=amplification(c,a,b,dx,0,fmax(1.05,b->props.raw_hk));
        *out=(AlulaBLTransition){0,0,u->s,amp};
    }
    return isfinite(amp)?ALULA_BL_OK:ALULA_BL_TRANSITION_ROOT;
}

AlulaBLStatus alula_bl_check_transition(const AlulaBLContext *c,const AlulaBLStation *upstream,
    const AlulaBLStation *downstream,double trip,AlulaBLTransition *out)
{
    if (!c || !out) return ALULA_BL_INVALID;
    Station a,b;
    AlulaBLStatus status=kinematics(c,upstream,&a);
    if (status!=ALULA_BL_OK) return status;
    status=kinematics(c,downstream,&b);
    if (status!=ALULA_BL_OK) return status;
    a.props.hk=fmax(a.props.raw_hk,1.05);
    return transition(c,&a,&b,trip,out);
}

static AlulaBLStatus mixed(const AlulaBLContext *c,const Station *a,Station *b,double trip,AlulaBLBlock *out)
{
    AlulaBLStatus status=transition(c,a,b,trip,&out->transition);
    if (status!=ALULA_BL_OK) return status;
    if (!out->transition.active) return ALULA_BL_NO_TRANSITION;
    double f=(out->transition.s-a->input.s)/(b->input.s-a->input.s);
    AlulaBLStation t=interpolate(&a->input,&b->input,f);
    t.s=out->transition.s;
    t.aux=c->ncrit;
    Station at;
    if (kinematics(c,&t,&at)!=ALULA_BL_OK || !closure(&at,ALULA_BL_LAMINAR)) return ALULA_BL_DOMAIN;
    double left[3];
    if (!differences(c,a,&at,ALULA_BL_LAMINAR,1,left) || !closure(&at,ALULA_BL_TURBULENT))
        return ALULA_BL_DOMAIN;
    at.input.aux=at.props.transition_shear;
    if (!closure(&at,ALULA_BL_TURBULENT) || !closure(b,ALULA_BL_TURBULENT) ||
        !differences(c,&at,b,ALULA_BL_TURBULENT,1,out->residual)) return ALULA_BL_DOMAIN;
    out->residual[1]+=left[1];
    out->residual[2]+=left[2];
    out->downstream=b->props;
    return ALULA_BL_OK;
}

AlulaBLStatus alula_bl_interval(const AlulaBLContext *c,const AlulaBLStation *upstream,
    const AlulaBLStation *downstream,AlulaBLRegime regime,double trip,double exponent,AlulaBLBlock *out)
{
    if (!c || !out || regime<ALULA_BL_SIMILARITY || regime>ALULA_BL_TRANSITION ||
        !isfinite(trip) || !isfinite(exponent)) return ALULA_BL_INVALID;
    Station a,b;
    AlulaBLStatus status=kinematics(c,downstream,&b);
    if (status!=ALULA_BL_OK) return status;
    if (regime!=ALULA_BL_SIMILARITY) {
        status=kinematics(c,upstream,&a);
        if (status!=ALULA_BL_OK) return status;
        if (!(downstream->s>upstream->s)) return ALULA_BL_INVALID;
    }
    *out=(AlulaBLBlock){0};
    if (regime==ALULA_BL_TRANSITION) {
        if (!closure(&a,ALULA_BL_LAMINAR)) return ALULA_BL_DOMAIN;
        return mixed(c,&a,&b,trip,out);
    }
    if (!closure(&b,regime)) return ALULA_BL_DOMAIN;
    if (regime==ALULA_BL_SIMILARITY) a=b;
    else if (!closure(&a,regime)) return ALULA_BL_DOMAIN;
    if (!differences(c,&a,&b,regime,exponent,out->residual)) return ALULA_BL_DOMAIN;
    out->downstream=b.props;
    return ALULA_BL_OK;
}

AlulaBLStatus alula_bl_trailing_edge(const AlulaBLContext *c,const AlulaBLStation *upper,
    const AlulaBLStation *lower,const AlulaBLStation *wake,double gap,AlulaBLBlock *out)
{
    if (!c || !out || !(gap>=0) || !isfinite(gap)) return ALULA_BL_INVALID;
    Station a,b,w;
    AlulaBLStatus status=kinematics(c,upper,&a);
    if (status!=ALULA_BL_OK) return status;
    status=kinematics(c,lower,&b);
    if (status!=ALULA_BL_OK) return status;
    status=kinematics(c,wake,&w);
    if (status!=ALULA_BL_OK) return status;
    if (!closure(&w,ALULA_BL_WAKE)) return ALULA_BL_DOMAIN;
    double theta=upper->theta+lower->theta;
    double aux=(upper->theta*upper->aux+lower->theta*lower->aux)/theta;
    *out=(AlulaBLBlock){.residual={wake->aux-aux,wake->theta-theta,
        wake->delta_star-(upper->delta_star+lower->delta_star+gap)},.downstream=w.props};
    return ALULA_BL_OK;
}
