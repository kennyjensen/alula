// SPDX-License-Identifier: GPL-2.0-or-later
#include "bl.h"

#include <stdio.h>
#include <stdlib.h>

static int station(AlulaBLStation *s)
{
    return scanf("%lf%lf%lf%lf%lf%lf",&s->s,&s->theta,&s->delta_star,&s->ue,&s->aux,&s->wake_gap)==6;
}

static void block(AlulaBLStatus status,const AlulaBLBlock *v)
{
    printf("%d",(int)status);
    if (status==ALULA_BL_OK) {
        const AlulaBLProperties *p=&v->downstream;
        const AlulaBLTransition *t=&v->transition;
        printf(" %.17g %.17g %.17g %d %d %.17g %.17g",v->residual[0],v->residual[1],v->residual[2],
            t->active,t->forced,t->s,t->amplification);
        printf(" %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g",
            p->rho,p->viscosity,p->mach_sq,p->h,p->raw_hk,p->hk,p->hstar,p->re_theta,p->cf,
            p->dissipation,p->equilibrium_shear,p->transition_shear);
    }
    putchar('\n');
}

static int march(const AlulaBLContext *c,int wake,double trip,double exponent)
{
    size_t count;
    AlulaBLStation upper,lower;
    if (wake && (!station(&upper) || !station(&lower))) return 0;
    if (scanf("%zu",&count)!=1 || count==0 || count>1000) return 0;
    AlulaBLEdge *edges=calloc(count,sizeof *edges);
    AlulaBLRegime *regimes=calloc(count,sizeof *regimes);
    AlulaBLStation *states=calloc(count,sizeof *states);
    int success=0;
    if (!edges || !regimes || !states) goto done;
    for (size_t i=0;i<count;++i)
        if (scanf("%lf%lf%lf",&edges[i].s,&edges[i].ue,&edges[i].wake_gap)!=3) goto done;
    AlulaBLMarchResult r=wake?alula_bl_march_wake(c,&upper,&lower,edges,count,states):
        alula_bl_march_surface(c,edges,count,trip,exponent,states,regimes);
    printf("%d %zu %d %d %.17g %d %d %.17g %.17g\n",(int)r.status,r.completed,r.iterations,r.backtracks,
        r.maximum_residual,r.transition.active,r.transition.forced,r.transition.s,r.transition.amplification);
    for (size_t i=0;i<r.completed;++i) {
        AlulaBLStation s=states[i];
        printf("%d %.17g %.17g %.17g %.17g %.17g %.17g\n",(int)(wake?ALULA_BL_WAKE:regimes[i]),
            s.s,s.theta,s.delta_star,s.ue,s.aux,s.wake_gap);
    }
    success=1;
done:
    free(states);free(regimes);free(edges);
    return success;
}

int main(void)
{
    int op;
    double re,mach,gamma,ncrit,trip,exponent;
    while (scanf("%d%lf%lf%lf%lf%lf%lf",&op,&re,&mach,&gamma,&ncrit,&trip,&exponent)==7) {
        AlulaBLContext c;
        AlulaBLStation a,b,w;
        if (!alula_bl_context(&c,re,mach,gamma,ncrit)) return 1;
        if (op==7 || op==8) {
            if (!march(&c,op==8,trip,exponent)) return 1;
            continue;
        }
        if (!station(&a) || !station(&b)) return 1;
        AlulaBLBlock result={0};
        AlulaBLStatus s;
        if (op==5) s=alula_bl_check_transition(&c,&a,&b,trip,&result.transition);
        else if (op==6) {
            if (!station(&w)) return 1;
            s=alula_bl_trailing_edge(&c,&a,&b,&w,trip,&result);
        } else s=alula_bl_interval(&c,&a,&b,(AlulaBLRegime)op,trip,exponent,&result);
        block(s,&result);
    }
    return feof(stdin)?0:1;
}
