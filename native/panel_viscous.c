// SPDX-License-Identifier: GPL-2.0-or-later
#include "panel_viscous.h"
#include "panel_internal.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

typedef struct {
    AlulaPanelSnapshot state;
    double chord[6];
    AlulaPanelProgress progress;
    void *context;
} Run;

static void loads(Run *run,const AlulaPanelBL *p,int iteration)
{
    AlulaPanelSnapshot *s=&run->state;
    const AlulaDisplacementData *d=alula_panel_displacement_data(s->outer);
    s->iteration=iteration;s->profile=p;
    s->equation_residual=fmax(p->interval_residual,p->outer_residual);
    double cx=0,cy=0,cm=0,drag=0,maximum=0;
    for (size_t e=0;e<d->bodies;++e) {
        size_t start=d->surface_start[e],end=d->surface_start[e+1]-1;
        for (size_t i=start;i<=end;++i) {
            size_t j=i==end?start:i+1;
            AlulaPoint a=d->points[i],b=d->points[j];
            double dx=b.x-a.x,dy=b.y-a.y;
            double cp_a=1-p->speed[i]*p->speed[i],cp_b=1-p->speed[j]*p->speed[j];
            double avg=.5*(cp_a+cp_b);
            /* Uniform base pressure integrates to the same force and moment
             * along every retained base segment as along its endpoint chord. */
            double change=i==end?0:cp_b-cp_a;
            cx-=avg*dy;cy+=avg*dx;
            cm-=avg*(dx*(.5*(a.x+b.x)-.25)+dy*.5*(a.y+b.y))+change*(dx*dx+dy*dy)/12;
            maximum=fmax(maximum,p->stations[i].delta_star/run->chord[e]);
        }
        size_t w=d->surface_stations+e*d->wake_count;
        for (size_t j=0;j<d->wake_count;++j)
            maximum=fmax(maximum,p->stations[w+j].delta_star/run->chord[e]);
        AlulaBLStation exit=p->stations[w+d->wake_count-1];
        double h=exit.delta_star/exit.theta;
        drag+=exit.wake_gap>0?NAN:2*exit.theta*pow(exit.ue,2+.5*(h+1));
    }
    double angle=s->alpha*ALULA_PI/180,cs=cos(angle),sn=sin(angle);
    s->cl=cy*cs-cx*sn;s->pressure_drag=cx*cs+cy*sn;
    s->cm=cm;s->wake_drag=drag;s->max_displacement=maximum;
}

static void accepted(void *context,const AlulaPanelBL *p,int iteration)
{
    Run *run=context;
    loads(run,p,iteration);
    run->progress(run->context,&run->state);
}

AlulaPanelViscousOptions alula_panel_viscous_defaults(void)
{
    return (AlulaPanelViscousOptions){.wake_count=24,.wake_length=2,.reynolds=1e6,.ncrit=9,
        .newton_limit=60,.wake_limit=40,.equation_tolerance=1e-8,.wake_tolerance=1e-6,.wake_relaxation=.8};
}

AlulaPanelViscousResult alula_panel_viscous_solve(const AlulaContour *contours,size_t count,
    const AlulaPanelViscousOptions *o,AlulaPanelProgress progress,void *context)
{
    AlulaPanelViscousResult r={.status=ALULA_VISCOUS_INPUT};
    if (!o || !count || count>6 || !contours || !(fabs(o->alpha)<=20) || !isfinite(o->reynolds) ||
        !(o->reynolds>0) || !(o->ncrit>=0 && isfinite(o->ncrit)) || o->wake_count<4 ||
        !(o->wake_length>0 && isfinite(o->wake_length)) || o->newton_limit<0 || o->wake_limit<0 ||
        !(o->equation_tolerance>0 && isfinite(o->equation_tolerance)) ||
        !(o->wake_tolerance>0 && isfinite(o->wake_tolerance)) || !(o->wake_relaxation>0 && o->wake_relaxation<=1))
        return r;
    for (size_t e=0;o->trips && e<count;++e) for (int side=0;side<2;++side)
        if (!(o->trips[e][side]>0 && o->trips[e][side]<=1)) return r;
    size_t nw=o->wake_count;
    AlulaPanelDisplacement *outer=alula_panel_displacement_create(contours,count,o->alpha,nw,o->wake_length,
        ALULA_WAKE_INVISCID,NULL);
    if (!outer) {r.status=ALULA_VISCOUS_GEOMETRY;return r;}
    const AlulaDisplacementData *d=alula_panel_displacement_data(outer);
    AlulaPanelBL *p=alula_panel_bl_initialize(outer,o->reynolds,o->ncrit,o->trips);
    Run run={.progress=progress,.context=context,.state={.outer=outer,.profile=p,.alpha=o->alpha,.target_alpha=o->alpha,
        .equation_residual=INFINITY,.wake_residual=INFINITY,.cl=NAN,.cm=NAN,.pressure_drag=NAN,.wake_drag=NAN}};
    r.status=ALULA_VISCOUS_SEED;
    if (!p) {r.bl_status=ALULA_BL_INVALID;r.final=run.state;return r;}
    size_t segments=count*(nw-1);
    double *lengths=malloc(segments*sizeof *lengths),*mass=malloc(d->stations*sizeof *mass);
    AlulaPoint *path=malloc(count*nw*sizeof *path),*vel=malloc(segments*sizeof *vel);
    if (!lengths || !mass || !path || !vel) goto done;
    for (size_t e=0;e<count;++e) {
        size_t start=d->surface_start[e],end=d->surface_start[e+1]-1,w=d->surface_stations+e*nw;
        AlulaPoint a=d->points[start],b=d->points[end],te={.5*(a.x+b.x),.5*(a.y+b.y)};
        for (size_t i=start;i<=end;++i) run.chord[e]=fmax(run.chord[e],hypot(d->points[i].x-te.x,d->points[i].y-te.y));
        for (size_t j=0;j+1<nw;++j) lengths[e*(nw-1)+j]=d->arc[w+j+1]-d->arc[w+j];
    }
    /* Large edge mismatch means the shared-field march can start over one
     * freestream speed from the displacement solution. Try the lower-load
     * seed on the very same wake, and ramp only when that seed improves it. */
    if (o->newton_limit>0 && fabs(o->alpha)>2 && (!p->complete || p->outer_residual>1)) {
        alula_panel_displacement_set_incidence(outer,0);
        AlulaPanelBL *lower=alula_panel_bl_initialize(outer,o->reynolds,o->ncrit,o->trips);
        if (lower && lower->complete && (!p->complete || lower->outer_residual<p->outer_residual)) {
            int ready=1;
            int stages=(int)ceil(fabs(o->alpha)/2);
            for (int stage=0;stage<stages;++stage) {
                double angle=o->alpha*((double)stage/stages);
                alula_panel_displacement_set_incidence(outer,angle);
                run.state.alpha=angle;
                AlulaPanelBLNewton startup=alula_panel_bl_solve(outer,lower,o->newton_limit,o->equation_tolerance,
                    progress?accepted:NULL,&run);
                if (!startup.converged) {ready=0;break;}
            }
            if (ready) {
                alula_panel_bl_free(p);
                p=lower;lower=NULL;
            }
        }
        alula_panel_bl_free(lower);
        alula_panel_displacement_set_incidence(outer,o->alpha);
        run.state.alpha=o->alpha;run.state.profile=p;
    }
    if (!p->complete) {r.bl_status=p->status;goto done;}
    for (int iter=0;;++iter) {
        run.state.wake_iteration=iter;run.state.wake_residual=INFINITY;
        AlulaPanelBLNewton step=alula_panel_bl_solve(outer,p,o->newton_limit,o->equation_tolerance,
            progress?accepted:NULL,&run);
        r.bl_status=step.status;
        loads(&run,p,step.iterations);
        run.state.equation_residual=step.residual;
        if (!step.converged) {r.status=ALULA_VISCOUS_NEWTON;break;}
        for (size_t i=0;i<d->stations;++i) mass[i]=p->speed[i]*p->stations[i].delta_star;
        if (!alula_panel_displacement_wake_velocity(outer,mass,d->stations,vel)) {r.status=ALULA_VISCOUS_GEOMETRY;break;}
        double residual=0;
        for (size_t e=0;e<count;++e) {
            const AlulaPoint *wake=d->points+d->surface_stations+e*nw;
            AlulaPoint *next=path+e*nw;
            next[0]=wake[0];
            for (size_t j=0;j+1<nw;++j) {
                size_t id=e*(nw-1)+j;
                AlulaPoint t=subtract(wake[j+1],wake[j]),v=vel[id];
                double h=hypot(t.x,t.y),speed=hypot(v.x,v.y);
                t.x/=h;t.y/=h;
                if (!(dot(v,t)>0)) {r.status=ALULA_VISCOUS_WAKE_REVERSE;goto done;}
                residual=fmax(residual,fabs(cross(t,v)));
                double tx=t.x+o->wake_relaxation*(v.x/speed-t.x),ty=t.y+o->wake_relaxation*(v.y/speed-t.y);
                double scale=lengths[id]/hypot(tx,ty);
                next[j+1]=(AlulaPoint){next[j].x+scale*tx,next[j].y+scale*ty};
            }
        }
        run.state.wake_residual=residual;
        if (progress) progress(context,&run.state);
        if (residual<=o->wake_tolerance) {
            r.status=ALULA_VISCOUS_OK;
            if (run.state.max_displacement>.1 || !isfinite(run.state.wake_drag)) {
                r.status=ALULA_VISCOUS_OUTSIDE_MODEL;
                run.state.cl=run.state.cm=run.state.pressure_drag=run.state.wake_drag=NAN;
            }
            break;
        }
        if (iter==o->wake_limit) {r.status=ALULA_VISCOUS_WAKE_LIMIT;break;}
        AlulaPanelDisplacement *updated=alula_panel_displacement_create(contours,count,o->alpha,nw,o->wake_length,
            ALULA_WAKE_INVISCID,path);
        if (!updated) {r.status=ALULA_VISCOUS_GEOMETRY;break;}
        alula_panel_displacement_free(outer);
        run.state.outer=outer=updated;d=alula_panel_displacement_data(outer);
        for (size_t e=0;e<count;++e) for (size_t j=0;j<nw;++j) {
            size_t i=d->surface_stations+e*nw+j;
            double gap=alula_bl_dead_air(p->fluid_gap[e],p->tangent_cross[e],d->arc[i]);
            p->stations[i].delta_star+=gap-p->stations[i].wake_gap;
            p->stations[i].wake_gap=gap;
        }
    }
done:
    free(path);free(vel);free(lengths);free(mass);r.final=run.state;
    return r;
}

void alula_panel_viscous_free(AlulaPanelViscousResult *r)
{
    if (!r) return;
    alula_panel_displacement_free((AlulaPanelDisplacement *)r->final.outer);
    alula_panel_bl_free((AlulaPanelBL *)r->final.profile);
    *r=(AlulaPanelViscousResult){0};
}
