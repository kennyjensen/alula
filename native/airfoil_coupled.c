// SPDX-License-Identifier: GPL-2.0-or-later
#include "airfoil_bl_internal.h"
#include "airfoil_grid_internal.h"

#include <stdlib.h>
#include <string.h>

static int displacement(const AlulaBody *b,const AlulaAirfoilBL *p,double *out,double *residual,double *edge)
{
    if (!alula_airfoil_bl_complete(p)) return 0;
    double chord=b->in.length_scale;
    *residual=0;
    *edge=0;
    for (int i=0;i<b->in.bodies;++i) {
        const AlulaAirfoilBLBranch *upper=alula_airfoil_bl_branch(p,i,ALULA_BL_UPPER);
        const AlulaAirfoilBLBranch *lower=alula_airfoil_bl_branch(p,i,ALULA_BL_LOWER);
        const AlulaAirfoilBLBranch *wake=alula_airfoil_bl_branch(p,i,ALULA_BL_MERGED_WAKE);
        double leading=.5*chord*(upper->stations[0].delta_star+lower->stations[0].delta_star);
        const AlulaAirfoilBLBranch *both[]={lower,upper};
        for (int side=0;side<2;++side) {
            const AlulaAirfoilBLBranch *branch=both[side];
            *out++=leading;
            for (size_t j=0;j<branch->count;++j) *out++=chord*branch->stations[j].delta_star;
            *residual=fmax(*residual,branch->result.maximum_residual);
            *edge=fmax(*edge,branch->result.edge_residual);
        }
        for (size_t j=1;j<wake->count;++j) *out++=chord*wake->stations[j].delta_star;
        *residual=fmax(*residual,wake->result.maximum_residual);
        *edge=fmax(*edge,wake->result.edge_residual);
    }
    return 1;
}

static double merit(double displacement_residual,double edge_residual)
{
    /* Balance speed and thickness changes for the line search. The final
     * convergence check tests each mismatch against its own threshold. */
    return fmax(displacement_residual,.1*edge_residual);
}

static double mismatch(const double *target,const double *applied,double *difference,size_t n,double chord)
{
    double r=0;
    for (size_t i=0;i<n;++i) {
        difference[i]=(target[i]-applied[i])/chord;
        r=fmax(r,fabs(difference[i]));
    }
    return r;
}

static AlulaAirfoilBL *march_interaction(AlulaAirfoilGrid *grid,const AlulaCoupledControls *controls,
    const double *applied,int *local)
{
    AlulaAirfoilBL *p=alula_airfoil_bl_create_interaction(grid,controls->reynolds,controls->ncrit,
        controls->trips,applied,*local);
    if (applied && !*local && !alula_airfoil_bl_complete(p)) {
        alula_airfoil_bl_free(p);
        *local=1;
        p=alula_airfoil_bl_create_interaction(grid,controls->reynolds,controls->ncrit,controls->trips,applied,1);
    }
    return p;
}

static double aitken(const double *last,const double *next,size_t count,double accepted)
{
    double numerator=0,denominator=0;
    for (size_t i=0;i<count;++i) {
        double change=next[i]-last[i];
        numerator+=last[i]*change;
        denominator+=change*change;
    }
    double proposal=denominator>0?-accepted*numerator/denominator:accepted;
    return isfinite(proposal) && proposal>0?fmax(.05,fmin(1.25,proposal)):fmin(1,accepted);
}

enum { ANDERSON_DEPTH=4, ANDERSON_STARTUP_UPDATES=5 };

static double dot_values(const double *a,const double *b,size_t n)
{
    double sum=0;
    for (size_t i=0;i<n;++i) sum+=a[i]*b[i];
    return sum;
}

static int anderson_direction(const double *defect,const double *dx,const double *df,
    double *basis,double *direction,size_t n,int history)
{
    /* Solve min ||defect - df * weights|| with a short QR factorization,
     * then apply the corresponding type-II Anderson correction. */
    double upper[ANDERSON_DEPTH][ANDERSON_DEPTH]={0},weights[ANDERSON_DEPTH]={0};
    double cutoff=1e-9*sqrt(dot_values(defect,defect,n));
    for (int j=0;j<history;++j) {
        double *q=basis+j*n;
        memcpy(q,df+j*n,n*sizeof *q);
        for (int k=0;k<j;++k) {
            const double *previous=basis+k*n;
            upper[k][j]=dot_values(q,previous,n);
            for (size_t i=0;i<n;++i) q[i]-=upper[k][j]*previous[i];
        }
        double norm=sqrt(dot_values(q,q,n));
        if (!(norm>cutoff)) return 0;
        upper[j][j]=norm;
        for (size_t i=0;i<n;++i) q[i]/=norm;
        weights[j]=dot_values(q,defect,n);
    }
    memcpy(direction,defect,n*sizeof *direction);
    for (int j=history-1;j>=0;--j) {
        for (int k=j+1;k<history;++k) weights[j]-=upper[j][k]*weights[k];
        weights[j]/=upper[j][j];
        for (size_t i=0;i<n;++i) direction[i]-=weights[j]*(dx[j*n+i]+df[j*n+i]);
    }
    return history>0;
}

AlulaCoupledResult alula_airfoil_coupled_solve(AlulaAirfoilGrid *g,const AlulaCoupledControls *c)
{
    AlulaCoupledResult r={.status=ALULA_COUPLED_INVALID,.displacement_residual=INFINITY,
        .bl_residual=INFINITY,.edge_residual=INFINITY,
        .euler={.status=ALULA_BODY_BAD_STATE,.norm=INFINITY,.max_mach=NAN}};
    if (!g || !g->converged || !c || !positive(c->reynolds) || !(c->ncrit>=0) || !isfinite(c->ncrit) ||
        !positive(c->euler_tolerance) || !positive(c->displacement_tolerance) ||
        !(c->displacement_scale>=0) || !isfinite(c->displacement_scale) ||
        !(c->relaxation>0 && c->relaxation<=1) || c->max_iterations<0 || c->euler_iterations<1) return r;
    const AlulaBody *b=g->body;
    size_t n=alula_body_coordinates_displacement_count(g->coordinates),rows=(size_t)b->size.rows;
    if (!n) return r;
    r.euler=g->last_solve;
    double chord=b->in.length_scale;
    double *work=malloc(((6+3*ANDERSON_DEPTH)*n+rows)*sizeof *work);
    AlulaPoint *checkpoint=malloc(alula_body_coordinates_checkpoint_count(g->coordinates)*sizeof *checkpoint);
    if (!work || !checkpoint) {
        r.status=ALULA_COUPLED_ALLOCATION;goto done;
    }
    double *applied=work,*target=applied+n,*trial=target+n,*old_defect=trial+n,*new_defect=old_defect+n;
    double *packed=new_defect+n;
    double *direction=packed+rows,*dx=direction+n,*df=dx+ANDERSON_DEPTH*n,*basis=df+ANDERSON_DEPTH*n;
    if (!alula_body_coordinates_get_displacement(g->coordinates,applied)) goto done;
    r.profiles=alula_airfoil_bl_create(g,c->reynolds,c->ncrit,c->trips);
    r.status=ALULA_COUPLED_BL_START;
    int interaction=!alula_airfoil_bl_complete(r.profiles);
    if (interaction) {
        int local=0;
        AlulaAirfoilBL *p=march_interaction(g,c,applied,&local);
        if (!alula_airfoil_bl_complete(p)) {alula_airfoil_bl_free(p);goto done;}
        alula_airfoil_bl_free(r.profiles);
        r.profiles=p;
        r.used_interaction=1;
        r.used_local_interaction=local;
    }
    if (!displacement(b,r.profiles,target,&r.bl_residual,&r.edge_residual)) goto done;
    r.displacement_residual=mismatch(target,applied,old_defect,n,chord);
    if (c->progress) c->progress(g,r.profiles,0,fmax(fmax(r.euler.norm,r.bl_residual),
        fmax(r.edge_residual,r.displacement_residual)),c->progress_context);
    double edge_tolerance=c->displacement_tolerance;
    double step=interaction?1:c->relaxation;
    int history=0;
    r.status=ALULA_COUPLED_ITERATION_LIMIT;
    for (int it=0;it<c->max_iterations && (r.displacement_residual>c->displacement_tolerance ||
        r.edge_residual>edge_tolerance);++it) {
        double euler_tolerance=c->euler_tolerance;
        /* Resolve the smaller flow changes that determine the final edge match. */
        if (c->displacement_scale>0 && r.edge_residual<1e-4)
            euler_tolerance=fmin(euler_tolerance,1e-12);
        alula_body_coordinates_save(g->coordinates,checkpoint);
        memcpy(packed,g->packed,rows*sizeof *packed);
        int accepted=0;
        /* Let the inverse seed settle before extrapolating the coupled iterates. */
        int accelerated=(c->displacement_scale==0 || r.edge_residual<.005) &&
            it>=ANDERSON_STARTUP_UPDATES && history &&
            anderson_direction(old_defect,dx,df,basis,direction,n,history);
        double trial_step=accelerated?1:step;
        for (int attempt=0;attempt<12;++attempt) {
            const double *update=accelerated?direction:old_defect;
            for (size_t i=0;i<n;++i) {
                double weight=c->displacement_scale>0 && !accelerated?
                    fmax(.2,fmin(8,c->displacement_scale/fmax(1e-15,applied[i]/chord))):1;
                trial[i]=applied[i]+trial_step*chord*update[i]*weight;
            }
            AlulaBodySolveResult e={.status=ALULA_BODY_BAD_STATE};
            if (alula_body_coordinates_set_displacement(g->coordinates,trial))
                e=alula_airfoil_grid_solve(g,c->euler_iterations,euler_tolerance);
            if (e.status==ALULA_BODY_INVALID_FINAL_GRID &&
                alula_airfoil_grid_repair(g,c->euler_iterations,c->euler_tolerance,0,&e))
                ++r.euler_grid_repairs;
            r.euler_updates+=e.iterations;
            int local=r.used_local_interaction;
            AlulaAirfoilBL *p=e.status==ALULA_BODY_SOLVED?march_interaction(g,c,interaction?trial:NULL,&local):NULL;
            double bl_residual,edge_residual;
            if (displacement(b,p,target,&bl_residual,&edge_residual)) {
                double norm=mismatch(target,trial,new_defect,n,chord);
                if (merit(norm,edge_residual)<merit(r.displacement_residual,r.edge_residual)*
                    (1-1e-4*fmin(1,trial_step)) || (norm<=c->displacement_tolerance && edge_residual<=edge_tolerance)) {
                    if (!accelerated) {
                        step=aitken(old_defect,new_defect,n,trial_step);
                        /* Aitken uses the unweighted defect. For a thickness-scaled
                         * direction, grow by at most one backtracking level. */
                        if (c->displacement_scale>0) step=fmin(step,2*trial_step);
                    }
                    if (local!=r.used_local_interaction) history=0;
                    else if (interaction) {
                        if (history==ANDERSON_DEPTH) {
                            --history;
                            memmove(dx,dx+n,history*n*sizeof *dx);
                            memmove(df,df+n,history*n*sizeof *df);
                        }
                        for (size_t i=0;i<n;++i) {
                            dx[history*n+i]=(trial[i]-applied[i])/chord;
                            df[history*n+i]=new_defect[i]-old_defect[i];
                        }
                        ++history;
                    }
                    memcpy(applied,trial,n*sizeof *applied);
                    memcpy(old_defect,new_defect,n*sizeof *old_defect);
                    alula_airfoil_bl_free(r.profiles);
                    r.profiles=p;
                    r.euler=e;
                    r.displacement_residual=norm;
                    r.bl_residual=bl_residual;
                    r.edge_residual=edge_residual;
                    r.used_local_interaction=local;
                    ++r.iterations;
                    if (c->progress) c->progress(g,p,r.iterations,fmax(fmax(e.norm,bl_residual),
                        fmax(edge_residual,norm)),c->progress_context);
                    accepted=1;
                    break;
                }
            }
            if (e.status!=ALULA_BODY_SOLVED) ++r.euler_rejections;
            else if (!alula_airfoil_bl_complete(p)) ++r.bl_rejections;
            else ++r.merit_rejections;
            alula_airfoil_bl_free(p);
            alula_body_coordinates_restore(g->coordinates,checkpoint);
            alula_body_coordinates_set_displacement(g->coordinates,applied);
            memcpy(g->packed,packed,rows*sizeof *packed);
            g->last_solve=r.euler;
            g->converged=1;
            ++r.backtracks;
            if (accelerated) {
                accelerated=0;
                trial_step=step;
            } else trial_step*=.5;
        }
        if (!accepted) {r.status=ALULA_COUPLED_LINE_SEARCH;break;}
    }
    if (r.displacement_residual<=c->displacement_tolerance && r.euler.norm<=c->euler_tolerance &&
        r.bl_residual<=c->euler_tolerance && r.edge_residual<=edge_tolerance)
        r.status=ALULA_COUPLED_OK;
    else if (r.status==ALULA_COUPLED_LINE_SEARCH && interaction && r.iterations<c->max_iterations)
        alula_airfoil_coupled_newton(g,c,&r,0);
done:
    free(checkpoint);
    free(work);
    return r;
}

AlulaCoupledResult alula_airfoil_coupled_solve_to_mach(AlulaAirfoilGrid *g,
    const AlulaCoupledControls *c,double target)
{
    if (!g || !(target>=g->mach && target<1))
        return (AlulaCoupledResult){.status=ALULA_COUPLED_INVALID,.displacement_residual=INFINITY,
            .bl_residual=INFINITY,.edge_residual=INFINITY,
            .euler={.status=ALULA_BODY_BAD_STATE,.norm=INFINITY,.max_mach=NAN}};
    AlulaCoupledResult result=alula_airfoil_coupled_solve(g,c);
    if (g->mach==target || result.status!=ALULA_COUPLED_OK) return result;
    size_t rows=g->body->size.rows,nd=alula_body_coordinates_displacement_count(g->coordinates);
    size_t nc=alula_body_coordinates_checkpoint_count(g->coordinates);
    double *packed=malloc((rows+nd)*sizeof(double));
    AlulaPoint *chart=malloc(nc*sizeof *chart);
    if (!packed || !chart) {
        result.status=ALULA_COUPLED_ALLOCATION;
        goto done;
    }
    double *thickness=packed+rows,step=.01;
    int joint=0,stages=0,retries=0;
    while (stages+retries<128 && g->mach<target) {
        double last=g->mach,next=fmin(target,last+step);
        alula_body_coordinates_save(g->coordinates,chart);
        alula_body_coordinates_get_displacement(g->coordinates,thickness);
        memcpy(packed,g->packed,rows*sizeof(double));
        if (!alula_body_set_mach(g->body,next)) break;
        g->mach=next;
        g->converged=0;
        AlulaCoupledResult trial={.status=ALULA_COUPLED_LINE_SEARCH};
        AlulaBodySolveResult flow=alula_airfoil_grid_solve(g,joint?0:c->euler_iterations,c->euler_tolerance);
        if (joint && (flow.status==ALULA_BODY_SOLVED || flow.status==ALULA_BODY_ITERATION_LIMIT)) {
            /* The previous profiles seed the new Mach. Newton replaces them
             * only after both the flow and the BL pass convergence. */
            trial.profiles=result.profiles;
            trial.euler=flow;
            alula_airfoil_coupled_newton(g,c,&trial,1);
        } else if (!joint && flow.status==ALULA_BODY_SOLVED) trial=alula_airfoil_coupled_solve(g,c);
        if (trial.status==ALULA_COUPLED_OK) {
            if (!joint) alula_airfoil_bl_free(result.profiles);
            result=trial;
            ++stages;
            step=fmin(.01,step*1.5);
            continue;
        }
        if (!joint) alula_airfoil_bl_free(trial.profiles);
        alula_body_coordinates_restore(g->coordinates,chart);
        alula_body_coordinates_set_displacement(g->coordinates,thickness);
        memcpy(g->packed,packed,rows*sizeof(double));
        alula_body_set_mach(g->body,last);
        g->mach=last;
        g->last_solve=result.euler;
        g->converged=1;
        ++retries;
        if (joint) step=.5*(next-last);
        joint=1;
        if (step<1e-4) break;
    }
    result.mach_stages=stages;
    result.mach_retries=retries;
    if (g->mach!=target)
        result.status=stages+retries==128?ALULA_COUPLED_ITERATION_LIMIT:ALULA_COUPLED_LINE_SEARCH;
done:
    free(packed);
    free(chart);
    return result;
}
