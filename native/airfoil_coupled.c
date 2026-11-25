// SPDX-License-Identifier: GPL-2.0-or-later
#include "airfoil_coupled.h"
#include "airfoil_grid_internal.h"
#include "body_internal.h"

#include <stdlib.h>
#include <string.h>

static int displacement(const AlulaBody *b,const AlulaAirfoilBL *p,double *out,double *residual)
{
    if (!alula_airfoil_bl_complete(p)) return 0;
    double chord=b->in.length_scale;
    *residual=0;
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
        }
        for (size_t j=1;j<wake->count;++j) *out++=chord*wake->stations[j].delta_star;
        *residual=fmax(*residual,wake->result.maximum_residual);
    }
    return 1;
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

AlulaCoupledResult alula_airfoil_coupled_solve(AlulaAirfoilGrid *g,const AlulaCoupledControls *c)
{
    AlulaCoupledResult r={.status=ALULA_COUPLED_INVALID,.displacement_residual=INFINITY,
        .bl_residual=INFINITY,.euler={.status=ALULA_BODY_BAD_STATE,.norm=INFINITY,.max_mach=NAN}};
    if (!g || !g->converged || !c || !positive(c->reynolds) || !(c->ncrit>=0) || !isfinite(c->ncrit) ||
        !positive(c->euler_tolerance) || !positive(c->displacement_tolerance) ||
        !(c->relaxation>0 && c->relaxation<=1) || c->max_iterations<0 || c->euler_iterations<1) return r;
    const AlulaBody *b=g->body;
    size_t n=alula_body_coordinates_displacement_count(g->coordinates),rows=(size_t)b->size.rows;
    if (!n) return r;
    r.euler=g->last_solve;
    double chord=b->in.length_scale;
    double *work=malloc((5*n+rows)*sizeof *work);
    AlulaPoint *checkpoint=malloc(alula_body_coordinates_checkpoint_count(g->coordinates)*sizeof *checkpoint);
    if (!work || !checkpoint) {
        r.status=ALULA_COUPLED_ALLOCATION;goto done;
    }
    double *applied=work,*target=applied+n,*trial=target+n,*old_defect=trial+n,*new_defect=old_defect+n;
    double *packed=new_defect+n;
    if (!alula_body_coordinates_get_displacement(g->coordinates,applied)) goto done;
    r.profiles=alula_airfoil_bl_create(g,c->reynolds,c->ncrit,c->trips);
    r.status=ALULA_COUPLED_BL_START;
    if (!displacement(b,r.profiles,target,&r.bl_residual)) goto done;
    r.displacement_residual=mismatch(target,applied,old_defect,n,chord);
    double step=c->relaxation;
    r.status=ALULA_COUPLED_ITERATION_LIMIT;
    for (int it=0;it<c->max_iterations && r.displacement_residual>c->displacement_tolerance;++it) {
        alula_body_coordinates_save(g->coordinates,checkpoint);
        memcpy(packed,g->packed,rows*sizeof *packed);
        int accepted=0;
        for (int attempt=0;attempt<12;++attempt,step*=.5) {
            for (size_t i=0;i<n;++i) trial[i]=applied[i]+step*chord*old_defect[i];
            AlulaBodySolveResult e={.status=ALULA_BODY_BAD_STATE};
            if (alula_body_coordinates_set_displacement(g->coordinates,trial))
                e=alula_airfoil_grid_solve(g,c->euler_iterations,c->euler_tolerance);
            r.euler_updates+=e.iterations;
            AlulaAirfoilBL *p=e.status==ALULA_BODY_SOLVED?
                alula_airfoil_bl_create(g,c->reynolds,c->ncrit,c->trips):NULL;
            double bl_residual;
            if (displacement(b,p,target,&bl_residual)) {
                double norm=mismatch(target,trial,new_defect,n,chord);
                if (norm<r.displacement_residual*(1-1e-4*fmin(1,step)) ||
                    norm<=c->displacement_tolerance) {
                    step=aitken(old_defect,new_defect,n,step);
                    memcpy(applied,trial,n*sizeof *applied);
                    memcpy(old_defect,new_defect,n*sizeof *old_defect);
                    alula_airfoil_bl_free(r.profiles);
                    r.profiles=p;
                    r.euler=e;
                    r.displacement_residual=norm;
                    r.bl_residual=bl_residual;
                    ++r.iterations;
                    accepted=1;
                    break;
                }
            }
            alula_airfoil_bl_free(p);
            alula_body_coordinates_restore(g->coordinates,checkpoint);
            alula_body_coordinates_set_displacement(g->coordinates,applied);
            memcpy(g->packed,packed,rows*sizeof *packed);
            g->last_solve=r.euler;
            g->converged=1;
            ++r.backtracks;
        }
        if (!accepted) {r.status=ALULA_COUPLED_LINE_SEARCH;break;}
    }
    if (r.displacement_residual<=c->displacement_tolerance && r.euler.norm<=c->euler_tolerance &&
        r.bl_residual<=c->euler_tolerance)
        r.status=ALULA_COUPLED_OK;
done:
    free(checkpoint);
    free(work);
    return r;
}
