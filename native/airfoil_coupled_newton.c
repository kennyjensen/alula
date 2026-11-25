// SPDX-License-Identifier: GPL-2.0-or-later
#include "airfoil_bl_internal.h"
#include "airfoil_grid_internal.h"
#include "dense.h"

#include <float.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Eliminate the Euler variables with the current sparse Jacobian. The
 * remaining unknowns per station are log(delta*), log(H-1), and either
 * amplification or log(sqrt(Ctau)). The first wake station is an exact
 * trailing-edge join and has no independent unknowns. */
typedef struct {size_t next; int color;} LocalColumn;

typedef struct {
    LocalColumn *local;
    size_t variables,rows,displacements;
    AlulaBLRegime *regimes;
    double *shear_scale,*state,*trial,*residual,*trial_residual,*jacobian,*direction;
    double *packed,*applied,*thickness,*trial_thickness;
    double *flow_source,*flow_direction,*response_rhs,*predictions;
    AlulaPoint *chart;
    AlulaAirfoilBL *base,*probe;
} Newton;

static int seed(Newton *w,const AlulaAirfoilGrid *g,const AlulaAirfoilBL *p)
{
    size_t k=0,offset=0;
    for (int e=0;e<p->bodies;++e) {
        size_t ns=p->branches[e][0].count,ns1=ns+1,nw=p->branches[e][ALULA_BL_MERGED_WAKE].count;
        size_t first_wake=k+2*ns;
        for (int side=0;side<3;++side) {
            const AlulaAirfoilBLBranch *b=&p->branches[e][side];
            for (size_t j=(side==ALULA_BL_MERGED_WAKE);j<b->count;++j,++k) {
                size_t index=side==ALULA_BL_MERGED_WAKE?2*ns1+j-1:(1-side)*ns1+j+1;
                double d=w->applied[offset+index]/g->body->in.length_scale;
                AlulaBLStation st=b->stations[j];
                AlulaBLRegime regime=b->regimes[j];
                double h=(st.delta_star-st.wake_gap)/st.theta;
                if (!(d>st.wake_gap && h>1 && st.theta>0) ||
                    (regime>=ALULA_BL_TURBULENT && !(st.aux>0))) return 0;
                w->state[3*k]=log(d);
                w->state[3*k+1]=log(h-1);
                w->state[3*k+2]=regime>=ALULA_BL_TURBULENT?log(st.aux):st.aux;
                w->regimes[k]=regime;
                /* A local change affects this block and its downstream
                 * neighbor. Both surface ends feed the first wake interval. */
                size_t next=j+1<b->count?k+1:(side!=ALULA_BL_MERGED_WAKE && nw>1?first_wake:k);
                w->local[k]=(LocalColumn){next,2*side+(int)(j%2)};

                w->shear_scale[k]=regime>=ALULA_BL_TURBULENT?1/(5*st.theta):1;
            }
        }
        offset+=2*ns1+nw-1;
    }
    return k*3==w->variables && offset==w->displacements;
}

static void displacement(const AlulaAirfoilGrid *g,const double *state,double *out)
{
    size_t k=0;
    double chord=g->body->in.length_scale;
    for (int e=0;e<g->elements;++e) {
        AlulaBodyCut cut=g->cuts[e];
        int ns=cut.trailing-cut.leading;
        out[0]=out[ns+1]=.5*chord*(exp(state[3*k])+exp(state[3*(k+ns)]));
        for (int side=0;side<2;++side) for (int j=0;j<ns;++j)
            out[(1-side)*(ns+1)+j+1]=chord*exp(state[3*k++]);
        out+=2*(ns+1);
        for (int j=cut.trailing;j<g->body->in.nx;++j) *out++=chord*exp(state[3*k++]);
    }
}

typedef struct {double rms,raw_maximum;} Norm;

static int evaluate(const Newton *w,AlulaAirfoilBL *p,const double *state,double *f,Norm *norm,
    int *transition_consistent)
{
    size_t k=0;
    double sum=0;
    *norm=(Norm){0};
    if (transition_consistent) *transition_consistent=1;
    p->complete=0;
    for (int e=0;e<p->bodies;++e) for (int side=0;side<3;++side) {
        AlulaAirfoilBLBranch *b=&p->branches[e][side];
        size_t start=(size_t)(b->stations-p->stations);
        const AlulaBLEdge *edges=p->edges+start;
        AlulaBLStation *stations=p->stations+start,virtual_station;
        const AlulaBLStation *previous=NULL;
        b->result=(AlulaBLMarchResult){0};
        int wake=side==ALULA_BL_MERGED_WAKE;
        double trip=wake?DBL_MAX:p->trip[e][side];
        if (wake) {
            const AlulaAirfoilBLBranch *u=&p->branches[e][0],*l=u+1;
            b->result=alula_bl_march_wake(&p->conditions,u->stations+u->count-1,
                l->stations+l->count-1,edges,1,stations);
            if (b->result.status!=ALULA_BL_OK) return 0;
            previous=stations;
            p->regimes[start]=ALULA_BL_WAKE;
        } else if (trip<=edges[0].s) {
            AlulaBLEdge first={.5*trip,.5*trip*edges[0].ue/edges[0].s,0};
            AlulaBLRegime regime;
            if (alula_bl_march_surface(&p->conditions,&first,1,DBL_MAX,1,&virtual_station,&regime).status
                !=ALULA_BL_OK) return 0;
            previous=&virtual_station;
        }
        for (size_t j=wake;j<b->count;++j,++k) {
            AlulaBLRegime regime=w->regimes[k];
            const double *x=state+3*k;
            double d=exp(x[0]),h=1+exp(x[1]),aux=regime>=ALULA_BL_TURBULENT?exp(x[2]):x[2];
            AlulaBLEdge edge=edges[j];
            stations[j]=(AlulaBLStation){edge.s,(d-edge.wake_gap)/h,d,edge.ue,aux,edge.wake_gap};
            p->regimes[start+j]=regime;
            AlulaBLBlock block;
            if ((previous==NULL)!=(regime==ALULA_BL_SIMILARITY) ||
                alula_bl_interval(&p->conditions,previous,stations+j,regime,trip,1,&block)!=ALULA_BL_OK) return 0;
            if (regime==ALULA_BL_TRANSITION) b->result.transition=block.transition;
            /* Perturbations use the accepted regime layout. Only the state
             * being checked for convergence needs a fresh transition check. */
            if (regime==ALULA_BL_LAMINAR && transition_consistent && *transition_consistent) {
                AlulaBLTransition event;
                if (alula_bl_check_transition(&p->conditions,previous,stations+j,trip,&event)!=ALULA_BL_OK ||
                    event.active) *transition_consistent=0;
            }
            for (int i=0;i<3;++i) {
                double v=block.residual[i];
                b->result.maximum_residual=fmax(b->result.maximum_residual,fabs(v));
                f[3*k+i]=i==0?v*w->shear_scale[k]:v;
                sum+=f[3*k+i]*f[3*k+i];
            }
            ++b->result.completed;
            previous=stations+j;
        }
        norm->raw_maximum=fmax(norm->raw_maximum,b->result.maximum_residual);
    }
    norm->rms=sqrt(sum/w->variables);
    p->complete=isfinite(norm->rms);
    return p->complete;
}

/* Shape and shear/amplification leave the Euler solution fixed. Alternating
 * stations on one branch have disjoint residual rows, so perturb an entire
 * color at once, including independently equivalent branches on other bodies. */
static int local_columns(Newton *w,int variable,int color)
{
    size_t m=w->variables,n=m/3;
    Norm norm;
    for (int sign=1;sign>=-1;sign-=2) {
        double delta=sign*2e-6;
        for (size_t k=0;k<n;++k) if (w->local[k].color==color)
            w->trial[3*k+variable]=w->state[3*k+variable]+delta;
        int valid=evaluate(w,w->base,w->trial,w->trial_residual,&norm,NULL);
        for (size_t k=0;k<n;++k) if (w->local[k].color==color) {
            size_t column=3*k+variable;
            w->trial[column]=w->state[column];
            if (!valid) continue;
            size_t rows[]={3*k,3*w->local[k].next};
            for (int block=0;block<2;++block) for (int j=0;j<3;++j) {
                size_t row=rows[block]+j;
                w->jacobian[row*m+column]=(w->trial_residual[row]-w->residual[row])/delta;
            }
        }
        if (valid) return 1;
    }
    return 0;
}

/* Lift the reduced BL step back into Euler coordinates. At a new Mach the
 * initial flow has a nonzero residual, so its particular solution also
 * contributes to the right-hand side of the reduced equations. */
static int solve_direction(Newton *w,AlulaAirfoilGrid *g,const AlulaCoupledControls *c)
{
    size_t m=w->variables;
    for (size_t i=0;i<m;++i) w->direction[i]=-w->residual[i];
    if (!w->flow_source) return alula_dense_solve(w->jacobian,w->direction,m,1);

    if (!alula_body_coordinates_set_displacement(g->coordinates,w->applied)) return 0;
    double largest=0;
    for (size_t i=0;i<w->rows;++i) largest=fmax(largest,fabs(w->flow_source[i]));
    double step=fmin(1,2e-6/fmax(DBL_MIN,largest));
    for (size_t i=0;i<w->rows;++i) g->packed[i]=w->packed[i]+step*w->flow_source[i];
    Norm norm;
    if (!alula_airfoil_bl_resample(w->probe,g,c->trips,1) ||
        !evaluate(w,w->probe,w->state,w->trial_residual,&norm,NULL)) return 0;
    for (size_t i=0;i<m;++i) w->direction[i]-=(w->trial_residual[i]-w->residual[i])/step;
    if (!alula_dense_solve(w->jacobian,w->direction,m,1)) return 0;

    largest=1;
    for (size_t i=0;i<m;i+=3) largest=fmax(largest,fabs(w->direction[i]));
    step=2e-6/largest;
    for (size_t i=0;i<m;++i) w->trial[i]=w->state[i]+step*w->direction[i];
    displacement(g,w->trial,w->trial_thickness);
    for (size_t i=0;i<w->displacements;++i)
        w->trial_thickness[i]=w->applied[i]+(w->trial_thickness[i]-w->thickness[i]);
    if (!alula_body_coordinates_set_displacement(g->coordinates,w->trial_thickness) ||
        !alula_body_newton_response(g->solver,w->packed,g->packed)) return 0;
    for (size_t i=0;i<w->rows;++i)
        w->flow_direction[i]=w->flow_source[i]+(g->packed[i]-w->packed[i])/step;
    return 1;
}

static int set_thickness_probe(Newton *w,AlulaAirfoilGrid *g,size_t column,double delta)
{
    w->trial[column]=w->state[column]+delta;
    displacement(g,w->trial,w->trial_thickness);
    w->trial[column]=w->state[column];
    /* The partitioned seed can have a different leading-edge average. */
    for (size_t j=0;j<w->displacements;++j)
        w->trial_thickness[j]+=w->applied[j]-w->thickness[j];
    return alula_body_coordinates_set_displacement(g->coordinates,w->trial_thickness);
}

static int thickness_columns(Newton *w,AlulaAirfoilGrid *g,const AlulaCoupledControls *c,
    size_t first,int count,double delta)
{
    for (int j=0;j<count;++j) {
        if (!set_thickness_probe(w,g,first+3*j,delta) ||
            !alula_body_newton_displacement_rhs(g->solver,w->packed,w->response_rhs+j*w->rows)) return 0;
    }
    if (!alula_body_newton_predict_batch(g->solver,w->packed,count,w->response_rhs,w->predictions)) return 0;
    for (int j=0;j<count;++j) {
        size_t column=first+3*j;
        if (count>1 && !set_thickness_probe(w,g,column,delta)) return 0;
        memcpy(g->packed,w->predictions+j*w->rows,w->rows*sizeof *g->packed);
        if (!alula_airfoil_bl_resample(w->probe,g,c->trips,1)) return 0;
        w->trial[column]=w->state[column]+delta;
        Norm norm;
        int valid=evaluate(w,w->probe,w->trial,w->trial_residual,&norm,NULL);
        w->trial[column]=w->state[column];
        if (!valid) return 0;
        for (size_t i=0;i<w->variables;++i)
            w->jacobian[i*w->variables+column]=(w->trial_residual[i]-w->residual[i])/delta;
    }
    return 1;
}

static int build_jacobian(Newton *w,AlulaAirfoilGrid *g,const AlulaCoupledControls *c)
{
    if (w->flow_source) {
        if (!alula_body_newton_correction(g->solver,w->packed,w->flow_source)) return 0;
    } else if (!alula_body_newton_linearize(g->solver,w->packed)) return 0;
    size_t m=w->variables;
    /* Dense LU overwrites even the structural zeros. */
    memset(w->jacobian,0,m*m*sizeof(double));
    memcpy(w->trial,w->state,m*sizeof(double));
    for (int variable=1;variable<=2;++variable) for (int color=0;color<6;++color)
        if (!local_columns(w,variable,color)) return 0;
    for (size_t column=0;column<m;column+=3*ALULA_SPARSE_BATCH) {
        int count=(int)((m-column)/3);
        if (count>ALULA_SPARSE_BATCH) count=ALULA_SPARSE_BATCH;
        if (thickness_columns(w,g,c,column,count,2e-6)) continue;
        /* A column near the edge of the admissible flow domain can require
         * a one-sided difference in the opposite direction. */
        for (int j=0;j<count;++j)
            if (!thickness_columns(w,g,c,column+3*j,1,2e-6) &&
                !thickness_columns(w,g,c,column+3*j,1,-2e-6)) return 0;
    }
    return solve_direction(w,g,c);
}

static int line_search(Newton *w,AlulaAirfoilGrid *g,const AlulaCoupledControls *c,
    AlulaCoupledResult *r,double rms)
{
    double step=1;
    for (size_t i=0;i<w->variables;++i)
        if (i%3!=2 || w->regimes[i/3]>=ALULA_BL_TURBULENT)
            step=fmin(step,1.5/fmax(DBL_MIN,fabs(w->direction[i])));
    for (int trial=0;trial<16;++trial,step*=.5) {
        alula_body_coordinates_restore(g->coordinates,w->chart);
        if (w->flow_direction)
            for (size_t i=0;i<w->rows;++i) g->packed[i]=w->packed[i]+step*w->flow_direction[i];
        else memcpy(g->packed,w->packed,w->rows*sizeof(double));
        for (size_t i=0;i<w->variables;++i) w->trial[i]=w->state[i]+step*w->direction[i];
        displacement(g,w->trial,w->trial_thickness);
        Norm norm;
        AlulaBodySolveResult er={.status=ALULA_BODY_BAD_STATE};
        if (alula_body_coordinates_set_displacement(g->coordinates,w->trial_thickness))
            er=alula_airfoil_grid_solve(g,w->flow_direction?0:c->euler_iterations,
                w->flow_direction?c->euler_tolerance:fmin(1e-10,c->euler_tolerance));
        r->euler_updates+=er.iterations;
        if (er.status!=ALULA_BODY_SOLVED && !(w->flow_direction && er.status==ALULA_BODY_ITERATION_LIMIT))
            ++r->euler_rejections;
        else if (!alula_airfoil_bl_resample(w->probe,g,c->trips,w->flow_direction!=NULL) ||
            !evaluate(w,w->probe,w->trial,w->trial_residual,&norm,NULL)) ++r->bl_rejections;
        else if ((w->flow_direction?hypot(norm.rms,er.norm):norm.rms)<rms &&
            (!w->flow_direction || !g->chart.rebase || g->chart.rebase(g->chart.context,g->packed))) {
            memcpy(w->state,w->trial,w->variables*sizeof(double));
            memcpy(w->applied,w->trial_thickness,w->displacements*sizeof(double));
            AlulaAirfoilBL *old=w->base;
            w->base=w->probe;
            w->probe=old;
            return 1;
        } else ++r->merit_rejections;
        ++r->backtracks;
    }
    return 0;
}

void alula_airfoil_coupled_newton(AlulaAirfoilGrid *g,const AlulaCoupledControls *c,AlulaCoupledResult *r,
    int coupled_flow)
{
    size_t nd=alula_body_coordinates_displacement_count(g->coordinates);
    /* Each element has two extra thickness values at the shared leading edge. */
    size_t n=nd-2*(size_t)g->elements,m=3*n,rows=(size_t)g->body->size.rows;
    size_t chart_count=alula_body_coordinates_checkpoint_count(g->coordinates);
    if (!n || m>SIZE_MAX/sizeof(double)/m) return;
    Newton w={.variables=m,.rows=rows,.displacements=nd};
    double *work=malloc((5*m+n+((coupled_flow?4:2)+2*ALULA_SPARSE_BATCH)*rows+4*nd)*sizeof(double));
    AlulaPoint *charts=NULL;
    w.regimes=malloc(n*sizeof *w.regimes);
    w.local=malloc(n*sizeof *w.local);
    if (!work || !w.regimes || !w.local) goto done;
    w.state=work;w.trial=w.state+m;w.direction=w.trial+m;w.residual=w.direction+m;
    w.trial_residual=w.residual+m;w.shear_scale=w.trial_residual+m;w.packed=w.shear_scale+n;
    double *initial_packed=w.packed+rows,*initial_d=initial_packed+rows;
    w.applied=initial_d+nd;w.thickness=w.applied+nd;w.trial_thickness=w.thickness+nd;
    w.response_rhs=w.trial_thickness+nd;
    w.predictions=w.response_rhs+ALULA_SPARSE_BATCH*rows;
    if (coupled_flow) {
        w.flow_source=w.predictions+ALULA_SPARSE_BATCH*rows;
        w.flow_direction=w.flow_source+rows;
    }
    if (!alula_body_coordinates_get_displacement(g->coordinates,w.applied) || !seed(&w,g,r->profiles)) goto done;
    w.jacobian=malloc(m*m*sizeof(double));
    charts=malloc(2*chart_count*sizeof *charts);
    if (!w.jacobian || !charts) goto done;
    w.base=alula_airfoil_bl_sample(g,c->reynolds,c->ncrit,c->trips,coupled_flow);
    w.probe=alula_airfoil_bl_sample(g,c->reynolds,c->ncrit,c->trips,coupled_flow);
    if (!w.base || !w.probe) goto done;
    w.chart=charts+chart_count;
    alula_body_coordinates_save(g->coordinates,charts);
    memcpy(initial_packed,g->packed,rows*sizeof(double));
    memcpy(initial_d,w.applied,nd*sizeof(double));
    AlulaBodySolveResult initial_euler=g->last_solve;
    for (int updates=0;;++updates) {
        /* The joint step rebases the moving grid after evaluating a trial.
         * Certify its Euler residual, convexity and wall samples on the
         * retained chart before publishing a converged solution. */
        if (coupled_flow && g->converged) {
            AlulaBodySolveResult verified=alula_airfoil_grid_solve(g,0,c->euler_tolerance);
            if (verified.status!=ALULA_BODY_SOLVED && verified.status!=ALULA_BODY_ITERATION_LIMIT) break;
            if (!alula_airfoil_bl_resample(w.base,g,c->trips,!g->converged)) break;
        }
        Norm norm;
        int transition_consistent;
        if (!evaluate(&w,w.base,w.state,w.residual,&norm,&transition_consistent)) break;
        displacement(g,w.state,w.thickness);
        double mismatch=0;
        for (size_t j=0;j<nd;++j) mismatch=fmax(mismatch,fabs(w.thickness[j]-w.applied[j])/g->body->in.length_scale);
        if (updates && c->progress) c->progress(g,w.base,r->iterations+updates,
            fmax(g->last_solve.norm,fmax(norm.raw_maximum,mismatch)),c->progress_context);
        if (transition_consistent && norm.raw_maximum<=c->euler_tolerance &&
            mismatch<=c->displacement_tolerance && g->converged && g->last_solve.norm<=c->euler_tolerance) {
            alula_airfoil_bl_free(r->profiles);
            r->profiles=w.base;
            w.base=NULL;
            r->bl_residual=norm.raw_maximum;
            r->edge_residual=0; /* Edge speeds are the Euler samples used in every interval. */
            r->displacement_residual=mismatch;
            r->euler=g->last_solve;
            r->iterations+=updates;
            r->status=ALULA_COUPLED_OK;
            goto done;
        }
        if (r->iterations+updates>=c->max_iterations) break;
        memcpy(w.packed,g->packed,rows*sizeof(double));
        alula_body_coordinates_save(g->coordinates,w.chart);
        double merit=coupled_flow?hypot(norm.rms,g->last_solve.norm):norm.rms;
        if (!build_jacobian(&w,g,c) || !line_search(&w,g,c,r,merit)) break;
        ++r->coupled_newton_updates;
    }
    alula_body_coordinates_restore(g->coordinates,charts);
    alula_body_coordinates_set_displacement(g->coordinates,initial_d);
    memcpy(g->packed,initial_packed,rows*sizeof(double));
    g->last_solve=initial_euler;
    g->converged=initial_euler.status==ALULA_BODY_SOLVED;
done:
    alula_airfoil_bl_free(w.base);
    alula_airfoil_bl_free(w.probe);
    free(w.local);free(w.regimes);free(w.jacobian);free(charts);free(work);
}
