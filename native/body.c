// SPDX-License-Identifier: GPL-2.0-or-later
#include "body_internal.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>

static int active(const AlulaBody *b,int g,int i)
{
    return i>=b->cuts[g].leading && i<=b->cuts[g].trailing;
}

enum { BLEND_ENTROPY, BLEND_MOMENTUM, BLEND_MIXED };

static int end_position(const AlulaBody *b,int g,int i,int j)
{
    /* A shared cut has one canonical lower passage position. An independent
     * wake has its upper bank equation in the wake-gap family. */
    return !(g>0 && j==0 && (i==0 || i==b->in.nx));
}

static void build_rows(AlulaBody *b)
{
    int nx=b->in.nx, bodies=b->in.bodies;
    Row *r=b->rows;
    for (int i=1;i<nx;++i) for (int g=0;g<=bodies;++g) for (int j=0;j<b->tubes[g];++j)
        *r++=(Row){ALULA_BODY_STREAMWISE,i,g,j};
    for (int g=0;g<=bodies;++g) for (int j=0;j<b->tubes[g];++j)
        *r++=(Row){ALULA_BODY_INLET_DENSITY,0,g,j};
    for (int i=1;i<nx;++i) {
        for (int g=0;g<=bodies;++g) for (int j=1;j<b->tubes[g];++j)
            *r++=(Row){ALULA_BODY_INTERNAL_PRESSURE,i,g,j};
        *r++=(Row){ALULA_BODY_FARFIELD_PRESSURE,i,0,0};
        *r++=(Row){ALULA_BODY_FARFIELD_PRESSURE,i,bodies,b->tubes[bodies]};
        for (int g=0;g<bodies;++g) if (!active(b,g,i))
            *r++=(Row){ALULA_BODY_CUT_PRESSURE,i,g,0};
    }
    for (int end=0;end<2;++end) {
        int i=end ? nx : 0;
        for (int g=0;g<=bodies;++g) for (int j=0;j<=b->tubes[g];++j)
            if (end_position(b,g,i,j)) *r++=(Row){ALULA_BODY_END_TANGENCY,i,g,j};
    }
    if (b->in.independent_wakes) for (int g=0;g<bodies;++g)
        for (int i=b->cuts[g].trailing+1;i<=nx;++i) *r++=(Row){ALULA_BODY_WAKE_GAP,i,g,0};
    for (int g=0;g<bodies;++g)
        *r++=(Row){ALULA_BODY_TRAILING_KUTTA,b->cuts[g].trailing,g,0};
    for (int g=0;g<bodies;++g) if (b->cuts[g].rounded)
        *r++=(Row){ALULA_BODY_LEADING_KUTTA,b->cuts[g].leading,g,0};
    for (int k=0;k<3;++k) *r++=(Row){ALULA_BODY_FARFIELD_MATCH,0,0,k};
    b->size.rows=(int)(r-b->rows);
}

int alula_body_set_mach(AlulaBody *b,double mach)
{
    if (!b || !(mach>0 && mach<1)) return 0;
    double m2=mach*mach,beta=sqrt(1-m2),beta3=beta*beta*beta;
    double gamma=b->in.gamma,gm=gamma-1;
    double p_inf=1/(gamma*m2),h0=1/(gm*m2)+.5,rho_total=pow(1+.5*gm*m2,1/gm);
    if (!positive(p_inf) || !positive(h0) || !positive(rho_total)) return 0;
    b->in.mach=mach;
    b->mach_sq=m2;
    b->beta=beta;
    b->pg_a=.25*((3-gamma)/beta+(gamma+1)/beta3);
    b->pg_b=(gamma+1)/16*(1/beta-1/beta3);
    b->pg_scale=m2/(TWO_PI*TWO_PI);
    b->p_inf=p_inf;
    b->h0=h0;
    b->rho_total=rho_total;
    ++b->residual_generation;
    return 1;
}

AlulaBody *alula_body_create(const AlulaBodyInput *in)
{
    if (!in || in->nx<4 || in->nx>50000 || in->bodies<1 || in->bodies>MAX_BODIES ||
        !in->tubes || !in->cuts || !positive(in->mach) || in->mach>=1 || !isfinite(in->alpha) ||
        !isfinite(in->gamma) || in->gamma<=1 || !positive(in->length_scale) ||
        !isfinite(in->correction) || in->correction<0 || !finite_point(in->center) ||
        in->mode<ALULA_BODY_MOMENTUM || in->mode>ALULA_BODY_HYBRID ||
        ((in->independent_wakes || in->bank_outlet) && !in->displaced) ||
        (in->mode==ALULA_BODY_ISENTROPIC && in->upwind) ||
        (in->mode==ALULA_BODY_HYBRID && (!in->upwind || in->ismom<1 || in->ismom>4 || !positive(in->epsilon))) ||
        (in->upwind && (!isfinite(in->mucon) || !isfinite(in->mcrit) || in->mcrit<0 || in->mcrit>1))) return NULL;
    int nt=0;
    for (int g=0;g<=in->bodies;++g) {
        if (in->tubes[g]<1 || in->tubes[g]>100000-nt) return NULL;
        nt+=in->tubes[g];
    }
    if (nt>100000/in->nx) return NULL;
    AlulaBody *b=calloc(1,sizeof *b);
    if (!b) return NULL;
    b->failed_initial_cell=-1;
    b->in=*in;
    b->in.tubes=b->tubes; b->in.cuts=b->cuts;
    int nx=in->nx;
    for (int g=0;g<=in->bodies;++g) {
        b->tubes[g]=in->tubes[g];
        b->column[g+1]=b->column[g]+b->tubes[g];
        b->node[g+1]=b->node[g]+(nx+1)*(b->tubes[g]+1);
        if (g==in->bodies) continue;
        AlulaBodyCut c=in->cuts[g];
        if (c.leading<1 || c.trailing<=c.leading || c.trailing>=nx) goto fail;
        if (c.through && (in->mode!=ALULA_BODY_HYBRID || in->ismom!=3 || c.through<c.leading ||
            c.through>=nx || c.lower_tube<0 || c.lower_tube>=b->column[g+1] ||
            c.upper_tube<=b->column[g+1] || c.upper_tube>nt)) goto fail;
        b->cuts[g]=c;
        b->wake[g+1]=b->wake[g]+(in->independent_wakes ? nx-c.trailing : 0);
    }
    b->size=(AlulaBodySize){.nodes=b->node[in->bodies+1],.density=nx*nt,.tubes=nt,.wake_widths=b->wake[in->bodies]};
    double angle=in->alpha*TWO_PI/360;
    b->freestream=(AlulaPoint){cos(angle),sin(angle)};
    if (!finite_point(b->freestream) || !alula_body_set_mach(b,in->mach)) goto fail;
    /* Each geometric position has at most one pressure/tangency/gap row.
     * Extra rows are the two edge conditions per body and three farfield modes. */
    size_t maxrows=(size_t)b->size.density+b->size.nodes+2*in->bodies+3;
    b->rows=malloc(maxrows*sizeof *b->rows);
    b->cells=malloc((size_t)(nx-1)*nt*sizeof *b->cells);
    b->sections=malloc((size_t)nx*nt*sizeof *b->sections);
    if (!b->rows || !b->cells || !b->sections) goto fail;
    if (in->mode==ALULA_BODY_HYBRID && in->ismom==4) {
        b->blend_regime=malloc((size_t)(nx-1)*nt*sizeof *b->blend_regime);
        if (!b->blend_regime) goto fail;
    }
    if (in->upwind && nx>2) {
        b->upwind_stencils=malloc((size_t)(nx-2)*nt*sizeof *b->upwind_stencils);
        if (!b->upwind_stencils) goto fail;
    }
    build_rows(b);
    return b;
fail:
    alula_body_free(b);
    return NULL;
}

AlulaBodySize alula_body_size(const AlulaBody *b) { return b ? b->size : (AlulaBodySize){0}; }

void alula_body_free(AlulaBody *b)
{
    if (!b) return;
    free(b->rows); free(b->cells); free(b->sections);
    free(b->upwind_stencils);
    free(b->blend_regime);
    free(b->jacobian_fields);
    free(b->section_tangents); free(b->cell_tangents); free(b);
}

int alula_body_farfield(const AlulaBody *b,const double *coeff,AlulaPoint p,Farfield *out)
{
    AlulaPoint d=sub(p,b->in.center), f=b->freestream;
    double x=dot(d,f)/b->beta,y=cross(f,d),r2=x*x+y*y;
    if (!positive(r2)) return 0;
    double r4=r2*r2,r6=r4*r2,logr=.5*log(r2);
    double cubic=x*x*x-3*x*y*y, a=b->pg_a,v=b->pg_b,w=b->pg_scale;
    AlulaPoint grad[]={
        {-y/(TWO_PI*r2),x/(TWO_PI*r2)},
        {x/(TWO_PI*r2),y/(TWO_PI*r2)},
        {(y*y-x*x)/(TWO_PI*r4),-2*x*y/(TWO_PI*r4)},
        {-2*x*y/(TWO_PI*r4),(x*x-y*y)/(TWO_PI*r4)},
        {w*(a*(x*x+logr*(y*y-x*x))/r4+v*((3*x*x-3*y*y)/r4-4*x*cubic/r6)),
         w*(a*x*y*(1-2*logr)/r4+v*(-6*x*y/r4-4*y*cubic/r6))}
    };
    AlulaPoint velocity=f;
    for (int k=0;k<5;++k) {
        double gx=grad[k].x, gy=grad[k].y;
        out->basis[k]=(AlulaPoint){gx*f.x/b->beta-gy*f.y,gx*f.y/b->beta+gy*f.x};
        velocity=add(velocity,scale(out->basis[k],coeff[k]));
    }
    double q2=dot(velocity,velocity), gm=b->in.gamma-1;
    double temp=1+.5*gm*b->mach_sq*(1-q2);
    if (!positive(temp) || !(b->mach_sq*q2/temp<1)) return 0;
    out->velocity=velocity;
    out->pressure=b->p_inf*pow(temp,b->in.gamma/gm);
    out->log_radius=logr;
    return positive(out->pressure);
}

static inline AlulaCellStatus initialize_cell(const AlulaBody *b,const AlulaBodyState *s,
    int i,int g,int j,AlulaCell *cell,const AlulaCell *previous)
{
    int k=index_section(b,g,i,j);
    AlulaCellInput p={.density={s->density[k],s->density[k+b->size.tubes]},
        .mass=s->mass[b->column[g]+j],.h0=b->h0,.gamma=b->in.gamma,
        .correction=b->in.correction,.positive_simple=b->in.positive_simple};
    for (int q=0;q<3;++q) {
        p.lower[q]=point(b,s,g,i+q,j);
        p.upper[q]=point(b,s,g,i+q,j+1);
    }
    return alula_cell_initialize(cell,&p,previous);
}

/* Probe only the wall-adjacent gas and geometry used for BL edge matching.
 * Keep the accepted Euler stencil intact while sampling a linear response. */
double alula_body_edge_speed(const AlulaBody *b,const AlulaBodyState *state,int g,int i,int j,
    AlulaEdgeCursor *cursor)
{
    int station=i<b->in.nx?i:b->in.nx-1;
    const AlulaCell *c=cell(b,g,station,j);
    if (state) {
        AlulaCell *current=&cursor->cells[station%2];
        if (cursor->station!=station) {
            const AlulaCell *previous=station>1 && cursor->station==station-1?&cursor->cells[(station-1)%2]:NULL;
            if (initialize_cell(b,state,station-1,g,j,current,previous)!=ALULA_CELL_OK) {
                cursor->station=0;
                return NAN;
            }
            cursor->station=station;
        }
        c=current;
    }
    const AlulaCellValue *v=&c->value;
    double q=.5*(v->state[0].speed+v->state[1].speed);
    double mach=.5*(v->state[0].mach_sq+v->state[1].mach_sq);
    return mach<1?q+.2*q*mach*(mach-1)*v->geometry.curvature:q;
}

static AlulaBodyStatus prepare(AlulaBody *b,const AlulaBodyState *s,double *max_mach)
{
    int nx=b->in.nx, nt=b->size.tubes, bodies=b->in.bodies, bad=0;
    /* Nearby line-search trials often fail on the same stencil. Retest
     * it first with the ordinary initializer and skip the full sweep
     * only when the new state also fails that exact check. */
    if (b->failed_initial_cell>=0) {
        int i=b->failed_initial_cell/nt, tube=b->failed_initial_cell%nt, g=0;
        while (tube>=b->column[g+1]) ++g;
        AlulaCell scratch;
        if (initialize_cell(b,s,i,g,tube-b->column[g],&scratch,NULL)!=ALULA_CELL_OK)
            return ALULA_BODY_CELL;
        b->failed_initial_cell=-1;
    }
    /* Consecutive stencils share a quadrilateral in each tube. Keep a
     * group of stations on one thread so only its first station has to
     * check both quadrilaterals. The six-point middle is always checked. */
    const int block=8;
    int failed_cell=INT_MAX;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(min:failed_cell) if(nx*nt>1024)
#endif
    for (int first=1;first<nx;first+=block) {
        int end=first+block<nx? first+block : nx;
        for (int i=first;i<end && failed_cell==INT_MAX;++i) {
            for (int g=0;g<=bodies && failed_cell==INT_MAX;++g) {
                for (int j=0;j<b->tubes[g];++j) {
                    int k=index_section(b,g,i-1,j);
                    AlulaCell *c=b->cells+k;
                    const AlulaCell *previous=i>first?c-nt:NULL;
                    if (initialize_cell(b,s,i-1,g,j,c,previous)!=ALULA_CELL_OK) {
                        failed_cell=k;
                        break;
                    }
                    for (int q=i==1?0:1;q<2;++q) {
                        AlulaSection gas=c->value.state[q];
                        b->sections[k+q*nt]=(Section){gas,c->value.geometry.length[q],gas.speed};
                    }
                }
            }
        }
    }
    if (failed_cell!=INT_MAX) {
        b->failed_initial_cell=failed_cell;
        return ALULA_BODY_CELL;
    }

    double max_sq=0;
    b->check_all_upwind_tangents=0;
    for (int k=0;k<nx;++k) for (int j=0;j<nt;++j) {
        Section *q=b->sections+k*nt+j;
        double m=q->gas.mach_sq;
        if (m>=1 && (!b->in.upwind || k==0 || k==nx-1))
            return k==0 || k==nx-1 ? ALULA_BODY_BOUNDARY : ALULA_BODY_CELL;
        max_sq=fmax(max_sq,m);
        if (!b->in.upwind || k<2) continue;
        const Section *a=q-2*nt,*z=q-nt;
        AlulaUpwindInput in={.speed={a->gas.speed,z->gas.speed,q->gas.speed},
            .mach_sq={z->gas.mach_sq,m},.spacing={.5*(a->length+z->length),.5*(z->length+q->length)},
            .mucon=b->in.mucon,.mcrit=b->in.mcrit,.gamma=b->in.gamma};
        AlulaUpwindValue value;
        AlulaUpwindStencil *stencil=b->upwind_stencils+(k-2)*nt+j;
        if (alula_upwind_prepare(&in,&value,stencil)!=ALULA_CELL_OK) return ALULA_BODY_CELL;
        if (stencil->kink || !isfinite(stencil->slope)) b->check_all_upwind_tangents=1;
        q->transport=value.speed;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(|:bad) if(nx*nt>1024)
#endif
    for (int k=0;k<(nx-1)*nt;++k) {
        if (bad) continue;
        AlulaCell *c=b->cells+k;
        if (b->in.upwind) {
            c->input.biased=1;
            c->input.transport[0]=b->sections[k].transport;
            c->input.transport[1]=b->sections[k+nt].transport;
        }
        if (alula_cell_flux(c)!=ALULA_CELL_OK) bad=1;
    }
    *max_mach=sqrt(max_sq);
    return bad ? ALULA_BODY_CELL : ALULA_BODY_OK;
}

int alula_body_streamwise(const AlulaBody *b,Row r,const AlulaCellValue *d,double *value)
{
    const AlulaCellValue *v=&cell(b,r.g,r.i,r.j)->value;
    int entropy=b->in.mode==ALULA_BODY_ISENTROPIC;
    if (b->in.mode==ALULA_BODY_HYBRID) switch(b->in.ismom) {
        case 2: entropy=1; break;
        case 3: {
            int tube=b->column[r.g]+r.j;
            for (int g=0;g<b->in.bodies;++g) {
                AlulaBodyCut c=b->cuts[g];
                int cut=b->column[g+1];
                int i=c.through ? c.through : c.leading+10;
                int first=c.through ? c.lower_tube : cut-4;
                int end=c.through ? c.upper_tube : cut+4;
                if (r.i<=i && tube>=first && tube<end) { entropy=1; break; }
            }
            break;
        }
        case 4: {
            int k=index_section(b,r.g,r.i-1,r.j);
            /* The primal blend has already selected a pure equation at
             * either endpoint. Only mixed cells need the blend derivative. */
            if (d && b->blend_regime[k]!=BLEND_MIXED) {
                *value=(b->blend_regime[k]==BLEND_ENTROPY?d->isentropic:d->momentum)/b->p_inf;
                return 1;
            }
            AlulaBlendValue blend, db;
            if (alula_streamtube_blend(v,b->in.epsilon,&blend,d,0,&db)!=ALULA_CELL_OK) return 0;
            if (!d) b->blend_regime[k]=blend.fraction==0?BLEND_ENTROPY:
                blend.fraction==1?BLEND_MOMENTUM:BLEND_MIXED;
            *value=(d ? db.residual : blend.residual)/b->p_inf;
            return 1;
        }
    }
    if (d) v=d;
    *value=(entropy ? v->isentropic : v->momentum)/b->p_inf;
    return 1;
}

static int tangency(const AlulaBody *b,const double *coeff,AlulaPoint a,AlulaPoint z,double *v)
{
    Farfield f,g;
    AlulaPoint edge=sub(z,a);
    double length=norm(edge);
    if (!positive(length) || !alula_body_farfield(b,coeff,a,&f) || !alula_body_farfield(b,coeff,z,&g)) return 0;
    *v=cross(edge,mean(f.velocity,g.velocity))/length;
    return isfinite(*v);
}

static int outer_match(const AlulaBody *b,const AlulaBodyState *s,const double *coeff,double *result)
{
    double match[3]={0},scale_sq[3]={0};
    int top=b->in.bodies;
    for (int side=0;side<2;++side) {
        int g=side?top:0,j=side?b->tubes[top]:0;
        for (int i=0;i<b->in.nx;++i) {
            AlulaPoint a=point(b,s,g,i,j),z=point(b,s,g,i+1,j),edge=sub(z,a);
            double length=norm(edge);
            Farfield f;
            if (!positive(length) || !alula_body_farfield(b,coeff,mean(a,z),&f)) return 0;
            double q=b->sections[index_section(b,g,i,side?j-1:0)].gas.speed;
            AlulaPoint u=scale(edge,q/length);
            double defect=cross(f.velocity,u);
            for (int k=0;k<3;++k) {
                double mode=cross(f.basis[k+1],u);
                match[k]+=length*defect*mode;
                scale_sq[k]+=length*mode*mode;
            }
        }
    }
    for (int k=0;k<3;++k) {
        if (!positive(scale_sq[k])) return 0;
        result[k]=match[k]/sqrt(scale_sq[k]);
    }
    return 1;
}

static AlulaBodyStatus gap_row(const AlulaBody *b,const AlulaBodyState *s,Row r,double *value)
{
    AlulaPoint l[3],u[3];
    int idx[]={r.i-1,r.i,r.i<b->in.nx?r.i+1:r.i};
    double scale_length=b->in.length_scale;
    for (int k=0;k<3;++k) {
        l[k]=point(b,s,r.g,idx[k],b->tubes[r.g]);
        u[k]=point(b,s,r.g+1,idx[k],0);
        scale_length=fmax(scale_length,fmax(fabs(l[k].x),fabs(l[k].y)));
        scale_length=fmax(scale_length,fmax(fabs(u[k].x),fabs(u[k].y)));
    }
    AlulaPoint edge=mean(sub(l[2],l[0]),sub(u[2],u[0]));
    double length=norm(edge);
    if (!positive(length)) return ALULA_BODY_GAP;
    double gap=cross(edge,sub(u[1],l[1]))/length;
    if (gap < -64*DBL_EPSILON*scale_length) return ALULA_BODY_GAP;
    int k=b->wake[r.g]+r.i-b->cuts[r.g].trailing-1;
    if (!isfinite(s->wake_width[k]) || s->wake_width[k]<0) return ALULA_BODY_INVALID;
    *value=(gap-s->wake_width[k])/b->in.length_scale;
    return ALULA_BODY_OK;
}

AlulaBodyStatus alula_body_residual(AlulaBody *b,const AlulaBodyState *s,
    double *r,AlulaBodyDiagnostics *diag)
{
    if (!b || !s || !s->nodes || !s->density || !s->mass || !r ||
        (b->size.wake_widths && !s->wake_width) || !isfinite(s->circulation) ||
        !isfinite(s->source) || !isfinite(s->doublet_x) || !isfinite(s->doublet_y))
        return ALULA_BODY_INVALID;
    AlulaBodyDiagnostics local={0};
    ++b->residual_generation;
    AlulaBodyStatus status=prepare(b,s,&local.max_mach);
    if (status!=ALULA_BODY_OK) return status;

    double coeff[]={-s->circulation,s->source,s->doublet_x,s->doublet_y,s->circulation*s->circulation};
    int nx=b->in.nx,bodies=b->in.bodies;
    double match[3],bank_mean[MAX_BODIES];
    if (!outer_match(b,s,coeff,match)) return ALULA_BODY_BOUNDARY;
    for (int g=0;g<bodies;++g) {
        double l,u;
        if (!tangency(b,coeff,point(b,s,g,nx-1,b->tubes[g]),point(b,s,g,nx,b->tubes[g]),&l)
            || !tangency(b,coeff,point(b,s,g+1,nx-1,0),point(b,s,g+1,nx,0),&u)) return ALULA_BODY_BOUNDARY;
        bank_mean[g]=.5*(l+u);
    }
    for (int k=0;k<b->size.rows;++k) {
        Row row=b->rows[k];
        int i=row.i,g=row.g,j=row.j;
        double value=0;
        switch (row.family) {
        case ALULA_BODY_STREAMWISE:
            if (!alula_body_streamwise(b,row,NULL,&value)) return ALULA_BODY_NONFINITE;
            break;
        case ALULA_BODY_INLET_DENSITY: {
            AlulaSection gas=b->sections[index_section(b,g,0,j)].gas;
            value=log(gas.density/b->rho_total)-log(gas.enthalpy/b->h0)/(b->in.gamma-1);
            break;
        }
        case ALULA_BODY_INTERNAL_PRESSURE:
            value=(cell(b,g,i,j-1)->value.pressure[1]-cell(b,g,i,j)->value.pressure[0])/b->p_inf;
            break;
        case ALULA_BODY_FARFIELD_PRESSURE: {
            Farfield f;
            if (!alula_body_farfield(b,coeff,point(b,s,g,i,j),&f)) return ALULA_BODY_BOUNDARY;
            int side=g==bodies;
            value=(cell(b,g,i,side?j-1:0)->value.pressure[side]-f.pressure)/b->p_inf;
            break;
        }
        case ALULA_BODY_CUT_PRESSURE:
        case ALULA_BODY_TRAILING_KUTTA:
        case ALULA_BODY_LEADING_KUTTA:
            value=(cell(b,g+1,i,0)->value.pressure[0]-cell(b,g,i,b->tubes[g]-1)->value.pressure[1])/b->p_inf;
            break;
        case ALULA_BODY_WAKE_GAP:
            status=gap_row(b,s,row,&value);
            if (status!=ALULA_BODY_OK) return status;
            break;
        case ALULA_BODY_END_TANGENCY: {
            int outlet=i==nx, cut=g<bodies && j==b->tubes[g];
            if (outlet && cut && b->in.bank_outlet) {
                value=bank_mean[g];
                break;
            }
            int first=outlet?nx-1:0;
            AlulaPoint a=point(b,s,g,first,j),z=point(b,s,g,first+1,j);
            if (outlet && cut && b->in.displaced) {
                a=mean(a,point(b,s,g+1,first,0));
                z=mean(z,point(b,s,g+1,first+1,0));
            }
            if (!tangency(b,coeff,a,z,&value)) return ALULA_BODY_BOUNDARY;
            break;
        }
        case ALULA_BODY_FARFIELD_MATCH:
            value=match[j];
            break;
        case ALULA_BODY_FAMILIES:
            return ALULA_BODY_INVALID;
        }
        if (!isfinite(value)) return ALULA_BODY_NONFINITE;
        r[k]=value;
        local.residual[row.family]=fmax(local.residual[row.family],fabs(value));
    }
    if (diag) *diag=local;
    return ALULA_BODY_OK;
}
