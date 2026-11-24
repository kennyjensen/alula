// SPDX-License-Identifier: GPL-2.0-or-later
#include "channel.h"
#include "streamtube_internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    AlulaSection state;
    AlulaPoint direction, gap;
    double length, area, transport, arc;
    AlulaUpwindInput filter;
} Section;

typedef struct { int column, station, bank; } Variable;

struct AlulaChannel {
    AlulaChannelInput in;
    int nx, nt, n, nd, base, reservoir;
    double height, max_mach_sq;
    double *fixed, *fraction, *y, *density, *mass, *work;
    AlulaCell *cells;
    Section *sections;
    AlulaSparse *matrix;
    int *slots;
};

static double norm_inf(const double *a, int n)
{
    double norm=0;
    for (int i=0; i<n; ++i) {
        if (!isfinite(a[i])) return INFINITY;
        norm=fmax(norm,fabs(a[i]));
    }
    return norm;
}

static int positive(double x) { return isfinite(x) && x > 0; }
static int node_index(const AlulaChannel *c, int i, int j) { return i*(c->nt+1)+j; }
static int position(const AlulaChannel *c, int i, int j) { return c->nd+i*(c->nt-1)+j-1; }
static int normal_row(const AlulaChannel *c, int i, int j) { return c->nd+(i-1)*(c->nt-1)+j; }

static const double *copy(double **next, const double *source, int n)
{
    double *dest=*next;
    if (source) memcpy(dest,source,(size_t)n*sizeof *dest);
    *next+=n;
    return dest;
}

static int variables(const AlulaChannel *c, int i, int j, Variable out[15])
{
    int first=c->reservoir && i>=3 ? i-3 : c->reservoir ? 0 : i-1, count=0;
    if (c->reservoir) out[count++]=(Variable){c->base+j,0,-2};
    for (int k=first; k<=i; ++k) out[count++]=(Variable){k*c->nt+j,k,-1};
    for (int bank=0; bank<2; ++bank) if (j+bank>0 && j+bank<c->nt)
        for (int k=first; k<=i+1; ++k) out[count++]=(Variable){position(c,k,j+bank),k,bank};
    return count;
}

enum { STREAM_ROW, INLET_ROW, UPPER_ROW, LOWER_ROW, OUTLET_ROW, CELL_ROWS };
static void rows(const AlulaChannel *c, int i, int j, int r[CELL_ROWS])
{
    r[STREAM_ROW]=(i-1)*c->nt+j;
    r[INLET_ROW]=i==1 ? (c->nx-1)*c->nt+j : -1;
    r[UPPER_ROW]=j<c->nt-1 ? normal_row(c,i,j) : -1;
    r[LOWER_ROW]=j>0 ? normal_row(c,i,j-1) : -1;
    r[OUTLET_ROW]=c->reservoir && i==c->nx-1 ? c->base+j : -1;
}

static int pattern(AlulaChannel *c)
{
    size_t count=0, cap=(size_t)(c->nx-1)*c->nt*15*CELL_ROWS+4*(size_t)c->nt;
    AlulaSparseEntry *entries=malloc(cap*sizeof *entries);
    if (!entries) return 0;
    for (int i=1; i<c->nx; ++i) for (int j=0; j<c->nt; ++j) {
        Variable vars[15];
        int r[CELL_ROWS], nv=variables(c,i,j,vars);
        rows(c,i,j,r);
        for (int k=0; k<nv; ++k)
            for (int q=0; q<CELL_ROWS; ++q)
                if (r[q]>=0) entries[count++]=(AlulaSparseEntry){r[q],vars[k].column};
    }
    int end=normal_row(c,c->nx,0);
    for (int j=1; j<c->nt; ++j,++end) {
        entries[count++]=(AlulaSparseEntry){end,position(c,0,j)};
        entries[count++]=(AlulaSparseEntry){end,position(c,1,j)};
    }
    for (int j=1; j<c->nt; ++j,++end) {
        if (c->reservoir) entries[count++]=(AlulaSparseEntry){end,position(c,0,j)};
        else {
            entries[count++]=(AlulaSparseEntry){end,position(c,c->nx-1,j)};
            entries[count++]=(AlulaSparseEntry){end,position(c,c->nx,j)};
        }
    }
    c->slots=malloc(count*sizeof *c->slots);
    if (c->slots) c->matrix=alula_sparse_create(c->n,entries,count,c->slots);
    free(entries);
    return c->matrix != NULL && end == c->base;
}

AlulaChannel *alula_channel_create(const AlulaChannelInput *in)
{
    if (!in || in->nx<2 || in->nt<1 || (size_t)in->nx*in->nt>100000
        || in->mode<ALULA_CHANNEL_MOMENTUM || in->mode>ALULA_CHANNEL_HYBRID
        || !in->x || !in->lower || !in->upper || !in->mass || !in->h0 || !in->rho_total
        || !isfinite(in->gamma) || in->gamma<=1 || !positive(in->rho_ref) || !positive(in->p_ref)
        || !isfinite(in->correction) || in->correction<0) return NULL;
    int nx=in->nx, nt=in->nt, reservoir=in->mode>=ALULA_CHANNEL_UPWIND;
    if (reservoir && (nx<3 || !in->outlet_pressure || in->outlet_slope
        || !isfinite(in->mucon) || !isfinite(in->mcrit) || in->mcrit<0 || in->mcrit>1
        || (in->mode==ALULA_CHANNEL_HYBRID && !positive(in->epsilon)))) return NULL;
    for (int i=0; i<=nx; ++i)
        if (!isfinite(in->x[i]) || (i && !positive(in->x[i]-in->x[i-1]))
            || !positive(in->upper[i]-in->lower[i])) return NULL;
    double total_mass=0;
    for (int j=0; j<nt; ++j) {
        if (!positive(in->mass[j]) || !positive(in->h0[j]) || !positive(in->rho_total[j])
            || (reservoir && !positive(in->outlet_pressure[j]))) return NULL;
        total_mass+=in->mass[j];
    }
    if (!positive(total_mass)) return NULL;
    for (int j=0; j<nt-1; ++j)
        if ((in->inlet_slope && !isfinite(in->inlet_slope[j])) ||
            (in->outlet_slope && !isfinite(in->outlet_slope[j]))) return NULL;
    AlulaChannel *c=calloc(1,sizeof *c);
    if (!c) return NULL;
    c->in=*in;
    c->nx=nx; c->nt=nt; c->nd=nx*nt; c->reservoir=reservoir;
    c->base=c->nd+(nx+1)*(nt-1);
    c->n=c->base+(reservoir ? nt : 0);
    c->height=in->upper[0]-in->lower[0];
    c->fixed=calloc(3*(size_t)(nx+1)+7*(size_t)nt,sizeof *c->fixed);
    c->y=malloc((size_t)(nx+1)*(nt+1)*sizeof *c->y);
    c->density=malloc((size_t)c->nd*sizeof *c->density);
    c->mass=malloc((size_t)nt*sizeof *c->mass);
    c->cells=malloc((size_t)(nx-1)*nt*sizeof *c->cells);
    c->sections=calloc((size_t)c->nd,sizeof *c->sections);
    c->work=malloc((size_t)5*c->n*sizeof *c->work);
    if (!c->fixed || !c->y || !c->density || !c->mass || !c->cells || !c->sections || !c->work) goto fail;
    double *next=c->fixed;
    c->in.x=copy(&next,in->x,nx+1);
    c->in.lower=copy(&next,in->lower,nx+1);
    c->in.upper=copy(&next,in->upper,nx+1);
    c->in.mass=copy(&next,in->mass,nt);
    c->in.h0=copy(&next,in->h0,nt);
    c->in.rho_total=copy(&next,in->rho_total,nt);
    c->in.inlet_slope=copy(&next,in->inlet_slope,nt-1);
    c->in.outlet_slope=copy(&next,in->outlet_slope,nt-1);
    c->in.outlet_pressure=copy(&next,reservoir ? in->outlet_pressure : NULL,nt);
    c->fraction=next;
    /* Fractions use the reference mass, even when the reservoir solve changes mass. */
    for (int j=1; j<nt; ++j) c->fraction[j]=c->fraction[j-1]+in->mass[j-1]/total_mass;
    if (!pattern(c)) goto fail;
    return c;
fail:
    alula_channel_free(c);
    return NULL;
}

void alula_channel_free(AlulaChannel *c)
{
    if (!c) return;
    alula_sparse_free(c->matrix);
    free(c->slots); free(c->fixed); free(c->y); free(c->density);
    free(c->mass); free(c->cells); free(c->sections); free(c->work);
    free(c);
}

int alula_channel_size(const AlulaChannel *c) { return c ? c->n : 0; }

static int prepare(AlulaChannel *c, const double *state)
{
    const AlulaChannelInput *in=&c->in;
    int nx=c->nx, nt=c->nt;
    for (int k=0; k<c->nd; ++k) {
        c->density[k]=in->rho_ref*exp(state[k]);
        if (!positive(c->density[k])) return 0;
    }
    for (int j=0; j<nt; ++j) {
        c->mass[j]=c->reservoir ? in->mass[j]*exp(state[c->base+j]) : in->mass[j];
        if (!positive(c->mass[j])) return 0;
    }
    for (int i=0; i<=nx; ++i) {
        double *y=c->y+node_index(c,i,0);
        y[0]=in->lower[i]; y[nt]=in->upper[i];
        for (int j=1; j<nt; ++j)
            y[j]=y[0]+c->fraction[j]*(y[nt]-y[0])+c->height*state[position(c,i,j)];
    }
    for (int i=1; i<nx; ++i) for (int j=0; j<nt; ++j) {
        AlulaCellInput p={.density={c->density[(i-1)*nt+j],c->density[i*nt+j]},
            .mass=c->mass[j],.h0=in->h0[j],.gamma=in->gamma,.correction=in->correction};
        for (int k=0; k<3; ++k) {
            int station=i-1+k, n=node_index(c,station,j);
            p.lower[k]=(AlulaPoint){in->x[station],c->y[n]};
            p.upper[k]=(AlulaPoint){in->x[station],c->y[n+1]};
        }
        AlulaCell *cell=c->cells+(i-1)*nt+j;
        if (alula_cell_initialize(cell,&p)!=ALULA_CELL_OK) return 0;
        const AlulaCellGeometry *g=&cell->value.geometry;
        for (int k=i==1 ? 0 : 1; k<2; ++k) {
            AlulaSection a=cell->value.state[k];
            c->sections[(i-1+k)*nt+j]=(Section){.state=a,.direction=g->direction[k],
                .gap=g->section[k],.length=g->length[k],.area=g->normal_area[k],.transport=a.speed};
        }
    }
    c->max_mach_sq=0;
    for (int j=0; j<nt; ++j) {
        double arc=0;
        for (int k=0; k<nx; ++k) {
            Section *s=c->sections+k*nt+j;
            double m=s->state.mach_sq;
            if (m>=1 && (!c->reservoir || k==0 || k==nx-1)) return 0;
            c->max_mach_sq=fmax(c->max_mach_sq,m);
            s->arc=arc+.5*s->length;
            arc+=s->length;
            if (!c->reservoir || k<2) continue;
            const Section *a=s-2*nt, *b=s-nt;
            s->filter=(AlulaUpwindInput){.speed={a->state.speed,b->state.speed,s->state.speed},
                .mach_sq={b->state.mach_sq,m},.spacing={b->arc-a->arc,s->arc-b->arc},
                .mucon=in->mucon,.mcrit=in->mcrit,.gamma=in->gamma};
            AlulaUpwindValue value;
            if (alula_streamtube_upwind(&s->filter,&value,NULL,NULL)!=ALULA_CELL_OK) return 0;
            s->transport=value.speed;
        }
    }
    for (int i=1; i<nx; ++i) for (int j=0; j<nt; ++j) {
        AlulaCell *cell=c->cells+(i-1)*nt+j;
        if (c->reservoir) {
            cell->input.biased=1;
            cell->input.transport[0]=c->sections[(i-1)*nt+j].transport;
            cell->input.transport[1]=c->sections[i*nt+j].transport;
        }
        if (alula_cell_flux(cell)!=ALULA_CELL_OK) return 0;
    }
    return 1;
}

static int residual(AlulaChannel *c, const double *state, double *r)
{
    if (!isfinite(norm_inf(state,c->n)) || !prepare(c,state)) return 0;
    const AlulaChannelInput *in=&c->in;
    int nt=c->nt, nx=c->nx;
    for (int i=1; i<nx; ++i) for (int j=0; j<nt; ++j) {
        int k=(i-1)*nt+j;
        const AlulaCellValue *v=&c->cells[k].value;
        double row=in->mode==ALULA_CHANNEL_ISENTROPIC ? v->isentropic : v->momentum;
        if (in->mode==ALULA_CHANNEL_HYBRID) {
            AlulaBlendValue blend;
            if (alula_streamtube_blend(v,in->epsilon,&blend,NULL,0,NULL)!=ALULA_CELL_OK) return 0;
            row=blend.residual;
        }
        r[k]=row/in->p_ref;
        if (i==1) {
            AlulaSection s=v->state[0];
            r[(nx-1)*nt+j]=log(s.density/in->rho_total[j])-log(s.enthalpy/in->h0[j])/(in->gamma-1);
        }
        if (j<nt-1) r[normal_row(c,i,j)]=(v->pressure[1]-c->cells[k+1].value.pressure[0])/in->p_ref;
        if (c->reservoir && i==nx-1)
            r[c->base+j]=(v->state[1].pressure-in->outlet_pressure[j])/in->p_ref;
    }
    int end=normal_row(c,nx,0);
    for (int j=1; j<nt; ++j)
        r[end++]=(c->y[node_index(c,1,j)]-c->y[node_index(c,0,j)]
            -in->inlet_slope[j-1]*(in->x[1]-in->x[0]))/c->height;
    for (int j=1; j<nt; ++j)
        r[end++]=c->reservoir ? state[position(c,0,j)] :
            (c->y[node_index(c,nx,j)]-c->y[node_index(c,nx-1,j)]
            -in->outlet_slope[j-1]*(in->x[nx]-in->x[nx-1]))/c->height;
    return isfinite(norm_inf(r,c->n));
}

int alula_channel_residual(AlulaChannel *c, const double *state, double *r)
{
    return c && state && r && state!=r && residual(c,state,r);
}

typedef struct { double speed, mach_sq, length; } SectionDirection;
static SectionDirection section_direction(const AlulaChannel *c, int k, int j, Variable var)
{
    const Section *s=c->sections+k*c->nt+j;
    AlulaSection v=s->state;
    double dlength=0, darea=0, density=0, mass=0;
    if (var.bank>=0 && (var.station==k || var.station==k+1)) {
        double dy=.5*c->height*(var.station==k ? -1 : 1);
        double dg=.5*c->height*(var.bank==0 ? -1 : 1);
        dlength=s->direction.y*dy;
        double dx=-s->direction.x*dlength/s->length;
        double dt=(dy-s->direction.y*dlength)/s->length;
        darea=dx*s->gap.y-dt*s->gap.x+s->direction.x*dg;
    } else if (var.bank==-1 && var.station==k) density=v.density;
    else if (var.bank==-2) mass=c->mass[j];
    double q=(-v.speed/v.density)*density+(v.speed/c->mass[j])*mass+(-v.speed/s->area)*darea;
    double h=-v.speed*q, gm=c->in.gamma-1;
    return (SectionDirection){q,(2*v.speed*q-v.mach_sq*gm*h)/(gm*v.enthalpy),dlength};
}

static int transport_direction(const AlulaChannel *c, int k, int j, Variable var, double *out)
{
    SectionDirection v=section_direction(c,k,j,var);
    if (k<2) { *out=v.speed; return 1; }
    SectionDirection a=section_direction(c,k-2,j,var), b=section_direction(c,k-1,j,var);
    AlulaUpwindInput d={.speed={a.speed,b.speed,v.speed},.mach_sq={b.mach_sq,v.mach_sq},
        .spacing={.5*(a.length+b.length),.5*(b.length+v.length)}};
    AlulaUpwindValue value, dv;
    if (alula_streamtube_upwind(&c->sections[k*c->nt+j].filter,&value,&d,&dv)!=ALULA_CELL_OK) return 0;
    *out=dv.speed;
    return 1;
}

static int jacobian(AlulaChannel *c)
{
    const AlulaChannelInput *in=&c->in;
    double *a=alula_sparse_values(c->matrix);
    memset(a,0,(size_t)alula_sparse_nonzeros(c->matrix)*sizeof *a);
    const int *slot=c->slots;
    for (int i=1; i<c->nx; ++i) for (int j=0; j<c->nt; ++j) {
        const AlulaCell *cell=c->cells+(i-1)*c->nt+j;
        Variable vars[15];
        int r[CELL_ROWS], nv=variables(c,i,j,vars);
        rows(c,i,j,r);
        for (int k=0; k<nv; ++k) {
            Variable var=vars[k];
            int local=var.station-(i-1);
            AlulaCellDirection d={0};
            if (var.bank==-2) d.mass=c->mass[j];
            else if (local>=0) {
                if (var.bank==-1) d.density[local]=cell->input.density[local];
                else (var.bank==0 ? d.lower : d.upper)[local].y=c->height;
            }
            if (c->reservoir)
                for (int q=0; q<2; ++q)
                    if (!transport_direction(c,i-1+q,j,var,&d.transport[q])) return 0;
            AlulaCellValue dv;
            if (alula_cell_apply(cell,&d,&dv)!=ALULA_CELL_OK) return 0;
            double value=in->mode==ALULA_CHANNEL_ISENTROPIC ? dv.isentropic : dv.momentum;
            if (in->mode==ALULA_CHANNEL_HYBRID) {
                AlulaBlendValue blend, db;
                if (alula_streamtube_blend(&cell->value,in->epsilon,&blend,&dv,0,&db)!=ALULA_CELL_OK) return 0;
                value=db.residual;
            }
            a[*slot++]+=value/in->p_ref;
            if (r[INLET_ROW]>=0) {
                AlulaSection s=cell->value.state[0], ds=dv.state[0];
                a[*slot++]+=ds.density/s.density-ds.enthalpy/((in->gamma-1)*s.enthalpy);
            }
            if (r[UPPER_ROW]>=0) a[*slot++]+=dv.pressure[1]/in->p_ref;
            if (r[LOWER_ROW]>=0) a[*slot++]-=dv.pressure[0]/in->p_ref;
            if (r[OUTLET_ROW]>=0) a[*slot++]+=dv.state[1].pressure/in->p_ref;
        }
    }
    for (int j=1; j<c->nt; ++j) { a[*slot++]-=1; a[*slot++]+=1; }
    for (int j=1; j<c->nt; ++j) {
        if (!c->reservoir) a[*slot++]-=1;
        a[*slot++]+=1;
    }
    return 1;
}

int alula_channel_jacobian_product(AlulaChannel *c, const double *state,
    const double *direction, double *product)
{
    if (!c || !state || !direction || !product || product==direction
        || !isfinite(norm_inf(direction,c->n)) || !residual(c,state,c->work) || !jacobian(c)) return 0;
    alula_sparse_product(c->matrix,direction,product);
    return isfinite(norm_inf(product,c->n));
}

AlulaChannelResult alula_channel_solve(AlulaChannel *c, double *state, int iterations, double tolerance)
{
    AlulaChannelResult out={.status=ALULA_CHANNEL_INVALID,.residual=INFINITY,.max_mach=NAN};
    if (!c || !state || iterations<0 || !positive(tolerance)) return out;
    int n=c->n;
    double *r=c->work, *rhs=r+n, *delta=rhs+n, *trial=delta+n, *rt=trial+n;
    if (!residual(c,state,r)) return out;
    out.residual=norm_inf(r,n);
    out.max_mach=sqrt(c->max_mach_sq);
    out.status=ALULA_CHANNEL_LIMIT;
    for (int it=0; it<iterations && out.residual>tolerance; ++it) {
        for (int i=0; i<n; ++i) rhs[i]=-r[i];
        if (!jacobian(c) || !alula_sparse_solve(c->matrix,rhs,delta)) {
            out.status=ALULA_CHANNEL_LINEAR_FAILURE;
            break;
        }
        int accepted=0;
        for (double step=1; step>=0x1p-20; step*=.5) {
            for (int i=0; i<n; ++i) trial[i]=state[i]+step*delta[i];
            if (residual(c,trial,rt)) {
                double norm=norm_inf(rt,n);
                if (norm <= (1-1e-4*step)*out.residual) {
                    memcpy(state,trial,(size_t)n*sizeof *state);
                    double *swap=r;
                    r=rt; rt=swap;
                    out.iterations++;
                    out.residual=norm;
                    out.max_mach=sqrt(c->max_mach_sq);
                    accepted=1;
                    break;
                }
            }
            out.backtracks++;
        }
        if (!accepted) { out.status=ALULA_CHANNEL_LINE_SEARCH_FAILURE; break; }
    }
    if (out.residual<=tolerance) out.status=ALULA_CHANNEL_CONVERGED;
    out.linear=alula_sparse_stats(c->matrix);
    return out;
}
