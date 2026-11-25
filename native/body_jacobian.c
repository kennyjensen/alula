// SPDX-License-Identifier: GPL-2.0-or-later
#include "body_internal.h"

#include <stdlib.h>
#include <string.h>

typedef struct { double value, tangent; } Dual;
typedef struct { Farfield *next; int fill; } FarfieldPass;
static Dual dadd(Dual a,Dual b) { return (Dual){a.value+b.value,a.tangent+b.tangent}; }
static Dual dsub(Dual a,Dual b) { return (Dual){a.value-b.value,a.tangent-b.tangent}; }
static Dual dscale(Dual a,double s) { return (Dual){a.value*s,a.tangent*s}; }
static Dual dmul(Dual a,Dual b)
{
    return (Dual){a.value*b.value,a.tangent*b.value+a.value*b.tangent};
}

static Dual ddiv(Dual a,Dual b)
{
    double q=a.value/b.value;
    return (Dual){q,(a.tangent-q*b.tangent)/b.value};
}

static double dcross(AlulaPoint a,AlulaPoint b,AlulaPoint da,AlulaPoint db)
{
    return cross(da,b)+cross(a,db);
}

static AlulaPoint point_direction(const AlulaBody *b,const AlulaBodyState *d,int g,int i,int j)
{
    return d->nodes ? point(b,d,g,i,j) : (AlulaPoint){0};
}

static double field_pressure_direction(const AlulaBody *b,const Farfield *f,AlulaPoint dv)
{
    double temp=1+.5*(b->in.gamma-1)*b->mach_sq*(1-dot(f->velocity,f->velocity));
    return -b->in.gamma*f->pressure*b->mach_sq/temp*dot(f->velocity,dv);
}

static int field_direction(const AlulaBody *b,const double *coeff,const double *dc,
    AlulaPoint p,AlulaPoint dp,FarfieldPass *pass,const Farfield **primal,Farfield *df)
{
    Farfield *f=pass->next++;
    if (pass->fill && !alula_body_farfield(b,coeff,p,f)) return 0;
    if (primal) *primal=f;
    if (dp.x==0 && dp.y==0) {
        *df=(Farfield){0};
        if (dc[0]==0 && dc[1]==0 && dc[2]==0 && dc[3]==0 && dc[4]==0) return 1;
        /* Strength-only colors have no geometric derivative. */
        for (int k=0;k<5;++k) df->velocity=add(df->velocity,scale(f->basis[k],dc[k]));
        df->pressure=field_pressure_direction(b,f,df->velocity);
        return finite_point(df->velocity) && isfinite(df->pressure);
    }
    AlulaPoint d=sub(p,b->in.center), fs=b->freestream;
    Dual x={dot(d,fs)/b->beta,dot(dp,fs)/b->beta}, y={cross(fs,d),cross(fs,dp)};
    Dual xx=dmul(x,x), yy=dmul(y,y), xy=dmul(x,y), r2=dadd(xx,yy);
    Dual r4=dmul(r2,r2), r6=dmul(r4,r2);
    Dual log_radius={f->log_radius,.5*r2.tangent/r2.value};
    Dual delta=dsub(xx,yy), cubic=dmul(x,dsub(xx,dscale(yy,3)));
    Dual gx[]={
        dscale(ddiv(y,r2),-1/TWO_PI),
        dscale(ddiv(x,r2),1/TWO_PI),
        dscale(ddiv(delta,r4),-1/TWO_PI),
        dscale(ddiv(xy,r4),-2/TWO_PI),
        dscale(dadd(dscale(ddiv(dsub(xx,dmul(log_radius,delta)),r4),b->pg_a),
            dscale(dsub(ddiv(dscale(delta,3),r4),ddiv(dscale(dmul(x,cubic),4),r6)),b->pg_b)),b->pg_scale)
    };
    Dual gy[]={
        dscale(ddiv(x,r2),1/TWO_PI),
        dscale(ddiv(y,r2),1/TWO_PI),
        dscale(ddiv(xy,r4),-2/TWO_PI),
        dscale(ddiv(delta,r4),1/TWO_PI),
        dscale(dadd(dscale(ddiv(dmul(xy,dsub((Dual){1,0},dscale(log_radius,2))),r4),b->pg_a),
            dscale(dsub(ddiv(dscale(xy,-6),r4),ddiv(dscale(dmul(y,cubic),4),r6)),b->pg_b)),b->pg_scale)
    };
    df->velocity=(AlulaPoint){0};
    for (int k=0;k<5;++k) {
        double dx=gx[k].tangent/b->beta,dy=gy[k].tangent;
        df->basis[k]=(AlulaPoint){dx*fs.x-dy*fs.y,dx*fs.y+dy*fs.x};
        df->velocity=add(df->velocity,add(scale(df->basis[k],coeff[k]),scale(f->basis[k],dc[k])));
    }
    df->pressure=field_pressure_direction(b,f,df->velocity);
    return finite_point(df->velocity) && isfinite(df->pressure);
}

static SectionTangent section_direction(const AlulaBody *b,const AlulaBodyState *d,
    int g,int k,int j)
{
    int column=index_section(b,g,k,j),q=k>0, i=k-q+1;
    const AlulaCell *c=cell(b,g,i,j);
    const AlulaCellGeometry *v=&c->value.geometry;
    AlulaPoint a=point_direction(b,d,g,k,j), z=point_direction(b,d,g,k+1,j);
    AlulaPoint u=point_direction(b,d,g,k,j+1), w=point_direction(b,d,g,k+1,j+1);
    double rho=d->density ? d->density[column] : 0;
    double mass=d->mass ? d->mass[b->column[g]+j] : 0;
    if (rho==0 && mass==0 && a.x==0 && a.y==0 && z.x==0 && z.y==0 &&
        u.x==0 && u.y==0 && w.x==0 && w.y==0) return (SectionTangent){0};
    AlulaPoint seg=mean(sub(z,a),sub(w,u)), gap=mean(sub(u,a),sub(w,z)), t=v->direction[q];
    double length_derivative=dot(t,seg);
    AlulaPoint dt=scale(sub(seg,scale(t,length_derivative)),1/v->length[q]);
    double area=dcross(t,v->section[q],dt,gap);
    AlulaSection gas=c->value.state[q];
    double speed=(-gas.speed/gas.density)*rho+(gas.speed/c->input.mass)*mass
        -(gas.speed/v->normal_area[q])*area;
    double h=-gas.speed*speed, gm=b->in.gamma-1;
    double m=(2*gas.speed*speed-gas.mach_sq*gm*h)/(gm*gas.enthalpy);
    return (SectionTangent){.speed=speed,.mach_sq=m,.enthalpy=h,.length=length_derivative,.density=rho,.transport=speed};
}

static int prepare_tangents(AlulaBody *b,const AlulaBodyState *d,int first,int end)
{
    int nx=b->in.nx,nt=b->size.tubes,bodies=b->in.bodies,bad=0;
    if (!b->section_tangents) b->section_tangents=malloc((size_t)nx*nt*sizeof *b->section_tangents);
    if (!b->cell_tangents) b->cell_tangents=malloc((size_t)(nx-1)*nt*sizeof *b->cell_tangents);
    if (!b->cell_tangents || !b->section_tangents) return 0;
    /* A cell at i uses sections i-1 and i. Upwind transport can also
     * carry a section perturbation two sections downstream. */
    int first_cell=first>0 ? first : 1;
    int end_cell=end+1+(b->in.upwind ? 2 : 0);
    if (end_cell>nx) end_cell=nx;
    if (first>0 || end<nx) {
        memset(b->section_tangents,0,(size_t)first*nt*sizeof *b->section_tangents);
        memset(b->section_tangents+end*nt,0,(size_t)(nx-end)*nt*sizeof *b->section_tangents);
        memset(b->cell_tangents,0,(size_t)(first_cell-1)*nt*sizeof *b->cell_tangents);
        memset(b->cell_tangents+(end_cell-1)*nt,0,(size_t)(nx-end_cell)*nt*sizeof *b->cell_tangents);
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(|:bad) if((end-first)*nt>1024)
#endif
    for (int k=first;k<end;++k) for (int g=0;g<=bodies;++g) for (int j=0;j<b->tubes[g];++j) {
        SectionTangent t=section_direction(b,d,g,k,j);
        if (!isfinite(t.speed) || !isfinite(t.mach_sq) || !isfinite(t.enthalpy) || !isfinite(t.length)) bad=1;
        b->section_tangents[index_section(b,g,k,j)]=t;
    }
    if (bad) return 0;
    if (b->in.upwind) {
        /* Filtered transport reads the current and two upstream sections.
         * An exceptional primal stencil must still be checked even where the
         * direction is zero, as in a full Jacobian product. */
        int start=first>2 ? first : 2;
        int stop=end+2<nx ? end+2 : nx;
        if (b->check_all_upwind_tangents) { start=2; stop=nx; }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(|:bad) if((stop-start)*nt>1024)
#endif
        for (int k=start*nt;k<stop*nt;++k) {
            SectionTangent *ds=b->section_tangents+k;
            const SectionTangent *da=ds-2*nt,*dz=ds-nt;
            AlulaUpwindInput dir={.speed={da->speed,dz->speed,ds->speed},
                .mach_sq={dz->mach_sq,ds->mach_sq},.spacing={.5*(da->length+dz->length),.5*(dz->length+ds->length)}};
            const AlulaUpwindStencil *stencil=b->upwind_stencils+(k-2*nt);
            double transport=alula_upwind_speed_direction(stencil,&dir);
            if (stencil->kink || !isfinite(transport)) bad=1;
            else ds->transport=transport;
        }
    }
    if (bad) return 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(|:bad) if((end_cell-first_cell)*nt>1024)
#endif
    for (int i=first_cell;i<end_cell;++i) for (int g=0;g<=bodies;++g) for (int j=0;j<b->tubes[g];++j) {
        int k=index_section(b,g,i-1,j);
        AlulaCellDirection dir={.density={b->section_tangents[k].density,b->section_tangents[k+nt].density},
            .mass=d->mass ? d->mass[b->column[g]+j] : 0};
        for (int q=0;q<3;++q) {
            dir.lower[q]=point_direction(b,d,g,i-1+q,j);
            dir.upper[q]=point_direction(b,d,g,i-1+q,j+1);
        }
        if (b->in.upwind) {
            dir.transport[0]=b->section_tangents[k].transport;
            dir.transport[1]=b->section_tangents[k+nt].transport;
        }
        CellTangent *t=b->cell_tangents+k;
        int active=dir.mass!=0 || dir.density[0]!=0 || dir.density[1]!=0 ||
            dir.transport[0]!=0 || dir.transport[1]!=0;
        for (int q=0;q<3;++q)
            active |= dir.lower[q].x!=0 || dir.lower[q].y!=0 || dir.upper[q].x!=0 || dir.upper[q].y!=0;
        if (!active) { *t=(CellTangent){0}; continue; }
        AlulaCellValue dv;
        if (alula_cell_apply(b->cells+k,&dir,&dv)!=ALULA_CELL_OK ||
            !alula_body_streamwise(b,(Row){ALULA_BODY_STREAMWISE,i,g,j},&dv,&t->streamwise)) {bad=1; continue;}
        t->pressure[0]=dv.pressure[0];
        t->pressure[1]=dv.pressure[1];
    }
    return !bad;
}

static int end_direction(const AlulaBody *b,const double *coeff,const double *dc,
    AlulaPoint a,AlulaPoint z,AlulaPoint da,AlulaPoint dz,FarfieldPass *pass,double *value)
{
    const Farfield *fa,*fz;
    Farfield dfa,dfz;
    AlulaPoint edge=sub(z,a), de=sub(dz,da);
    double length=norm(edge);
    if (!positive(length) || !field_direction(b,coeff,dc,a,da,pass,&fa,&dfa) ||
        !field_direction(b,coeff,dc,z,dz,pass,&fz,&dfz)) return 0;
    AlulaPoint velocity=mean(fa->velocity,fz->velocity), dv=mean(dfa.velocity,dfz.velocity);
    double dn=dot(edge,de)/length;
    *value=(dcross(edge,velocity,de,dv)-cross(edge,velocity)*dn/length)/length;
    return isfinite(*value);
}

static int boundary_direction(const AlulaBody *b,const AlulaBodyState *s,const AlulaBodyState *d,
    const double *coeff,const double *dc,Row row,FarfieldPass *pass,double *value)
{
    int g=row.g,j=row.j,outlet=row.i==b->in.nx,cut=g<b->in.bodies && j==b->tubes[g];
    int i=outlet?b->in.nx-1:0;
    AlulaPoint a=point(b,s,g,i,j),z=point(b,s,g,i+1,j);
    AlulaPoint da=point_direction(b,d,g,i,j),dz=point_direction(b,d,g,i+1,j);
    if (outlet && cut && b->in.displaced) {
        AlulaPoint u=point(b,s,g+1,i,0),w=point(b,s,g+1,i+1,0);
        AlulaPoint du=point_direction(b,d,g+1,i,0),dw=point_direction(b,d,g+1,i+1,0);
        if (b->in.bank_outlet) {
            double lower,upper;
            if (!end_direction(b,coeff,dc,a,z,da,dz,pass,&lower) ||
                !end_direction(b,coeff,dc,u,w,du,dw,pass,&upper)) return 0;
            *value=.5*(lower+upper);
            return 1;
        }
        a=mean(a,u); z=mean(z,w); da=mean(da,du); dz=mean(dz,dw);
    }
    return end_direction(b,coeff,dc,a,z,da,dz,pass,value);
}

static double gap_direction(const AlulaBody *b,const AlulaBodyState *s,const AlulaBodyState *d,Row r)
{
    int j=b->tubes[r.g], i=r.i, next=i<b->in.nx ? i+1 : i;
    AlulaPoint edge=mean(sub(point(b,s,r.g,next,j),point(b,s,r.g,i-1,j)),
        sub(point(b,s,r.g+1,next,0),point(b,s,r.g+1,i-1,0)));
    AlulaPoint de=mean(sub(point_direction(b,d,r.g,next,j),point_direction(b,d,r.g,i-1,j)),
        sub(point_direction(b,d,r.g+1,next,0),point_direction(b,d,r.g+1,i-1,0)));
    AlulaPoint span=sub(point(b,s,r.g+1,i,0),point(b,s,r.g,i,j));
    AlulaPoint ds=sub(point_direction(b,d,r.g+1,i,0),point_direction(b,d,r.g,i,j));
    double length=norm(edge), dn=dot(edge,de)/length;
    double derivative=(dcross(edge,span,de,ds)-cross(edge,span)*dn/length)/length;
    int k=b->wake[r.g]+i-b->cuts[r.g].trailing-1;
    return (derivative-(d->wake_width ? d->wake_width[k] : 0))/b->in.length_scale;
}

static int match_direction(const AlulaBody *b,const AlulaBodyState *s,const AlulaBodyState *d,
    const double *coeff,const double *dc,FarfieldPass *pass,double out[3],double *edges)
{
    double match[3]={0},sq[3]={0},dm[3]={0},ds[3]={0};
    int top=b->in.bodies;
    for (int side=0;side<2;++side) {
        int g=side?top:0,j=side?b->tubes[top]:0,tube=side?j-1:0;
        for (int i=0;i<b->in.nx;++i) {
            AlulaPoint a=point(b,s,g,i,j),z=point(b,s,g,i+1,j),edge=sub(z,a);
            AlulaPoint da=point_direction(b,d,g,i,j),dz=point_direction(b,d,g,i+1,j),de=sub(dz,da);
            double length=norm(edge), dn=dot(edge,de)/length;
            const Farfield *f;
            Farfield df;
            if (!field_direction(b,coeff,dc,mean(a,z),mean(da,dz),pass,&f,&df)) return 0;
            int k=index_section(b,g,i,tube);
            double q=b->sections[k].gas.speed,dq=b->section_tangents[k].speed;
            AlulaPoint t=scale(edge,1/length),u=scale(t,q);
            AlulaPoint du=add(scale(t,dq),scale(sub(de,scale(t,dn)),q/length));
            double defect=cross(f->velocity,u), ddefect=dcross(f->velocity,u,df.velocity,du);
            for (int l=0;l<3;++l) {
                double mode=cross(f->basis[l+1],u), dmode=dcross(f->basis[l+1],u,df.basis[l+1],du);
                match[l]+=length*defect*mode;
                sq[l]+=length*mode*mode;
                double a=dn*defect*mode+length*(ddefect*mode+defect*dmode);
                double z=dn*mode*mode+2*length*mode*dmode;
                dm[l]+=a;
                ds[l]+=z;
                if (edges) {
                    size_t k=6*((size_t)side*b->in.nx+i)+2*l;
                    edges[k]=a; edges[k+1]=z;
                }
            }
        }
    }
    for (int k=0;k<3;++k) {
        if (!positive(sq[k])) return 0;
        out[k]=(dm[k]-.5*match[k]*(ds[k]/sq[k]))/sqrt(sq[k]);
        if (!isfinite(out[k])) return 0;
        if (edges) for (int i=0;i<2*b->in.nx;++i) {
            size_t j=6*(size_t)i+2*k;
            edges[j]=(edges[j]-.5*match[k]*(edges[j+1]/sq[k]))/sqrt(sq[k]);
        }
    }
    return 1;
}

AlulaBodyStatus alula_body_apply_prepared(AlulaBody *b,const AlulaBodyState *s,const AlulaBodyState *d,
    int first,int end,double *product,double *edges)
{
    if (!prepare_tangents(b,d,first,end)) return ALULA_BODY_NONFINITE;
    if (!b->jacobian_fields) {
        /* Two exterior midpoint sequences, one field per pressure row,
         * and at most four endpoints per tangency row. */
        size_t count=2*(size_t)b->in.nx;
        for (int k=0;k<b->size.rows;++k) {
            if (b->rows[k].family==ALULA_BODY_FARFIELD_PRESSURE) ++count;
            if (b->rows[k].family==ALULA_BODY_END_TANGENCY) count+=4;
        }
        b->jacobian_fields=malloc(count*sizeof *b->jacobian_fields);
        if (!b->jacobian_fields) return ALULA_BODY_NONFINITE;
    }
    /* Every color visits these points in the same topology-defined order. */
    FarfieldPass pass={b->jacobian_fields,b->jacobian_field_generation!=b->residual_generation};

    double coeff[]={-s->circulation,s->source,s->doublet_x,s->doublet_y,s->circulation*s->circulation};
    double dc[]={-d->circulation,d->source,d->doublet_x,d->doublet_y,2*s->circulation*d->circulation};
    double match[3];
    if (!match_direction(b,s,d,coeff,dc,&pass,match,edges)) return ALULA_BODY_NONFINITE;
    for (int k=0;k<b->size.rows;++k) {
        Row row=b->rows[k];
        int i=row.i,g=row.g,j=row.j;
        double v=0;
        switch (row.family) {
        case ALULA_BODY_STREAMWISE:
            v=b->cell_tangents[index_section(b,g,i-1,j)].streamwise;
            break;
        case ALULA_BODY_INLET_DENSITY: {
            int l=index_section(b,g,0,j);
            AlulaSection a=b->sections[l].gas;
            SectionTangent t=b->section_tangents[l];
            v=t.density/a.density-t.enthalpy/((b->in.gamma-1)*a.enthalpy);
            break;
        }
        case ALULA_BODY_INTERNAL_PRESSURE: {
            int l=index_section(b,g,i-1,j);
            v=(b->cell_tangents[l-1].pressure[1]-b->cell_tangents[l].pressure[0])/b->p_inf;
            break;
        }
        case ALULA_BODY_FARFIELD_PRESSURE: {
            Farfield df;
            if (!field_direction(b,coeff,dc,point(b,s,g,i,j),point_direction(b,d,g,i,j),&pass,NULL,&df))
                return ALULA_BODY_NONFINITE;
            int side=g==b->in.bodies, l=index_section(b,g,i-1,side?j-1:0);
            v=(b->cell_tangents[l].pressure[side]-df.pressure)/b->p_inf;
            break;
        }
        case ALULA_BODY_CUT_PRESSURE:
        case ALULA_BODY_TRAILING_KUTTA:
        case ALULA_BODY_LEADING_KUTTA: {
            int lower=index_section(b,g,i-1,b->tubes[g]-1),upper=index_section(b,g+1,i-1,0);
            v=(b->cell_tangents[upper].pressure[0]-b->cell_tangents[lower].pressure[1])/b->p_inf;
            break;
        }
        case ALULA_BODY_END_TANGENCY:
            if (!boundary_direction(b,s,d,coeff,dc,row,&pass,&v)) return ALULA_BODY_NONFINITE;
            break;
        case ALULA_BODY_WAKE_GAP:
            v=gap_direction(b,s,d,row);
            break;
        case ALULA_BODY_FARFIELD_MATCH:
            v=match[j];
            break;
        case ALULA_BODY_FAMILIES:
            return ALULA_BODY_INVALID;
        }
        if (!isfinite(v)) return ALULA_BODY_NONFINITE;
        product[k]=v;
    }
    b->jacobian_field_generation=b->residual_generation;
    return ALULA_BODY_OK;
}


AlulaBodyStatus alula_body_linearize(AlulaBody *b,const AlulaBodyState *s,const AlulaBodyState *d,
    double *residual,double *product,AlulaBodyDiagnostics *diag)
{
    if (!d || !product || product==residual || !isfinite(d->circulation) || !isfinite(d->source) ||
        !isfinite(d->doublet_x) || !isfinite(d->doublet_y)) return ALULA_BODY_INVALID;
    AlulaBodyStatus status=alula_body_residual(b,s,residual,diag);
    return status==ALULA_BODY_OK ? alula_body_apply_prepared(b,s,d,0,b->in.nx,product,NULL) : status;
}
