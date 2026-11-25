// SPDX-License-Identifier: GPL-2.0-or-later
#include "streamtube_internal.h"
#include <math.h>

static AlulaPoint add(AlulaPoint a, AlulaPoint b) { return (AlulaPoint){a.x+b.x, a.y+b.y}; }
static AlulaPoint sub(AlulaPoint a, AlulaPoint b) { return (AlulaPoint){a.x-b.x, a.y-b.y}; }
static AlulaPoint scale(AlulaPoint a, double k) { return (AlulaPoint){a.x*k, a.y*k}; }
static AlulaPoint mean(AlulaPoint a, AlulaPoint b) { return scale(add(a, b), .5); }
static double dot(AlulaPoint a, AlulaPoint b) { return a.x*b.x+a.y*b.y; }
static double cross(AlulaPoint a, AlulaPoint b) { return a.x*b.y-a.y*b.x; }
static double dcross(AlulaPoint a, AlulaPoint b, AlulaPoint da, AlulaPoint db)
{
    return cross(da, b)+cross(a, db);
}

static double ddot(AlulaPoint a, AlulaPoint b, AlulaPoint da, AlulaPoint db)
{
    return dot(da, b)+dot(a, db);
}

static int finite_point(AlulaPoint a) { return isfinite(a.x) && isfinite(a.y); }
static int positive(double x) { return isfinite(x) && x > 0; }
static int nonnegative(double x) { return isfinite(x) && x >= 0; }

static int between(AlulaPoint a, AlulaPoint b, AlulaPoint c)
{
    return c.x >= fmin(a.x,b.x) && c.x <= fmax(a.x,b.x)
        && c.y >= fmin(a.y,b.y) && c.y <= fmax(a.y,b.y);
}

static int opposite(double a, double b) { return (a < 0 && b > 0) || (a > 0 && b < 0); }
static double orient(AlulaPoint a, AlulaPoint b, AlulaPoint c) { return cross(sub(b,a),sub(c,a)); }

static int intersects(AlulaPoint a, AlulaPoint b, AlulaPoint c, AlulaPoint d)
{
    double u=orient(a,b,c), v=orient(a,b,d), w=orient(c,d,a), z=orient(c,d,b);
    return (opposite(u,v) && opposite(w,z)) || (u == 0 && between(a,b,c))
        || (v == 0 && between(a,b,d)) || (w == 0 && between(c,d,a))
        || (z == 0 && between(c,d,b));
}

static int positive_simple(const AlulaPoint *vertices, int n)
{
    AlulaPoint p[6];
    for (int i=0; i<n; ++i) {
        p[i]=sub(vertices[i],vertices[0]);
        if (!finite_point(p[i])) return 0;
    }
    double area=0;
    for (int i=0; i<n; ++i) {
        AlulaPoint a=p[i], b=p[(i+1)%n], c=p[(i+2)%n];
        AlulaPoint ab=sub(b,a), bc=sub(c,b);
        if ((ab.x == 0 && ab.y == 0) || (orient(a,b,c) == 0 && dot(ab,bc) <= 0)) return 0;
        if (i>0 && i<n-1) area+=cross(p[i],p[i+1]);
        for (int j=i+2; j<n; ++j)
            if (!(i==0 && j==n-1) && intersects(a,b,p[j],p[(j+1)%n])) return 0;
    }
    return positive(area);
}

static int valid_geometry(const AlulaCellInput *in)
{
    const AlulaPoint *l=in->lower, *u=in->upper;
    AlulaPoint ml[2], mu[2];
    for (int i=0; i<3; ++i) if (!finite_point(l[i]) || !finite_point(u[i])) return 0;
    for (int i=0; i<2; ++i) {
        AlulaPoint p[]={l[i],l[i+1],u[i+1],u[i]};
        if (in->positive_simple) {
            if (!positive_simple(p,4)) return 0;
            ml[i]=mean(l[i],l[i+1]);
            mu[i]=mean(u[i],u[i+1]);
            AlulaPoint left[]={l[i],ml[i],mu[i],u[i]};
            AlulaPoint right[]={ml[i],l[i+1],u[i+1],mu[i]};
            if (!positive_simple(left,4) || !positive_simple(right,4)) return 0;
        } else {
            for (int j=0; j<4; ++j)
                if (!(cross(sub(p[(j+1)%4],p[j]),sub(p[(j+2)%4],p[(j+1)%4])) > 0)) return 0;
        }
    }
    if (!in->positive_simple) return 1;
    AlulaPoint middle[]={ml[0],l[1],ml[1],mu[1],u[1],mu[0]};
    return positive_simple(middle,6);
}

static int geometry(AlulaCell *cell)
{
    const AlulaCellInput *in=&cell->input;
    const AlulaPoint *l=in->lower, *u=in->upper;
    AlulaCellGeometry *g=&cell->value.geometry;
    for (int i=0; i<3; ++i) cell->gap[i]=sub(u[i],l[i]);
    for (int i=0; i<2; ++i) {
        cell->segment[i]=mean(sub(l[i+1],l[i]),sub(u[i+1],u[i]));
        g->length[i]=hypot(cell->segment[i].x,cell->segment[i].y);
        if (!positive(g->length[i])) return 0;
        g->direction[i]=scale(cell->segment[i],1/g->length[i]);
        g->section[i]=mean(cell->gap[i],cell->gap[i+1]);
        g->normal_area[i]=cross(g->direction[i],g->section[i]);
        if (!positive(g->normal_area[i])) return 0;
    }
    g->side[0]=scale(sub(l[2],l[0]),.5);
    g->side[1]=scale(sub(u[2],u[0]),.5);
    g->streamwise=mean(g->side[0],g->side[1]);
    g->transverse=mean(g->section[0],g->section[1]);
    g->area=cross(g->streamwise,g->transverse);
    if (!positive(g->area)) return 0;
    double bend=cross(cell->segment[0],sub(cell->gap[1],cell->gap[2]))
        +cross(sub(cell->gap[0],cell->gap[1]),cell->segment[1]);
    g->curvature=bend/(2*g->area);
    return isfinite(g->curvature);
}

static int finite_section(AlulaSection s)
{
    return isfinite(s.density) && isfinite(s.speed) && isfinite(s.pressure)
        && isfinite(s.enthalpy) && isfinite(s.mach_sq);
}

int alula_cell_geometry_only(AlulaCell *cell,const AlulaCellInput *in)
{
    if (!valid_geometry(in)) return 0;
    cell->input=*in;
    return geometry(cell);
}

AlulaCellStatus alula_cell_initialize(AlulaCell *cell, const AlulaCellInput *in)
{
    if (!cell || !in || !positive(in->mass) || !positive(in->h0) || !positive(in->density[0])
        || !positive(in->density[1]) || !isfinite(in->gamma) || in->gamma <= 1
        || !nonnegative(in->correction)) return ALULA_CELL_INPUT;
    if (!valid_geometry(in)) return ALULA_CELL_GEOMETRY;
    cell->input=*in;
    AlulaCellValue *v=&cell->value;
    *v=(AlulaCellValue){0};
    if (!geometry(cell)) return ALULA_CELL_GEOMETRY;
    const AlulaCellGeometry *g=&v->geometry;
    double gm=in->gamma-1, k=gm/in->gamma;
    for (int i=0; i<2; ++i) {
        double rho=in->density[i], q=in->mass/(rho*g->normal_area[i]);
        double h=in->h0-.5*q*q;
        if (!positive(h)) return ALULA_CELL_ENTHALPY;
        v->state[i]=(AlulaSection){rho,q,k*rho*h,h,q*q/(gm*h)};
        if (!positive(v->state[i].pressure) || !finite_section(v->state[i])) return ALULA_CELL_NONFINITE;
    }
    return ALULA_CELL_OK;
}

AlulaCellStatus alula_cell_flux(AlulaCell *cell)
{
    AlulaCellInput *in=&cell->input;
    AlulaCellValue *v=&cell->value;
    const AlulaCellGeometry *g=&v->geometry;
    double gm=in->gamma-1, ratio[2]={0,0};
    for (int i=0; i<2; ++i) {
        double q=v->state[i].speed;
        v->transport[i]=in->biased ? in->transport[i] : q;
        double qt=v->transport[i];
        if (!nonnegative(qt)) return ALULA_CELL_INPUT;
        v->artificial_h[i]=in->h0-.5*qt*qt;
        if (!positive(v->artificial_h[i])) return ALULA_CELL_TRANSPORT;
        if (in->biased) ratio[i]=log1p(.5*(q-qt)*(q+qt)/v->state[i].enthalpy);
    }
    AlulaSection a=v->state[0], b=v->state[1];
    double pmean=.5*(a.pressure+b.pressure), mmean=.5*(a.mach_sq+b.mach_sq);
    double qa=v->transport[0], qb=v->transport[1];
    double pc=mmean < 1 ? in->correction*in->gamma*pmean*mmean*(1-mmean)*g->curvature : 0;
    v->correction=pc;
    double normal=in->mass*(qa*dot(g->direction[0],g->transverse)-qb*dot(g->direction[1],g->transverse))/g->area;
    double diff=normal+pc*cross(g->section[0],g->section[1])/g->area;
    double sum=a.pressure+b.pressure+2*pc;
    v->pressure[0]=.5*(sum-diff);
    v->pressure[1]=.5*(sum+diff);
    if (!positive(v->pressure[0]) || !positive(v->pressure[1])) return ALULA_CELL_PRESSURE;
    v->momentum=in->mass*(qb*dot(g->direction[1],g->streamwise)-qa*dot(g->direction[0],g->streamwise))/g->area
        +b.pressure-a.pressure-pc*cross(g->side[0],g->side[1])/g->area;
    v->entropy=log1p((b.enthalpy-a.enthalpy)/a.enthalpy)/gm-log1p((b.density-a.density)/a.density);
    v->artificial_entropy=v->entropy+in->gamma/gm*(ratio[1]-ratio[0]);
    v->isentropic=-pmean*v->artificial_entropy;
    return isfinite(v->momentum) && isfinite(v->entropy) && isfinite(v->isentropic) ? ALULA_CELL_OK : ALULA_CELL_NONFINITE;
}

AlulaCellStatus alula_cell_prepare(AlulaCell *cell, const AlulaCellInput *in)
{
    AlulaCellStatus status=alula_cell_initialize(cell,in);
    return status == ALULA_CELL_OK ? alula_cell_flux(cell) : status;
}

static void geometry_direction(const AlulaCell *cell, const AlulaCellDirection *d, AlulaCellGeometry *dg)
{
    const AlulaCellGeometry *g=&cell->value.geometry;
    AlulaPoint gap[3], segment[2];
    for (int i=0; i<3; ++i) gap[i]=sub(d->upper[i],d->lower[i]);
    for (int i=0; i<2; ++i) {
        segment[i]=mean(sub(d->lower[i+1],d->lower[i]),sub(d->upper[i+1],d->upper[i]));
        dg->length[i]=dot(g->direction[i],segment[i]);
        dg->direction[i]=scale(sub(segment[i],scale(g->direction[i],dg->length[i])),1/g->length[i]);
        dg->section[i]=mean(gap[i],gap[i+1]);
        dg->normal_area[i]=dcross(g->direction[i],g->section[i],dg->direction[i],dg->section[i]);
    }
    dg->side[0]=scale(sub(d->lower[2],d->lower[0]),.5);
    dg->side[1]=scale(sub(d->upper[2],d->upper[0]),.5);
    dg->streamwise=mean(dg->side[0],dg->side[1]);
    dg->transverse=mean(dg->section[0],dg->section[1]);
    dg->area=dcross(g->streamwise,g->transverse,dg->streamwise,dg->transverse);
    double bend=dcross(cell->segment[0],sub(cell->gap[1],cell->gap[2]),segment[0],sub(gap[1],gap[2]))
        +dcross(sub(cell->gap[0],cell->gap[1]),cell->segment[1],sub(gap[0],gap[1]),segment[1]);
    dg->curvature=(.5*bend-g->curvature*dg->area)/g->area;
}

static double quotient_direction(double numerator, double derivative, const AlulaCellGeometry *g,
    const AlulaCellGeometry *dg)
{
    return (derivative-numerator*dg->area/g->area)/g->area;
}

AlulaCellStatus alula_cell_apply(const AlulaCell *cell, const AlulaCellDirection *d, AlulaCellValue *dv)
{
    if (!cell || !d || !dv || dv == &cell->value || !isfinite(d->mass) || !isfinite(d->h0)
        || !isfinite(d->correction)) return ALULA_CELL_INPUT;
    const AlulaCellInput *in=&cell->input;
    for (int i=0; i<3; ++i)
        if (!finite_point(d->lower[i]) || !finite_point(d->upper[i])) return ALULA_CELL_INPUT;
    for (int i=0; i<2; ++i)
        if (!isfinite(d->density[i]) || !isfinite(d->transport[i]) || (!in->biased && d->transport[i] != 0))
            return ALULA_CELL_INPUT;

    const AlulaCellValue *v=&cell->value;
    const AlulaCellGeometry *g=&v->geometry;
    AlulaCellGeometry *dg=&dv->geometry;
    geometry_direction(cell,d,dg);
    double gm=in->gamma-1, k=gm/in->gamma;
    for (int i=0; i<2; ++i) {
        AlulaSection s=v->state[i];
        double area=g->normal_area[i], da=dg->normal_area[i], rho=d->density[i];
        double q=(-s.speed/s.density)*rho+(s.speed/in->mass)*d->mass+(-s.speed/area)*da;
        double p=k*(in->h0+.5*s.speed*s.speed)*rho-k*s.speed/area*d->mass
            +k*s.density*s.speed*s.speed/area*da+k*s.density*d->h0;
        double h=d->h0-s.speed*q;
        dv->state[i]=(AlulaSection){rho,q,p,h,(2*s.speed*q-s.mach_sq*gm*h)/(gm*s.enthalpy)};
        if (!finite_section(dv->state[i])) return ALULA_CELL_NONFINITE;
        dv->transport[i]=in->biased ? d->transport[i] : q;
        dv->artificial_h[i]=d->h0-v->transport[i]*dv->transport[i];
    }

    AlulaSection a=v->state[0], b=v->state[1], da=dv->state[0], db=dv->state[1];
    double pmean=.5*(a.pressure+b.pressure), dpmean=.5*(da.pressure+db.pressure);
    double mmean=.5*(a.mach_sq+b.mach_sq), dmmean=.5*(da.mach_sq+db.mach_sq);
    double pc=v->correction, dpc=0;
    if (mmean < 1) {
        double mach=mmean*(1-mmean);
        dpc=in->gamma*(d->correction*pmean*mach*g->curvature
            +in->correction*((dpmean*mach+pmean*(1-2*mmean)*dmmean)*g->curvature
            +pmean*mach*dg->curvature));
    }
    dv->correction=dpc;
    double norm_term[2], dnorm_term[2], stream_term[2], dstream_term[2];
    for (int i=0; i<2; ++i) {
        AlulaPoint t=g->direction[i], dt=dg->direction[i];
        double q=v->transport[i], dq=dv->transport[i];
        norm_term[i]=q*dot(t,g->transverse);
        dnorm_term[i]=dq*dot(t,g->transverse)+q*ddot(t,g->transverse,dt,dg->transverse);
        stream_term[i]=q*dot(t,g->streamwise);
        dstream_term[i]=dq*dot(t,g->streamwise)+q*ddot(t,g->streamwise,dt,dg->streamwise);
    }
    double c=cross(g->section[0],g->section[1]);
    double dc=dcross(g->section[0],g->section[1],dg->section[0],dg->section[1]);
    double normal=norm_term[0]-norm_term[1];
    double dn=quotient_direction(in->mass*normal,
        d->mass*normal+in->mass*(dnorm_term[0]-dnorm_term[1]),g,dg);
    double diff=dn+quotient_direction(pc*c,dpc*c+pc*dc,g,dg);
    double sum=da.pressure+db.pressure+2*dpc;
    dv->pressure[0]=.5*(sum-diff);
    dv->pressure[1]=.5*(sum+diff);

    c=cross(g->side[0],g->side[1]);
    dc=dcross(g->side[0],g->side[1],dg->side[0],dg->side[1]);
    double stream=stream_term[1]-stream_term[0];
    dv->momentum=quotient_direction(in->mass*stream,
        d->mass*stream+in->mass*(dstream_term[1]-dstream_term[0]),g,dg)
        +db.pressure-da.pressure-quotient_direction(pc*c,dpc*c+pc*dc,g,dg);
    dv->entropy=(db.enthalpy/b.enthalpy-da.enthalpy/a.enthalpy)/gm-db.density/b.density+da.density/a.density;
    dv->artificial_entropy=in->biased ? in->gamma/gm*(dv->artificial_h[1]/v->artificial_h[1]
        -dv->artificial_h[0]/v->artificial_h[0])-db.pressure/b.pressure+da.pressure/a.pressure : dv->entropy;
    dv->isentropic=-dpmean*v->artificial_entropy-pmean*dv->artificial_entropy;
    return isfinite(dv->momentum) && isfinite(dv->entropy) && isfinite(dv->isentropic)
        && isfinite(dpc) && isfinite(dv->pressure[0]) && isfinite(dv->pressure[1])
        ? ALULA_CELL_OK : ALULA_CELL_NONFINITE;
}

static int finite_upwind(AlulaUpwindValue v)
{
    return isfinite(v.speed) && isfinite(v.coefficient) && isfinite(v.correction);
}

AlulaCellStatus alula_streamtube_upwind(const AlulaUpwindInput *in, AlulaUpwindValue *value,
    const AlulaUpwindInput *d, AlulaUpwindValue *dv)
{
    if (!in || !value || (d && !dv) || !isfinite(in->mucon) || !nonnegative(in->mcrit)
        || in->mcrit > 1 || !isfinite(in->gamma) || in->gamma <= 1) return ALULA_CELL_INPUT;
    for (int i=0; i<3; ++i)
        if (!nonnegative(in->speed[i]) || (d && !isfinite(d->speed[i]))) return ALULA_CELL_INPUT;
    for (int i=0; i<2; ++i)
        if (!nonnegative(in->mach_sq[i]) || !positive(in->spacing[i]) ||
            (d && (!isfinite(d->mach_sq[i]) || !isfinite(d->spacing[i])))) return ALULA_CELL_INPUT;

    double x=.5*in->mach_sq[0]+.5*in->mach_sq[1], e=1-in->mcrit;
    double base=0, dx=0, de=0;
    if (x > 0) {
        double a=1-1/x;
        if (e == 0) {
            base=fmax(0,a);
            if (x > 1) dx=1/x/x;
        } else {
            double z=a/e, az=fabs(z), tail=exp(-az);
            double sigmoid=(z >= 0 ? 1 : tail)/(1+tail);
            double log_tail=log1p(tail);
            base=fmax(a,0)+e*log_tail;
            dx=sigmoid == 0 ? 0 : sigmoid/x/x;
            de=log_tail+(isfinite(az) ? az*tail/(1+tail) : 0);
        }
    }
    int second=in->mucon >= 0;
    double ratio=in->spacing[1]/in->spacing[0], scale=fabs(in->mucon)/in->gamma;
    double delta=in->speed[1]-in->speed[0];
    double difference=-(in->speed[2]-in->speed[1])+(second ? ratio*delta : 0);
    value->coefficient=scale*base;
    value->correction=value->coefficient*difference;
    value->speed=in->speed[2]+value->correction;
    if (!isfinite(ratio) || !finite_upwind(*value)) return ALULA_CELL_NONFINITE;
    if (!d) return ALULA_CELL_OK;

    if (!isfinite(d->mucon) || !isfinite(d->mcrit) || !isfinite(d->gamma)
        || (e == 0 && (x == 1 || d->mcrit != 0)) || (in->mucon == 0 && d->mucon != 0))
        return ALULA_CELL_INPUT;
    double dscale=(in->mucon < 0 ? -d->mucon : d->mucon)/in->gamma
        -scale*(d->gamma/in->gamma);
    double dr=(d->spacing[1]-ratio*d->spacing[0])/in->spacing[0];
    double dd=-(d->speed[2]-d->speed[1])+(second ? ratio*(d->speed[1]-d->speed[0])+delta*dr : 0);
    dv->coefficient=dscale*base+scale*(dx*(.5*d->mach_sq[0]+.5*d->mach_sq[1])-de*d->mcrit);
    dv->correction=dv->coefficient*difference+value->coefficient*dd;
    dv->speed=d->speed[2]+dv->correction;
    return finite_upwind(*dv) ? ALULA_CELL_OK : ALULA_CELL_NONFINITE;
}

typedef struct { double value, derivative; } Step;

static Step smooth_step(double t)
{
    if (t <= 0) return (Step){0,0};
    if (t >= 1) return (Step){1,0};
    double x=t <= .5 ? t : 1-t, s=x*x*x*(10+x*(-15+6*x));
    return (Step){t <= .5 ? s : 1-s,30*t*t*(1-t)*(1-t)};
}

AlulaCellStatus alula_streamtube_blend(const AlulaCellValue *v, double epsilon,
    AlulaBlendValue *value, const AlulaCellValue *d, double de, AlulaBlendValue *dv)
{
    if (!v || !value || (d && !dv) || !positive(epsilon) || !isfinite(v->momentum)
        || !isfinite(v->isentropic) || (d && !isfinite(de))) return ALULA_CELL_INPUT;
    double phi[2], bias[2], dphi[2];
    for (int i=0; i<2; ++i) {
        AlulaSection s=v->state[i];
        if (!positive(s.density) || !positive(s.speed) || !positive(s.pressure)
            || !nonnegative(v->transport[i])) return ALULA_CELL_INPUT;
        phi[i]=s.density*(s.speed/s.pressure);
        bias[i]=v->transport[i]-s.speed;
        if (d) {
            AlulaSection ds=d->state[i];
            dphi[i]=phi[i]*(ds.density/s.density+ds.speed/s.speed-ds.pressure/s.pressure);
        }
    }
    double delta=phi[1]-phi[0], qa=v->state[0].speed, qb=v->state[1].speed;
    double relative=(qa-qb)/qb, scale=sqrt(epsilon);
    value->compression=fabs(relative)<.5 ? log1p(relative) : log(qa)-log(qb);
    Step gate=smooth_step(value->compression/scale), loss[2];
    for (int i=0; i<2; ++i) {
        value->loss[i]=-bias[i]*delta;
        if (!isfinite(value->loss[i])) return ALULA_CELL_NONFINITE;
        loss[i]=smooth_step(value->loss[i]/epsilon);
    }
    double either=loss[0].value+(1-loss[0].value)*loss[1].value;
    double f=gate.value*either;
    value->fraction=f;
    value->residual=f==0 ? v->isentropic : f==1 ? v->momentum : f*v->momentum+(1-f)*v->isentropic;
    if (!isfinite(value->compression) || !isfinite(value->residual)) return ALULA_CELL_NONFINITE;
    if (!d) return ALULA_CELL_OK;

    double ds[2], ddelta=dphi[1]-dphi[0];
    for (int i=0; i<2; ++i) {
        dv->loss[i]=-(d->transport[i]-d->state[i].speed)*delta-bias[i]*ddelta;
        ds[i]=loss[i].derivative==0 ? 0 : loss[i].derivative*
            (dv->loss[i]-value->loss[i]/epsilon*de)/epsilon;
    }
    dv->compression=d->state[0].speed/qa-d->state[1].speed/qb;
    double dg=gate.derivative==0 ? 0 : gate.derivative*
        (dv->compression-.5*value->compression*(de/epsilon))/scale;
    double df=f==0 || f==1 ? 0 : dg*either+gate.value*
        ((1-loss[1].value)*ds[0]+(1-loss[0].value)*ds[1]);
    dv->fraction=df;
    dv->residual=f==0 ? d->isentropic : f==1 ? d->momentum :
        f*d->momentum+(1-f)*d->isentropic+df*(v->momentum-v->isentropic);
    return isfinite(df) && isfinite(dv->residual) && isfinite(dv->compression)
        && isfinite(dv->loss[0]) && isfinite(dv->loss[1]) ? ALULA_CELL_OK : ALULA_CELL_NONFINITE;
}
