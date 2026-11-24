// SPDX-License-Identifier: GPL-2.0-or-later
#include "streamtube.h"
#include <stdio.h>

static int pair(AlulaPoint *p)
{
    return scanf("%lf%lf",&p->x,&p->y) == 2;
}

static void emit(const AlulaCellValue *v)
{
    const AlulaCellGeometry *g=&v->geometry;
    printf(" %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g",
        g->area,g->curvature,g->normal_area[0],g->normal_area[1],g->length[0],g->length[1],
        v->pressure[0],v->pressure[1],v->momentum,v->entropy,v->isentropic,v->correction);
    for (int i=0; i<2; ++i) {
        AlulaSection s=v->state[i];
        printf(" %.17g %.17g %.17g %.17g %.17g",s.density,s.speed,s.pressure,s.enthalpy,s.mach_sq);
    }
}

static int upwind_input(AlulaUpwindInput *in)
{
    return scanf("%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
        &in->speed[0],&in->speed[1],&in->speed[2],&in->mach_sq[0],&in->mach_sq[1],
        &in->spacing[0],&in->spacing[1],&in->mucon,&in->mcrit,&in->gamma)==10;
}

static void emit_blend(const AlulaBlendValue *v)
{
    printf(" %.17g %.17g %.17g %.17g %.17g",v->fraction,v->residual,v->loss[0],v->loss[1],v->compression);
}

int main(void)
{
    AlulaCellInput in;
    AlulaCellDirection d;
    char kind;
    int fields;
    while ((fields=scanf(" %c",&kind)) == 1) {
        if (kind == 'U') {
            AlulaUpwindInput in, d;
            AlulaUpwindValue v, dv;
            if (!upwind_input(&in) || !upwind_input(&d)) return 2;
            AlulaCellStatus status=alula_streamtube_upwind(&in,&v,&d,&dv);
            printf("%d",(int)status);
            if (status == ALULA_CELL_OK) printf(" %.17g %.17g %.17g %.17g %.17g %.17g",
                v.speed,v.coefficient,v.correction,dv.speed,dv.coefficient,dv.correction);
            putchar('\n');
            continue;
        }
        double epsilon=0, depsilon=0;
        if (kind == 'B') {
            if (scanf("%lf%lf",&epsilon,&depsilon)!=2) return 2;
        } else if (kind != 'C') return 2;
        if (scanf("%d%d%lf%lf%lf%lf%lf%lf%lf%lf",
        &in.positive_simple,&in.biased,&in.mass,&in.h0,&in.gamma,&in.correction,
        &in.density[0],&in.density[1],&in.transport[0],&in.transport[1]) != 10) return 2;
        for (int i=0; i<3; ++i)
            if (!pair(&in.lower[i]) || !pair(&in.upper[i])) return 2;
        if (scanf("%lf%lf%lf%lf%lf%lf%lf",&d.mass,&d.h0,&d.correction,
            &d.density[0],&d.density[1],&d.transport[0],&d.transport[1]) != 7) return 2;
        for (int i=0; i<3; ++i)
            if (!pair(&d.lower[i]) || !pair(&d.upper[i])) return 2;
        AlulaCell cell;
        AlulaCellStatus status=alula_cell_prepare(&cell,&in);
        printf("%d",(int)status);
        if (status == ALULA_CELL_OK) {
            emit(&cell.value);
            AlulaCellValue dv;
            status=alula_cell_apply(&cell,&d,&dv);
            printf(" %d",(int)status);
            if (status == ALULA_CELL_OK) {
                emit(&dv);
                if (kind == 'B') {
                    AlulaBlendValue v, db;
                    status=alula_streamtube_blend(&cell.value,epsilon,&v,&dv,depsilon,&db);
                    printf(" %d",(int)status);
                    if (status==ALULA_CELL_OK) { emit_blend(&v); emit_blend(&db); }
                }
            }
        }
        putchar('\n');
    }
    return fields == EOF ? 0 : 2;
}
