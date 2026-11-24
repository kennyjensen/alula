// SPDX-License-Identifier: GPL-2.0-or-later
#define _POSIX_C_SOURCE 200809L
#include "channel.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static int read_vector(double *v,int n)
{
    for (int i=0; i<n; ++i) if (scanf("%lf",v+i)!=1) return 0;
    return 1;
}

static void emit(const double *v,int n)
{
    for (int i=0; i<n; ++i) printf("%.17g%c",v[i],i==n-1 ? '\n' : ' ');
}

static double milliseconds(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC,&t);
    return 1000.*t.tv_sec+1e-6*t.tv_nsec;
}

int main(void)
{
    int nx,nt,mode,iterations,diagnostics;
    AlulaChannelInput in={0};
    double tolerance;
    if (scanf("%d%d%d%d%d%lf",&nx,&nt,&mode,&iterations,&diagnostics,&tolerance)!=6 ||
        nx<2 || nt<1 || (size_t)nx*nt>100000 ||
        scanf("%lf%lf%lf%lf%lf%lf%lf",&in.gamma,&in.rho_ref,&in.p_ref,&in.correction,
            &in.mucon,&in.mcrit,&in.epsilon)!=7) return 2;
    in.nx=nx; in.nt=nt; in.mode=(AlulaChannelMode)mode;
    int n=nx*nt+(nx+1)*(nt-1)+(mode>=ALULA_CHANNEL_UPWIND ? nt : 0);
    double *fixed=calloc(3*(size_t)(nx+1)+6*(size_t)nt,sizeof *fixed);
    double *vectors=calloc(3*(size_t)n,sizeof *vectors);
    AlulaChannel *c=NULL;
    int status=2;
    if (!fixed || !vectors) goto done;
    double *cursor=fixed;
    const int lengths[]={nx+1,nx+1,nx+1,nt,nt,nt,nt-1,nt-1,nt};
    const double **inputs[]={&in.x,&in.lower,&in.upper,&in.mass,&in.h0,&in.rho_total,
        &in.inlet_slope,&in.outlet_slope,&in.outlet_pressure};
    for (size_t i=0; i<sizeof lengths/sizeof *lengths; ++i) {
        if (!read_vector(cursor,lengths[i])) goto done;
        *inputs[i]=cursor;
        cursor+=lengths[i];
    }
    if (mode>=ALULA_CHANNEL_UPWIND) in.outlet_slope=NULL;
    if (!read_vector(vectors,n)) goto done;
    c=alula_channel_create(&in);
    if (!c || alula_channel_size(c)!=n) goto done;
    if (diagnostics) {
        if (!read_vector(vectors+n,n) || !alula_channel_residual(c,vectors,vectors+2*n)) goto done;
        emit(vectors+2*n,n);
        if (!alula_channel_jacobian_product(c,vectors,vectors+n,vectors+2*n)) goto done;
        emit(vectors+2*n,n);
    }
    double begin=milliseconds();
    AlulaChannelResult r=alula_channel_solve(c,vectors,iterations,tolerance);
    double elapsed=milliseconds()-begin;
    printf("%d %d %d %.17g %.17g %d %d %d %d %.17g %.17g\n",
        (int)r.status,r.iterations,r.backtracks,r.residual,r.max_mach,
        r.linear.nonzeros,r.linear.factor_nonzeros,r.linear.factor_count,r.linear.refactor_count,
        r.linear.max_relative_residual,elapsed);
    emit(vectors,n);
    status=0;
done:
    alula_channel_free(c);
    free(fixed);
    free(vectors);
    return status;
}
