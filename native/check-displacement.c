// SPDX-License-Identifier: GPL-2.0-or-later
#define _POSIX_C_SOURCE 200809L
#include "cli_geometry.h"
#include "panel_bl.h"
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    AlulaContour contours[6];
    AlulaPoint points[ALULA_PANEL_MAX_TOTAL+6];
    size_t bodies,used,nw;
    double alpha,extent;
    int mode,full;
    if (!alula_cli_contours(1,contours,points,&bodies,&used) ||
        scanf("%lf%zu%lf%d%d",&alpha,&nw,&extent,&mode,&full)!=5 || nw>1000) return 1;
    AlulaPoint *path=calloc(bodies*nw,sizeof *path);
    if (!path) return 1;
    if (mode==2) for (size_t i=0;i<bodies*nw;++i)
        if (scanf("%lf%lf",&path[i].x,&path[i].y)!=2) return 1;
    struct timespec start,finish;
    int timing=getenv("ALULA_PANEL_BENCH") && clock_gettime(CLOCK_MONOTONIC,&start)==0;
    AlulaPanelDisplacement *o=alula_panel_displacement_create(contours,bodies,alpha,nw,extent,
        mode==1?ALULA_WAKE_INVISCID:ALULA_WAKE_STRAIGHT,mode==2?path:NULL);
    if (timing && clock_gettime(CLOCK_MONOTONIC,&finish)==0)
        fprintf(stderr,"native assembly %.2f ms\n",1000.0*(finish.tv_sec-start.tv_sec)+1e-6*(finish.tv_nsec-start.tv_nsec));
    free(path);
    if (!o) return 2;
    const AlulaDisplacementData *d=alula_panel_displacement_data(o);
    size_t n=d->stations;
    double *mass=malloc(n*sizeof *mass),*speed=malloc(n*sizeof *speed);
    AlulaPoint *velocity=malloc(bodies*(nw-1)*sizeof *velocity);
    if (!mass || !speed || !velocity) return 1;
    for (size_t i=0;i<n;++i) if (scanf("%lf",mass+i)!=1) return 1;
    if (!alula_panel_displacement_evaluate(o,mass,n,speed) ||
        !alula_panel_displacement_wake_velocity(o,mass,n,velocity)) return 1;
    for (size_t i=0;i<n;++i)
        printf("%.17g %.17g %.17g %.17g %.17g\n",d->points[i].x,d->points[i].y,d->arc[i],d->inviscid[i],speed[i]);
    for (size_t i=0;i<bodies*(nw-1);++i) printf("%.17g %.17g\n",velocity[i].x,velocity[i].y);
    if (full) for (size_t i=0;i<n;++i) {
        for (size_t j=0;j<n;++j) printf(" %.17g",d->influence[i*n+j]);
        putchar('\n');
    }
    double re,ncrit,trips[6][2];
    if (scanf("%lf%lf",&re,&ncrit)==2) {
        for (size_t e=0;e<bodies;++e) if (scanf("%lf%lf",&trips[e][0],&trips[e][1])!=2) return 1;
        AlulaPanelBL *p=alula_panel_bl_initialize(o,re,ncrit,(const double (*)[2])trips);
        if (!p || !p->complete) {
            fprintf(stderr,"BL seed failed with status %d at body %zu, station %zu\n",
                p?(int)p->status:-1,p?p->failed_body:0,p?p->failed_station:0);
            return 3;
        }
        int iterations;
        if (scanf("%d",&iterations)==1) {
            AlulaPanelBLNewton r=alula_panel_bl_solve(o,p,iterations,1e-9,NULL,NULL);
            fprintf(stderr,"Newton status %d, %d iterations, %d backtracks, residual %.12g\n",
                r.status,r.iterations,r.backtracks,r.residual);
            if (!r.converged) return 4;
        }
        printf("%.17g %.17g\n",p->interval_residual,p->outer_residual);
        for (size_t e=0;e<bodies;++e) printf("%zu %zu %.17g %.17g %.17g\n",
            p->upper_first[e],p->lower_first[e],p->stagnation[e],p->trip_arc[e][0],p->trip_arc[e][1]);
        for (size_t i=0;i<n;++i) {
            AlulaBLStation a=p->stations[i];
            printf("%d %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                p->regimes[i],p->speed[i],a.s,a.ue,a.theta,a.delta_star,a.aux,a.wake_gap);
        }
        alula_panel_bl_free(p);
    }
    mass[n/2]=NAN;
    if (alula_panel_displacement_evaluate(o,mass,n,speed) ||
        alula_panel_displacement_evaluate(o,mass,n-1,speed)) return 1;
    free(mass);free(speed);free(velocity);alula_panel_displacement_free(o);
    return 0;
}
