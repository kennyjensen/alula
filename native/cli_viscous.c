// SPDX-License-Identifier: GPL-2.0-or-later
#define _POSIX_C_SOURCE 200809L
#include "cli_geometry.h"
#include "panel_viscous.h"
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void number(FILE *out,double value)
{
    if (isfinite(value)) fprintf(out,"%.17g",value);
    else fputs("null",out);
}

static void field(FILE *out,const char *name,double value)
{
    fprintf(out,",\"%s\":",name);number(out,value);
}

static void plots(FILE *out,const AlulaPanelSnapshot *s)
{
    const AlulaDisplacementData *d=alula_panel_displacement_data(s->outer);
    const AlulaPanelBL *p=s->profile;
    if (!d || !p || !p->complete) return;
    fputs(",\"cp\":[",out);
    for (size_t e=0;e<d->bodies;++e) {
        size_t start=d->surface_start[e],end=d->surface_start[e+1]-1,nose=start;
        AlulaPoint a=d->points[start],b=d->points[end],te={.5*(a.x+b.x),.5*(a.y+b.y)};
        double max_dist=0;
        for (size_t i=start;i<=end;++i) {
            double dist=hypot(d->points[i].x-te.x,d->points[i].y-te.y);
            if (dist>max_dist) {max_dist=dist;nose=i;}
        }
        for (size_t i=start;i<=end;++i)
            fprintf(out,"%s{\"contour\":%zu,\"surface\":\"%s\",\"x\":%.12g,\"cp\":%.12g}",
                i?",":"",e,i<=nose?"upper":"lower",d->points[i].x,1-p->speed[i]*p->speed[i]);
    }
    fputs("],\"boundaryLayerCurves\":[",out);
    for (size_t e=0;e<d->bodies;++e) for (int side=0;side<2;++side) {
        fprintf(out,"%s{\"contour\":%zu,\"surface\":\"%s\",\"samples\":[",e||side?",":"",e,side?"lower":"upper");
        size_t first=side?p->lower_first[e]:p->upper_first[e],last=side?d->surface_start[e+1]-1:d->surface_start[e];
        for (size_t i=first;;side?++i:--i) {
            fprintf(out,"%s[%.12g,%.12g]",i==first?"":",",p->stations[i].s,p->stations[i].delta_star);
            if (i==last) break;
        }
        fputs("]}",out);
    }
    fputs("],\"streamlines\":[",out);
    for (size_t e=0;e<d->bodies;++e) {
        fputs(e?",[":"[",out);
        for (size_t j=0;j<d->wake_count;++j) {
            AlulaPoint v=d->points[d->surface_stations+e*d->wake_count+j];
            fprintf(out,"%s[%.12g,%.12g]",j?",":"",v.x,v.y);
        }
        fputc(']',out);
    }
    fputc(']',out);
}

static void fields(FILE *out,const AlulaPanelSnapshot *s)
{
    int wake=isfinite(s->wake_residual);
    fprintf(out,",\"mach\":0,\"stage\":\"%s\",\"iteration\":%d,\"wakeIteration\":%d",
        wake?"wake":"panel/BL",wake?s->wake_iteration:s->iteration,s->wake_iteration);
    field(out,"actualAlpha",s->alpha);field(out,"targetAlpha",s->target_alpha);
    field(out,"equationResidual",s->equation_residual);field(out,"wakeResidual",s->wake_residual);
    field(out,"residual",wake?s->wake_residual:s->equation_residual);
    field(out,"cl",s->cl);field(out,"cm",s->cm);field(out,"pressureDrag",s->pressure_drag);
    field(out,"viscousWakeDrag",s->wake_drag);
}

typedef struct {int sampled,has_frame;double next_plot;} Output;
static void progress(void *context,const AlulaPanelSnapshot *s)
{
    Output *output=context;
    struct timespec ts;
    int clock_ok=clock_gettime(CLOCK_MONOTONIC,&ts)==0;
    double now=clock_ok?ts.tv_sec+ts.tv_nsec*1e-9:0;
    fputs("{\"type\":\"iterate\"",stderr);fields(stderr,s);
    if (!output->sampled || !output->has_frame || !clock_ok || now>=output->next_plot) {
        plots(stderr,s);output->has_frame=1;output->next_plot=now+.05;
    } else fputs(",\"plotsDeferred\":true",stderr);
    fputs("}\n",stderr);fflush(stderr);
}

static int parse(const char *s,double *v)
{
    char *end;errno=0;*v=strtod(s,&end);
    return !errno && end!=s && !*end && isfinite(*v);
}

static int integer(const char *s,int min,int max,int *v)
{
    char *end;errno=0;
    long val=strtol(s,&end,10);
    if (errno || end==s || *end || val<min || val>max) return 0;
    *v=(int)val;return 1;
}

int main(int argc,char **argv)
{
    AlulaPanelViscousOptions o=alula_panel_viscous_defaults();
    double trips[6][2];
    size_t trip_count=0;
    int first=1,topology=0,emit=0;
    Output output={0};
    while (first<argc-3) {
        const char *flag=argv[first++];
        if (!strcmp(flag,"--topology")) topology=1;
        else if (!strcmp(flag,"--progress") || !strcmp(flag,"--live-progress")) {
            emit=1;output.sampled=!strcmp(flag,"--live-progress");
        } else if (!strcmp(flag,"--trip") && first+2<=argc-3 && trip_count<6) {
            for (int k=0;k<2;++k) if (!parse(argv[first++],&trips[trip_count][k]) ||
                !(trips[trip_count][k]>0 && trips[trip_count][k]<=1)) goto usage;
            ++trip_count;
        } else if (!strcmp(flag,"--wake-count") && first<argc-3) {
            int n;
            if (!integer(argv[first++],4,256,&n)) goto usage;
            o.wake_count=(size_t)n;
        } else if (!strcmp(flag,"--wake-limit") && first<argc-3) {
            if (!integer(argv[first++],0,200,&o.wake_limit)) goto usage;
        } else if (!strcmp(flag,"--iterations") && first<argc-3) {
            if (!integer(argv[first++],0,200,&o.newton_limit)) goto usage;
        } else goto usage;
    }
    if (argc!=first+3 || !parse(argv[first],&o.alpha) || !parse(argv[first+1],&o.reynolds) ||
        !parse(argv[first+2],&o.ncrit) || fabs(o.alpha)>20 || o.reynolds<5e4 || o.reynolds>5e7 ||
        o.ncrit<0 || o.ncrit>14) goto usage;
    AlulaContour contours[6];
    AlulaPoint points[ALULA_PANEL_MAX_TOTAL+6];
    size_t count,used;
    if (!alula_cli_contours(topology,contours,points,&count,&used)) {fputs("Invalid contours.\n",stderr);return 2;}
    if (trip_count && trip_count!=count) goto usage;
    o.trips=trip_count?(const double (*)[2])trips:NULL;
    AlulaPanelViscousResult r=alula_panel_viscous_solve(contours,count,&o,emit?progress:NULL,&output);
    static const char *reason[]={"converged","invalid input","invalid panel or wake geometry","boundary-layer initialization failed",
        "panel/BL Newton solve stopped","wake flow reversed","wake iteration limit","outside thin-layer model"};
    int solved=r.status==ALULA_VISCOUS_OK;
    printf("{\"model\":\"panel-bl\",\"converged\":%s,\"reason\":\"%s\",\"blStatus\":%d,\"alpha\":%.17g",
        solved?"true":"false",reason[r.status],r.bl_status,o.alpha);
    if (r.final.outer && r.final.profile && r.final.profile->complete) {
        fields(stdout,&r.final);plots(stdout,&r.final);
    }
    puts("}");alula_panel_viscous_free(&r);
    return solved?0:1;
usage:
    fprintf(stderr,"Usage: %s [--topology] [--progress | --live-progress] [--trip TOP BOTTOM]... "
        "[--wake-count N] [--wake-limit N] [--iterations N] ALPHA RE NCRIT < contours.txt\n",argv[0]);
    return 2;
}
