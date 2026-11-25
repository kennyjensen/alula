// SPDX-License-Identifier: GPL-2.0-or-later
#define _POSIX_C_SOURCE 200809L
#include "airfoil_grid_internal.h"
#include "airfoil_bl.h"
#include "airfoil_coupled.h"
#include "cli_geometry.h"
#include "body_internal.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

enum { MAX_INTERVALS=512, MAX_TUBES=256 };
typedef struct { int intervals,tubes; } CoarseStage;

static int number(const char *s,double *value)
{
    char *end;
    errno=0;
    double parsed=strtod(s,&end);
    if (errno || end==s || *end || !isfinite(parsed)) return 0;
    *value=parsed;
    return 1;
}

static int resolution(const char *s,int min,int max,int *value)
{
    char *end;
    errno=0;
    long parsed=strtol(s,&end,10);
    if (errno || end==s || *end || parsed<min || parsed>max) return 0;
    *value=(int)parsed;
    return 1;
}

static int coarse_stages(const char *s,CoarseStage *stages)
{
    int count=0,previous=3;
    for (;;) {
        char *end;
        errno=0;
        long intervals=strtol(s,&end,10);
        if (errno || end==s || intervals<=previous || intervals>=MAX_INTERVALS) return 0;
        long tubes=0;
        if (*end=='x') {
            s=end+1;
            tubes=strtol(s,&end,10);
            if (errno || end==s || tubes<1 || tubes>MAX_TUBES) return 0;
        }
        stages[count++]=(CoarseStage){(int)intervals,(int)tubes};
        previous=(int)intervals;
        if (!*end) return count;
        if (*end!=',') return 0;
        s=end+1;
    }
}

static int usage(const char *program)
{
    fprintf(stderr,"Usage: %s [--topology] [--potential-grid] [--coarse-start N[xT][,N[xT]...]]\n"
        "       [--inlet-intervals N] [--wake-intervals N] [--domain XMIN XMAX YMIN YMAX]\n"
        "       [--surface-uniformity FRACTION]\n"
        "       [--mesh-only | --bl RE NCRIT | --coupled RE NCRIT] [--refine-coupled] [--mach-start MACH]\n"
        "       [--cp] [--progress | --live-progress] [--trip UPPER LOWER]...\n"
        "       ALPHA MACH [SURFACE_INTERVALS [TUBES]] < contours.txt\n",program);
    return 2;
}

static void write_scalar(FILE *stream,const char *label,double value)
{
    if (isfinite(value)) fprintf(stream,",\"%s\":%.17g",label,value);
    else fprintf(stream,",\"%s\":null",label);
}

static void scalar(const char *label,double value) { write_scalar(stdout,label,value); }
typedef struct { FILE *stream; int started; } PressureWriter;

static void pressure_sample(void *context,int contour,AlulaAirfoilPressurePart part,AlulaPoint p,double cp)
{
    PressureWriter *writer=context;
    FILE *stream=writer->stream;
    const char *parts[]={"upper","lower","base"};
    fputs(writer->started?",":",\"cp\":[",stream);
    writer->started=1;
    fprintf(stream,"{\"contour\":%d,\"surface\":\"%s\"",contour,parts[part]);
    write_scalar(stream,"x",p.x);
    write_scalar(stream,"y",p.y);
    write_scalar(stream,"cp",cp);
    fputc('}',stream);
}

static void family_residuals(const double *residuals)
{
    static const char *const names[ALULA_BODY_FAMILIES]={
        [ALULA_BODY_STREAMWISE]="streamwise", [ALULA_BODY_INLET_DENSITY]="inletDensity",
        [ALULA_BODY_INTERNAL_PRESSURE]="internalPressure", [ALULA_BODY_FARFIELD_PRESSURE]="farfieldPressure",
        [ALULA_BODY_CUT_PRESSURE]="cutPressure", [ALULA_BODY_END_TANGENCY]="endTangency",
        [ALULA_BODY_WAKE_GAP]="wakeGap", [ALULA_BODY_TRAILING_KUTTA]="trailingKutta",
        [ALULA_BODY_LEADING_KUTTA]="leadingKutta", [ALULA_BODY_FARFIELD_MATCH]="farfieldMatch",
    };
    fputs(",\"residualByFamily\":{",stdout);
    for (int k=0;k<ALULA_BODY_FAMILIES;++k)
        printf("%s\"%s\":%.17g",k?",":"",names[k],residuals[k]);
    putchar('}');
}

static int boundary_layer(const AlulaAirfoilBL *profiles,size_t count,const char *edge_source)
{
    int complete=alula_airfoil_bl_complete(profiles);
    printf(",\"boundaryLayer\":{\"edgeSource\":\"%s\",\"complete\":%s,\"profiles\":[",
        edge_source,complete?"true":"false");
    if (profiles) for (size_t e=0;e<count;++e) for (int side=0;side<3;++side) {
        const AlulaAirfoilBLBranch *p=alula_airfoil_bl_branch(profiles,(int)e,(AlulaBLBranch)side);
        const char *names[]={"upper","lower","wake"};
        printf("%s{\"body\":%zu,\"side\":\"%s\",\"status\":%d,\"stations\":%zu,\"completed\":%zu",
            e||side?",":"",e,names[side],(int)p->result.status,p->count,p->result.completed);
        if (p->result.transition.active) scalar("transitionS",p->result.transition.s);
        if (p->result.completed) {
            const AlulaBLStation *s=p->stations+p->result.completed-1;
            scalar("momentumThickness",s->theta);
            scalar("displacementThickness",s->delta_star);
            scalar("maximumResidual",p->result.maximum_residual);
        }
        if (p->result.edge_residual>0) {
            size_t index=p->result.edge_mismatch_index;
            const AlulaBLStation *s=p->stations+index;
            printf(",\"edgeMismatch\":{\"station\":%zu",index);
            scalar("s",s->s);
            scalar("residual",p->result.edge_residual);
            scalar("eulerSpeed",p->result.euler_edge_speed);
            scalar("blSpeed",s->ue);
            putchar('}');
        }
        putchar('}');
    }
    fputs("]}",stdout);
    return complete;
}

typedef struct {
    AlulaAirfoilGrid *grid;
    AlulaAirfoilGridReport report;
    AlulaBodySolveResult result;
    int transferred,wall_following,repaired;
    double uniformity;
    int coarse;
} EulerStart;

static int progress_enabled,live_progress;

/* The desktop can repaint at most once per polling interval. Keep all scalar
 * samples while avoiding repeated serialization of paths between repaints. */
static int plot_due(const char *stage,double mach)
{
    static const char *last_stage;
    static double last_mach,last_plot;
    struct timespec now;
    if (!live_progress || clock_gettime(CLOCK_MONOTONIC,&now)!=0) return 1;
    double seconds=(double)now.tv_sec+1e-9*now.tv_nsec;
    if (last_stage && strcmp(last_stage,stage)==0 && mach==last_mach && seconds-last_plot<.05)
        return 0;
    last_plot=seconds;
    last_stage=stage;
    last_mach=mach;
    return 1;
}

/* Sample tube centers from the same evaluated cell geometry as the pressure
 * integral. Bound the plotting traffic independently of the mesh resolution. */
static void write_streamlines(const AlulaAirfoilGrid *g)
{
    const AlulaBody *b=g->body;
    int first=1, stride=(b->in.nx+62)/64;
    fputs(",\"streamlines\":[",stderr);
    for (int group=0;group<=b->in.bodies;++group) {
        int tubes=b->tubes[group],count=tubes<8?tubes:8;
        for (int k=0;k<count;++k) {
            int tube=count==1?0:k*(tubes-1)/(count-1);
            fputs(first?"[":",[",stderr);
            first=0;
            for (int i=1;i<b->in.nx;i+=stride) {
                const AlulaCellInput *input=&cell(b,group,i,tube)->input;
                AlulaPoint center=mean(input->lower[1],input->upper[1]);
                fprintf(stderr,"%s[%.7g,%.7g]",i==1?"":",",center.x,center.y);
            }
            fputc(']',stderr);
        }
    }
    fputc(']',stderr);
}

/* Surface arc from stagnation and displacement thickness, both in
 * reference chords. Partial marches publish only completed stations. */
static void write_bl_curves(const AlulaAirfoilGrid *g,const AlulaAirfoilBL *profiles)
{
    if (!profiles) return;
    fputs(",\"boundaryLayerCurves\":[",stderr);
    int written=0;
    for (int e=0;e<g->elements;++e) for (int side=0;side<2;++side) {
        const AlulaAirfoilBLBranch *branch=alula_airfoil_bl_branch(profiles,e,(AlulaBLBranch)side);
        if (!branch || !branch->result.completed) continue;
        fprintf(stderr,"%s{\"contour\":%d,\"surface\":\"%s\",\"samples\":[",written?",":"",
            g->contour_index[e],side==ALULA_BL_UPPER?"upper":"lower");
        for (size_t i=0;i<branch->result.completed;++i) {
            const AlulaBLStation *p=branch->stations+i;
            fprintf(stderr,"%s[%.9g,%.9g]",i?",":"",p->s,p->delta_star);
        }
        fputs("]}",stderr);
        written=1;
    }
    fputc(']',stderr);
}

static void progress_event(AlulaAirfoilGrid *g,const AlulaAirfoilBL *profiles,int iteration,
    double residual,const char *stage)
{
    double mach=alula_airfoil_grid_mach(g);
    int plot=plot_due(stage,mach);
    fprintf(stderr,"{\"type\":\"iterate\",\"stage\":\"%s\",\"iteration\":%d",stage,iteration);
    write_scalar(stderr,"mach",mach);
    write_scalar(stderr,"residual",residual);
    PressureWriter writer={.stream=stderr};
    AlulaAirfoilPressureLoads loads;
    int valid=alula_airfoil_pressure_prepared(g,&loads,plot?pressure_sample:NULL,&writer);
    if (writer.started) fputc(']',stderr);
    if (valid) {
        write_scalar(stderr,"cl",loads.cl);
        write_scalar(stderr,"cm",loads.cm);
        write_scalar(stderr,"pressureDrag",loads.pressure_drag);
        if (plot) write_streamlines(g);
        else fputs(",\"plotsDeferred\":true",stderr);
    }
    if (plot) write_bl_curves(g,profiles);
    double drag;
    if (alula_airfoil_bl_wake_drag(profiles,&drag)) write_scalar(stderr,"viscousWakeDrag",drag);
    fputs("}\n",stderr);
    fflush(stderr);
}

static void euler_progress(void *context,int iteration,double residual)
{
    EulerStart *start=context;
    progress_event(start->grid,NULL,iteration,residual,start->coarse?"Coarse Euler":"Euler");
}

static void transfer_progress(void *grid,int iteration,double residual)
{
    progress_event(grid,NULL,iteration,residual,"Displacement transfer");
}

static void coupled_progress(AlulaAirfoilGrid *g,const AlulaAirfoilBL *profiles,int iteration,
    double residual,void *stage)
{
    progress_event(g,profiles,iteration,residual,stage);
}

static void transfer_start(EulerStart *start,AlulaAirfoilGrid *coarse)
{
    start->transferred=alula_airfoil_grid_transfer(start->grid,coarse,.75);
    if (start->transferred) start->report.flow_seed=ALULA_BODY_OK;
}

static void solve_start(EulerStart *start)
{
    start->grid->progress=progress_enabled?euler_progress:NULL;
    start->grid->progress_context=start;
    start->result=alula_airfoil_grid_solve(start->grid,40,1e-9);
    if (start->wall_following || (start->transferred && start->result.status==ALULA_BODY_INVALID_FINAL_GRID)) {
        start->repaired=alula_airfoil_grid_repair(start->grid,40,1e-9,0,&start->result);
        if (!start->repaired && start->result.status==ALULA_BODY_INVALID_FINAL_GRID)
            start->repaired=alula_airfoil_grid_repair(start->grid,40,1e-9,1,&start->result);
    }
    start->grid->progress=NULL;
    start->grid->progress_context=NULL;
}

static int stalled_grid(const EulerStart *start)
{
    return start->result.status==ALULA_BODY_LINE_SEARCH_FAILURE ||
        start->result.status==ALULA_BODY_INVALID_FINAL_GRID;
}

static void retry_start(EulerStart *start,AlulaAirfoilGridSeed *seed,AlulaAirfoilGridInput mesh,
    AlulaAirfoilGrid *coarse,int wall_following)
{
    EulerStart candidate={.wall_following=wall_following,.uniformity=mesh.surface_uniformity,.coarse=start->coarse};
    mesh.crosslines=ALULA_CROSSLINE_POTENTIAL;
    candidate.grid=alula_airfoil_grid_create_seeded(seed,&mesh,wall_following,NULL,&candidate.report);
    if (!candidate.grid) return;
    transfer_start(&candidate,coarse);
    solve_start(&candidate);
    if (candidate.result.status==ALULA_BODY_SOLVED || candidate.result.norm<start->result.norm) {
        alula_airfoil_grid_free(start->grid);
        *start=candidate;
    } else alula_airfoil_grid_free(candidate.grid);
}

static AlulaAirfoilGrid *coarse_solution(AlulaAirfoilGridSeed *seed,size_t count,
    AlulaAirfoilGridInput input,const CoarseStage *stages,int stage_count,int *solved)
{
    int tubes[7],default_tubes=input.tubes[0]>4?4:input.tubes[0];
    input.crosslines=ALULA_CROSSLINE_POTENTIAL;
    input.tubes=tubes;
    input.surface_uniformity=0;
    AlulaAirfoilGrid *coarse=NULL;
    for (int i=0;i<stage_count;++i) {
        int nt=stages[i].tubes?stages[i].tubes:default_tubes;
        for (size_t j=0;j<=count;++j) tubes[j]=nt;
        input.surface_intervals=stages[i].intervals;
        input.inlet_intervals=input.outlet_intervals=input.surface_intervals/2+4;
        EulerStart stage={.coarse=1};
        stage.grid=alula_airfoil_grid_create_seeded(seed,&input,0,NULL,&stage.report);
        if (!stage.grid) continue;
        transfer_start(&stage,coarse);
        solve_start(&stage);
        if (count>1 && stalled_grid(&stage))
            retry_start(&stage,seed,input,coarse,1);
        if (stage.result.status==ALULA_BODY_SOLVED) {
            alula_airfoil_grid_free(coarse);
            coarse=stage.grid;
            ++*solved;
        } else alula_airfoil_grid_free(stage.grid);
    }
    return coarse;
}

static const double (*order_trips(AlulaAirfoilGrid *grid,size_t count,
    double (*trips)[2],double (*ordered)[2]))[2]
{
    if (!grid || !count) return NULL;
    for (size_t i=0;i<count;++i)
        memcpy(ordered[i],trips[alula_airfoil_grid_contour_index(grid,(int)i)],sizeof ordered[i]);
    return (const double (*)[2])ordered;
}

int main(int argc,char **argv)
{
    int first=1,topology=0,mesh_only=0,bl=0,potential_grid=0,stage_count=0;
    int inlet=0,wake=0,fixed_domain=0,pressure=0,refine_coupled=0;
    CoarseStage stages[MAX_INTERVALS-4];
    double reynolds=0,ncrit=0,uniformity=-1,start_mach=0;
    double trips[6][2],domain[4];
    size_t trip_count=0;
    while (first<argc) {
        if (strcmp(argv[first],"--topology")==0 && !topology) topology=1;
        else if (strcmp(argv[first],"--mesh-only")==0 && !mesh_only) mesh_only=1;
        else if (strcmp(argv[first],"--cp")==0 && !pressure) pressure=1;
        else if (strcmp(argv[first],"--progress")==0 && !progress_enabled) progress_enabled=1;
        else if (strcmp(argv[first],"--live-progress")==0 && !progress_enabled) progress_enabled=live_progress=1;
        else if (strcmp(argv[first],"--refine-coupled")==0 && !refine_coupled) refine_coupled=1;
        else if (strcmp(argv[first],"--potential-grid")==0 && !potential_grid) potential_grid=1;
        else if (!strcmp(argv[first],"--mach-start") && !start_mach && first+1<argc) {
            if (!number(argv[++first],&start_mach) || !(start_mach>0 && start_mach<1)) return usage(argv[0]);
        }
        else if (!strcmp(argv[first],"--surface-uniformity") && uniformity<0 && first+1<argc) {
            if (!number(argv[++first],&uniformity) || uniformity<0 || uniformity>1) return usage(argv[0]);
        }
        else if (!strcmp(argv[first],"--inlet-intervals") && !inlet && first+1<argc) {
            if (!resolution(argv[++first],2,MAX_INTERVALS,&inlet)) return usage(argv[0]);
        }
        else if (!strcmp(argv[first],"--wake-intervals") && !wake && first+1<argc) {
            if (!resolution(argv[++first],2,MAX_INTERVALS,&wake)) return usage(argv[0]);
        }
        else if (!strcmp(argv[first],"--domain") && !fixed_domain && first+4<argc) {
            for (int i=0;i<4;++i) if (!number(argv[++first],domain+i)) return usage(argv[0]);
            if (!(domain[0]<domain[1] && domain[2]<domain[3])) return usage(argv[0]);
            fixed_domain=1;
        }
        else if (strcmp(argv[first],"--coarse-start")==0 && !stage_count && first+1<argc) {
            stage_count=coarse_stages(argv[++first],stages);
            if (!stage_count) return usage(argv[0]);
        }
        else if ((!strcmp(argv[first],"--bl") || !strcmp(argv[first],"--coupled")) && !bl && first+2<argc) {
            if (!number(argv[first+1],&reynolds) || !(reynolds>0) ||
                !number(argv[first+2],&ncrit) || ncrit<0) return usage(argv[0]);
            bl=strcmp(argv[first],"--bl")==0?1:2;
            first+=2;
        }
        else if (!strcmp(argv[first],"--trip") && trip_count<6 && first+2<argc) {
            for (int side=0;side<2;++side) {
                double *fraction=&trips[trip_count][side];
                if (!number(argv[first+1+side],fraction) || !(*fraction>0 && *fraction<=1))
                    return usage(argv[0]);
            }
            ++trip_count;
            first+=2;
        }
        else break;
        ++first;
    }
    if (argc-first<2 || argc-first>4 || (mesh_only && (bl || stage_count || pressure)) ||
        (trip_count && !bl) || (refine_coupled && (bl!=2 || !stage_count))) return usage(argv[0]);
    AlulaAirfoilGridInput in={.surface_uniformity=fmax(0,uniformity),.chord=1,.moment={.25,0},.surface_intervals=6,
        .mode=ALULA_BODY_HYBRID,.ismom=4,
        .crosslines=potential_grid?ALULA_CROSSLINE_POTENTIAL:ALULA_CROSSLINE_AUTO};
    int tubes=4,counts[7];
    if (!number(argv[first],&in.alpha) || fabs(in.alpha)>90 ||
        !number(argv[first+1],&in.mach) || !(in.mach>0 && in.mach<1) ||
        (argc>first+2 && !resolution(argv[first+2],4,MAX_INTERVALS,&in.surface_intervals)) ||
        (argc>first+3 && !resolution(argv[first+3],1,MAX_TUBES,&tubes))) return usage(argv[0]);
    if (stage_count && stages[stage_count-1].intervals>=in.surface_intervals) return usage(argv[0]);
    double target_mach=in.mach;
    if (start_mach) {
        if (bl!=2 || start_mach>target_mach) return usage(argv[0]);
        in.mach=start_mach;
    }
    for (int i=0;i<7;++i) counts[i]=tubes;
    in.tubes=counts;
    in.inlet_intervals=inlet?inlet:in.surface_intervals/2+4;
    in.outlet_intervals=wake?wake:in.surface_intervals/2+4;
    AlulaContour contours[6];
    AlulaPoint points[ALULA_PANEL_MAX_TOTAL+6];
    size_t elements,used;
    if (!alula_cli_contours(topology,contours,points,&elements,&used)) {
        fputs("Invalid geometry input. Use closed counterclockwise contours.\n",stderr);
        return 2;
    }
    if (trip_count && trip_count!=elements) return usage(argv[0]);
    in.inlet_x=in.lower_y=INFINITY;
    in.outlet_x=in.upper_y=-INFINITY;
    for (size_t i=0;i<used;++i) {
        in.inlet_x=fmin(in.inlet_x,points[i].x);
        in.outlet_x=fmax(in.outlet_x,points[i].x);
        in.lower_y=fmin(in.lower_y,points[i].y);
        in.upper_y=fmax(in.upper_y,points[i].y);
    }
    if (fixed_domain) {
        if (!(domain[0]<in.inlet_x && domain[1]>in.outlet_x &&
            domain[2]<in.lower_y && domain[3]>in.upper_y)) {
            fputs("Domain must strictly enclose all contours.\n",stderr);
            return 2;
        }
        in.inlet_x=domain[0];in.outlet_x=domain[1];in.lower_y=domain[2];in.upper_y=domain[3];
    } else {
        in.inlet_x-=3;in.outlet_x+=3;in.lower_y-=3;in.upper_y+=3;
    }
    AlulaAirfoilGridStatus status;
    EulerStart start={.uniformity=in.surface_uniformity};
    if (progress_enabled) setvbuf(stderr,NULL,_IOFBF,BUFSIZ);
    AlulaCoupledControls controls={.reynolds=reynolds,.ncrit=ncrit,.displacement_tolerance=1e-8,
        .euler_tolerance=1e-9,.relaxation=.65,.max_iterations=32,.euler_iterations=24,
        .progress=progress_enabled?coupled_progress:NULL,.progress_context="Euler/BL"};
    AlulaAirfoilGridSeed *seed=alula_airfoil_grid_prepare(contours,elements,&in,!mesh_only,&status);
    int coarse_solved=0;
    AlulaAirfoilGrid *coarse=seed?coarse_solution(seed,elements,in,stages,stage_count,&coarse_solved):NULL;
    start.grid=seed?alula_airfoil_grid_create_seeded(seed,&in,0,&status,&start.report):NULL;
    if (start.grid && !mesh_only) {
        int algebraic=start.report.crosslines==ALULA_CROSSLINE_ALGEBRAIC;
        int retry_spacing=uniformity<0 && elements==1 && in.surface_intervals>=32;
        int can_retry=algebraic || elements>1 || retry_spacing;
        transfer_start(&start,coarse);
        if (!can_retry && !refine_coupled) {
            alula_airfoil_grid_free(coarse);
            alula_airfoil_grid_seed_free(seed);
            coarse=NULL;
            seed=NULL;
        }
        solve_start(&start);
        if (start.result.status!=ALULA_BODY_SOLVED && algebraic)
            retry_start(&start,seed,in,coarse,0);
        if (elements>1 && stalled_grid(&start))
            retry_start(&start,seed,in,coarse,1);
        if (elements>1 && uniformity<0 && stalled_grid(&start)) {
            /* Rebuild clustered inlet and nose stations on the requested mesh
             * when the moving cut still folds with wall-following coordinates. */
            AlulaAirfoilGridInput retry=in;
            retry.surface_uniformity=.7;
            retry_start(&start,seed,retry,NULL,1);
        }
        if (retry_spacing && start.result.status!=ALULA_BODY_SOLVED) {
            if (!coarse) {
                CoarseStage stage={8,4};
                coarse=coarse_solution(seed,elements,in,&stage,1,&coarse_solved);
            }
            const double spacing[]={1,.7};
            for (size_t k=0;k<2 && start.result.status!=ALULA_BODY_SOLVED;++k) {
                AlulaAirfoilGridInput retry=in;
                retry.surface_uniformity=spacing[k];
                retry_start(&start,seed,retry,coarse,0);
            }
        }
    }
    double refinement_fraction=0;
    if (refine_coupled && coarse && start.result.status==ALULA_BODY_SOLVED && start.grid) {
        double coarse_trips[6][2];
        AlulaCoupledControls coarse_controls=controls;
        coarse_controls.progress_context="Coarse Euler/BL";
        coarse_controls.trips=order_trips(coarse,trip_count,trips,coarse_trips);
        AlulaCoupledResult result=alula_airfoil_coupled_solve(coarse,&coarse_controls);
        if (result.status==ALULA_COUPLED_OK) {
            start.grid->progress=progress_enabled?transfer_progress:NULL;
            start.grid->progress_context=start.grid;
            refinement_fraction=alula_airfoil_grid_continue_displacement(start.grid,coarse);
            start.grid->progress=NULL;
            start.grid->progress_context=NULL;
        }
        alula_airfoil_bl_free(result.profiles);
        start.result=start.grid->last_solve;
    }
    alula_airfoil_grid_free(coarse);
    alula_airfoil_grid_seed_free(seed);
    AlulaAirfoilGrid *g=start.grid;
    AlulaAirfoilGridReport report=start.report;
    double grid_trips[6][2];
    const double (*material_trips)[2]=order_trips(g,trip_count,trips,grid_trips);
    controls.trips=material_trips;
    if (refinement_fraction>0) {
        controls.displacement_scale=.0005;
        controls.max_iterations=800;
    }
    printf("{\"alpha\":%.17g,\"mach\":%.17g,\"gridStatus\":%d,\"intervals\":%d,\"invalidCells\":%d",
        in.alpha,target_mach,(int)status,report.nx,report.invalid_cells);
    if (start_mach) scalar("startMach",start_mach);
    if (report.nx) scalar("minimumCornerSine",report.minimum_corner_sine);
    if (report.crosslines!=ALULA_CROSSLINE_AUTO)
        printf(",\"crosslines\":\"%s\"",report.crosslines==ALULA_CROSSLINE_POTENTIAL?"potential":"algebraic");
    if (start.uniformity>0 || uniformity>=0) scalar("surfaceUniformity",start.uniformity);
    if (stage_count || coarse_solved)
        printf(",\"coarseStartApplied\":%s,\"coarseStagesSolved\":%d",start.transferred?"true":"false",coarse_solved);
    if (refine_coupled) scalar("coupledRefinementFraction",refinement_fraction);
    if (start.wall_following) printf(",\"startupRecovery\":\"wall-following\"");
    if (start.wall_following || start.repaired)
        printf(",\"startupGridRepaired\":%s",start.repaired?"true":"false");
    if (report.invalid_cells)
        printf(",\"firstInvalidCell\":[%d,%d,%d]",report.first_group,report.first_interval,report.first_tube);
    int converged=0,profile_ok=!bl;
    if (g) {
        AlulaBodySize size=alula_airfoil_grid_size(g);
        printf(",\"flowSeedStatus\":%d,\"nodes\":%d,\"unknowns\":%d,\"wakeWidths\":%d",
            (int)report.flow_seed,size.nodes,size.rows,size.wake_widths);
        if (!mesh_only) {
            AlulaBodySolveResult r=start.result;
            converged=r.status==ALULA_BODY_SOLVED;
            AlulaAirfoilBL *profiles=NULL;
            if (bl==2 && converged) {
                AlulaCoupledResult result=start_mach?alula_airfoil_coupled_solve_to_mach(g,&controls,target_mach):
                    alula_airfoil_coupled_solve(g,&controls);
                r=result.euler;
                profiles=result.profiles;
                converged=result.status==ALULA_COUPLED_OK;
                printf(",\"coupled\":{\"method\":\"%s\",\"coupledNewtonUpdates\":%d,\"status\":%d,\"updates\":%d,"
                    "\"backtracks\":%d,\"eulerUpdates\":%d,\"eulerGridRepairs\":%d,"
                    "\"usedInteraction\":%s,\"usedLocalInteraction\":%s,"
                    "\"rejections\":{\"euler\":%d,\"boundaryLayer\":%d,\"merit\":%d}",
                    converged && result.coupled_newton_updates?"coupled-newton":"partitioned-displacement",
                    result.coupled_newton_updates,(int)result.status,result.iterations,result.backtracks,result.euler_updates,result.euler_grid_repairs,
                    result.used_interaction?"true":"false",result.used_local_interaction?"true":"false",
                    result.euler_rejections,result.bl_rejections,result.merit_rejections);
                scalar("displacementResidual",result.displacement_residual);
                scalar("blResidual",result.bl_residual);
                scalar("edgeResidual",result.edge_residual);
                if (start_mach) printf(",\"machStages\":%d,\"machRetries\":%d",result.mach_stages,result.mach_retries);
                putchar('}');
                const char *edge_source=result.iterations?"displaced-euler":"converged-euler";
                if (!converged && result.used_interaction) edge_source="interaction-iterate";
                profile_ok=boundary_layer(profiles,elements,edge_source);
            }
            if (start_mach) scalar("actualMach",alula_airfoil_grid_mach(g));
            printf(",\"solveStatus\":%d,\"converged\":%s,\"iterations\":%d,\"backtracks\":%d",
                (int)r.status,converged?"true":"false",r.iterations,r.backtracks);
            printf(",\"linear\":{\"factors\":%d,\"refactors\":%d,\"factorReuses\":%d,\"factorNonzeros\":%d",
                r.linear.factor_count,r.linear.refactor_count,r.linear.factor_reuses,r.linear.factor_nonzeros);
            if (r.status==ALULA_BODY_LINEAR_FAILURE && r.linear.failed_relative_residual>0)
                scalar("rejectedRelativeResidual",r.linear.failed_relative_residual);
            putchar('}');
            scalar("residual",r.norm);
            scalar("maxMach",r.max_mach);
            if (!converged && isfinite(r.norm)) family_residuals(r.residual_by_family);
            if (r.backtracks) {
                printf(",\"stepRejections\":{\"chart\":%d,\"flowDomain\":%d,\"residual\":%d",
                    r.rejected_chart,r.rejected_flow,r.rejected_residual);
                if (r.rejected_flow) printf(",\"lastFlowStatus\":%d",(int)r.last_rejected_flow);
                putchar('}');
            }
            if (!converged) {
                AlulaAirfoilGridReport final;
                if (alula_airfoil_grid_current_quality(g,&final)) {
                    printf(",\"finalGrid\":{\"invalidCells\":%d",final.invalid_cells);
                    scalar("minimumCornerSine",final.minimum_corner_sine);
                    if (final.invalid_cells)
                        printf(",\"firstInvalidCell\":[%d,%d,%d]",final.first_group,final.first_interval,final.first_tube);
                    putchar('}');
                }
            }
            AlulaBodyState state;
            if (alula_airfoil_grid_state(g,&state)) scalar("circulation",state.circulation);
            if (bl==1 && converged) {
                profiles=alula_airfoil_bl_create(g,reynolds,ncrit,material_trips);
                profile_ok=boundary_layer(profiles,elements,"converged-euler");
            }
            AlulaAirfoilPressureLoads loads;
            PressureWriter writer={.stream=stdout};
            int have_loads=alula_airfoil_grid_pressure_results(g,&loads,pressure?pressure_sample:NULL,&writer);
            if (writer.started) {
                putchar(']');
                if (!converged || !have_loads) fputs(",\"cpProvisional\":true",stdout);
            }
            if (have_loads) {
                scalar("cl",loads.cl);
                scalar("cm",loads.cm);
                scalar("pressureDrag",loads.pressure_drag);
                if (progress_enabled) {
                    fputs("{\"type\":\"retained\"",stderr);
                    write_streamlines(g);
                    write_bl_curves(g,profiles);
                    fputs("}\n",stderr);
                    fflush(stderr);
                }
            }
            double wake_drag;
            if (alula_airfoil_bl_wake_drag(profiles,&wake_drag)) scalar("viscousWakeDrag",wake_drag);
            alula_airfoil_bl_free(profiles);
        }
    }
    puts("}");
    alula_airfoil_grid_free(g);
    return status==ALULA_GRID_OK && (mesh_only || converged) && profile_ok ? 0 : 1;
}
