// SPDX-License-Identifier: GPL-2.0-or-later
#include "airfoil_grid.h"
#include "airfoil_bl.h"
#include "airfoil_coupled.h"
#include "cli_geometry.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

static int usage(const char *program)
{
    fprintf(stderr,"Usage: %s [--topology] [--mesh-only | --bl RE NCRIT | --coupled RE NCRIT] ALPHA MACH [SURFACE_INTERVALS [TUBES]] < contours.txt\n",program);
    return 2;
}

static void scalar(const char *label,double value)
{
    if (isfinite(value)) printf(",\"%s\":%.17g",label,value);
    else printf(",\"%s\":null",label);
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
        putchar('}');
    }
    fputs("]}",stdout);
    return complete;
}

int main(int argc,char **argv)
{
    int first=1,topology=0,mesh_only=0,bl=0;
    double reynolds=0,ncrit=0;
    while (first<argc) {
        if (strcmp(argv[first],"--topology")==0 && !topology) topology=1;
        else if (strcmp(argv[first],"--mesh-only")==0 && !mesh_only) mesh_only=1;
        else if ((!strcmp(argv[first],"--bl") || !strcmp(argv[first],"--coupled")) && !bl && first+2<argc) {
            if (!number(argv[first+1],&reynolds) || !(reynolds>0) ||
                !number(argv[first+2],&ncrit) || ncrit<0) return usage(argv[0]);
            bl=strcmp(argv[first],"--bl")==0?1:2;
            first+=2;
        }
        else break;
        ++first;
    }
    if (argc-first<2 || argc-first>4 || (mesh_only && bl)) return usage(argv[0]);
    AlulaAirfoilGridInput in={.chord=1,.moment={.25,0},.surface_intervals=6,
        .mode=ALULA_BODY_HYBRID,.ismom=4};
    int tubes=4,counts[7];
    if (!number(argv[first],&in.alpha) || fabs(in.alpha)>90 ||
        !number(argv[first+1],&in.mach) || !(in.mach>0 && in.mach<1) ||
        (argc>first+2 && !resolution(argv[first+2],4,512,&in.surface_intervals)) ||
        (argc>first+3 && !resolution(argv[first+3],1,256,&tubes))) return usage(argv[0]);
    for (int i=0;i<7;++i) counts[i]=tubes;
    in.tubes=counts;
    in.inlet_intervals=in.outlet_intervals=in.surface_intervals/2+4;
    AlulaContour contours[6];
    AlulaPoint points[ALULA_PANEL_MAX_TOTAL+6];
    size_t elements,used;
    if (!alula_cli_contours(topology,contours,points,&elements,&used)) {
        fputs("Invalid geometry input. Use closed counterclockwise contours.\n",stderr);
        return 2;
    }
    in.inlet_x=in.lower_y=INFINITY;
    in.outlet_x=in.upper_y=-INFINITY;
    for (size_t i=0;i<used;++i) {
        in.inlet_x=fmin(in.inlet_x,points[i].x);
        in.outlet_x=fmax(in.outlet_x,points[i].x);
        in.lower_y=fmin(in.lower_y,points[i].y);
        in.upper_y=fmax(in.upper_y,points[i].y);
    }
    in.inlet_x-=3;in.outlet_x+=3;in.lower_y-=3;in.upper_y+=3;
    AlulaAirfoilGridStatus status;
    AlulaAirfoilGridReport report;
    AlulaAirfoilGrid *g=alula_airfoil_grid_create(contours,elements,&in,&status,&report);
    printf("{\"alpha\":%.17g,\"mach\":%.17g,\"gridStatus\":%d,\"intervals\":%d,\"invalidCells\":%d",
        in.alpha,in.mach,(int)status,report.nx,report.invalid_cells);
    if (report.nx) scalar("minimumCornerSine",report.minimum_corner_sine);
    if (report.invalid_cells)
        printf(",\"firstInvalidCell\":[%d,%d,%d]",report.first_group,report.first_interval,report.first_tube);
    int converged=0,profile_ok=!bl;
    if (g) {
        AlulaBodySize size=alula_airfoil_grid_size(g);
        printf(",\"flowSeedStatus\":%d,\"nodes\":%d,\"unknowns\":%d,\"wakeWidths\":%d",
            (int)report.flow_seed,size.nodes,size.rows,size.wake_widths);
        if (!mesh_only) {
            AlulaBodySolveResult r=alula_airfoil_grid_solve(g,40,1e-9);
            converged=r.status==ALULA_BODY_SOLVED;
            AlulaAirfoilBL *profiles=NULL;
            if (bl==2 && converged) {
                AlulaCoupledControls controls={.reynolds=reynolds,.ncrit=ncrit,
                    .displacement_tolerance=1e-8,.euler_tolerance=1e-9,.relaxation=.65,
                    .max_iterations=32,.euler_iterations=24};
                AlulaCoupledResult result=alula_airfoil_coupled_solve(g,&controls);
                r=result.euler;
                profiles=result.profiles;
                converged=result.status==ALULA_COUPLED_OK;
                printf(",\"coupled\":{\"method\":\"partitioned-displacement\",\"status\":%d,\"updates\":%d,"
                    "\"backtracks\":%d,\"eulerUpdates\":%d",(int)result.status,result.iterations,
                    result.backtracks,result.euler_updates);
                scalar("displacementResidual",result.displacement_residual);
                scalar("blResidual",result.bl_residual);
                putchar('}');
                profile_ok=boundary_layer(profiles,elements,result.iterations?"displaced-euler":"converged-euler");
            }
            printf(",\"solveStatus\":%d,\"converged\":%s,\"iterations\":%d,\"backtracks\":%d",
                (int)r.status,converged?"true":"false",r.iterations,r.backtracks);
            scalar("residual",r.norm);
            scalar("maxMach",r.max_mach);
            AlulaBodyState state;
            if (alula_airfoil_grid_state(g,&state)) scalar("circulation",state.circulation);
            if (bl==1 && converged) {
                profiles=alula_airfoil_bl_create(g,reynolds,ncrit,NULL);
                profile_ok=boundary_layer(profiles,elements,"converged-euler");
            }
            AlulaAirfoilPressureLoads loads;
            if (alula_airfoil_grid_pressure_loads(g,&loads)) {
                scalar("cl",loads.cl);
                scalar("cm",loads.cm);
                scalar("pressureDrag",loads.pressure_drag);
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
