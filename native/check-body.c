// SPDX-License-Identifier: GPL-2.0-or-later
#define _POSIX_C_SOURCE 200809L
#include "body_newton.h"
#include "body_chart.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int read_vector(double *v,int n)
{
    for (int i=0;i<n;++i) if (scanf("%lf",v+i)!=1) return 0;
    return 1;
}

static double milliseconds(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC,&t);
    return 1000.*t.tv_sec+1e-6*t.tv_nsec;
}

typedef struct {
    AlulaBodyPhysicalLayout physical;
    const int *start,*column;
    const double *slope,*base,*initial;
} TestChart;

static int decode(void *context,const double *packed,double *out,double *derivatives)
{
    TestChart *c=context;
    for (int i=0;i<c->physical.size;++i) {
        double dx=0;
        int log_density=i>=c->physical.density && i<c->physical.mass;
        for (int k=c->start[i];k<c->start[i+1];++k)
            dx+=c->slope[k]*(packed[c->column[k]]-c->initial[c->column[k]]);
        out[i]=log_density ? c->base[i]*exp(dx) : c->base[i]+dx;
        if (derivatives) for (int k=c->start[i];k<c->start[i+1];++k)
            derivatives[k]=c->slope[k]*(log_density ? out[i] : 1);
    }
    return 1;
}

static int run_newton(AlulaBody *b,const AlulaBodyState *state)
{
    AlulaBodySize s=alula_body_size(b);
    AlulaBodyPhysicalLayout p=alula_body_physical_layout(s);
    int nnz,iterations;
    double tolerance;
    if (scanf("%d%d%lf",&nnz,&iterations,&tolerance)!=3 || nnz<1 || nnz>20*p.size) return 2;
    int *indices=calloc((size_t)p.size+1+nnz+s.rows,sizeof *indices);
    double *data=calloc((size_t)p.size+nnz+4*s.rows,sizeof *data);
    AlulaBodyNewton *solver=NULL;
    int ret=2;
    if (!indices || !data) goto done;
    int *start=indices,*column=start+p.size+1,*order=column+nnz;
    double *base=data,*slope=base+p.size,*x0=slope+nnz,*direction=x0+s.rows;
    double *product=direction+s.rows,*x=product+s.rows;
    int k=0;
    for (int i=0;i<p.size;++i) {
        int count;
        if (scanf("%d",&count)!=1 || count<0 || count>nnz-k) goto done;
        for (int j=0;j<count;++j,++k)
            if (scanf("%d%lf",column+k,slope+k)!=2) goto done;
        start[i+1]=k;
    }
    if (k!=nnz) goto done;
    for (int i=0;i<s.rows;++i) if (scanf("%d",order+i)!=1) goto done;
    if (!read_vector(x0,s.rows) || !read_vector(direction,s.rows)) goto done;
    for (int i=0;i<s.nodes;++i) {base[2*i]=state->nodes[i].x;base[2*i+1]=state->nodes[i].y;}
    for (int i=0;i<s.density;++i) base[p.density+i]=state->density[i];
    for (int i=0;i<s.tubes;++i) base[p.mass+i]=state->mass[i];
    for (int i=0;i<s.wake_widths;++i) base[p.wake+i]=state->wake_width[i];
    base[p.strength]=state->circulation;
    base[p.strength+1]=state->source;
    base[p.strength+2]=state->doublet_x;
    base[p.strength+3]=state->doublet_y;
    for (int i=0;i<s.rows;++i) x[i]=x0[i];
    TestChart context={p,start,column,slope,base,x0};
    AlulaBodyChart chart={.starts=start,.columns=column,.decode=decode,.context=&context,.equation_order=order};
    solver=alula_body_newton_create(b,&chart);
    if (!solver || !alula_body_newton_product(solver,x,direction,product)) goto done;
    for (int i=0;i<s.rows;++i) printf("%.17g%c",product[i],i==s.rows-1?'\n':' ');
    double t=milliseconds();
    AlulaBodySolveResult result=alula_body_newton_solve(solver,x,iterations,tolerance);
    double elapsed=milliseconds()-t;
    printf("%d %d %d %d %.17g %.17g %.17g %d %d %.17g %d\n",result.status,result.iterations,result.backtracks,
        result.colors,result.norm,result.max_mach,elapsed,result.linear.factor_count,result.linear.refactor_count,
        result.linear.max_relative_residual,result.linear.factor_reuses);
    for (int i=0;i<s.rows;++i) printf("%.17g%c",x[i],i==s.rows-1?'\n':' ');
    ret=0;
done:
    alula_body_newton_free(solver);
    free(indices);free(data);
    return ret;
}

static void emit_vector(const double *v,int n)
{
    for(int i=0;i<n;++i) printf("%.17g%c",v[i],i==n-1?'\n':' ');
}

static int run_coordinates(AlulaBody *b,const AlulaBodyInput *in)
{
    AlulaBodySize s=alula_body_size(b);
    AlulaBodyPhysicalLayout p=alula_body_physical_layout(s);
    AlulaBodySurface surfaces[6]={0};
    AlulaPoint *curve[6]={0};
    double *shape[6]={0};
    AlulaPoint *raw=malloc((size_t)s.nodes*sizeof *raw);
    double *weights=malloc((size_t)s.tubes*sizeof *weights),levels[8];
    double *work=malloc(((size_t)2*p.size+3*(size_t)s.rows)*sizeof *work);
    AlulaBodyCoordinates *coordinates=NULL;
    AlulaBodyNewton *solver=NULL;
    AlulaBodyCoordinateInput input={.raw_nodes=raw,.mass_weights=weights,.capture_levels=levels,.surfaces=surfaces};
    int ret=2,iterations;
    double tolerance;
    if(!raw||!weights||!work || scanf("%d%d%d%d%d%lf",&input.primary_body,&input.walls_only,
        &input.body_station_normals,&input.te_center_wake,&iterations,&tolerance)!=6) goto done;
    for(int i=0;i<s.nodes;++i) if(scanf("%lf%lf",&raw[i].x,&raw[i].y)!=2) goto done;
    if(!read_vector(weights,s.tubes) || !read_vector(levels,in->bodies+2)) goto done;
    for(int g=0;g<in->bodies;++g) {
        AlulaBodySurface *surface=surfaces+g;
        int stations=in->cuts[g].trailing-in->cuts[g].leading+1, wakes=in->nx-in->cuts[g].trailing;
        if(scanf("%d%lf%d",&surface->count,&surface->stagnation,&surface->material_wake)!=3 ||
            surface->count<4 || surface->count>4097) goto done;
        curve[g]=malloc((size_t)surface->count*sizeof **curve);
        shape[g]=malloc((size_t)(4*stations+wakes)*sizeof **shape);
        if(!curve[g]||!shape[g]) goto done;
        for(int i=0;i<surface->count;++i)
            if(scanf("%lf%lf",&curve[g][i].x,&curve[g][i].y)!=2) goto done;
        surface->surface=curve[g];
        surface->lower_fraction=shape[g];
        surface->upper_fraction=shape[g]+stations;
        if(!read_vector(shape[g],2*stations)) goto done;
        if(in->displaced) {
            surface->lower_thickness=shape[g]+2*stations;
            surface->upper_thickness=shape[g]+3*stations;
            surface->wake_width=shape[g]+4*stations;
            if(!read_vector(shape[g]+2*stations,2*stations+wakes)) goto done;
        }
    }
    double *before=work,*after=before+p.size,*x=after+p.size,*direction=x+s.rows,*product=direction+s.rows;
    if(!read_vector(x,s.rows)||!read_vector(direction,s.rows)) goto done;
    coordinates=alula_body_coordinates_create(b,&input);
    if(!coordinates) goto done;
    AlulaBodyChart chart=alula_body_coordinates_chart(coordinates);
    solver=alula_body_newton_create(b,&chart);
    if(!solver || !chart.decode(chart.context,x,before,NULL) ||
        !alula_body_newton_product(solver,x,direction,product) ||
        !chart.rebase(chart.context,x) || !chart.decode(chart.context,x,after,NULL)) goto done;
    double error=0;
    for(int i=0;i<p.size;++i) error=fmax(error,fabs(before[i]-after[i]));
    emit_vector(before,p.size);
    emit_vector(product,s.rows);
    double begin=milliseconds();
    AlulaBodySolveResult r=alula_body_newton_solve(solver,x,iterations,tolerance);
    double elapsed=milliseconds()-begin;
    if(!chart.decode(chart.context,x,after,NULL)) goto done;
    printf("%d %d %d %d %.17g %.17g %.17g\n",r.status,r.iterations,r.backtracks,r.colors,r.norm,error,elapsed);
    emit_vector(after,p.size);
    emit_vector(x,s.rows);
    ret=0;
done:
    alula_body_newton_free(solver);alula_body_coordinates_free(coordinates);
    free(raw);free(weights);free(work);
    for(int g=0;g<6;++g){free(curve[g]);free(shape[g]);}
    return ret;
}

static int check_sparse_cancellation(void)
{
    if (LDBL_MANT_DIG<=DBL_MANT_DIG) return 0;
    const AlulaSparseEntry entries[]={{0,0},{1,0},{0,1},{1,1}};
    const double values[]={0x1.cc35e34ep+0,0x1.465a9254p+0,0x1.cc35de09be825p+0,0x1.465a995bf76f5p+0};
    const double expected[]={0x1.7d784p+26,-0x1.7d783ff721e2ep+26};
    const double rhs[]={0x1.fa3ff0752d4d9p+4,-0x1.4ddb5844e5a4bp+5};
    AlulaSparse *m=alula_sparse_create(2,entries,4,NULL);
    if (!m) return 1;
    memcpy(alula_sparse_values(m),values,sizeof values);
    double x[2];
    int ok=alula_sparse_solve(m,rhs,x);
    if (ok) for (int i=0;i<2;++i) {
        long double residual=(long double)rhs[i]-(long double)values[i]*x[0]-(long double)values[i+2]*x[1];
        if (fabsl(residual)>1e-10*fabs(rhs[i]) || fabs(x[i]-expected[i])>.01) ok=0;
    }
    alula_sparse_free(m);
    if (!ok) fputs("Sparse cancellation regression failed.\n",stderr);
    return ok?0:1;
}

static int check_sparse_reuse(void)
{
    const AlulaSparseEntry entries[]={{0,0},{1,1}};
    const double rhs[]={2,4};
    AlulaSparse *m=alula_sparse_create(2,entries,2,NULL);
    if (!m) return 1;
    double *a=alula_sparse_values(m),x[2];
    a[0]=2;a[1]=4;
    int ok=alula_sparse_solve(m,rhs,x);
    a[0]=2.002;
    ok=ok && alula_sparse_solve(m,rhs,x) && fabs(a[0]*x[0]-rhs[0])<=4e-10 && x[1]==1;
    AlulaSparseStats refined=alula_sparse_stats(m);
    a[0]=8;
    ok=ok && alula_sparse_solve(m,rhs,x) && x[0]==.25 && x[1]==1;
    AlulaSparseStats refactored=alula_sparse_stats(m);
    ok=ok && refined.factor_reuses==1 && refined.refactor_count==0 &&
        refactored.factor_reuses==1 && refactored.refactor_count==1;
    alula_sparse_free(m);
    if (!ok) fputs("Sparse factor reuse regression failed.\n",stderr);
    return ok?0:1;
}

int main(int argc,char **argv)
{
    if (argc==2 && !strcmp(argv[1],"--check-sparse"))
        return check_sparse_cancellation() || check_sparse_reuse();
    AlulaBodyInput in={0};
    int mode, repeat, jacobian, tubes[7];
    AlulaBodyCut cuts[6];
    if (scanf("%d%d%d%d%d",&in.nx,&in.bodies,&mode,&repeat,&jacobian)!=5 ||
        in.bodies<1 || in.bodies>6 || repeat<1 ||
        scanf("%lf%lf%lf%lf%lf%lf%lf",&in.gamma,&in.mach,&in.alpha,&in.length_scale,
            &in.correction,&in.center.x,&in.center.y)!=7 ||
        scanf("%d%d%d%d%d%d%lf%lf%lf",&in.upwind,&in.ismom,&in.positive_simple,
            &in.displaced,&in.independent_wakes,&in.bank_outlet,&in.mucon,&in.mcrit,&in.epsilon)!=9) return 2;
    for (int g=0;g<=in.bodies;++g) if (scanf("%d",tubes+g)!=1) return 2;
    for (int g=0;g<in.bodies;++g) {
        AlulaBodyCut *c=cuts+g;
        if (scanf("%d%d%d%d%d%d",&c->leading,&c->trailing,&c->rounded,
            &c->through,&c->lower_tube,&c->upper_tube)!=6) return 2;
    }
    in.mode=(AlulaBodyMode)mode;
    in.tubes=tubes; in.cuts=cuts;
    AlulaBody *b=alula_body_create(&in);
    if (!b) return 2;
    if (jacobian==3) {
        int result=run_coordinates(b,&in);
        alula_body_free(b);
        return result;
    }
    AlulaBodySize size=alula_body_size(b);
    AlulaPoint *nodes=malloc((size_t)size.nodes*sizeof *nodes);
    size_t width_values=(size_t)size.density+size.tubes+size.wake_widths+size.rows;
    double *data=malloc(width_values*sizeof *data);
    double *directions=jacobian ? malloc(width_values*sizeof *directions) : NULL;
    AlulaPoint *node_directions=jacobian ? malloc((size_t)size.nodes*sizeof *node_directions) : NULL;
    AlulaBodyState state={0}, direction={0};
    int result=2;
    if (!nodes || !data || (jacobian && (!directions || !node_directions))) goto done;
    double *mass=data+size.density,*width=mass+size.tubes,*res=width+size.wake_widths;
    for (int i=0;i<size.nodes;++i)
        if (scanf("%lf%lf",&nodes[i].x,&nodes[i].y)!=2) goto done;
    if (!read_vector(data,size.density) || !read_vector(mass,size.tubes) || !read_vector(width,size.wake_widths)
        || scanf("%lf%lf%lf%lf",&state.circulation,&state.source,&state.doublet_x,&state.doublet_y)!=4) goto done;
    state.nodes=nodes; state.density=data; state.mass=mass; state.wake_width=width;
    if (jacobian==2) { result=run_newton(b,&state); goto done; }
    double *product=NULL;
    if (jacobian) {
        double *dm=directions+size.density,*dw=dm+size.tubes;
        product=dw+size.wake_widths;
        for (int i=0;i<size.nodes;++i)
            if (scanf("%lf%lf",&node_directions[i].x,&node_directions[i].y)!=2) goto done;
        if (!read_vector(directions,size.density) || !read_vector(dm,size.tubes) || !read_vector(dw,size.wake_widths)
            || scanf("%lf%lf%lf%lf",&direction.circulation,&direction.source,&direction.doublet_x,&direction.doublet_y)!=4) goto done;
        direction.nodes=node_directions; direction.density=directions; direction.mass=dm; direction.wake_width=dw;
    }
    AlulaBodyDiagnostics diag;
    AlulaBodyStatus status=ALULA_BODY_INVALID;
    double start=milliseconds();
    for (int i=0;i<repeat;++i) {
        status=jacobian ? alula_body_linearize(b,&state,&direction,res,product,&diag) :
            alula_body_residual(b,&state,res,&diag);
        if (status!=ALULA_BODY_OK) break;
    }
    double time=(milliseconds()-start)/repeat;
    printf("%d %d %.17g\n",status,size.rows,time);
    if (status==ALULA_BODY_OK) {
        printf("%.17g",diag.max_mach);
        for (int k=0;k<ALULA_BODY_FAMILIES;++k) printf(" %.17g",diag.residual[k]);
        putchar('\n');
        for (int i=0;i<size.rows;++i) printf("%.17g%c",res[i],i==size.rows-1?'\n':' ');
        if (jacobian) for (int i=0;i<size.rows;++i) printf("%.17g%c",product[i],i==size.rows-1?'\n':' ');
    }
    result=0;
done:
    free(data); free(nodes); free(directions); free(node_directions);
    alula_body_free(b);
    return result;
}
