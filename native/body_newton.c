// SPDX-License-Identifier: GPL-2.0-or-later
#include "body_newton.h"
#include "body_internal.h"

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct { int row,column; } Dependency;
typedef struct { int slot,virtual_row; } Scatter;
typedef struct { int physical,derivative; } Seed;

struct AlulaBodyNewton {
    AlulaBody *body;
    AlulaBodyPhysicalLayout physical;
    AlulaBodyDecode decode;
    int (*rebase)(void *,double *);
    void *context;
    AlulaSparse *matrix;
    int colors, *order, *scatter_start, *seed_start;
    Scatter *scatter;
    Seed *seed;
    double *physical_values,*physical_direction,*derivatives,*edge_product,*work;
    AlulaPoint *nodes,*node_direction;
};

typedef struct {
    AlulaBodyNewton *solver;
    const AlulaBodyChart *map;
    int row,*seen,*start;
    size_t count,capacity;
    Dependency *edges;
    int failed;
} Builder;

AlulaBodyPhysicalLayout alula_body_physical_layout(AlulaBodySize size)
{
    AlulaBodyPhysicalLayout p;
    p.density=2*size.nodes;
    p.mass=p.density+size.density;
    p.strength=p.mass+size.tubes;
    p.wake=p.strength+4;
    p.size=p.wake+size.wake_widths;
    return p;
}

static void dependency(Builder *b,int physical)
{
    if (b->failed) return;
    for (int k=b->map->starts[physical];k<b->map->starts[physical+1];++k) {
        int col=b->map->columns[k];
        if (b->seen[col]==b->row) continue;
        b->seen[col]=b->row;
        if (b->count==b->capacity) {
            if (b->capacity>INT_MAX/2) { b->failed=1; return; }
            size_t cap=b->capacity ? 2*b->capacity : 1024;
            Dependency *edges=realloc(b->edges,cap*sizeof *edges);
            if (!edges) { b->failed=1; return; }
            b->edges=edges; b->capacity=cap;
        }
        b->edges[b->count++]=(Dependency){b->row,col};
    }
}

static void dep_node(Builder *b,int g,int i,int j)
{
    const AlulaBody *body=b->solver->body;
    int k=2*(body->node[g]+i*(body->tubes[g]+1)+j);
    dependency(b,k);
    dependency(b,k+1);
}

static void dep_section(Builder *b,int g,int i,int j)
{
    const AlulaBody *body=b->solver->body;
    dependency(b,b->solver->physical.density+index_section(body,g,i,j));
    dependency(b,b->solver->physical.mass+body->column[g]+j);
    for (int k=i;k<=i+1;++k) {
        dep_node(b,g,k,j); dep_node(b,g,k,j+1);
    }
}

static void dep_cell(Builder *b,int g,int i,int j)
{
    const AlulaBody *body=b->solver->body;
    int start=body->in.upwind ? (i>=3 ? i-3 : 0) : i-1;
    for (int k=start;k<=i;++k) dep_section(b,g,k,j);
}

static void dep_strengths(Builder *b)
{
    int start=b->solver->physical.strength;
    for (int k=0;k<4;++k) dependency(b,start+k);
}

static void row_pattern(Builder *b)
{
    const AlulaBody *body=b->solver->body;
    int nx=body->in.nx,top=body->in.bodies,n=body->size.rows;
    /* The three global matching rows are accumulated from virtual edge
     * rows. Coloring the global rows themselves would couple the entire
     * exterior and destroy the bounded stencil count. */
    for (b->row=0;b->row<n-3;++b->row) {
        b->start[b->row]=(int)b->count;
        Row r=body->rows[b->row];
        int i=r.i,g=r.g,j=r.j;
        switch (r.family) {
        case ALULA_BODY_STREAMWISE:
            dep_cell(b,g,i,j); break;
        case ALULA_BODY_INLET_DENSITY:
            dep_section(b,g,0,j); break;
        case ALULA_BODY_INTERNAL_PRESSURE:
            dep_cell(b,g,i,j-1); dep_cell(b,g,i,j); break;
        case ALULA_BODY_FARFIELD_PRESSURE:
            dep_cell(b,g,i,g==top?j-1:0);
            dep_node(b,g,i,j); dep_strengths(b); break;
        case ALULA_BODY_CUT_PRESSURE:
        case ALULA_BODY_TRAILING_KUTTA:
        case ALULA_BODY_LEADING_KUTTA:
            dep_cell(b,g,i,body->tubes[g]-1); dep_cell(b,g+1,i,0); break;
        case ALULA_BODY_WAKE_GAP:
            for (int k=i-1;k<=i+1 && k<=nx;++k) {
                dep_node(b,g,k,body->tubes[g]); dep_node(b,g+1,k,0);
            }
            dependency(b,b->solver->physical.wake+body->wake[g]+i-body->cuts[g].trailing-1);
            break;
        case ALULA_BODY_END_TANGENCY: {
            int start=i==nx ? nx-1 : 0;
            for (int k=start;k<=start+1;++k) {
                dep_node(b,g,k,j);
                if (i==nx && g<top && j==body->tubes[g] && body->in.displaced) dep_node(b,g+1,k,0);
            }
            dep_strengths(b);
            break;
        }
        default: b->failed=1; return;
        }
    }
    for (int side=0;side<2;++side) {
        int g=side?top:0,j=side?body->tubes[top]-1:0;
        for (int i=0;i<nx;++i,++b->row) {
            b->start[b->row]=(int)b->count;
            dep_section(b,g,i,j); dep_strengths(b);
        }
    }
    b->start[b->row]=(int)b->count;
}

static int *coloring(const Builder *b,int n,int *total)
{
    int *start=calloc((size_t)n+1,sizeof *start);
    int *next=malloc((size_t)n*sizeof *next);
    int *incidence=malloc(b->count*sizeof *incidence);
    int *colors=malloc((size_t)n*sizeof *colors);
    int *seen=calloc((size_t)n,sizeof *seen);
    if (!start || !next || !incidence || !colors || !seen) {
        free(colors); colors=NULL; goto done;
    }
    for (size_t k=0;k<b->count;++k) start[b->edges[k].column+1]++;
    for (int i=0;i<n;++i) start[i+1]+=start[i];
    memcpy(next,start,(size_t)n*sizeof *next);
    for (size_t k=0;k<b->count;++k) incidence[next[b->edges[k].column]++]=b->edges[k].row;
    *total=0;
    for (int col=0;col<n;++col) {
        for (int k=start[col];k<start[col+1];++k) {
            int row=incidence[k];
            for (int t=b->start[row];t<b->start[row+1];++t) {
                int neighbor=b->edges[t].column;
                if (neighbor<col) seen[colors[neighbor]]=col+1;
            }
        }
        int c=0;
        while (seen[c]==col+1) ++c;
        colors[col]=c;
        if (c>=*total) *total=c+1;
    }
done:
    free(start); free(next); free(incidence); free(seen);
    return colors;
}

static int build(AlulaBodyNewton *s,const AlulaBodyChart *map)
{
    const AlulaBody *body=s->body;
    int n=body->size.rows, match=n-3, vrows=match+2*body->in.nx;
    Builder b={.solver=s,.map=map};
    b.start=malloc(((size_t)vrows+1)*sizeof *b.start);
    b.seen=malloc((size_t)n*sizeof *b.seen);
    if (!b.start || !b.seen) {free(b.start);free(b.seen);return 0;}
    for (int i=0;i<n;++i) b.seen[i]=-1;
    row_pattern(&b);
    free(b.seen);
    int ok=0,*colors=NULL,*slots=NULL,*next=NULL;
    AlulaSparseEntry *entries=NULL;
    if (b.failed || b.count==0) goto done;
    colors=coloring(&b,n,&s->colors);
    if (!colors) goto done;
    int nc=s->colors;
    int derivatives=map->starts[s->physical.size];
    s->scatter_start=calloc((size_t)nc+1,sizeof *s->scatter_start);
    s->seed_start=calloc((size_t)nc+1,sizeof *s->seed_start);
    s->seed=malloc((size_t)derivatives*sizeof *s->seed);
    size_t count=(size_t)b.start[match]+3*(b.count-b.start[match]);
    if (count>INT_MAX) goto done;
    s->scatter=malloc(count*sizeof *s->scatter);
    entries=malloc(count*sizeof *entries);
    slots=malloc(count*sizeof *slots);
    next=malloc((size_t)nc*sizeof *next);
    if (!s->scatter_start || !s->seed_start || !s->seed || !s->scatter || !entries || !slots || !next) goto done;
    size_t k=0;
    for (size_t j=0;j<b.count;++j) {
        Dependency dep=b.edges[j];
        int c=colors[dep.column];
        if (dep.row<match) {
            entries[k++]=(AlulaSparseEntry){s->order[dep.row],dep.column};
            s->scatter_start[c+1]++;
        } else for (int q=0;q<3;++q) {
            entries[k++]=(AlulaSparseEntry){s->order[match+q],dep.column};
            s->scatter_start[c+1]++;
        }
    }
    s->matrix=alula_sparse_create(n,entries,count,slots);
    if (!s->matrix) goto done;
    for (int c=0;c<nc;++c) s->scatter_start[c+1]+=s->scatter_start[c];
    memcpy(next,s->scatter_start,(size_t)nc*sizeof *next);
    k=0;
    for (size_t j=0;j<b.count;++j) {
        Dependency dep=b.edges[j];
        int c=colors[dep.column];
        if (dep.row<match) s->scatter[next[c]++]=(Scatter){slots[k++],dep.row};
        else for (int q=0;q<3;++q)
            s->scatter[next[c]++]=(Scatter){slots[k++],n+6*(dep.row-match)+2*q};
    }
    for (int k=0;k<derivatives;++k) s->seed_start[colors[map->columns[k]]+1]++;
    for (int c=0;c<nc;++c) s->seed_start[c+1]+=s->seed_start[c];
    memcpy(next,s->seed_start,(size_t)nc*sizeof *next);
    for (int i=0;i<s->physical.size;++i) for (int k=map->starts[i];k<map->starts[i+1];++k)
        s->seed[next[colors[map->columns[k]]]++]=(Seed){i,k};
    ok=1;
done:
    free(b.start);free(b.edges);free(colors);free(slots);free(entries);free(next);
    return ok;
}

AlulaBodyNewton *alula_body_newton_create(AlulaBody *body,const AlulaBodyChart *chart)
{
    if (!body || !chart || !chart->decode || !chart->starts || !chart->columns || chart->starts[0]!=0)
        return NULL;
    AlulaBodyPhysicalLayout p=alula_body_physical_layout(body->size);
    int n=body->size.rows;
    for (int i=0;i<p.size;++i) {
        int start=chart->starts[i],end=chart->starts[i+1];
        if (end<start) return NULL;
        int previous=-1;
        for (int k=start;k<end;++k) {
            int col=chart->columns[k];
            if (col<=previous || col>=n) return NULL;
            previous=col;
        }
    }
    int nd=chart->starts[p.size];
    if (nd==0) return NULL;
    AlulaBodyNewton *s=calloc(1,sizeof *s);
    if (!s) return NULL;
    s->body=body; s->physical=p; s->decode=chart->decode; s->context=chart->context; s->rebase=chart->rebase;
    s->order=malloc((size_t)n*sizeof *s->order);
    s->physical_values=malloc((size_t)p.size*sizeof *s->physical_values);
    s->physical_direction=calloc((size_t)p.size,sizeof *s->physical_direction);
    s->derivatives=malloc((size_t)nd*sizeof *s->derivatives);
    s->nodes=malloc((size_t)body->size.nodes*sizeof *s->nodes);
    s->node_direction=malloc((size_t)body->size.nodes*sizeof *s->node_direction);
    s->work=malloc(((size_t)6*n+12*(size_t)body->in.nx)*sizeof *s->work);
    unsigned char *seen=calloc((size_t)n,sizeof *seen);
    int ok=s->order && s->physical_values && s->physical_direction && s->derivatives &&
        s->nodes && s->node_direction && s->work && seen;
    if (ok) for (int i=0;i<n;++i) {
        int row=chart->equation_order ? chart->equation_order[i] : i;
        if (row<0 || row>=n || seen[row]) {ok=0;break;}
        seen[row]=1; s->order[i]=row;
    }
    free(seen);
    if (ok && build(s,chart)) {
        s->edge_product=s->work+6*(size_t)n;
        return s;
    }
    alula_body_newton_free(s);
    return NULL;
}

void alula_body_newton_free(AlulaBodyNewton *s)
{
    if (!s) return;
    alula_sparse_free(s->matrix);
    free(s->order); free(s->scatter_start); free(s->seed_start);
    free(s->scatter); free(s->seed); free(s->physical_values); free(s->physical_direction);
    free(s->derivatives); free(s->work); free(s->nodes); free(s->node_direction); free(s);
}

static AlulaBodyState view(const AlulaBodyNewton *s,const double *a,AlulaPoint *nodes)
{
    AlulaBodyPhysicalLayout p=s->physical;
    for (int i=0;i<s->body->size.nodes;++i) nodes[i]=(AlulaPoint){a[2*i],a[2*i+1]};
    return (AlulaBodyState){.nodes=nodes,.density=a+p.density,.mass=a+p.mass,
        .circulation=a[p.strength],.source=a[p.strength+1],
        .doublet_x=a[p.strength+2],.doublet_y=a[p.strength+3],.wake_width=a+p.wake};
}

static int evaluate(AlulaBodyNewton *s,const double *x,double *r,AlulaBodyDiagnostics *diag,int deriv)
{
    for (int i=0;i<s->body->size.rows;++i) if (!isfinite(x[i])) return 0;
    if (!s->decode(s->context,x,s->physical_values,deriv?s->derivatives:NULL)) return 0;
    AlulaBodyState state=view(s,s->physical_values,s->nodes);
    return alula_body_residual(s->body,&state,r,diag)==ALULA_BODY_OK;
}

static int assemble(AlulaBodyNewton *s,const double *x,double *r)
{
    if (!evaluate(s,x,r,NULL,1)) return 0;
    AlulaBodyState state=view(s,s->physical_values,s->nodes);
    int n=s->body->size.rows;
    double *a=alula_sparse_values(s->matrix), *jv=s->work+5*(size_t)n;
    memset(a,0,(size_t)alula_sparse_nonzeros(s->matrix)*sizeof *a);
    for (int color=0;color<s->colors;++color) {
        int start=s->seed_start[color], end=s->seed_start[color+1];
        if (start==end) continue;
        for (int k=start;k<end;++k) {
            Seed q=s->seed[k];
            s->physical_direction[q.physical]+=s->derivatives[q.derivative];
        }
        AlulaBodyState d=view(s,s->physical_direction,s->node_direction);
        AlulaBodyStatus status=alula_body_apply_prepared(s->body,&state,&d,jv,s->edge_product);
        for (int k=start;k<end;++k) s->physical_direction[s->seed[k].physical]=0;
        if (status!=ALULA_BODY_OK) return 0;
        for (int k=s->scatter_start[color];k<s->scatter_start[color+1];++k) {
            Scatter q=s->scatter[k];
            a[q.slot]+=jv[q.virtual_row];
        }
    }
    return 1;
}

int alula_body_newton_product(AlulaBodyNewton *s,const double *x,const double *d,double *product)
{
    if (!s || !x || !d || !product) return 0;
    int n=s->body->size.rows;
    for (int i=0;i<n;++i) if (!isfinite(d[i])) return 0;
    if (!assemble(s,x,s->work)) return 0;
    double *temp=s->work+n;
    alula_sparse_product(s->matrix,d,temp);
    for (int i=0;i<n;++i) if (!isfinite(temp[i])) return 0;
    for (int i=0;i<n;++i) product[i]=temp[s->order[i]];
    return 1;
}

static double max_abs(const double *x,int n)
{
    double m=0;
    for (int i=0;i<n;++i) m=fmax(m,fabs(x[i]));
    return m;
}

static int final_convex(const AlulaBodyNewton *s)
{
    const AlulaBody *b=s->body;
    for (int g=0;g<=b->in.bodies;++g)
        for (int i=0;i<b->in.nx;++i)
            for (int j=0;j<b->tubes[g];++j) {
                int k=b->node[g]+i*(b->tubes[g]+1)+j,m=k+b->tubes[g]+1;
                AlulaPoint p[]={s->nodes[k],s->nodes[m],s->nodes[m+1],s->nodes[k+1]};
                for (int q=0;q<4;++q)
                    if (!(cross(sub(p[(q+1)%4],p[q]),sub(p[(q+2)%4],p[(q+1)%4]))>0)) return 0;
            }
    return 1;
}

AlulaBodySolveResult alula_body_newton_solve(AlulaBodyNewton *s,double *x,int iterations,double tolerance)
{
    AlulaBodySolveResult result={.status=ALULA_BODY_BAD_STATE,.norm=INFINITY,.max_mach=NAN};
    if (!s || !x || iterations<0 || !positive(tolerance)) return result;
    int n=s->body->size.rows;
    double *r=s->work,*rhs=r+n,*delta=rhs+n,*trial=delta+n,*rt=trial+n;
    AlulaBodyDiagnostics diag;
    if (!evaluate(s,x,r,&diag,0)) return result;
    result.norm=max_abs(r,n);
    result.max_mach=diag.max_mach;
    result.status=ALULA_BODY_ITERATION_LIMIT;
    result.colors=s->colors;
    for (int it=0;it<iterations && result.norm>tolerance;++it) {
        if (!assemble(s,x,r)) {result.status=ALULA_BODY_LINEAR_FAILURE;break;}
        for (int i=0;i<n;++i) rhs[s->order[i]]=-r[i];
        if (!alula_sparse_solve(s->matrix,rhs,delta)) {result.status=ALULA_BODY_LINEAR_FAILURE;break;}
        int accepted=0;
        for (double step=1;step>=0x1p-20;step*=.5) {
            for (int i=0;i<n;++i) trial[i]=x[i]+step*delta[i];
            if (evaluate(s,trial,rt,&diag,0)) {
                double norm=max_abs(rt,n);
                if (norm<=(1-1e-4*step)*result.norm) {
                    if (s->rebase && !s->rebase(s->context,trial)) {
                        result.backtracks++;
                        continue;
                    }
                    memcpy(x,trial,(size_t)n*sizeof *x);
                    result.norm=norm;
                    result.max_mach=diag.max_mach;
                    result.iterations++;
                    accepted=1;
                    break;
                }
            }
            result.backtracks++;
        }
        if (!accepted) {result.status=ALULA_BODY_LINE_SEARCH_FAILURE;break;}
    }
    if (result.norm<=tolerance)
        result.status=final_convex(s) ? ALULA_BODY_SOLVED : ALULA_BODY_INVALID_FINAL_GRID;
    result.linear=alula_sparse_stats(s->matrix);
    return result;
}
