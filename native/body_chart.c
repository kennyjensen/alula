// SPDX-License-Identifier: GPL-2.0-or-later
#include "body_chart.h"
#include "body_internal.h"
#include "contour.h"

#include <stdlib.h>
#include <string.h>

#define NODE_COLUMNS 8

typedef struct {
    AlulaPoint p, slope[NODE_COLUMNS];
    int count, col[NODE_COLUMNS];
} Node;

typedef struct {
    int group, station, bank, cut, wall, column, canonical;
} NodeSpec;

typedef struct {
    int stations, material;
    AlulaContourSpline curve;
    double *fraction[2], *thickness[2], *width, *arc_lengths;
    double stagnation;
} Surface;

struct AlulaBodyCoordinates {
    AlulaBody *body;
    AlulaBodyPhysicalLayout physical;
    Surface surfaces[MAX_BODIES];
    int walls_only, body_stations, te_center;
    int *starts,*columns,*equations;
    int offset, strength[4], stag[MAX_BODIES], capture[MAX_BODIES];
    NodeSpec *spec;
    AlulaPoint *origin,*normal,*next_normal;
    Node *raw,*decoded,*wall,*wall_anchor;
    double *mass_fraction,*eta;
    double level[MAX_BODIES+2], scale, span;
};

static int id(const AlulaBodyCoordinates *c,int g,int i,int j)
{
    return c->body->node[g]+i*(c->body->tubes[g]+1)+j;
}

static int wall_id(const AlulaBodyCoordinates *c,int b,int i,int side)
{
    return side ? id(c,b+1,i,0) : id(c,b,i,c->body->tubes[b]);
}

static inline int on_body(const AlulaBody *b,int body,int i)
{
    return body>=0 && i>=b->cuts[body].leading && i<=b->cuts[body].trailing;
}

static void term(Node *n,int col,AlulaPoint p)
{
    if (col<0) return;
    for(int k=0;k<n->count;++k) if(n->col[k]==col) {
        n->slope[k]=add(n->slope[k],p); return;
    }
    /* At most three free stencil positions and two adjacent body parameters
     * contribute in any node chart. */
    n->col[n->count]=col;
    n->slope[n->count]=p;
    n->count++;
}

static void combine(Node *a,const Node *b,double factor)
{
    a->p=add(a->p,scale(b->p,factor));
    for(int k=0;k<b->count;++k) term(a,b->col[k],scale(b->slope[k],factor));
}

static int normal(AlulaPoint v,AlulaPoint *t,double *length)
{
    *length=norm(v);
    if (!positive(*length)) return 0;
    *t=scale(v,1/ *length);
    return 1;
}

static AlulaPoint rotate(AlulaPoint p) { return (AlulaPoint){-p.y,p.x}; }
static AlulaPoint unit_derivative(AlulaPoint v,AlulaPoint t,double length)
{
    return scale(sub(v,scale(t,dot(t,v))),1/length);
}

static void surface_free(Surface *s)
{
    alula_contour_spline_free(&s->curve);free(s->fraction[0]);free(s->fraction[1]);
    free(s->thickness[0]);free(s->thickness[1]);free(s->width);free(s->arc_lengths);
}

static int surface_create(Surface *s,const AlulaBodySurface *in,AlulaBodyCut cut,int nx,int displaced)
{
    int n=in->count, stations=cut.trailing-cut.leading+1;
    if(n<4 || n>4097 || !in->surface || !in->lower_fraction || !in->upper_fraction ||
        (displaced && (!in->lower_thickness || !in->upper_thickness || !in->wake_width))) return 0;
    s->stations=stations;s->stagnation=in->stagnation;s->material=in->material_wake;
    if (!alula_contour_spline_create(&s->curve,in->surface,n) ||
        !(s->stagnation>0 && s->stagnation<s->curve.arc[n-1])) return 0;
    if (displaced) {
        s->arc_lengths=malloc((size_t)n*sizeof *s->arc_lengths);
        if (!s->arc_lengths) return 0;
        alula_contour_arc_build(&s->curve,s->arc_lengths);
    }

    for(int side=0;side<2;++side) {
        const double *f=side?in->upper_fraction:in->lower_fraction;
        const double *d=side?in->upper_thickness:in->lower_thickness;
        if(f[0]!=0 || f[stations-1]!=1) return 0;
        s->fraction[side]=malloc((size_t)stations*sizeof *f);
        s->thickness[side]=calloc((size_t)stations,sizeof *d);
        if(!s->fraction[side] || !s->thickness[side]) return 0;
        for(int k=0;k<stations;++k) {
            if(k && !(f[k]>f[k-1])) return 0;
            if(displaced && (!isfinite(d[k]) || d[k]<0)) return 0;
            s->fraction[side][k]=f[k];
            if(displaced) s->thickness[side][k]=d[k];
        }
    }
    if(displaced && in->lower_thickness[0]!=in->upper_thickness[0]) return 0;
    if(displaced) {
        int count=nx-cut.trailing;
        s->width=malloc((size_t)count*sizeof *s->width);
        if(!s->width) return 0;
        for(int i=0;i<count;++i) {
            if (!isfinite(in->wake_width[i]) || in->wake_width[i]<0) return 0;
            s->width[i]=in->wake_width[i];
        }
    }
    return 1;
}

static int wall_node(const AlulaBodyCoordinates *c,int b,int side,int station,double stag,
    Node *raw,Node *displaced,int deriv)
{
    const Surface *s=c->surfaces+b;
    double end=s->curve.arc[s->curve.count-1],f=s->fraction[side][station];
    double arc=side ? (1-f)*stag : stag+f*(end-stag);
    if(f==1) arc=side?0:end;
    AlulaPoint p,t,dd;
    alula_contour_spline_value(&s->curve,arc,&p,&t,&dd);
    *raw=(Node){.p=p};
    double slope=(1-f)*end;
    if(deriv) term(raw,c->stag[b],scale(t,slope));
    *displaced=*raw;
    if(!c->body->in.displaced) return 1;
    AlulaPoint tangent;
    double length;
    if(!normal(t,&tangent,&length)) return 0;
    double thickness=s->thickness[side][station];
    displaced->p=add(p,scale(rotate(tangent),-thickness));
    if(deriv) term(displaced,c->stag[b],scale(rotate(unit_derivative(dd,tangent,length)),-thickness*slope));
    return 1;
}

static void build_layout(AlulaBodyCoordinates *c,int primary)
{
    AlulaBody *b=c->body;
    int nx=b->in.nx, top=b->in.bodies, next=b->size.density;
    for(int g=0;g<=top;++g) for(int i=0;i<=nx;++i) for(int j=0;j<=b->tubes[g];++j) {
        int k=id(c,g,i,j);
        int cut=j==0 && g>0 ? g-1 : j==b->tubes[g] && g<top ? g : -1;
        NodeSpec v={g,i,j,cut,on_body(b,cut,i),-1,k};
        if(!v.wall) {
            if(cut>=0 && g==cut+1 && !(b->in.independent_wakes && i>b->cuts[cut].trailing)) {
                v.canonical=wall_id(c,cut,i,0);
                v.column=c->spec[v.canonical].column;
            } else v.column=next++;
        }
        c->spec[k]=v;
    }
    c->offset=next;
    c->strength[0]=next++;
    for(int g=0;g<top;++g) c->stag[g]=b->cuts[g].rounded ? next++ : -1;
    for(int g=0;g<top;++g) c->capture[g]=g!=primary ? next++ : -1;
    for(int k=1;k<4;++k) c->strength[k]=next++;
    if(next!=b->size.rows) c->offset=-1;

    for(int k=0;k<b->size.rows;++k) {
        Row r=b->rows[k];
        int i=r.i,g=r.g,j=r.j;
        switch(r.family) {
        case ALULA_BODY_STREAMWISE: case ALULA_BODY_INLET_DENSITY:
            c->equations[k]=index_section(b,g,i,j);break;
        case ALULA_BODY_INTERNAL_PRESSURE:case ALULA_BODY_FARFIELD_PRESSURE:case ALULA_BODY_END_TANGENCY:
            c->equations[k]=c->spec[id(c,g,i,j)].column;break;
        case ALULA_BODY_CUT_PRESSURE:
            c->equations[k]=c->spec[wall_id(c,g,i,0)].column;break;
        case ALULA_BODY_WAKE_GAP:
            c->equations[k]=c->spec[wall_id(c,g,i,1)].column;break;
        case ALULA_BODY_TRAILING_KUTTA:
            c->equations[k]=g==primary?c->strength[0]:c->capture[g];break;
        case ALULA_BODY_LEADING_KUTTA:
            c->equations[k]=c->stag[g];break;
        case ALULA_BODY_FARFIELD_MATCH:
            c->equations[k]=c->strength[j+1];break;
        default:c->equations[k]=-1;break;
        }
    }
}

static Node shift(const AlulaBodyCoordinates *c,int b,int i,int side)
{
    if(b<0 || i>c->body->cuts[b].trailing || c->stag[b]<0) return (Node){0};
    AlulaBodyCut cut=c->body->cuts[b];
    int station=i<cut.leading?cut.leading:i;
    int k=wall_id(c,b,station,side);
    Node d=c->raw[k];
    d.p=sub(d.p,c->wall_anchor[k].p);
    if(i<cut.leading) {
        double w=(double)i/cut.leading;
        d.p=scale(d.p,w);
        for(int j=0;j<d.count;++j) d.slope[j]=scale(d.slope[j],w);
    }
    return d;
}

static int coordinates(const AlulaBodyCoordinates *c,const double *x,int derivatives,int edges_only)
{
    AlulaBody *b=c->body;
    int nx=b->in.nx,top=b->in.bodies;
    double stag[MAX_BODIES];
    for(int g=0;g<top;++g) {
        const Surface *s=c->surfaces+g;
        double length=s->curve.arc[s->curve.count-1];
        stag[g]=s->stagnation+(c->stag[g]>=0 ? length*x[c->stag[g]] : 0);
        if(!(stag[g]>0 && stag[g]<length)) return 0;
        for(int side=0;side<2;++side)
            for(int i=b->cuts[g].leading;i<=b->cuts[g].trailing;++i) {
                int k=wall_id(c,g,i,side);
                if(!wall_node(c,g,side,i-b->cuts[g].leading,stag[g],c->raw+k,c->wall+k,derivatives)) return 0;
            }
    }
    for (int g=0;g<=top;++g) for (int i=0;i<=nx;++i) {
        /* All interior nodes on the same crossline interpolate the same
         * two material-wall shifts. Compute each shift once per passage. */
        Node lower={0},upper={0};
        if (!c->walls_only) {
            lower=shift(c,g-1,i,1);
            upper=shift(c,g<top?g:-1,i,0);
        }
        int start=id(c,g,i,0);
        int tubes=b->tubes[g];
        for (int j=0;j<=tubes;++j) {
            if (edges_only && j==2 && tubes>3) j=tubes-1;
            int k=start+j;
            NodeSpec v=c->spec[k];
            if (v.wall) {c->decoded[k]=c->wall[k];continue;}
            if (v.canonical!=k) {
                c->decoded[k]=c->raw[k]=c->raw[v.canonical];
                continue;
            }
            AlulaPoint direction=scale(c->normal[v.column-b->size.density],c->scale);
            const Node *cut_shift=v.cut==g?&upper:&lower;
            if (!derivatives) {
                AlulaPoint p=c->origin[k];
                if (!c->walls_only) {
                    if (v.cut>=0) p=add(p,cut_shift->p);
                    else {
                        p=add(p,scale(lower.p,1-c->eta[k]));
                        p=add(p,scale(upper.p,c->eta[k]));
                    }
                }
                c->decoded[k].p=c->raw[k].p=add(p,scale(direction,x[v.column]));
                c->decoded[k].count=c->raw[k].count=0;
            } else {
                Node p={.p=c->origin[k]};
                if (!c->walls_only) {
                    if (v.cut>=0) combine(&p,cut_shift,1);
                    else {
                        combine(&p,&lower,1-c->eta[k]);
                        combine(&p,&upper,c->eta[k]);
                    }
                }
                p.p=add(p.p,scale(direction,x[v.column]));
                term(&p,v.column,direction);
                c->decoded[k]=c->raw[k]=p;
            }
        }
    }
    if(!b->in.displaced) return 1;
    for(int g=0;g<top;++g) {
        int te=b->cuts[g].trailing,a=wall_id(c,g,te,0),z=wall_id(c,g,te,1);
        const Surface *s=c->surfaces+g;
        AlulaPoint ref=s->material ? sub(s->curve.points[0],s->curve.points[s->curve.count-1]) : (AlulaPoint){0};
        Node center={0};
        combine(&center,c->wall+a,.5);combine(&center,c->wall+z,.5);
        AlulaPoint translation=sub(center.p,mean(c->raw[a].p,c->raw[z].p));
        for(int i=te+1;i<=nx;++i) {
            int k0=wall_id(c,g,i,0),k1=wall_id(c,g,i,1);
            if(b->in.independent_wakes) {
                if(c->te_center) {
                    c->decoded[k0].p=add(c->raw[k0].p,translation);
                    c->decoded[k1].p=add(c->raw[k1].p,translation);
                }
                continue;
            }
            const Node *left=i==te+1 ? &center : c->raw+wall_id(c,g,i-1,0);
            const Node *right=c->raw+wall_id(c,g,i<nx?i+1:i,0);
            Node edge={0};
            combine(&edge,right,1);combine(&edge,left,-1);
            AlulaPoint t;
            double length;
            if(!normal(edge.p,&t,&length)) return 0;
            double width=s->width[i-te-1],along=dot(ref,t);
            AlulaPoint offset=add(scale(rotate(t),width),scale(t,along));
            for(int side=0;side<2;++side) {
                int k=side?k1:k0;
                double sign=side?.5:-.5;
                Node n=c->raw[k0];
                n.p=add(n.p,scale(offset,sign));
                for(int q=0;q<edge.count;++q) {
                    AlulaPoint dt=unit_derivative(edge.slope[q],t,length);
                    AlulaPoint ds=add(scale(rotate(dt),width),add(scale(t,dot(ref,dt)),scale(dt,along)));
                    term(&n,edge.col[q],scale(ds,sign));
                }
                c->decoded[k]=n;
            }
        }
    }
    return 1;
}

static AlulaPoint normal_point(const AlulaBodyCoordinates *c, NodeSpec v,int i)
{
    AlulaBody *b=c->body;
    if(v.cut<0 || (b->in.independent_wakes && v.station>b->cuts[v.cut].trailing))
        return c->decoded[id(c,v.group,i,v.bank)].p;
    return mean(c->decoded[wall_id(c,v.cut,i,0)].p,c->decoded[wall_id(c,v.cut,i,1)].p);
}

static int normals(AlulaBodyCoordinates *c,AlulaPoint *result)
{
    AlulaBody *b=c->body;
    int top=b->in.bodies,nx=b->in.nx;
    for(int k=0;k<b->size.nodes;++k) {
        NodeSpec v=c->spec[k];
        if(v.wall || v.canonical!=k) continue;
        int i=v.station,left=i>0?i-1:i,right=i<nx?i+1:i;
        if(c->body_stations && v.cut<0) {
            int leading=0,trailing=0;
            for(int g=v.group-1;g<=v.group;++g) if(g>=0 && g<top) {
                leading |= b->cuts[g].leading==i;
                trailing |= b->cuts[g].trailing==i;
            }
            if(leading && trailing) return 0;
            if(leading) right=i;
            if(trailing) left=i;
        }
        double length;
        AlulaPoint t;
        if(!normal(sub(normal_point(c,v,right),normal_point(c,v,left)),&t,&length)) return 0;
        result[v.column-b->size.density]=rotate(t);
    }
    return 1;
}

static int map_create(AlulaBodyCoordinates *c)
{
    int n=c->physical.size;
    c->starts=calloc((size_t)n+1,sizeof *c->starts);
    c->columns=malloc((size_t)n*NODE_COLUMNS*sizeof *c->columns);
    if(!c->starts || !c->columns) return 0;
    int k=0;
    for(int i=0;i<n;++i) {
        int list[NODE_COLUMNS],count=0;
        if(i<c->physical.density) {
            const Node *p=c->decoded+i/2;
            count=p->count;
            for(int j=0;j<count;++j) list[j]=p->col[j];
        } else if(i<c->physical.mass) {list[0]=i-c->physical.density;count=1;}
        else if(i<c->physical.strength) {
            int tube=i-c->physical.mass,g=0;
            while(tube>=c->body->column[g+1]) ++g;
            if(g>0 && c->capture[g-1]>=0) list[count++]=c->capture[g-1];
            if(g<c->body->in.bodies && c->capture[g]>=0) list[count++]=c->capture[g];
        } else if(i<c->physical.wake) {
            list[0]=c->strength[i-c->physical.strength];count=1;
        }
        for(int j=1;j<count;++j) {
            int v=list[j],p=j;
            while(p>0 && list[p-1]>v) {list[p]=list[p-1];--p;}
            list[p]=v;
        }
        for(int j=0;j<count;++j) c->columns[k++]=list[j];
        c->starts[i+1]=k;
    }
    return 1;
}

static double node_slope(const Node *p,int col,int axis)
{
    for(int k=0;k<p->count;++k) if(p->col[k]==col) return axis?p->slope[k].y:p->slope[k].x;
    return 0;
}

static int decode_masses(const AlulaBodyCoordinates *c,const double *x,double *mass,double *derivatives)
{
    int top=c->body->in.bodies;
    double cap[MAX_BODIES+2];
    cap[0]=c->level[0];cap[top+1]=c->level[top+1];
    for(int g=0;g<top;++g)
        cap[g+1]=c->level[g+1]+(c->capture[g]>=0 ? c->span*x[c->capture[g]]:0);
    for(int g=0;g<=top;++g) {
        double amount=cap[g+1]-cap[g];
        if(!positive(amount)) return 0;
        for(int j=c->body->column[g];j<c->body->column[g+1];++j) {
            mass[j]=c->mass_fraction[j]*amount;
            if(derivatives) for(int k=c->starts[c->physical.mass+j];k<c->starts[c->physical.mass+j+1];++k)
                derivatives[k]=(g>0 && c->columns[k]==c->capture[g-1] ? -1 : 1)*c->span*c->mass_fraction[j];
        }
    }
    return 1;
}

int alula_body_coordinates_edge_state(AlulaBodyCoordinates *c,const double *x,AlulaPoint *nodes,
    double *gas,AlulaBodyState *state)
{
    if (!c || !x || !nodes || !gas || !state || !coordinates(c,x,0,1)) return 0;
    const AlulaBody *b=c->body;
    double *mass=gas+b->size.density;
    if (!decode_masses(c,x,mass,NULL)) return 0;
    for (int g=0;g<=b->in.bodies;++g) {
        int tubes=b->tubes[g];
        for (int i=0;i<=b->in.nx;++i) {
            int k=id(c,g,i,0);
            nodes[k]=c->decoded[k].p;
            nodes[k+1]=c->decoded[k+1].p;
            nodes[k+tubes-1]=c->decoded[k+tubes-1].p;
            nodes[k+tubes]=c->decoded[k+tubes].p;
            if (i<b->in.nx) {
                int first=index_section(b,g,i,0),last=first+tubes-1;
                gas[first]=exp(x[first]);
                gas[last]=exp(x[last]);
            }
        }
    }
    *state=(AlulaBodyState){.nodes=nodes,.density=gas,.mass=mass};
    return 1;
}

static int decode(void *context,const double *x,double *physical,double *derivatives)
{
    AlulaBodyCoordinates *c=context;
    AlulaBody *b=c->body;
    AlulaBodyPhysicalLayout p=c->physical;
    int top=b->in.bodies;
    if(!coordinates(c,x,derivatives!=NULL,0)) return 0;
    for(int i=0;i<b->size.nodes;++i) {
        const Node *v=c->decoded+i;
        physical[2*i]=v->p.x;
        physical[2*i+1]=v->p.y;
        if(derivatives) for(int axis=0;axis<2;++axis) {
            int scalar=2*i+axis;
            for(int k=c->starts[scalar];k<c->starts[scalar+1];++k)
                derivatives[k]=node_slope(v,c->columns[k],axis);
        }
    }
    for(int i=0;i<b->size.density;++i) {
        int scalar=p.density+i;
        double rho=exp(x[i]);
        physical[scalar]=rho;
        if(derivatives) derivatives[c->starts[scalar]]=rho;
    }
    if (!decode_masses(c,x,physical+p.mass,derivatives)) return 0;
    for(int k=0;k<4;++k) {
        int scalar=p.strength+k;
        double scale_length=k<2?c->scale:c->scale*c->scale;
        physical[scalar]=scale_length*x[c->strength[k]];
        if(derivatives) derivatives[c->starts[scalar]]=scale_length;
    }
    if(b->in.independent_wakes) for(int g=0;g<top;++g)
        memcpy(physical+p.wake+b->wake[g],c->surfaces[g].width,(size_t)(b->wake[g+1]-b->wake[g])*sizeof *physical);
    return 1;
}

static int rebase_decoded(void *context,double *x)
{
    AlulaBodyCoordinates *c=context;
    AlulaBody *b=c->body;
    if(!normals(c,c->next_normal)) return 0;
    for(int k=0;k<b->size.nodes;++k) {
        NodeSpec v=c->spec[k];
        if(v.wall || v.canonical!=k) continue;
        int normal_index=v.column-b->size.density;
        c->origin[k]=add(c->origin[k],scale(c->normal[normal_index],c->scale*x[v.column]));
        x[v.column]=0;
        c->normal[normal_index]=c->next_normal[normal_index];
    }
    return 1;
}

static int rebase(void *context,double *x)
{
    return coordinates(context,x,0,0) && rebase_decoded(context,x);
}

AlulaBodyCoordinates *alula_body_coordinates_create(AlulaBody *b,const AlulaBodyCoordinateInput *in)
{
    if(!b || !in || !in->raw_nodes || !in->mass_weights || !in->capture_levels || !in->surfaces ||
        in->primary_body<0 || in->primary_body>=b->in.bodies ||
        (in->te_center_wake && !b->in.independent_wakes)) return NULL;
    AlulaBodyCoordinates *c=calloc(1,sizeof *c);
    if(!c) return NULL;
    c->body=b;c->physical=alula_body_physical_layout(b->size);
    c->walls_only=in->walls_only;c->body_stations=in->body_station_normals;c->te_center=in->te_center_wake;
    c->scale=b->in.length_scale;
    int size=b->size.nodes,top=b->in.bodies,nx=b->in.nx;
    c->spec=malloc((size_t)size*sizeof *c->spec);
    c->equations=malloc((size_t)b->size.rows*sizeof *c->equations);
    c->origin=malloc((size_t)size*sizeof *c->origin);
    c->normal=calloc((size_t)size,sizeof *c->normal);
    c->next_normal=calloc((size_t)size,sizeof *c->next_normal);
    c->raw=calloc((size_t)size,sizeof *c->raw);
    c->decoded=calloc((size_t)size,sizeof *c->decoded);
    c->wall=calloc((size_t)size,sizeof *c->wall);
    c->wall_anchor=calloc((size_t)size,sizeof *c->wall_anchor);
    c->mass_fraction=malloc((size_t)b->size.tubes*sizeof *c->mass_fraction);
    c->eta=malloc((size_t)size*sizeof *c->eta);
    double *x=calloc((size_t)b->size.rows,sizeof *x);
    if(!c->spec || !c->equations || !c->origin || !c->normal || !c->next_normal || !c->raw ||
        !c->decoded || !c->wall || !c->wall_anchor || !c->mass_fraction || !c->eta || !x) goto fail;
    build_layout(c,in->primary_body);
    if(c->offset<0) goto fail;
    for(int i=0;i<size;++i) {
        AlulaPoint p=in->raw_nodes[i];
        if(!finite_point(p)) goto fail;
        NodeSpec v=c->spec[i];
        if(v.canonical!=i && norm(sub(p,in->raw_nodes[v.canonical]))>1e-10*c->scale) goto fail;
        c->origin[i]=p;
    }
    for(int g=0;g<top+2;++g) {
        double level=in->capture_levels[g];
        if(!isfinite(level) || (g && !(level>c->level[g-1]))) goto fail;
        c->level[g]=level;
    }
    c->span=c->level[top+1]-c->level[0];
    if(!positive(c->span)) goto fail;
    for(int g=0;g<=top;++g) {
        double sum=0;
        int t=b->tubes[g],offset=b->column[g];
        for(int j=0;j<t;++j) {
            double w=in->mass_weights[offset+j];
            if(!positive(w)) goto fail;
            sum+=w;
        }
        if(!positive(sum)) goto fail;
        double eta=0;
        for(int j=0;j<=t;++j) {
            for(int i=0;i<=nx;++i) c->eta[id(c,g,i,j)]=j==t ? 1 : eta;
            if(j<t) {
                double fraction=in->mass_weights[offset+j]/sum;
                c->mass_fraction[offset+j]=fraction;
                eta+=fraction;
            }
        }
        if(g==top) continue;
        Surface *s=c->surfaces+g;
        if(!surface_create(s,in->surfaces+g,b->cuts[g],nx,b->in.displaced)) goto fail;
        for(int side=0;side<2;++side)
            for(int i=b->cuts[g].leading;i<=b->cuts[g].trailing;++i) {
                int k=wall_id(c,g,i,side);
                Node displaced;
                if(!wall_node(c,g,side,i-b->cuts[g].leading,s->stagnation,c->wall_anchor+k,&displaced,0) ||
                    norm(sub(c->wall_anchor[k].p,c->origin[k]))>1e-10*c->scale) goto fail;
            }
    }
    if(!coordinates(c,x,1,0) || !normals(c,c->normal) || !map_create(c)) goto fail;
    free(x);return c;
fail:
    free(x);alula_body_coordinates_free(c);return NULL;
}

AlulaBodyChart alula_body_coordinates_chart(AlulaBodyCoordinates *c)
{
    return c ? (AlulaBodyChart){.starts=c->starts,.columns=c->columns,.decode=decode,.context=c,
        .equation_order=c->equations,.rebase=rebase,.rebase_decoded=rebase_decoded} : (AlulaBodyChart){0};
}

int alula_body_coordinates_stagnation_column(const AlulaBodyCoordinates *c,int body)
{
    return c && body>=0 && body<c->body->in.bodies ? c->stag[body] : -1;
}

void alula_body_coordinates_seed_circulation(const AlulaBodyCoordinates *c,double *x,double circulation)
{
    if (c && x && isfinite(circulation)) x[c->strength[0]]=circulation/c->scale;
}

int alula_body_coordinates_seed_state(AlulaBodyCoordinates *c,double *x,
    const double *stagnation,const double *passage_mass,const AlulaBodyState *source)
{
    double level=c->level[0];
    int bodies=c->body->in.bodies;
    for (int g=0;g<bodies;++g) {
        level+=passage_mass[g];
        if (c->capture[g]>=0) x[c->capture[g]]=(level-c->level[g+1])/c->span;
        else if (fabs(level-c->level[g+1])>1e-12*c->span) return 0;
        if (c->stag[g]>=0) {
            const Surface *s=c->surfaces+g;
            x[c->stag[g]]=(stagnation[g]-s->stagnation)/s->curve.arc[s->curve.count-1];
        }
    }
    if (fabs(level+passage_mass[bodies]-c->level[bodies+1])>1e-12*c->span) return 0;
    const double scales[]={c->scale,c->scale,c->scale*c->scale,c->scale*c->scale};
    const double coefficients[]={source->circulation,source->source,source->doublet_x,source->doublet_y};
    for (int k=0;k<4;++k) x[c->strength[k]]=coefficients[k]/scales[k];
    return 1;
}

int alula_body_coordinates_seed_nodes(AlulaBodyCoordinates *c,double *x,const AlulaPoint *nodes)
{
    if (!coordinates(c,x,0,0)) return 0;
    for (int k=0;k<c->body->size.nodes;++k) {
        NodeSpec v=c->spec[k];
        if (v.wall || v.canonical!=k) continue;
        c->origin[k]=add(c->origin[k],sub(nodes[k],c->decoded[k].p));
    }
    return rebase(c,x);
}

int alula_body_coordinates_surface(const AlulaBodyCoordinates *c,const double *x,int element,
    int upper,AlulaBodySurfaceView *view)
{
    if (!c || !x || !view || element<0 || element>=c->body->in.bodies || (upper!=0 && upper!=1)) return 0;
    const Surface *s=c->surfaces+element;
    double arc=s->curve.arc[s->curve.count-1];
    double position=s->stagnation+(c->stag[element]>=0?arc*x[c->stag[element]]:0);
    if (!(position>0 && position<arc)) return 0;
    *view=(AlulaBodySurfaceView){&s->curve,s->fraction[upper],position,s->stagnation,s->arc_lengths};
    return 1;
}

size_t alula_body_coordinates_displacement_count(const AlulaBodyCoordinates *c)
{
    if (!c || !c->body->in.displaced) return 0;
    const AlulaBody *b=c->body;
    size_t n=0;
    for (int i=0;i<b->in.bodies;++i)
        n+=(size_t)(2*c->surfaces[i].stations+b->in.nx-b->cuts[i].trailing);
    return n;
}

int alula_body_coordinates_get_displacement(const AlulaBodyCoordinates *c,double *value)
{
    if (!c || !value || !c->body->in.displaced) return 0;
    const AlulaBody *b=c->body;
    for (int i=0;i<b->in.bodies;++i) {
        const Surface *s=c->surfaces+i;
        int nw=b->in.nx-b->cuts[i].trailing;
        for (int side=0;side<2;++side) {
            memcpy(value,s->thickness[side],(size_t)s->stations*sizeof *value);
            value+=s->stations;
        }
        memcpy(value,s->width,(size_t)nw*sizeof *value);
        value+=nw;
    }
    return 1;
}

int alula_body_coordinates_set_displacement(AlulaBodyCoordinates *c,const double *values)
{
    if (!c || !values || !c->body->in.displaced) return 0;
    const AlulaBody *b=c->body;
    const double *check=values;
    for (int i=0;i<b->in.bodies;++i) {
        int ns=c->surfaces[i].stations,n=2*ns+b->in.nx-b->cuts[i].trailing;
        if (check[0]!=check[ns]) return 0;
        for (int k=0;k<n;++k) if (!(check[k]>=0) || !isfinite(check[k])) return 0;
        check+=n;
    }
    for (int i=0;i<b->in.bodies;++i) {
        Surface *s=c->surfaces+i;
        int nw=b->in.nx-b->cuts[i].trailing;
        for (int side=0;side<2;++side) {
            memcpy(s->thickness[side],values,(size_t)s->stations*sizeof *values);
            values+=s->stations;
        }
        memcpy(s->width,values,(size_t)nw*sizeof *values);
        values+=nw;
    }
    return 1;
}

size_t alula_body_coordinates_checkpoint_count(const AlulaBodyCoordinates *c)
{
    return c ? (size_t)c->body->size.nodes+(size_t)(c->offset-c->body->size.density) : 0;
}

int alula_body_coordinates_save(const AlulaBodyCoordinates *c,AlulaPoint *checkpoint)
{
    if (!c || !checkpoint) return 0;
    int n=c->body->size.nodes,normals=c->offset-c->body->size.density;
    memcpy(checkpoint,c->origin,(size_t)n*sizeof *checkpoint);
    memcpy(checkpoint+n,c->normal,(size_t)normals*sizeof *checkpoint);
    return 1;
}

int alula_body_coordinates_restore(AlulaBodyCoordinates *c,const AlulaPoint *checkpoint)
{
    if (!c || !checkpoint) return 0;
    int n=c->body->size.nodes,normals=c->offset-c->body->size.density;
    memcpy(c->origin,checkpoint,(size_t)n*sizeof *checkpoint);
    memcpy(c->normal,checkpoint+n,(size_t)normals*sizeof *checkpoint);
    return 1;
}

void alula_body_coordinates_free(AlulaBodyCoordinates *c)
{
    if(!c) return;
    for(int g=0;g<MAX_BODIES;++g) surface_free(c->surfaces+g);
    free(c->starts);free(c->columns);free(c->equations);free(c->spec);free(c->origin);
    free(c->normal);free(c->next_normal);free(c->raw);free(c->decoded);free(c->wall);
    free(c->wall_anchor);free(c->mass_fraction);free(c->eta);free(c);
}
