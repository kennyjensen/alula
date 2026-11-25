// SPDX-License-Identifier: GPL-2.0-or-later
#include "airfoil_grid_internal.h"
#include "grid_relax.h"
#include "panel_guide.h"
#include "panel_internal.h"
#include "streamtube_internal.h"

#include <stdlib.h>

static AlulaPoint mix(AlulaPoint a, AlulaPoint b,double w)
{
    return (AlulaPoint){a.x+w*(b.x-a.x),a.y+w*(b.y-a.y)};
}

static double norm(AlulaPoint v) { return hypot(v.x,v.y); }

static int ascending(const void *a,const void *b)
{
    double x=*(const double *)a,y=*(const double *)b;
    return (x>y)-(x<y);
}

static int save_base(AlulaMaterialBase *base,const AlulaContour *source)
{
    if (!source->finite_base) return 1;
    size_t unique=source->count-1;
    base->count=(source->upper_index+unique-source->lower_index)%unique+1;
    base->points=malloc(base->count*sizeof *base->points);
    if (!base->points) return 0;
    for (size_t i=0;i<base->count;++i)
        base->points[i]=source->points[(source->lower_index+i)%unique];
    return 1;
}

static void extend(double *x,int count,double anchor,double boundary,double first)
{
    double extent=fabs(boundary-anchor),sign=copysign(1,boundary-anchor);
    first=fmin(first,extent/count);
    double lower=1,upper=pow(extent/first,1.0/(count-1));
    for (int i=0;i<55;++i) {
        double ratio=.5*(lower+upper),sum=first,step=first;
        for (int j=1;j<count;++j) {step*=ratio;sum+=step;}
        if (sum<extent) lower=ratio;
        else upper=ratio;
    }
    double offset=0,step=first,ratio=.5*(lower+upper);
    for (int j=0;j<count-1;++j) {offset+=step;x[j]=anchor+sign*offset;step*=ratio;}
    x[count-1]=boundary;
}

static int station_index(const double *x,int count,double value,double tolerance)
{
    for (int i=0;i<count;++i) if (fabs(x[i]-value)<=tolerance) return i;
    return -1;
}

static int stations(const AlulaAirfoilGridInput *in,const AlulaPanelBodyGuide *guide,int bodies,
    AlulaBodyCut *cuts,double *x)
{
    int count=0,ni=in->surface_intervals;
    double first=INFINITY,last=-INFINITY,leading_step=0,trailing_step=0;
    for (int g=0;g<bodies;++g) {
        double a=guide[g].stagnation_point.x,b=guide[g].trailing_edge.x;
        if (!(b>a)) return 0;
        double step=.5*(1-cos(ALULA_PI/ni))*(b-a);
        if (a<first) {first=a;leading_step=fmax(step,2*guide[g].seed_distance);}
        if (b>last) {last=b;trailing_step=fmax(step,2*guide[g].seed_distance);}
        for (int i=0;i<=ni;++i)
            x[count++]=i==ni?b:a+(b-a)*.5*(1-cos(ALULA_PI*i/ni));
    }
    if (!(in->inlet_x<first && in->outlet_x>last)) return 0;
    extend(x+count,in->inlet_intervals,first,in->inlet_x,leading_step);
    count+=in->inlet_intervals;
    extend(x+count,in->outlet_intervals,last,in->outlet_x,trailing_step);
    count+=in->outlet_intervals;
    qsort(x,(size_t)count,sizeof *x,ascending);
    double tolerance=1e-9*in->chord;
    int kept=1;
    for (int i=1;i<count;++i)
        if (x[i]-x[kept-1]>tolerance) x[kept++]=x[i];
    for (int g=0;g<bodies;++g) {
        int le=station_index(x,kept,guide[g].stagnation_point.x,tolerance);
        int te=station_index(x,kept,guide[g].trailing_edge.x,tolerance);
        if (le<=0 || te<=le || te>=kept-1) return 0;
        cuts[g]=(AlulaBodyCut){.leading=le,.trailing=te,.rounded=1};
    }
    return kept;
}

static int guide_point(const AlulaPanelGuides *guides,AlulaGuidePath path,
    AlulaPoint edge,double x,AlulaPoint *point)
{
    AlulaPoint seed=path.points[0].position;
    /* The seeds sit a short distance off the body. Connect only that missing
     * interval, leaving all resolved guide intersections untouched. */
    if (fmin(seed.x,edge.x)<x && x<fmax(seed.x,edge.x)) {
        *point=mix(seed,edge,(x-seed.x)/(edge.x-seed.x));
        return 1;
    }
    AlulaTracePoint hit;
    if (alula_panel_guides_sample_x(guides,path,x,&hit)!=ALULA_GUIDE_OK) return 0;
    *point=hit.position;
    return 1;
}

static int boundary_points(const AlulaPanelField *field,const AlulaPanelGuides *guides,
    const AlulaPanelBodyGuide *profiles,AlulaAirfoilGrid *grid,const double *x,int nx,int count,
    AlulaBodySurface *surfaces,double *fractions,AlulaPoint *boundaries)
{
    int rows=nx+1;
    for (int side=0;side<2;++side) {
        AlulaGuidePath path=alula_panel_guides_farfield(guides,side);
        for (int i=0;i<=nx;++i) {
            AlulaTracePoint hit;
            if (alula_panel_guides_sample_x(guides,path,x[i],&hit)!=ALULA_GUIDE_OK) return 0;
            boundaries[(side?2*count+1:0)*rows+i]=hit.position;
        }
    }
    for (int g=0;g<count;++g) {
        const AlulaPanelBodyGuide *profile=profiles+g;
        const AlulaContourSpline *curve=profile->curve;
        AlulaBodyCut cut=grid->cuts[g];
        AlulaPoint *lower=boundaries+(2*g+1)*rows,*upper=lower+rows;
        double end=curve->arc[curve->count-1],stag=profile->stagnation;
        double *f=fractions+(size_t)g*4*rows;
        double *width=f+2*rows;
        surfaces[g]=(AlulaBodySurface){.surface=curve->points,.count=curve->count,.stagnation=stag,
            .lower_fraction=f,.upper_fraction=f,.lower_thickness=f+rows,.upper_thickness=f+rows,
            .wake_width=width,.material_wake=1};
        AlulaPoint a,b,d;
        alula_contour_spline_value(curve,0,NULL,&a,NULL);
        alula_contour_spline_value(curve,end,NULL,&b,NULL);
        d=subtract(b,a);
        double magnitude=norm(d);
        if (!(magnitude>0) || !isfinite(magnitude)) return 0;
        AlulaPoint tangent={d.x/magnitude,d.y/magnitude}, vector=profile->trailing_vector;
        double along=dot(tangent,vector),normal=cross(tangent,vector);
        if (normal<0) return 0;
        grid->base_width[g]=normal;
        for (int i=0;i<=nx;++i) {
            if (i<cut.leading) {
                if (!guide_point(guides,profile->incoming,profile->stagnation_point,x[i],lower+i)) return 0;
                upper[i]=lower[i];
            } else if (i<=cut.trailing) {
                f[i-cut.leading]=(x[i]-x[cut.leading])/(x[cut.trailing]-x[cut.leading]);
                double t=f[i-cut.leading];
                alula_contour_spline_value(curve,(1-t)*stag,upper+i,NULL,NULL);
                alula_contour_spline_value(curve,stag+t*(end-stag),lower+i,NULL,NULL);
            } else {
                AlulaPoint p;
                if (!guide_point(guides,profile->wake,profile->trailing_edge,x[i],&p)) return 0;
                if (normal>0) {
                    AlulaPoint q;
                    if (!alula_panel_field_sample(field,p,&q,NULL)) return 0;
                    magnitude=norm(q);
                    if (!(magnitude>0) || !isfinite(magnitude)) return 0;
                    tangent=(AlulaPoint){q.x/magnitude,q.y/magnitude};
                }
                AlulaPoint gap={along*tangent.x-normal*tangent.y,along*tangent.y+normal*tangent.x};
                lower[i]=(AlulaPoint){p.x-.5*gap.x,p.y-.5*gap.y};
                upper[i]=(AlulaPoint){p.x+.5*gap.x,p.y+.5*gap.y};
                width[i-cut.trailing-1]=normal;
            }
        }
    }
    return 1;
}

static AlulaPoint crossline_direction(const AlulaPanelBodyGuide *profile,AlulaBodyCut cut,
    const double *fraction,const double *x,int i,int upper_surface,AlulaPoint across,double gap)
{
    double chord=profile->trailing_edge.x-profile->stagnation_point.x;
    int station=i<cut.leading?cut.leading:i>cut.trailing?cut.trailing:i;
    double t=fraction[station-cut.leading],end=profile->curve->arc[profile->curve->count-1];
    double parameter=upper_surface?(1-t)*profile->stagnation:profile->stagnation+t*(end-profile->stagnation);
    AlulaPoint arc;
    alula_contour_spline_value(profile->curve,parameter,NULL,&arc,NULL);
    double speed=norm(arc), sign=upper_surface?1:-1;
    if (!(speed>0)) return across;
    AlulaPoint normal={sign*gap*arc.y/speed,-sign*gap*arc.x/speed};
    double weight=exp(-fabs(x[station]-x[i])/(.12*chord));
    return mix(across,normal,weight);
}

static void fill_passages(AlulaAirfoilGrid *grid,const AlulaPanelBodyGuide *profiles,
    const AlulaPoint *boundaries,const double *fractions,const double *x,double *weights,int nx,int count)
{
    size_t offset=0;
    for (int g=0;g<=count;++g) {
        int n=grid->tubes[g],rows=nx+1;
        const AlulaPoint *lower=boundaries+2*g*rows,*upper=lower+rows;
        double growth=pow(24,1.0/fmax(1,n-1)),sum=0;
        for (int j=0;j<n;++j) {
            weights[j]=g==0?pow(growth,n-1-j):g==count?pow(growth,j):sin(ALULA_PI*(j+.5)/n);
            sum+=weights[j];
        }
        double lower_decay=g==0?0:.045*(profiles[g-1].trailing_edge.x-profiles[g-1].stagnation_point.x);
        double upper_decay=g==count?0:.045*(profiles[g].trailing_edge.x-profiles[g].stagnation_point.x);
        for (int i=0;i<=nx;++i) {
            AlulaPoint across=subtract(upper[i],lower[i]);
            double gap=norm(across);
            AlulaPoint dl={0},du={0};
            if (g>0) dl=subtract(crossline_direction(profiles+g-1,grid->cuts[g-1],
                fractions+(size_t)(g-1)*4*rows,x,i,1,across,gap),across);
            if (g<count) du=subtract(crossline_direction(profiles+g,grid->cuts[g],
                fractions+(size_t)g*4*rows,x,i,0,across,gap),across);
            double eta=0;
            for (int j=0;j<=n;++j) {
                AlulaPoint p=j==n?upper[i]:mix(lower[i],upper[i],eta);
                double from_lower=lower_decay>0?eta*pow(1-eta,1+gap/lower_decay):0;
                double from_upper=upper_decay>0?(1-eta)*pow(eta,1+gap/upper_decay):0;
                grid->nodes[offset+i*(n+1)+j]=(AlulaPoint){
                    p.x+from_lower*dl.x-from_upper*du.x,p.y+from_lower*dl.y-from_upper*du.y};
                if (j<n) eta+=weights[j]/sum;
            }
        }
        offset+=(size_t)rows*(n+1);
        weights+=n;
    }
}

static void check_quality(const AlulaAirfoilGrid *grid,AlulaAirfoilGridReport *report)
{
    int nx=report->nx, offset=0;
    report->minimum_corner_sine=1;
    for (int g=0;g<=report->elements;++g) {
        int width=grid->tubes[g]+1;
        for (int i=0;i<nx;++i) for (int j=0;j<width-1;++j) {
            int k=offset+i*width+j,m=k+width;
            AlulaPoint p[]={grid->nodes[k],grid->nodes[m],grid->nodes[m+1],grid->nodes[k+1]};
            int valid=1;
            for (int q=0;q<4;++q) {
                AlulaPoint a=subtract(p[(q+1)%4],p[q]),b=subtract(p[(q+2)%4],p[(q+1)%4]);
                double denominator=norm(a)*norm(b);
                double sine=denominator>0?cross(a,b)/denominator:-1;
                if (!(sine>0)) valid=0;
                report->minimum_corner_sine=fmin(report->minimum_corner_sine,sine);
            }
            if (!valid && report->invalid_cells++==0) {
                report->first_group=g;report->first_interval=i;report->first_tube=j;
            }
        }
        offset+=(nx+1)*width;
    }
}

static double seed_density(double flux,double specific,double stag,double sonic,double capacity)
{
    /* The isentropic branch folds at choking. Give an over-capacity
     * provisional section a positive enthalpy and let Newton update its
     * area and density together. The normal Euler residual still applies. */
    if (flux>=capacity) return flux/(.995*sonic);
    double low=0,high=sonic,q=fmin(flux,.9*sonic);
    for (int i=0;i<48;++i) {
        double temp=stag-specific*q*q,rho=pow(temp,2.5);
        double error=rho*q-flux;
        if (fabs(error)<=2e-14*fmax(1,flux)) return rho;
        if (error<0) low=q;
        else high=q;
        double derivative=rho*(1-5*specific*q*q/temp);
        double next=q-error/derivative;
        q=next>low && next<high?next:.5*(low+high);
    }
    return pow(stag-specific*q*q,2.5);
}

static int seed_gas(AlulaAirfoilGrid *g,int nx,int elements,double mach)
{
    AlulaBodyState state;
    if (!alula_airfoil_grid_state(g,&state)) return 0;
    double specific=.2*mach*mach,stag=1+specific;
    double sonic=sqrt(stag/(1.2*mach*mach));
    double capacity=sonic*pow(stag-specific*sonic*sonic,2.5);
    int offset=0,tube=0,total=alula_body_size(g->body).tubes;
    for (int k=0;k<=elements;++k) {
        int width=g->tubes[k]+1;
        for (int i=0;i<nx;++i) for (int j=0;j<width-1;++j) {
            int id=offset+i*width+j,next=id+width;
            const AlulaPoint *p=state.nodes;
            AlulaPoint d=mix(subtract(p[next],p[id]),subtract(p[next+1],p[id+1]),.5);
            AlulaPoint s=mix(subtract(p[id+1],p[id]),subtract(p[next+1],p[next]),.5);
            double area=cross(d,s)/norm(d);
            if (!(area>0) || !isfinite(area)) return 0;
            double rho=seed_density(state.mass[tube+j]/area,specific,stag,sonic,capacity);
            if (!(rho>0) || !isfinite(rho)) return 0;
            g->packed[i*total+tube+j]=log(rho);
        }
        offset+=(nx+1)*width;
        tube+=width-1;
    }
    return 1;
}

/* The selected Euler residual is the final admissibility check. This
 * fallback ports the subsonic uniform-density pressure bound used by the
 * web initializer to avoid starting from a nonpositive interface pressure. */
static int pressure_domain_seed(AlulaAirfoilGrid *g,int nx,int elements,double mach)
{
    const double gamma=1.4,gm=gamma-1,h0=1/(gm*mach*mach)+.5;
    const double k=gm/gamma, sonic_sq=2*gm*h0/(gamma+1), correction=.1;
    const double rho_total=pow(1+.5*gm*mach*mach,1/gm);
    double needed=rho_total*rho_total;
    AlulaBodyState state;
    if (!alula_airfoil_grid_state(g,&state)) return 0;
    int offset=0,tube=0;
    for (int p=0;p<=elements;++p) {
        int width=g->tubes[p]+1;
        for (int i=1;i<nx;++i) for (int j=0;j<width-1;++j) {
            int idx=offset+(i-1)*width+j;
            AlulaCellInput cell={.positive_simple=1};
            for (int m=0;m<3;++m) {
                cell.lower[m]=state.nodes[idx+m*width];
                cell.upper[m]=state.nodes[idx+m*width+1];
            }
            AlulaCell scratch;
            if (!alula_cell_geometry_only(&scratch,&cell)) return 0;
            const AlulaCellGeometry *c=&scratch.value.geometry;
            double mass=state.mass[tube+j],a=mass/c->normal_area[0],b=mass/c->normal_area[1];
            double energy=.5*(a*a+b*b);
            double normal=mass*(a*dot(c->direction[0],c->transverse)
                -b*dot(c->direction[1],c->transverse))/c->area;
            double area=cross(c->section[0],c->section[1])/c->area;
            double pressure=correction*.5*(gamma+1)*energy*fabs(c->curvature);
            needed=fmax(needed,(k*energy+fabs(normal)+pressure*(2+fabs(area)))/(2*k*h0));
            needed=fmax(needed,fmax(a*a,b*b)/sonic_sq);
        }
        offset+=(nx+1)*width;
        tube+=width-1;
    }
    if (!isfinite(needed)) return 0;
    double log_density=.5*log(needed/.95);
    for (int i=0;i<nx*tube;++i) g->packed[i]=log_density;
    return 1;
}

static AlulaBodyStatus evaluate_seed(AlulaAirfoilGrid *g,double *residual)
{
    AlulaBodyState state;
    if (!alula_airfoil_grid_state(g,&state)) return ALULA_BODY_INVALID;
    return alula_body_residual(g->body,&state,residual,NULL);
}

AlulaAirfoilGrid *alula_airfoil_grid_create(const AlulaContour *contours,size_t count,
    const AlulaAirfoilGridInput *in,AlulaAirfoilGridStatus *status,AlulaAirfoilGridReport *diagnostics)
{
    AlulaAirfoilGridStatus code=ALULA_GRID_INVALID;
    AlulaAirfoilGridReport report={.first_group=-1,.first_interval=-1,.first_tube=-1,
        .flow_seed=ALULA_BODY_INVALID};
    AlulaAirfoilGrid *g=NULL;
    AlulaPanelSystem *panel=NULL;
    AlulaPanelField *field=NULL;
    AlulaPanelGuides *guides=NULL;
    double *x=NULL,*weights=NULL,*fractions=NULL;
    AlulaPoint *banks=NULL;
    if (!in || !contours || count<1 || count>6 || !in->tubes ||
        in->surface_intervals<4 || in->surface_intervals>512 ||
        in->inlet_intervals<2 || in->inlet_intervals>512 ||
        in->outlet_intervals<2 || in->outlet_intervals>512 ||
        !(in->mach>0 && in->mach<1) || !(in->chord>0) || !isfinite(in->chord)) goto done;
    int total=0;
    for (size_t i=0;i<=count;++i) {
        if (in->tubes[i]<1 || in->tubes[i]>100000-total) goto done;
        total+=in->tubes[i];
    }
    code=ALULA_GRID_PANEL;
    panel=alula_panel_prepare(contours,count,in->chord,in->moment);
    field=alula_panel_field_create(panel,in->alpha);
    AlulaPanelResult flow;
    if (!field || !alula_panel_solve(panel,in->alpha,&flow,NULL)) goto done;
    code=ALULA_GRID_GUIDE;
    AlulaPanelGuideControls controls={in->inlet_x,in->outlet_x,in->lower_y,in->upper_y,
        in->chord,2e-9,.002,NULL};
    guides=alula_panel_guides_create(field,&controls,NULL);
    if (!guides) goto done;
    code=ALULA_GRID_ALLOCATION;
    g=calloc(1,sizeof *g);
    size_t cap=count*(in->surface_intervals+1)+in->inlet_intervals+in->outlet_intervals;
    x=malloc(cap*sizeof *x);
    if (!g || !x) goto done;
    const AlulaPanelBodyGuide *profiles=alula_panel_guides_bodies(guides,NULL);
    code=ALULA_GRID_STATIONS;
    int rows=stations(in,profiles,(int)count,g->cuts,x);
    if (!rows || total>100000/(rows-1)) goto done;
    report.nx=rows-1;
    report.elements=(int)count;
    int primary=0,finite=0;
    for (size_t i=0;i<=count;++i) g->tubes[i]=in->tubes[i];
    code=ALULA_GRID_ALLOCATION;
    for (size_t i=0;i<count;++i) {
        const AlulaContour *source=contours+profiles[i].element;
        finite|=source->finite_base!=0;
        if (!save_base(g->bases+i,source)) goto done;
        if (profiles[i].curve->arc[profiles[i].curve->count-1]>
            profiles[primary].curve->arc[profiles[primary].curve->count-1]) primary=(int)i;
    }
    AlulaBodyInput body={.nx=rows-1,.bodies=(int)count,.cuts=g->cuts,.tubes=g->tubes,
        .gamma=1.4,.mach=in->mach,.alpha=in->alpha,.length_scale=in->chord,.center=in->moment,.correction=.1,
        .mode=in->mode,.positive_simple=1,.displaced=1,.independent_wakes=finite,.bank_outlet=finite,
        .upwind=in->mode!=ALULA_BODY_ISENTROPIC,.ismom=in->ismom,.mucon=1,.mcrit=.99,.epsilon=1e-5};
    code=ALULA_GRID_ALLOCATION;
    g->body=alula_body_create(&body);
    if (!g->body) goto done;
    AlulaBodySize size=alula_body_size(g->body);
    g->nodes=malloc((size_t)size.nodes*sizeof *g->nodes);
    fractions=calloc(4*(size_t)rows*count,sizeof *fractions);
    banks=malloc(2*(count+1)*rows*sizeof *banks);
    weights=malloc((size_t)size.tubes*sizeof *weights);
    g->packed=calloc((size_t)size.rows,sizeof *g->packed);
    g->physical=malloc((size_t)alula_body_physical_layout(size).size*sizeof *g->physical);
    if (!g->nodes || !fractions || !banks || !weights || !g->packed || !g->physical) goto done;
    AlulaBodySurface surface[6];
    code=ALULA_GRID_STATIONS;
    if (!boundary_points(field,guides,profiles,g,x,rows-1,(int)count,surface,fractions,banks)) goto done;
    fill_passages(g,profiles,banks,fractions,x,weights,rows-1,(int)count);
    check_quality(g,&report);
    if (report.invalid_cells) {
        AlulaPoint *passage=g->nodes;
        for (size_t i=0;i<=count;++i) {
            alula_repair_passage(passage,rows-1,g->tubes[i]);
            passage+=(size_t)rows*(g->tubes[i]+1);
        }
        report.invalid_cells=0;
        report.first_group=report.first_interval=report.first_tube=-1;
        check_quality(g,&report);
    }
    code=ALULA_GRID_FOLDED;
    if (report.invalid_cells) goto done;
    code=ALULA_GRID_STATIONS;
    AlulaBodyCoordinateInput coordinate={.raw_nodes=g->nodes,.mass_weights=weights,
        .capture_levels=alula_panel_guides_levels(guides),.surfaces=surface,.primary_body=primary,
        .body_station_normals=1,.te_center_wake=finite,.walls_only=1};
    g->coordinates=alula_body_coordinates_create(g->body,&coordinate);
    if (!g->coordinates) goto done;
    g->chart=alula_body_coordinates_chart(g->coordinates);
    alula_body_coordinates_seed_circulation(g->coordinates,g->packed,.5*in->chord*flow.circulation_lift);
    if (!seed_gas(g,rows-1,(int)count,in->mach)) goto done;
    double *residual=malloc((size_t)size.rows*sizeof *residual);
    if (!residual) {code=ALULA_GRID_ALLOCATION;goto done;}
    report.flow_seed=evaluate_seed(g,residual);
    if (report.flow_seed!=ALULA_BODY_OK && pressure_domain_seed(g,rows-1,(int)count,in->mach))
        report.flow_seed=evaluate_seed(g,residual);
    free(residual);
    code=ALULA_GRID_OK;
done:
    free(x);free(weights);free(fractions);free(banks);
    alula_panel_guides_free(guides);alula_panel_field_free(field);alula_panel_free(panel);
    if (status) *status=code;
    if (diagnostics) *diagnostics=report;
    if (code==ALULA_GRID_OK) return g;
    alula_airfoil_grid_free(g);
    return NULL;
}

AlulaBodySize alula_airfoil_grid_size(const AlulaAirfoilGrid *g)
{
    return g?alula_body_size(g->body):(AlulaBodySize){0};
}

const AlulaBodyCut *alula_airfoil_grid_cuts(const AlulaAirfoilGrid *g) { return g?g->cuts:NULL; }
const int *alula_airfoil_grid_tubes(const AlulaAirfoilGrid *g) { return g?g->tubes:NULL; }

int alula_airfoil_grid_state(AlulaAirfoilGrid *g,AlulaBodyState *state)
{
    if (!g || !state || !g->chart.decode(g->chart.context,g->packed,g->physical,NULL)) return 0;
    AlulaBodySize size=alula_body_size(g->body);
    AlulaBodyPhysicalLayout p=alula_body_physical_layout(size);
    for (int i=0;i<size.nodes;++i)
        g->nodes[i]=(AlulaPoint){g->physical[2*i],g->physical[2*i+1]};
    *state=(AlulaBodyState){.nodes=g->nodes,.density=g->physical+p.density,.mass=g->physical+p.mass,
        .circulation=g->physical[p.strength],.source=g->physical[p.strength+1],
        .doublet_x=g->physical[p.strength+2],.doublet_y=g->physical[p.strength+3],
        .wake_width=size.wake_widths?g->physical+p.wake:NULL};
    return 1;
}

AlulaBodySolveResult alula_airfoil_grid_solve(AlulaAirfoilGrid *g,int iterations,double tolerance)
{
    if (!g) return (AlulaBodySolveResult){.status=ALULA_BODY_BAD_STATE,.norm=INFINITY,.max_mach=NAN};
    if (!g->solver) g->solver=alula_body_newton_create(g->body,&g->chart);
    AlulaBodySolveResult r=alula_body_newton_solve(g->solver,g->packed,iterations,tolerance);
    g->converged=r.status==ALULA_BODY_SOLVED;
    g->last_solve=r;
    return r;
}

void alula_airfoil_grid_free(AlulaAirfoilGrid *g)
{
    if (!g) return;
    for (int i=0;i<6;++i) free(g->bases[i].points);
    alula_body_newton_free(g->solver);alula_body_coordinates_free(g->coordinates);alula_body_free(g->body);
    free(g->physical);free(g->packed);free(g->nodes);free(g);
}
