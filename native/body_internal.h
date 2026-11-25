// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_BODY_INTERNAL_H
#define ALULA_BODY_INTERNAL_H

#include "body.h"
#include "streamtube_internal.h"
#include <math.h>
#include <stdint.h>

#define MAX_BODIES 6
#define TWO_PI 6.283185307179586476925286766559

typedef struct {
    AlulaBodyFamily family;
    int i, g, j;
} Row;

typedef struct {
    AlulaPoint velocity, basis[5];
    double pressure, log_radius;
} Farfield;

typedef struct {
    AlulaSection gas;
    double length, transport;
} Section;

typedef struct { double speed, mach_sq, enthalpy, length, density, transport; } SectionTangent;
typedef struct { double streamwise, pressure[2]; } CellTangent;

struct AlulaBody {
    AlulaBodyInput in;
    AlulaBodySize size;
    uint64_t residual_generation;
    AlulaBodyCut cuts[MAX_BODIES];
    int tubes[MAX_BODIES+1], column[MAX_BODIES+2], node[MAX_BODIES+2];
    int wake[MAX_BODIES+1];
    AlulaPoint freestream;
    double mach_sq, beta, pg_a, pg_b, pg_scale, p_inf, h0, rho_total;
    Row *rows;
    Section *sections;
    AlulaUpwindStencil *upwind_stencils;
    int check_all_upwind_tangents;
    AlulaCell *cells;
    unsigned char *blend_regime;
    Farfield *jacobian_fields;
    uint64_t jacobian_field_generation;
    SectionTangent *section_tangents;
    CellTangent *cell_tangents;
    int failed_initial_cell;
};

/* Update the freestream and invalidate physical stencils from the previous Mach. */
int alula_body_set_mach(AlulaBody *body,double mach);

static inline AlulaPoint add(AlulaPoint a, AlulaPoint b) { return (AlulaPoint){a.x+b.x,a.y+b.y}; }
static inline AlulaPoint sub(AlulaPoint a, AlulaPoint b) { return (AlulaPoint){a.x-b.x,a.y-b.y}; }
static inline AlulaPoint scale(AlulaPoint a,double b) { return (AlulaPoint){a.x*b,a.y*b}; }
static inline AlulaPoint mean(AlulaPoint a, AlulaPoint b) { return scale(add(a,b),.5); }
static inline double cross(AlulaPoint a, AlulaPoint b) { return a.x*b.y-a.y*b.x; }
static inline double dot(AlulaPoint a, AlulaPoint b) { return a.x*b.x+a.y*b.y; }
static inline double norm(AlulaPoint a) { return hypot(a.x,a.y); }
static inline int positive(double x) { return isfinite(x) && x>0; }
static inline int finite_point(AlulaPoint p) { return isfinite(p.x) && isfinite(p.y); }

static inline int index_section(const AlulaBody *b,int g,int i,int j)
{
    return i*b->size.tubes+b->column[g]+j;
}

static inline AlulaPoint point(const AlulaBody *b,const AlulaBodyState *s,int g,int i,int j)
{
    return s->nodes[b->node[g]+i*(b->tubes[g]+1)+j];
}

static inline AlulaCell *cell(const AlulaBody *b,int g,int i,int j)
{
    return b->cells+index_section(b,g,i-1,j);
}

/* One cursor per passage bank and probe state. Initialize station to zero. */
typedef struct {int station; AlulaCell cells[2];} AlulaEdgeCursor;
double alula_body_edge_speed(const AlulaBody *b,const AlulaBodyState *probe,int g,int i,int j,
    AlulaEdgeCursor *cursor);
int alula_body_farfield(const AlulaBody *body,const double *coeff,AlulaPoint p,Farfield *out);
int alula_body_streamwise(const AlulaBody *body,Row row,const AlulaCellValue *direction,double *value);

/* The state must match the latest successful residual evaluation. The
 * half-open section range must cover all geometry, mass and density changes.
 * Edges optionally receives three normalized matching contributions per
 * outer edge in the even slots of a 6 * 2 * nx workspace. */
AlulaBodyStatus alula_body_apply_prepared(AlulaBody *body,const AlulaBodyState *state,
    const AlulaBodyState *direction,int first_section,int end_section,double *product,double *edges);

#endif
