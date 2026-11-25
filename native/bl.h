// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_BL_H
#define ALULA_BL_H

#include <stddef.h>

typedef enum {
    ALULA_BL_SIMILARITY, ALULA_BL_LAMINAR, ALULA_BL_TURBULENT, ALULA_BL_WAKE, ALULA_BL_TRANSITION
} AlulaBLRegime;
typedef enum {
    ALULA_BL_OK, ALULA_BL_INVALID, ALULA_BL_DOMAIN, ALULA_BL_NO_TRANSITION,
    ALULA_BL_TRANSITION_ROOT, ALULA_BL_LINEAR, ALULA_BL_ITERATION_LIMIT, ALULA_BL_LINE_SEARCH
} AlulaBLStatus;

typedef struct {
    double reynolds,mach,gamma,ncrit;
    double gm1,rho_total,enthalpy_inverse,re_total;
} AlulaBLContext;

/* Lengths in reference chords; velocity normalized by freestream. aux
 * denotes N for laminar flow and square root of the shear coefficient
 * for turbulent flow or wake. delta_star includes the dead-air wake gap. */
typedef struct { double s,theta,delta_star,ue,aux,wake_gap; } AlulaBLStation;

typedef struct {
    double rho,viscosity,mach_sq,h,raw_hk,hk,hstar,re_theta,cf,dissipation;
    double equilibrium_shear,transition_shear;
} AlulaBLProperties;

typedef struct {
    int active,forced;
    double s,amplification;
} AlulaBLTransition;

typedef struct {
    double residual[3]; /* amplification/shear, momentum, kinetic energy */
    AlulaBLProperties downstream;
    AlulaBLTransition transition;
} AlulaBLBlock;

/* XFOIL dead-air thickness downstream of a finite trailing edge. */
double alula_bl_dead_air(double width,double tangent_cross,double distance);
/* Initialized contexts are immutable and may be shared across threads. */
int alula_bl_context(AlulaBLContext *ctx,double reynolds,double mach,double gamma,double ncrit);
AlulaBLStatus alula_bl_station(const AlulaBLContext *ctx,const AlulaBLStation *station,
    AlulaBLRegime regime,AlulaBLProperties *out);
AlulaBLStatus alula_bl_check_transition(const AlulaBLContext *ctx,const AlulaBLStation *upstream,
    const AlulaBLStation *downstream,double trip_s,AlulaBLTransition *out);
/* Similarity uses only downstream. Pass a finite trip beyond the stations
 * (e.g. DBL_MAX) to allow natural transition without a forced trip. The
 * transition block handles mixed intervals downstream of similarity. */
AlulaBLStatus alula_bl_interval(const AlulaBLContext *ctx,const AlulaBLStation *upstream,
    const AlulaBLStation *downstream,AlulaBLRegime regime,double trip_s,
    double similarity_exponent,AlulaBLBlock *out);
AlulaBLStatus alula_bl_trailing_edge(const AlulaBLContext *ctx,const AlulaBLStation *upper,
    const AlulaBLStation *lower,const AlulaBLStation *wake,double gap,AlulaBLBlock *out);

typedef struct { double s,ue,wake_gap; } AlulaBLEdge;
typedef struct {
    AlulaBLStatus status;
    size_t completed;
    int iterations,backtracks;
    double maximum_residual,edge_residual;
    AlulaBLTransition transition;
    size_t edge_mismatch_index; /* Defined when edge_residual is positive. */
    double euler_edge_speed;
} AlulaBLMarchResult;

/* Direct profile initialization at fixed edge velocities. The result owns
 * no storage. On failure, only the first completed entries are published.
 * For an untripped march, use DBL_MAX as the trip position. */
AlulaBLMarchResult alula_bl_march_surface(const AlulaBLContext *ctx,
    const AlulaBLEdge *edges,size_t count,double trip_s,double similarity_exponent,
    AlulaBLStation *profile,AlulaBLRegime *regimes);
/* Pass count nonnegative displacements in reference chords, or NULL for a
 * fixed-edge march. Positive entries use a local interaction law. Zero
 * uses the direct march with an inverse fallback for separation. The caller
 * must check the edge mismatch separately. */
AlulaBLMarchResult alula_bl_march_surface_interaction(const AlulaBLContext *ctx,
    const AlulaBLEdge *edges,size_t count,double trip_s,double similarity_exponent,
    const double *displacement,AlulaBLStation *profile,AlulaBLRegime *regimes);
/* Local mode keeps Euler edge speeds upstream of the first failed
 * interval that needs interaction, then enables the law downstream. */
AlulaBLMarchResult alula_bl_march_surface_controlled(const AlulaBLContext *ctx,
    const AlulaBLEdge *edges,size_t count,double trip_s,double similarity_exponent,
    const double *displacement,int local,AlulaBLStation *profile,AlulaBLRegime *regimes);
/* The two incoming TE stations use positive turbulent shear in aux.
 * The first wake edge supplies the material gap at the junction. */
AlulaBLMarchResult alula_bl_march_wake(const AlulaBLContext *ctx,
    const AlulaBLStation *upper,const AlulaBLStation *lower,
    const AlulaBLEdge *edges,size_t count,AlulaBLStation *profile);
/* Allow local recovery in the wake. Displacement has count-1 entries in
 * reference chords, including each station's dead-air gap. */
AlulaBLMarchResult alula_bl_march_wake_interaction(const AlulaBLContext *ctx,
    const AlulaBLStation *upper,const AlulaBLStation *lower,
    const AlulaBLEdge *edges,size_t count,const double *displacement,AlulaBLStation *profile);

#endif
