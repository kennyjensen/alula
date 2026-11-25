// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_BL_INTERNAL_H
#define ALULA_BL_INTERNAL_H

#include "bl.h"

typedef struct {
    AlulaBLStation input;
    AlulaBLProperties props;
    double d,hc,us,thickness;
} AlulaBLPreparedStation;

/* The upstream station stays fixed during a local Newton solve. A mixed
 * transition interval uses the laminar closure at its upstream end. */
AlulaBLStatus alula_bl_prepare_upstream(const AlulaBLContext *c,const AlulaBLStation *upstream,
    AlulaBLRegime regime,AlulaBLPreparedStation *prepared);
AlulaBLStatus alula_bl_prepared_interval(const AlulaBLContext *c,const AlulaBLPreparedStation *upstream,
    const AlulaBLStation *downstream,AlulaBLRegime regime,double trip,double exponent,AlulaBLBlock *out);

AlulaBLMarchResult alula_bl_march_wake_join(const AlulaBLContext *c,
    const AlulaBLStation *upper,const AlulaBLStation *lower,const AlulaBLEdge *edges,size_t count,
    double join_gap,const double *displacement,AlulaBLStation *profile);
#endif
