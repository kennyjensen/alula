// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_PANEL_SOURCES_H
#define ALULA_PANEL_SOURCES_H
#include "panel_internal.h"

/* Influence for the source densities 1, xi and xi^2, with xi from 0 to 1.
 * The body uses only degree zero. The wake also uses the two higher moments. */
void alula_source_psi(const Panel *p,AlulaPoint target,double heading,double chart,int degree,double psi[3]);
void alula_source_velocity(const Panel *p,AlulaPoint target,int degree,int principal,
    int midpoint,AlulaPoint velocity[3]);

#endif
