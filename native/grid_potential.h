// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_GRID_POTENTIAL_H
#define ALULA_GRID_POTENTIAL_H

#include "airfoil_grid_internal.h"
#include "panel_guide.h"

/* Adjust the potential samples from guide crossings in place. Body guides
 * cover the inlet through stagnation and the relative wake. */
int alula_grid_guide_potential(const AlulaPanelField *field,
    const AlulaPanelBodyGuide *profiles,const AlulaAirfoilGrid *grid,const double *x,int nx,int count,double *phi);

/* Place interior nodes on traced streamfunction levels and keep both
 * supplied passage banks fixed. The caller owns the streamline caches,
 * with one slot per tube in each passage. On failure nodes are partial. */
int alula_fill_potential_grid(const AlulaPanelField *field,const AlulaPanelGuides *guides,
    const AlulaPanelBodyGuide *profiles,AlulaAirfoilGrid *grid,const AlulaPoint *boundaries,
    const double *guide_phi,const double *fractions,const double *weights,int nx,int count,double chord,
    AlulaPanelStreamline **const *lines);

#endif
