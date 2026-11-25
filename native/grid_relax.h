// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_GRID_RELAX_H
#define ALULA_GRID_RELAX_H

#include "panel.h"

/* Bounded geometry repair of one structured passage. Boundaries stay fixed.
 * Storage is station major with (tubes+1) points per station. normal_weight
 * adds a penalty on squared motion across the original streamline. */
int alula_repair_passage(AlulaPoint *nodes,int nx,int tubes,double normal_weight);

#endif
