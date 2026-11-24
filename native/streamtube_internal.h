// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_STREAMTUBE_INTERNAL_H
#define ALULA_STREAMTUBE_INTERNAL_H

#include "streamtube.h"

/* Whole-grid assembly uses the physical sections before applying the shared
 * speed filter. No undissipated interface-pressure gate applies in between. */
AlulaCellStatus alula_cell_initialize(AlulaCell *cell, const AlulaCellInput *input);
AlulaCellStatus alula_cell_flux(AlulaCell *cell);

#endif
