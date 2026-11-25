// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef ALULA_CLI_GEOMETRY_H
#define ALULA_CLI_GEOMETRY_H

#include "panel.h"

/* Read the shared alula-panel/alula-euler geometry format from stdin. */
int alula_cli_contours(int topology,AlulaContour *contours,AlulaPoint *points,
    size_t *element_count,size_t *point_count);

#endif
