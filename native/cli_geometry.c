// SPDX-License-Identifier: GPL-2.0-or-later
#include "cli_geometry.h"

#include <stdio.h>

int alula_cli_contours(int topology,AlulaContour *contours,AlulaPoint *points,
    size_t *element_count,size_t *point_count)
{
    size_t count;
    if (scanf("%zu", &count) != 1 || count < 1 || count > 6) return 0;
    size_t used = 0;
    for (size_t e = 0; e < count; ++e) {
        size_t length;
        if (scanf("%zu", &length) != 1 || length < 4 || length > ALULA_PANEL_MAX_TOTAL+6-used) return 0;
        contours[e] = (AlulaContour){.points=points+used, .count=length};
        if (topology) {
            char kind;
            if (scanf(" %c", &kind) != 1 || (kind != 'S' && kind != 'F')) return 0;
            if (kind == 'F') {
                contours[e].finite_base = 1;
                if (scanf("%zu%zu", &contours[e].upper_index, &contours[e].lower_index) != 2) return 0;
            }
        }
        for (size_t i = 0; i < length; ++i)
            if (scanf("%lf%lf", &points[used+i].x, &points[used+i].y) != 2) return 0;
        used += length;
    }
    *element_count=count;
    *point_count=used;
    return 1;
}
