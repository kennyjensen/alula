// SPDX-License-Identifier: GPL-2.0-or-later
#include "cli_geometry.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
    int pressure = 0, topology = 0, first = 1;
    while (first < argc) {
        if (strcmp(argv[first], "--cp") == 0 && !pressure) pressure = 1;
        else if (strcmp(argv[first], "--topology") == 0 && !topology) topology = 1;
        else break;
        ++first;
    }
    if (argc <= first) {
        fprintf(stderr, "Usage: %s [--cp] [--topology] ALPHA [ALPHA ...] < contours.txt\n", argv[0]);
        return 2;
    }

    size_t count, used;
    AlulaContour contours[6];
    AlulaPoint points[ALULA_PANEL_MAX_TOTAL+6];
    AlulaPanelSystem *system = NULL;
    int status = 1;
    if (!alula_cli_contours(topology,contours,points,&count,&used)) goto invalid;
    system = alula_panel_prepare(contours, count, 1, (AlulaPoint){.25, 0});
    if (!system) goto invalid;
    for (int i = first; i < argc; ++i) {
        char *end;
        errno = 0;
        double alpha = strtod(argv[i], &end);
        double cp[ALULA_PANEL_MAX_TOTAL];
        AlulaPanelResult result;
        if (errno || end == argv[i] || *end || !alula_panel_solve(system, alpha, &result, pressure ? cp : NULL)) {
            fprintf(stderr, "Invalid angle or nonfinite solution: %s\n", argv[i]);
            goto done;
        }
        printf("{\"alpha\":%.17g,\"cl\":%.17g,\"cm\":%.17g,\"pressureDrag\":%.17g,"
            "\"circulationLift\":%.17g,\"residualBound\":%.17g",
            alpha, result.cl, result.cm, result.pressure_drag,
            result.circulation_lift, result.residual_bound);
        if (pressure) {
            printf(",\"cp\":[");
            for (size_t j = 0; j < used-count; ++j)
                printf("%s%.17g", j ? "," : "", cp[j]);
            putchar(']');
        }
        puts("}");
    }
    status = 0;
    goto done;
invalid:
    fprintf(stderr, "Invalid geometry or singular panel system. Use prepared CCW contours with explicit finite-base corners.\n");
done:
    alula_panel_free(system);
    return status;
}
