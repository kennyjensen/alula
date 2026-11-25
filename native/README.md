# Alula native C solver

This directory contains the incompressible linear vortex panel solver for
Ubuntu. It handles one to six elements with up to 700 wetted panels.
Assemblies with a finite trailing edge use the nodal streamfunction system,
with source and vortex influences from each retained base segment.

The native Euler verification solver assembles moving streamlines, the
conservative cell equations, inlet and outlet conditions, and a sparse
Jacobian. It supports subsonic channels and transonic nozzle shocks. The
airfoil body evaluator accepts prepared multi-element grids, including
separate finite trailing-edge wake banks. It provides analytic Jacobian
products and a KLU Newton solve using native coordinates. The native Euler
command also builds a startup grid from contours. Finer and closely spaced
multi-element grids still need startup improvements.

```sh
make -C native
```

GCC or Clang and the C math library suffice for this build. For threaded
influence assembly with GCC and LAPACK factorization through OpenBLAS, install
Ubuntu's `build-essential` and `libopenblas-dev`, then run:

```sh
make -C native clean
make -C native OPENMP=1 LAPACK_LIBS=-lopenblas
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 make -C native check
```

The parity check requires Node.js and compares surface pressure and integrated
loads with the JavaScript panel solver. It also checks intrinsic Euler residuals
and analytic directional derivatives for geometry, density, and shock
transport against the JavaScript kernels. The C executable itself does not
use Node.js.

Pass the number of elements, then for each element its point count and its
whitespace-separated coordinate pairs on stdin. Each contour must be disjoint,
counterclockwise, and closed by repeating its first point. A sharp contour
starts at the trailing edge.
Coordinates use a unit reference chord and the moment reference is (0.25, 0).
The command emits one JSON line per angle in degrees:

```sh
native/alula-panel --cp -4 0 4 8 < contours.txt
```

With `--topology`, append `S` to each sharp element's point count, or
`F UPPER LOWER` to a finite element's point count. The two zero-based indices
identify the upper and lower trailing-edge vertices in the supplied contour.
The closing vertex counts toward the point count but cannot be used as a
corner index. The total limit including base segments is 4,096 panels.

Omit `--cp` for integrated loads only. Pressure drag is the integrated
inviscid pressure force projected along the freestream. The solver factors the
influence matrix once per geometry and caches two freestream responses.
Integrated loads then take constant work per angle. Returning surface pressure
takes one pass over the panels.

`alula_panel_field_create` prepares the velocity and streamfunction field for
one angle using the existing factorization responses. The field supports
sharp and finite trailing edges, including the source branch jumps behind
finite bases. `trace.h` advances this field in velocity potential or physical
arc length. It can stop on a plane, reject steps that intersect the solid using
caller-supplied geometry checks, and stream accepted points to the grid
initializer without allocating a path.

`panel_guide.h` prepares ordered inlet streamfunction levels and traces the
incoming dividing curves, wake centerlines, and outer streamlines. It uses
the collocation-speed stagnation estimate and rejects steps through the solid
contours. Interior streamlines start from a specified inlet streamfunction
level. Guide intersections with a supplied x-plane must be unique in the
resolved path.

`alula-euler` reads the same contours and trailing-edge metadata as
`alula-panel`. For a compressible Euler run, pass angle in degrees, Mach,
and optionally the surface interval and tube counts. The defaults of six
surface intervals and four tubes per passage are a coarse startup setting.
The command emits a JSON grid and solve report and exits with status 1 on
grid failure or nonconvergence.

```sh
native/alula-euler --topology 2.68 0.74 6 4 < contours.txt
native/alula-euler --topology --mesh-only 2.68 0.74 6 4 < contours.txt
```

`airfoil_grid.h` exposes grid setup, the moving-grid solve, and the current
nodes, density, passage mass flow, and farfield strengths. The startup grid
uses contour splines and panel-guided boundaries with interpolated interior
cross-lines. Bounded vertex repair preserves the solid and guide boundaries.
The contour-to-Euler checks include coarse NACA 2412 and stacked NACA 0012
cases, a coarse transonic RAE 2822 solve, and the starting grid for NLR 7301
with flap.

`alula-euler --bl RE NCRIT` marches the integral boundary layer and merged
wake on the converged Euler edge speeds. The report includes completed
stations, transition locations, and terminal thicknesses for each branch.
`--coupled RE NCRIT` alternates between the boundary-layer march and a warm
Euler solve on the displaced grid. It checks the Euler and boundary-layer
residuals and requires the displacement thickness to agree within 1e-8
reference chords. The report gives lift and moment from surface pressure,
along with separate pressure and viscous wake drag estimates.

```sh
native/alula-euler --topology --bl 1000000 9 0 0.1 8 4 < contours.txt
native/alula-euler --topology --coupled 1000000 9 0 0.1 8 4 < contours.txt
```

Coupled checks cover coarse subsonic sharp and finite-base NACA 2412 and
stacked NACA 0012. The direct march can stop in adverse pressure gradients,
including during coupled startup for the coarse transonic RAE 2822 case.

`bl.h` exposes similarity, laminar, mixed transition, turbulent and wake
residuals and fixed-edge profile initialization. The C equations are
derived from Mark Drela's GPL XFOIL and Alula's adapted Vibefoil port.
See [source provenance](../src/viscous/xfoil/PROVENANCE.json) and
[the vendored XFOIL sources](../third_party/Xfoil).

`libalula.a` exposes the panel, Euler cell, channel, and body interfaces in
`panel.h`, `streamtube.h`, `channel.h`, and `body.h`. An Euler cell can be prepared
once and reused for several Jacobian columns. The channel solves density and
streamline positions together, with free tube masses for reservoir inflow
and static-pressure outflow. It uses the existing
[SuiteSparse KLU sources](../third_party/klu/README.md), reuses symbolic
analysis and numerical factor storage, and verifies each linear solution
against the original matrix.

The body API evaluates the same pressure matching, edge conditions,
farfield modes, wake gaps, and ISMOM choices as the compressible web solver
on supplied physical nodes, densities, and tube masses. It reuses work
arrays across evaluations. The OpenMP build parallelizes cell initialization
and flux evaluation on larger grids. Use `ALULA_BODY_BENCH=1 make -C native check`
to time the residual parity cases.

`alula_body_linearize` accepts a direction in the same physical variables
and returns the residual and its directional derivative. It propagates shared
upwind speeds, farfield normalization, and independent wake widths. Callers
using packed coordinates supply the node and gas directions from their
coordinate chart.

`body_newton.h` accepts a native coordinate-chart callback and its derivative
pattern. It assembles the square sparse Newton matrix with graph coloring
and reuses the KLU ordering across accepted steps. Each trial checks the
physical Euler equations and the returned root must pass a convex-grid check.
The tests converge single and multi-element grids on a fixed chart and
compare sparse matrix products for rounded bodies and independent wakes.

`body_chart.h` builds the coordinates from a prepared grid, contour splines,
and tube mass weights. It supports moving rounded leading-edge stagnation
points and prescribed displacement, and updates grid normals after accepted
Newton steps while preserving the accepted physical geometry.
