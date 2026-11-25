# Alula native C solver

## Desktop interface

From the repository root on an Ubuntu desktop, run:

```sh
sudo apt install build-essential python3-tk
make -C native gui
```

Select NACA 0012, RAE 2822, NLR 7301 with flap, or open a file in the
native contour format described below. Choose panel flow or Euler, each available with integral boundary layers,
and click **Run**. Inviscid panel flow
uses the input contour panels and angle of attack. Its direct solve returns
surface pressure, loads and a residual bound. For Euler, set the flow
conditions and mesh to follow the iterative solve. **Cancel**
stops the active solver process. The airfoil view shows streamtube centers
from the current Euler mesh or the current wake paths in viscous panel mode. The pressure plot separates upper and
lower surfaces. During coupled iterations, select **Displacement thickness**
above the pressure plot to see δ*/c on both surfaces of each element,
against distance from stagnation. Pressure drag is available during Euler iterations and wake
drag appears when a complete boundary-layer profile is available. The
status shows the actual Mach during continuation. A converged result
replaces the live iterate when the solver finishes.

Run `make -C native check-gui` on the desktop to exercise the live plots,
cancellation and all three solver executables. On a headless Ubuntu installation
with `xvfb` and `xauth`, use `xvfb-run -a make -C native check-gui`.

For other viewers, `alula-euler --progress` emits newline-delimited JSON
on standard error for the Euler seed, accepted precursor steps, coarse-to-fine
displacement transfer and coupled iterates. Each record includes its stage,
local iteration, Mach, residual and available streamlines, pressure samples
and loads. `boundaryLayerCurves` carries the completed surface stations
as pairs of normalized arc distance and displacement thickness. A final
`retained` record supplies the available curves and streamlines for the reported state.
The desktop uses `--live-progress` to sample the plot geometry every 50 ms
and at stage or Mach changes. It receives coefficients and residuals at every step.
Coupled residuals report the largest Euler, boundary-layer, edge or
displacement residual. Iteration
numbers restart across grids and continuation attempts. Standard output
contains the usual final JSON report.

## Native solvers

This directory contains the incompressible linear vortex panel solver for
Ubuntu. It handles one to six elements with up to 700 wetted panels.
Assemblies with a finite trailing edge use the nodal streamfunction system,
with source and vortex influences from each retained base segment.

Use `alula-viscous` for a Mach-zero viscous analysis. Supply the angle in degrees,
Reynolds number and Ncrit. The input geometry follows the same contour format
as `alula-panel`. For a single finite-base element with upper and lower trips
at 8% and 13% chord, run

```sh
native/alula-viscous --topology --trip .08 .13 --progress 4 1000000 9 < contours.txt
```

The C solver couples all surface and wake boundary layers to one panel
displacement field. Each Newton step eliminates three local unknowns per
station and checks the reduced solution against the full linearization.
The two freestream responses on each fixed wake allow incidence continuation
without rebuilding the panel matrix. For difficult starts the solver tries a
lower-incidence seed, then advances to the requested incidence in steps of at
most 2°. Progress includes both actual and requested incidence.

Each wake advances in its local mean flow direction while preserving segment
lengths and its trailing-edge origin. Success requires a maximum equation
residual of `1e-8` and a maximum wake-normal speed of `1e-6 U∞`.
Pressure loads integrate nodal Cp around each element. Wake drag uses
Squire–Young extrapolation. Thickness above 10% of the local chord or a wake
exit inside the dead-air region withholds final loads. `--iterations` limits
Newton updates per fixed geometry and incidence. `--wake-limit` limits wake
updates. `--wake-count` sets the number of stations per wake (default 24).
The library interface is in [panel_viscous.h](panel_viscous.h).
After building `check-displacement` and `alula-viscous`, run
`ALULA_PANEL_BENCH=1 node native/check-displacement.mjs` from the repository
root to time operator assembly and the complete 96-panel NACA 2412 viscous
case. The end-to-end comparison uses the native CLI and the web solver in
Node, checks agreement of the loads and reports five-run medians.

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

GCC or Clang and the C math library suffice for this build. With GCC,
`make -C native OPENMP=1` parallelizes independent guide and interior
streamline traces, panel influence assembly, and larger Euler kernels.
Set `OMP_NUM_THREADS` to control the thread count at runtime. For LAPACK
factorization through OpenBLAS, install
Ubuntu's `build-essential` and `libopenblas-dev`, then run:

```sh
make -C native clean
make -C native OPENMP=1 LAPACK_LIBS=-lopenblas
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 make -C native OPENMP=1 LAPACK_LIBS=-lopenblas check
```

The parity check requires Node.js and compares surface pressure and integrated
loads with the JavaScript panel solver. It also checks intrinsic Euler residuals
and analytic directional derivatives for geometry, density, and shock
transport against the JavaScript kernels. The C executable itself does not
use Node.js.

OpenBLAS accelerates the dense panel factorization during geometry setup.
The Euler Newton steps use sparse KLU. For a large single-block matrix,
the solver compares AMD ordering with and without KLU's preliminary row
matching. It uses the alternative when the symbolic estimates predict
fewer operations and no extra fill, and retains the matched ordering in
case the numerical factorization or residual check fails.

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

Distant velocity samples use 15 precomputed multipole moments for groups of
16 wetted panels. Within 16 group radii the solver evaluates each panel
directly. Finite trailing-edge bases retain direct evaluation throughout.

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
For provisional solves, `finalGrid` measures the last retained grid,
including the first folded cell if present. `stepRejections` breaks down
Euler backtracks by chart, flow domain and residual. A linear solve that
fails its accuracy check reports `linear.rejectedRelativeResidual` when
the measured error is available. `residualByFamily` shows which Euler
equations dominate the last accepted state.

```sh
native/alula-euler --topology 2.68 0.74 6 4 < contours.txt
native/alula-euler --topology --mesh-only 2.68 0.74 6 4 < contours.txt
```

`--cp` adds pressure samples at the material coordinates used to integrate
the loads. Each sample includes its input contour index, surface, coordinates
and pressure coefficient. A finite trailing-edge base includes its original
points and uses the same interpolated pressure as the load integral. When
Euler has converged and boundary-layer coupling is still provisional, the
report marks those samples with `cpProvisional`.

`--inlet-intervals N` and `--wake-intervals N` set the requested grid's
upstream and downstream interval counts independently. Intermediate coarse
grids retain their resolution-based counts. `--domain XMIN XMAX YMIN YMAX`
sets the outer bounds for every stage in the contour coordinate frame. The
domain must enclose all contours. For example, use the MIT RAE 2822 domain
with 16 intervals on each end to build a fine grid:

```sh
native/alula-euler --topology --potential-grid --mesh-only \
  --domain -2 3 -3 3.5 --inlet-intervals 16 --wake-intervals 16 \
  2.68 .74 64 24 < rae-contour.txt
```

`airfoil_grid.h` exposes grid setup, the moving-grid solve, and the current
nodes, density, passage mass flow, and farfield strengths. The startup grid
uses contour splines and panel-guided boundaries with interpolated interior
crosslines. If they fold on a rounded nose, a second construction traces
interior streamlines and aligns them by relative velocity potential.
Guide crossings supply both position and potential. Retries with gentler
surface spacing reuse those samples and the interior streamline traces.
`--potential-grid` selects this construction immediately. The default CLI
also retries it if the Euler solve stalls on an algebraic grid. The report
identifies the active `crosslines` method. Bounded vertex repair preserves
the solid and guide boundaries. If an isentropic seed has a nonpositive
interface pressure, the initializer first raises density locally and keeps
the rest of the seed. Grid tests cover NLR 7301 with flap and RAE 2822 with
64 surface intervals and up to 24 tubes per passage. If a single-element
Euler solve stalls at 32 or more surface intervals, the CLI retries with
less clustering at the leading and trailing edges and a coarse Euler seed.
This converges the 64×24 RAE 2822 Euler case on the MIT domain at Mach 0.74
and 2.68 degrees.
To specify the station distribution, pass `--surface-uniformity F`, where
zero selects cosine spacing, one selects uniform spacing, and intermediate
values blend the two. A spacing retry reports the chosen `surfaceUniformity`
in the JSON output.

When a multi-element Euler startup stalls, the CLI can retry with a
coordinate chart that carries stagnation-point motion into neighboring
streamlines. At the same resolution, the retry reuses the initial potential
mesh and creates a fresh coordinate chart and solver. If cells remain
concave, it can reposition interior nodes and solve again. If this repair
fails on a grid that met the Euler residual tolerance, it retries from the
same state, favoring motion along the streamlines. Each retry must meet
the original residual and convexity criteria.
The JSON report identifies this route with `startupRecovery` and records
successful repairs in `startupGridRepaired`. This converges the coarse
NLR 7301 with flap Euler case at Mach 0.185 and 6 degrees on both six-
and eight-interval meshes with four tubes per passage. The Euler case
also converges at 12 intervals and six tubes with `--surface-uniformity .7`.

At 32 surface intervals and 24 tubes per passage, the NLR 7301 Euler solve
converges with automatic spacing. If the original station distribution
stalls, the solver retries on a fresh wall-following mesh with 70% uniform
spacing. If a repair halves the flow residual but Newton folds a cell again,
the next repair favors movement along the local streamline. The solver
checks convexity and convergence on the requested resolution before
returning loads. An explicit `--surface-uniformity` takes precedence.
At the same resolution with natural transition, Re 2.51 million and
Ncrit 9, the coupled solve transfers the full coarse displacement, then
stalls with an edge-speed mismatch of about 0.0205.

Use `--coarse-start N` to solve an inviscid grid at a lower surface
resolution with at most four tubes per passage, then seed the requested
grid with its geometry, density, stagnation locations, and farfield
strengths. The target keeps its own wall spline and resolution. The JSON
field `coarseStartApplied` gives the transfer status, and `converged`
describes the requested grid. For RAE 2822, the 16-interval Euler solve
converges with an eight-interval seed:

```sh
native/alula-euler --topology --potential-grid --coarse-start 8 2.68 .74 16 4 < rae-contour.txt
```

Use comma-separated stages to refine both grid directions. Each stage is
`N` surface intervals or `NxT` for `T` tubes per passage. The surface
interval counts increase at each stage and stay below the requested mesh.
Each stage starts from the last converged stage when a transfer is possible.
When interpolation folds a cell, the transfer can repair interior nodes
with the walls and passage boundaries fixed. It accepts the repaired
seed only when its decoded grid is convex and the Euler evaluator admits it.
The stages share one panel factorization and set of guide paths. A passage
also reuses its interior streamline traces when its tube count stays the same.
`coarseStagesSolved` reports how many intermediate grids converged.
The 24-interval transonic RAE 2822 Euler case uses this sequence:

```sh
native/alula-euler --topology --potential-grid --coarse-start 8,16x4,20x12 2.68 .74 24 12 < rae-contour.txt
```

The NLR 7301 with flap converges at Mach 0.185 and 6 degrees on a
10-interval Euler mesh with four tubes per passage using the same transfer:

```sh
native/alula-euler --topology --potential-grid --coarse-start 8x4 6 .185 10 4 < nlr7301-contours.txt
```

`alula-euler --bl RE NCRIT` marches the integral boundary layer and merged
wake on the converged Euler edge speeds. The report includes completed
stations, transition locations, and terminal thicknesses for each branch.
During each station's Newton solve, the march reuses the upstream fluid
properties and closure values across downstream trials. It recomputes them
when the interval switches to the transition equations. Sampling the wall
uses cumulative arc lengths from the fixed material contour.

For the 64×24 coupled RAE 2822 case, start at Mach 0.64 and advance the
converged state on the same grid to Mach 0.74.

```sh
awk 'NR == 1 {print "1\n129 S"} NR > 2' third_party/airfoils/rae2822-mses/blade.rae | \
  native/alula-euler --topology --potential-grid --coarse-start 8 \
  --domain -2 3 -3 3.5 --inlet-intervals 16 --wake-intervals 16 \
  --surface-uniformity .7 --mach-start .64 --coupled 2700000 4 \
  --trip .03 .07 2.68 .74 64 24
```

`--mach-start` first tries alternating flow and boundary-layer solves at
each higher Mach. If a stage fails, it restores the last converged state
and advances both sets of equations in the same Newton step. It reduces
the Mach increment if needed. When assembling the coupled Newton matrix,
the solver groups independent boundary-layer displacement perturbations
into batches. On a single-block Euler matrix it applies each triangular
factor entry to a contiguous row of right-hand sides, so the CPU can reuse
the entry and vectorize the updates. KLU handles smaller batches and
block-triangular systems. Each column has its own residual check and
receives iterative refinement as needed.
Boundary-layer probes decode the two wall-adjacent tubes in each passage.
Consecutive samples along a wall or wake share the geometry and gas state
of their common section. The final Euler residual still uses the entire grid.
The report includes the requested `mach`,
`actualMach`, and the counts of accepted stages and retries. The residuals
and iteration statistics describe the last retained stage. If continuation
stops short, the process exits with a failure status and reports
`converged: false`. In C, use `alula_airfoil_coupled_solve_to_mach` and
`alula_airfoil_grid_mach` to run and inspect the same calculation.
Set `ALULA_FINE_COUPLED=1` when running `make -C native check` to include
this longer case in the test suite.

`--coupled RE NCRIT` alternates between the boundary-layer march and a warm
Euler solve on the displaced grid. It checks the Euler and boundary-layer
residuals and requires the displacement thickness to agree within 1e-8
reference chords and the edge speeds within 1e-8 times freestream speed.
If the fixed-edge start stalls, the coupled mode can seed a separated
laminar or turbulent surface by solving for its edge speed at a prescribed
shape, then iterate with a displacement-based interaction law and Anderson
acceleration. When
that march stalls, the solver can retain Euler edge speeds upstream of the
first interval that needs interaction and activate it from that interval
onward. The wake can also retry a failed interval with interaction, and its
edge-speed mismatch counts toward convergence. The JSON report counts
rejected steps by Euler failure, boundary-layer failure, and insufficient
reduction in the coupling error. Each affected profile reports the
station, arc length, and both speeds at its largest `edgeMismatch`.
The report also gives lift and moment from surface pressure and separate
pressure and viscous wake drag estimates.
If an inner Euler solve meets the residual tolerance on a concave grid,
the solver can repair the grid and solve again before trying a smaller
coupling step. `eulerGridRepairs` records successful inner repairs.

`--refine-coupled` requires `--coupled` and `--coarse-start`. After the
requested Euler grid converges, it couples the last converged coarse grid
and gradually applies its displacement to the requested grid. The report
shows the fraction of that displacement reached in
`coupledRefinementFraction`. After a thickness-scaled displacement step is
accepted, the next relaxation factor can grow to at most twice its accepted
value. This limits repeated Euler solves on oversized trial steps.
Tighter inner Euler solves resolve the final edge-speed mismatch. This path
can take more iterations than a direct solve. For NLR 7301 with flap:

```sh
native/alula-euler --topology --potential-grid --coarse-start 8x4 --refine-coupled \
  --coupled 2510000 9 6 .185 10 4 < nlr7301-contours.txt
```

```sh
native/alula-euler --topology --bl 1000000 9 0 0.1 8 4 < contours.txt
native/alula-euler --topology --coupled 1000000 9 0 0.1 8 4 < contours.txt
```

Use one `--trip UPPER LOWER` pair per contour, in input order, to prescribe
material surface fractions. A fraction of one selects natural transition
with terminal wake initialization. Coupled checks cover coarse subsonic
sharp and finite-base NACA 2412, stacked NACA 0012, the built-in main
element and flap with 20% material trips, and transonic RAE 2822 with trips
at six and eight surface intervals. This example checks convergence on a
startup grid:

```sh
native/alula-euler --topology --potential-grid --coupled 2700000 4 --trip .03 .07 2.68 .74 8 4 < rae-contour.txt
```

With 16 surface intervals and `--coarse-start 8`, the tripped RAE 2822
also converges in coupled mode. If displacement iteration stalls on a
complete interaction profile, the native solver can solve the BL station
equations together. It reuses one sparse Euler Jacobian for the effect of
each station's displacement and groups independent shape and shear changes.
For a local wall displacement, the solver differentiates the cells touching
the moved crosslines and the downstream cells reached by the upwind stencil.
It still evaluates the wake, farfield and other boundary equations over the full grid.
The reduced dense system uses LAPACK when enabled. Each accepted trial
solves the nonlinear Euler equations. The final check includes the original
BL equations, transition selection, displacement, and Euler residual.
A failed recovery restores the previous grid and profiles.

The NLR 7301 with flap also converges with natural transition on a coarse
grid at Mach 0.185, Reynolds number 2.51 million, and 6 degrees.

```sh
native/alula-euler --topology --potential-grid --coupled 2510000 9 6 .185 8 4 < nlr7301-contours.txt
```

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
against the current matrix. If the previous factors produce contracting
corrections, it can finish a linear solve without refactorization. The
Euler JSON report gives the factor, refactor, and reuse counts. Refinement
uses extended-precision residuals when double-precision products cannot
meet the linear tolerance. Where hardware FMA is available and
`long double` uses binary128, refinement can form corrections with
compensated double-precision sums before verifying the solution in binary128.

The body API evaluates the same pressure matching, edge conditions,
farfield modes, wake gaps, and ISMOM choices as the compressible web solver
on supplied physical nodes, densities, and tube masses. It reuses work
arrays across evaluations. Consecutive stencils reuse the shared interval's
validated geometry and physical gas state. Boundary-layer marching and
pressure integration reuse the prepared stencils while the decoded Euler
state matches exactly. During line search, a previously failed cell can
reject a trial if it fails the same initialization check. The OpenMP build
parallelizes cell initialization and flux evaluation on larger grids.
Use `ALULA_BODY_BENCH=1 make -C native check`
to time the residual parity cases.

`alula_body_linearize` accepts a direction in the same physical variables
and returns the residual and its directional derivative. It caches the
shared upwind stencil coefficients for each linearization and propagates
farfield normalization and independent wake widths. Callers
using packed coordinates supply the node and gas directions from their
coordinate chart.

`body_newton.h` accepts a native coordinate-chart callback and its derivative
pattern. It assembles the square sparse Newton matrix with graph coloring
and reuses the KLU ordering across accepted steps. An accepted trial's
residual and prepared cell data carry into the next Jacobian assembly
when the chart decodes to the same physical state. The colors in each
Jacobian assembly share the base farfield samples. Each trial checks the
physical Euler equations and the returned root must pass a convex-grid check.
The tests converge single and multi-element grids on a fixed chart and
compare sparse matrix products for rounded bodies and independent wakes.

`body_chart.h` builds the coordinates from a prepared grid, contour splines,
and tube mass weights. It supports moving rounded leading-edge stagnation
points and prescribed displacement, and updates grid normals after accepted
Newton steps while preserving the accepted physical geometry.
