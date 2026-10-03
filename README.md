# alula

A multi-element and transonic airfoil solver inspired by
[MSES](https://web.mit.edu/drela/Public/web/mses/). Alula runs aerodynamic
analysis in the browser, with interactive geometry, flow visualization and
integral boundary-layer calculations.

[Open Alula](https://alula.vibefoil.com) ·
[Report a bug or request a feature](https://github.com/kennyjensen/alula/issues)

## Main features

- **Build airfoil assemblies.** Start with a single element, main element and
  flap, or slat–main–flap configuration. Edit element geometry and placement,
  add elements, or import DAT and MSES blade files. RAE 2822 and NLR 7301 with
  flap are included as benchmark presets.
- **Choose a flow model.** Use incompressible panel flow, panel flow with
  boundary layers, compressible Euler flow, or coupled Euler and boundary
  layers. The coupled Euler model is the default.
- **Control the analysis.** Set angle of attack, Mach number, Reynolds number,
  transition amplification threshold (Ncrit), and upper/lower trip limits.
  Euler modes provide surface, inlet, wake and streamtube resolution controls
  with optional initial grid smoothing.
- **Inspect the solution.** View pressure coefficient, speed-colored
  streamlines, Mach contours, mesh geometry and boundary-layer development.
  The airfoil overlay shows conditions, lift, drag, pitching moment and L/D;
  the Solution Status panel tracks convergence and numerical diagnostics.
- **Export and debug.** Export coordinates and analysis results, or copy a
  diagnostic record containing the case settings and solver failure details.
  The interface uses swipeable Settings and Airfoil Vis tabs on phones and
  two columns on desktop.

## Solvers

| Mode | Main calculation | Entry point |
| --- | --- | --- |
| Panel flow | Incompressible linear-vortex panel solution for interacting elements | [`solveInviscid`](src/inviscid/linear-vortex.js) |
| Panel flow + boundary layers | Incompressible viscous interaction with integral boundary layers and element wakes | [`solveViscousAssembly`](src/viscous/result.js) |
| Quad Euler | Compressible flow on a moving quadrilateral streamtube grid | [`solveStreamtubeAssembly`](src/euler/streamtube-result.js) |
| Quad Euler + boundary layers | Coupled Euler, boundary-layer, displacement and wake equations, with fixed-trip or automatic eᴺ transition | [`solveCoupledStreamtubeAssembly`](src/euler/streamtube-coupled-assembly.js) |

The browser sends cases to a Web Worker so calculations do not block the
interface. The worker reports mesh, pressure, coefficients and convergence
progress as they become available. Numerical code uses JavaScript double
precision, with WebAssembly helpers including SuiteSparse KLU for sparse
linear solves. Geometry and flow calculations run locally; no solver backend
or geometry upload is required.

Alula is under active development. Transonic and multi-element viscous results
remain research estimates; numerical convergence alone does not establish
physical accuracy. Check mesh sensitivity and suitable reference data before
relying on a prediction. Alula is an independent project, not an official MSES
release.

## Codebase guide

| Location | Responsibility |
| --- | --- |
| [`src/ui/`](src/ui/) | Application controls, plots, mobile layout, result display and exports |
| [`src/worker/solver.js`](src/worker/solver.js) | Solver dispatch and progress/result messages between numerical code and the interface |
| [`src/geometry/`](src/geometry/) | Airfoil generation and parsing, contour validation, benchmark coordinates, streamtube meshing and smoothing |
| [`src/inviscid/`](src/inviscid/) | Panel influence calculations, inviscid flow, displacement interaction and streamline tracing |
| [`src/viscous/`](src/viscous/) | Integral boundary layers, transition, wakes, viscous coupling and XFOIL-derived numerical routines |
| [`src/euler/`](src/euler/) | Streamtube Euler equations, Jacobians, nonlinear iteration, coupled startup/recovery and refinement; also a separate finite-volume reference solver |
| [`src/numerics/`](src/numerics/) | Linear algebra, sparse systems and nonlinear-solver utilities |
| [`tests/`](tests/) and [`scripts/`](scripts/) | Unit and browser tests, numerical validation, profiling, development server and static build |
| [`third_party/`](third_party/) | Upstream code, reference material, benchmark source data, licenses and provenance |

For the main application flow, start with [`src/ui/app.js`](src/ui/app.js),
then [`src/worker/solver.js`](src/worker/solver.js) and the solver entry points
above. [`euler.html`](euler.html) is a separate finite-volume verification
interface.

## Run locally

Requires Node.js 22 or newer and npm.

```sh
npm ci
npm run dev
```

Open [localhost:5173](http://localhost:5173). To create the static deployment:

```sh
npm run build
```

The build writes the application, runtime assets and required third-party
notices to `dist/`. It can be served by a static web host.

## Tests and validation

```sh
npm test                    # Node.js test suite
npm run validate            # Numerical validation scripts
npx playwright install chromium
npm run test:browser        # Browser tests
```

Some regression tests and investigation scripts require historical checkpoint
fixtures from a separate local `docs/` archive. That archive and generated
results are ignored by Git and excluded from the build, so the full test suite
is not self-contained in a fresh checkout. CI runs a selected self-contained
test set, geometry checks and browser solves; test sources and fixtures under
`tests/fixtures/` remain versioned.

## License and attribution

Alula's project code is licensed under [GPL-2.0-or-later](LICENSE). Third-party
code and coordinate data retain their own licenses and notices. See
[NOTICE.md](NOTICE.md) and [third_party/README.md](third_party/README.md) for
attribution, provenance and component-specific terms.
