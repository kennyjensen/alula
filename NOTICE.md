# Third-party notices

alula is an independent implementation project. It is not an official MSES release.

New project code is offered under GPL-2.0-or-later. The full GPL version 2 text is in [LICENSE](LICENSE). There is no warranty.

The finite-volume Euler modules added in version 0.3 are new implementations of perfect-gas conservation and HLLC equations. They do not contain MSES source code or incorporate another CFD implementation.

The browser meshing bundle `third_party/cdt2d/cdt2d.js` contains cdt2d by
Mikola Lysenko and its dependencies. Their individual licenses are preserved
in `third_party/cdt2d/cdt2d-LICENSES.txt`; source versions and hashes are
recorded in `third_party/cdt2d/cdt2d-PROVENANCE.json`. These MIT notices
remain applicable to the bundled code.

The subcritical coupled solver uses SuiteSparse KLU in WebAssembly for
pivoted sparse LU. KLU and BTF are copyright (C) 2004–2013 University of
Florida, by Timothy A. Davis and Ekanathan Palamadai, LGPL-2.1-or-later;
AMD and COLAMD use BSD 3-clause licenses. The unmodified selected sources,
complete license texts, ABI wrapper, build instructions, linked runtime
notices, and provenance are included in
[`third_party/klu/`](third_party/klu/README.md), including in
the static distribution. The flow equations remain JavaScript.

`src/boundary-layer/tests/closures.js` translates the HKIN, HSL, CFL and DIL routines from the supplied `third_party/Xfoil/src/xblsys.f`, copyright (C) 2000 Mark Drela, GPL v2 or later. The source header records the translation. The earlier panel and dense numerical kernels are new implementations.

The supplied Vibefoil repository is credited to Kenny Jensen by its HTML metadata. Its bundled license is GPL version 2, and its JavaScript XFOIL port carries GPL v2-or-later notices. The interface takes visual inspiration from Vibefoil. Its `favicon.ico` is copied unchanged to `public/favicon.ico` (SHA-256 `8aac731ffd412e0ff42b87764c01efaa39f617b0110f6a85a5b6a1af52f9c21f`). Version 0.2 vendors ten of its XFOIL numerical modules into `src/viscous/xfoil/`, preserving their GPL notices. These modules run in the production worker; its app and analytics are not loaded. `src/viscous/context.js` also adapts the panel allocation in its solver worker. Original-file hashes and local modifications are recorded in `src/viscous/xfoil/PROVENANCE.json`. The modified numerical modules derive from XFOIL, copyright (C) 2000 Mark Drela. Tests additionally import its original laminar closures as a secondary translation cross-check. The bundled license text supplied with Vibefoil is reproduced in `LICENSE`.

The supplied MSES User's Guide 3.05 by Mark Drela (July 2007) is a technical reference. It is not MSES source code. The GPL designation for new project code does not relicense this PDF or other third-party materials. Original files under `third_party/` are kept as supplied and retain their own notices. The static build includes the ISES and GRAPE reference directories; it excludes the other supplied reference collections, except the benchmark coordinate files and attribution documents listed below.

The selected unmodified GRAPE files in [`third_party/grape/`](third_party/grape/README.md)
come from the Public Domain Aeronautical Software distribution of NASA Ames
COSMIC ARC-11379. The original description, distributor's readme and file hashes
are retained. `src/geometry/grape-poisson-drift.js` ports its sign-biased Poisson
drift stencil, extending uniform computational intervals to nonuniform ones.
The browser does not execute the Fortran reference. Native validation covers
the extracted coefficient loop and smooth-boundary P/Q source blocks, including
a documented Q=0 specialization. It does not establish full GRAPE or MSET parity.

Before distributing the supplied reference collection, review the notices in those files separately. The browser application and its new sources can be built with `npm run build`.

## Benchmark coordinate data

The GPL license for the application does not relicense coordinate data. Source
receipts, attribution and available rights information are retained with each dataset:

| Dataset | Attribution | Source documentation |
| --- | --- | --- |
| RAE 2822, MSES sample | Mark Drela, MIT; public `blade.rae` sample | [README](third_party/airfoils/rae2822-mses/README.md) |
| RAE 2822, measured coordinates | Cook, McDonald and Firmin, AGARD AR 138 (1979), distributed through NASA's NPARC archive | [README](third_party/airfoils/rae2822/README.md) |
| NLR 7301 | SAAB Military Aircraft / FLOWNET, via ERCOFTAC UFR 3-01; case description by Jan Vos, CFS Engineering | [README](third_party/airfoils/nlr7301/README.md) |
| 30P30N | UIUC Applied Aerodynamics Group coordinate database; contributed by Brent Pomeroy (2016) | [README](third_party/airfoils/30p30n/README.md) |

These upstream coordinate downloads do not include explicit data-license grants;
credit alone does not establish permission to relicense them. The individual
READMEs distinguish coordinate sources from separately supplied papers and CAD files.
The static distribution includes the coordinate files used by the app and their
attribution documents, but not the supplied papers, grid archives or CAD files.

## Directory convention

Keep upstream originals and reference material under `third_party/`. Runtime
ports and vendor bundles may remain under `src/.../vendor/` or `src/viscous/xfoil/`
so their imports and build assets stay together. Each component must retain its
license, copyright notices, source/version provenance and a record of local
modifications, and be listed here. See [third_party/README.md](third_party/README.md).
