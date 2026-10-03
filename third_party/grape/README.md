# GRAPE reference source

These four files are unmodified selections from the GRAPE archive distributed
by [Public Domain Aeronautical Software](https://www.pdas.com/grape.html).
The program was developed by NASA Ames and released through COSMIC as
ARC-11379. Original notices and the distributor's readme are retained.
[PROVENANCE.json](PROVENANCE.json) records the archive URL, retrieval date
and exact file hashes.

- `original.src`: original COSMIC source, including CDC-specific statements.
- `grape.f90`: the distributor's modern Fortran version.
- `arc11379.txt`: original COSMIC program description.
- `readme.txt`: distributor's build and usage information.

The relevant original `RELAX` loop is at lines 3123–3131. It uses a forward
first difference for positive P+R and a backward first difference for negative
P+R, adding the magnitude to the negative diagonal. This is the direct source
for `src/geometry/grape-poisson-drift.js`. The original has uniform index
increments; the JavaScript extension uses the corresponding left/right
interval on a nonuniform coordinate. Its source derivative is first order.

`node scripts/generate-grape-drift-reference.js` extracts that loop without
editing it, compiles a small Float64 Fortran harness, and saves 22 coefficient
cases. It uses `$FC`, an available local reference toolchain, or `gfortran`.
This is a stencil reference, not execution of the complete GRAPE program.

`node scripts/generate-grape-boundary-reference.js` additionally extracts the
unmodified inner/outer boundary P/Q factor and source blocks. Forty Float64
cases check the Q=0 smooth-boundary specialization against native equations.
Nonuniform rows are independently resampled from a cubic Hermite interpolant
before calling the original uniform-step routine. This does not execute
the complete GRAPE iteration. [Scope and results](../../docs/GRAPE_BOUNDARY_LINEARIZATION.md).

GRAPE is a geometric Poisson grid generator and is not ISET, MSET or MSES.
The browser runtime remains JavaScript/WASM. Our harmonic-mass restriction,
angle-only boundary reconstruction, source extension and iteration differ
from the complete original algorithm; see
[the validation record](../../docs/GRAPE_SOURCE_STENCIL.md).
