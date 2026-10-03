# Grid smoothing kernel

`smoothing.ts` is the AssemblyScript/Float64 implementation of the paired
boundary Newton-SLOR and fixed-boundary scalar Giles SLOR kernels. The web
app loads the committed `smoothing.wasm`; no compiler runs in the browser.

Rebuild with `npm ci` followed by `npm run build:smoothing`. The pinned
AssemblyScript compiler runs under Node, including on ARM64. `BUILD.json`
records compiler options and source/binary hashes; a test prevents stale
binaries from passing validation. The compiler runtime's MIT license is
included in `ASSEMBLYSCRIPT-LICENSE`.

The adapter in `../wasm-smoothing.js` packs coordinates, prescribed controls,
and difference stencils once into each instance. All line assembly, boundary
reconstruction, metric evaluation, row-pivoted block Thomas solves and paired
Armijo/convexity trials run inside WASM. Buffers are reused between sweeps.
There are no host callbacks in the numerical loops. Progress observers remain
JavaScript, called between sweeps; instances are independent for nested solves.

The original JS equations remain in `../elliptic-streamtube-grid.js`.
Both smoothing drivers accept `{ backend: 'javascript' }` for reference tests.
WASM is the default for paired fixed-boundary full-metric Giles smoothing and
fixed-boundary scalar Giles smoothing (harmonic, fixed F, or fixed P). Moving
farfield curves, quadratic discretization, mapped controls, and single-row
implicit angle variants retain their JS implementations. Scalar whole-grid
convex backtracking remains in the JS driver; its expensive sweeps are WASM.

A failed WASM sweep is replayed from the original input in JS to preserve
existing detailed errors and typed termination. Partially updated WASM nodes
are never committed on failure. Results expose non-enumerable `backend` and
`referenceReplays` fields; parity tests require zero reference replays.
Ordinary runtime/bounds errors are not swallowed.

The tolerance and line search are unchanged. AssemblyScript and V8 `hypot`
can differ by a few ulps in reported cell quality; test physical coordinates,
residuals, accepted step fractions and convergence, rather than requiring
bitwise equality of those diagnostic angles across runtimes.

Validation commands:

```sh
node --test tests/wasm-smoothing.test.js
npx playwright test tests/browser/wasm-smoothing.spec.js
node scripts/validation/benchmark-smoothing-kernel.js 128 24 20
node scripts/validation/benchmark-grid-smoothing.js 64 11
```
