# SuiteSparse KLU in WebAssembly

`klu.wasm` is the real, 32-bit-index KLU sparse LU solver from SuiteSparse
v5.13.0, with AMD/COLAMD ordering and BTF. It uses IEEE double precision and
row scaling with numerical pivoting. No fast-math flags are used. The module
has no host imports and runs in a browser worker or Node without a server,
filesystem access, Emscripten runtime, or network access during a solve.

The selected upstream C sources and headers are included unchanged in
`source/`. Their hashes and the original archive URL/hash are in
`PROVENANCE.json`. `bridge.c` is the project's small ABI wrapper; the flow
equations and convergence checks remain JavaScript. `BUILD.json` records the
compiler, flags, wrapper, script, and generated module hash.

The JavaScript adapter certifies each answer against the original CSR matrix
and right-hand side. Its default `ordering: 'auto'` tries AMD first and can
retry COLAMD if the certificate or a singular factorization fails. A caller
can choose `preferredOrdering` or force one ordering. If these attempts fail,
the adapter retries with maximum-magnitude partial pivoting (`Common.tol=1`)
instead of KLU's default diagonal preference (`0.001`). The same original
matrix/RHS and residual limit apply to every attempt. No diagonal shift or
regularization is used. `pivotTolerance` selects an explicit threshold;
`pivotFallback: false` disables recovery for diagnostic controls. The new
bridge entry point accepts the threshold; the original ABI remains available.

Each attempt has bounded iterative refinement; failed factors are freed
before retrying. Plain coupled Newton remembers the successful ordering;
the ISES path remembers the successful ordering/pivot threshold and records
each attempt and residual. Accepted checkpoints preserve these preferences. Neither
path changes flow equations or convergence tolerances. Recovery can cost
additional factorizations and memory; it runs only when normal attempts fail.

KLU and BTF are copyright (C) 2004–2013 University of Florida, by Timothy
A. Davis and Ekanathan Palamadai, under LGPL-2.1-or-later. AMD and COLAMD use
BSD 3-clause licenses. The notices are in the corresponding
`source/*/Doc/License.txt` files; the full LGPL 2.1 text is in
`COPYING.LESSER`. SuiteSparse_config has no
licensing restrictions, as stated in its source. Notices for the linked
WASI libc/libm and LLVM compiler builtins are included alongside this file.
The wrapper and build script are GPL-2.0-or-later.

To rebuild or replace the library, install a WASI-capable clang/LLVM,
wasm-ld, WASI libc, and wasm32 compiler builtins, then run from the project
root:

```sh
WASM_CC=clang WASI_SYSROOT=/path/to/wasi/sysroot npm run build:klu
node --test tests/klu.test.js
```

Put the matching wasm-ld executable on `PATH`. The shipped module was built
with clang 19.1.7 and Ubuntu's wasi-libc
`0.0~git20250726.3f7eb4c-4`; a different compiler/libc can produce a different
binary hash. Run the numerical tests after rebuilding. The build writes a
temporary file and replaces the module only after successful compilation
and verification of the upstream source hashes and absence of host imports.

The static application build includes this entire directory: the complete
library source needed to rebuild, the wrapper, build script, provenance,
and license notices remain available with every shipped binary. Replace
`klu.wasm` with a rebuilt ABI-compatible module to use a modified library.
