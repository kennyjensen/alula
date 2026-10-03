// SPDX-License-Identifier: GPL-2.0-or-later
// Validation-only KLU with explicit P/Q. Uses an isolated WASM build; the
// original CSR and RHS are never permuted or modified by JavaScript.
const url = new URL('../../docs/linear-ordering-prototype/klu-given.wasm', import.meta.url);
const bytes = await (await import('node:fs/promises')).readFile(url);
const module = await WebAssembly.compile(bytes);
if (WebAssembly.Module.imports(module).length) throw new Error('Given-order KLU must not require host imports.');
const wasm = (await WebAssembly.instantiate(module, {})).exports;
wasm._initialize();

function norm(values) {
  let scale = 0, sum = 0;
  for (const value of values) {
    if (!Number.isFinite(value)) return Infinity;
    const a = Math.abs(value);
    if (a > scale) { sum = 1 + sum * (scale / a) ** 2; scale = a; }
    else if (a) sum += (a / scale) ** 2;
  }
  return scale * Math.sqrt(sum);
}
function permutation(input, n, name) {
  if (!(Array.isArray(input) || ArrayBuffer.isView(input)) || input.length !== n)
    throw new Error(`${name} must be a bijection of 0..n-1.`);
  const seen = new Uint8Array(n), result = new Int32Array(n);
  for (let k = 0; k < n; k++) {
    const i = input[k];
    if (!Number.isInteger(i) || i < 0 || i >= n || seen[i]) throw new Error(`${name} must be a bijection of 0..n-1.`);
    seen[i] = 1; result[k] = i;
  }
  return result;
}
// Same original-CSR conversion convention as the production adapter: structural
// zeros are omitted from numerical CSC, while residual checks retain all entries.
function compressedColumns(a) {
  const { n, rowPtr, colIndex, values } = a;
  if (!Number.isInteger(n) || n < 1 || n >= 2 ** 31 || !rowPtr || !colIndex || !values
    || rowPtr.length !== n + 1 || rowPtr[0] !== 0 || rowPtr[n] !== values.length || values.length !== colIndex.length)
    throw new Error('Invalid sparse matrix dimensions.');
  const ap = new Int32Array(n + 1);
  for (let i = 0; i < n; i++) {
    if (!Number.isInteger(rowPtr[i + 1]) || rowPtr[i + 1] < rowPtr[i]) throw new Error('Invalid sparse row pointers.');
    let previous = -1;
    for (let k = rowPtr[i]; k < rowPtr[i + 1]; k++) {
      const j = colIndex[k];
      if (!Number.isInteger(j) || j <= previous || j >= n || !Number.isFinite(values[k])) throw new Error('Invalid sparse matrix entry.');
      previous = j; if (values[k] !== 0) ap[j + 1]++;
    }
  }
  for (let j = 0; j < n; j++) ap[j + 1] += ap[j];
  const ai = new Int32Array(ap[n]), ax = new Float64Array(ap[n]), next = ap.slice();
  for (let i = 0; i < n; i++) for (let k = rowPtr[i]; k < rowPtr[i + 1]; k++) if (values[k] !== 0) {
    const p = next[colIndex[k]]++; ai[p] = i; ax[p] = values[k];
  }
  return { ap, ai, ax };
}

export function solveSparseDirectGiven(a, b, options = {}) {
  const { rowPermutation, columnPermutation, pivotTolerance = 1, tolerance = 1e-10, maxRefinements = 4, btf = false } = options;
  const allowed = ['rowPermutation', 'columnPermutation', 'pivotTolerance', 'tolerance', 'maxRefinements', 'btf'];
  if (Object.keys(options).some(k => !allowed.includes(k)) || !a || !(Array.isArray(b) || ArrayBuffer.isView(b))
    || b.length !== a.n || !b.every(Number.isFinite) || !Number.isFinite(tolerance) || tolerance <= 0
    || !Number.isFinite(pivotTolerance) || pivotTolerance <= 0 || pivotTolerance > 1 || typeof btf !== 'boolean'
    || !Number.isInteger(maxRefinements) || maxRefinements < 0 || maxRefinements > 4)
    throw new Error('Invalid given-order sparse solve controls or right-hand side.');
  const { ap, ai, ax } = compressedColumns(a);
  const p = permutation(rowPermutation, a.n, 'rowPermutation'), q = permutation(columnPermutation, a.n, 'columnPermutation');
  const bnorm = norm(b), allocations = []; let factor = 0;
  const copy = data => {
    const pointer = wasm.malloc(Math.max(1, data.byteLength));
    if (!pointer) throw new Error('WebAssembly given-order sparse solver ran out of memory.');
    allocations.push(pointer);
    new data.constructor(wasm.memory.buffer, pointer, data.length).set(data); return pointer;
  };
  const attempts = [];
  try {
    const apPointer = copy(ap), aiPointer = copy(ai), axPointer = copy(ax), pPointer = copy(p), qPointer = copy(q), rhsPointer = copy(Float64Array.from(b));
    factor = wasm.mses_klu_factor_given(a.n, apPointer, aiPointer, axPointer, pPointer, qPointer, pivotTolerance, btf ? 1 : 0);
    if (!factor) {
      const status = wasm.mses_klu_status(); attempts.push({ ordering: 'given', pivotTolerance, btf, status });
      const error = new Error(`Given-order sparse LU factorization failed (KLU status ${status}).`);
      error.code = 'KLU_GIVEN_FACTORIZATION'; error.status = status; error.attempts = attempts; throw error;
    }
    const actualBtf = wasm.mses_klu_btf(factor) !== 0, symbolicBlocks = wasm.mses_klu_blocks(factor);
    if (actualBtf !== btf) throw new Error('KLU did not retain the requested BTF control.');
    const solve = rhs => {
      // Factorization/solves may grow WASM memory; always obtain a fresh view.
      new Float64Array(wasm.memory.buffer, rhsPointer, a.n).set(rhs);
      if (!wasm.mses_klu_solve(factor, rhsPointer)) {
        const error = new Error(`Given-order sparse LU solve failed (KLU status ${wasm.mses_klu_status()}).`);
        error.code = 'KLU_GIVEN_SOLVE'; error.attempts = attempts; throw error;
      }
      return new Float64Array(wasm.memory.buffer, rhsPointer, a.n).slice();
    };
    const x = solve(b); let relativeResidual = Infinity, refinements = 0;
    for (; refinements <= maxRefinements; refinements++) {
      const residual = Float64Array.from(b);
      for (let i = 0; i < a.n; i++) for (let k = a.rowPtr[i]; k < a.rowPtr[i + 1]; k++) residual[i] -= a.values[k] * x[a.colIndex[k]];
      const rnorm = norm(residual); relativeResidual = bnorm ? rnorm / bnorm : rnorm;
      if (x.every(Number.isFinite) && relativeResidual <= tolerance) break;
      if (!Number.isFinite(relativeResidual) || refinements === maxRefinements) break;
      const correction = solve(residual); for (let i = 0; i < a.n; i++) x[i] += correction[i];
    }
    const factorNonzeros = wasm.mses_klu_nnz(factor);
    attempts.push({ ordering: 'given', pivotTolerance, btf: actualBtf, relativeResidual, refinements, factorNonzeros, symbolicBlocks });
    if (!x.every(Number.isFinite) || relativeResidual > tolerance) {
      const error = new Error(`Given-order sparse LU missed the requested relative residual (${relativeResidual} > ${tolerance}).`);
      error.code = 'KLU_RESIDUAL_LIMIT'; error.relativeResidual = relativeResidual; error.attempts = attempts; throw error;
    }
    return { x, relativeResidual, refinements, factorNonzeros, symbolicBlocks, backend: 'klu-given-wasm', ordering: 'given',
      pivotTolerance, btf: actualBtf, attempts, permutationConvention: 'new-position-to-original-index; Aordered[k,l]=A[P[k],Q[l]]',
      factorizationCalls: 1, originalSystemResidual: true };
  } finally {
    if (factor) wasm.mses_klu_free(factor);
    for (const pointer of allocations) wasm.free(pointer);
  }
}
