import test from 'node:test';
import assert from 'node:assert/strict';
import { solveStreamwiseScalarLines } from '../src/numerics/streamwise-scalar-lines.js';
import { sparseMatrix, sparseAdd, sparseDense, sparseProduct } from '../src/numerics/sparse.js';
import { solveLinear } from '../src/numerics/linear.js';

const close = (a, b, tolerance = 2e-11) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);
const makeSystem = (nx = 5, nt = 4) => {
  const n = (nx - 1) * (nt - 1), id = (i, j) => (i - 1) * (nt - 1) + j - 1;
  const matrix = sparseMatrix(Array.from({ length: n }, (_, row) => {
    const i = 1 + Math.floor(row / (nt - 1)), j = 1 + row % (nt - 1), cols = [];
    for (let a = Math.max(1, i - 1); a <= Math.min(nx - 1, i + 1); a++)
      for (let b = Math.max(1, j - 1); b <= Math.min(nt - 1, j + 1); b++) cols.push(id(a, b));
    return cols;
  }));
  // Symmetric positive definite by strict diagonal dominance; includes
  // mixed-derivative-like diagonal-neighbor couplings of both signs.
  for (let row = 0; row < n; row++) {
    let sum = 0;
    for (let p = matrix.rowPtr[row]; p < matrix.rowPtr[row + 1]; p++) {
      const col = matrix.colIndex[p];
      if (col !== row) {
        const value = (row + col) % 3 === 0 ? .17 : -.31;
        matrix.values[p] = value; sum += Math.abs(value);
      }
    }
    sparseAdd(matrix, row, row, .8 + .03 * row + sum);
  }
  const rhs = Float64Array.from({ length: n }, (_, k) => Math.sin(.7 * (k + 1)) + .1 * k);
  return { matrix, rhs, nx, nt };
};

test('scalar streamwise SLOR agrees with independent dense solve and original residual', () => {
  const { matrix, rhs, nx, nt } = makeSystem(), before = structuredClone({ matrix, rhs });
  const expected = solveLinear(sparseDense(matrix), rhs), result = solveStreamwiseScalarLines(matrix, rhs, { nx, nt, tolerance: 1e-12 });
  assert.equal(result.converged, true); assert.ok(result.sweeps > 1 && result.sweeps < 100);
  result.x.forEach((value, i) => close(value, expected[i]));
  const residual = sparseProduct(matrix, result.x).map((value, i) => value - rhs[i]);
  const norm = a => Math.sqrt(a.reduce((sum, value) => sum + value * value, 0));
  close(result.relativeResidual, norm(residual) / norm(rhs), 1e-15);
  assert.deepEqual({ matrix, rhs }, before);
  for (const scale of [1e-120, 1e120]) {
    const scaled = { ...matrix, values: matrix.values.map(value => scale * value) };
    const actual = solveStreamwiseScalarLines(scaled, rhs.map(value => scale * value), { nx, nt, tolerance: 1e-12 });
    assert.equal(actual.converged, true); actual.x.forEach((value, i) => close(value, expected[i]));
  }
});

test('one complete overrelaxed sweep agrees with independently assembled dense line solves', () => {
  const { matrix, rhs, nx, nt } = makeSystem(), dense = sparseDense(matrix), expected = new Float64Array(matrix.n), omega = 1.47;
  for (let j = 1; j < nt; j++) {
    const ids = Array.from({ length: nx - 1 }, (_, k) => k * (nt - 1) + j - 1);
    const a = ids.flatMap(row => ids.map(col => dense[row * matrix.n + col]));
    const b = ids.map(row => rhs[row] - expected.reduce((sum, value, col) => sum + dense[row * matrix.n + col] * value, 0));
    const delta = solveLinear(a, b);
    ids.forEach((row, k) => { expected[row] += omega * delta[k]; });
  }
  const actual = solveStreamwiseScalarLines(matrix, rhs, { nx, nt, omega, maxSweeps: 1, tolerance: 1e-14 });
  assert.equal(actual.sweeps, 1); assert.equal(actual.converged, false);
  actual.x.forEach((value, k) => close(value, expected[k], 1e-14));
});

test('scalar lines report limits, zero right-hand sides and indefinite failures without repairing matrices', () => {
  const { matrix, rhs, nx, nt } = makeSystem();
  const stopped = solveStreamwiseScalarLines(matrix, rhs, { nx, nt, maxSweeps: 0 });
  assert.equal(stopped.converged, false); assert.equal(stopped.sweeps, 0); assert.equal(stopped.relativeResidual, 1);
  assert.ok(stopped.x.every(value => value === 0));
  const zero = solveStreamwiseScalarLines(matrix, new Float64Array(matrix.n), { nx, nt });
  assert.equal(zero.converged, true); assert.equal(zero.sweeps, 0); assert.equal(zero.relativeResidual, 0);
  // Each eta line has a positive pivot, but the full matrix is indefinite.
  // The original residual must expose failure instead of certifying updates.
  const indefinite = sparseMatrix([[0, 1], [0, 1]]);
  indefinite.values.set([1, 2, 2, 1]); const before = structuredClone(indefinite);
  const failed = solveStreamwiseScalarLines(indefinite, [1, 0], { nx: 2, nt: 3, omega: 1, maxSweeps: 5 });
  assert.equal(failed.converged, false); assert.ok(failed.relativeResidual > 1);
  assert.deepEqual(indefinite, before);
  const negative = sparseMatrix([[0]]); negative.values[0] = -1;
  assert.throws(() => solveStreamwiseScalarLines(negative, [1], { nx: 2, nt: 2 }), /pivot/);
  const singular = sparseMatrix([[0, 1], [0, 1]]); singular.values.fill(1);
  assert.throws(() => solveStreamwiseScalarLines(singular, [1, 1], { nx: 3, nt: 2 }), /pivot/);
});

test('scalar lines reject invalid controls, nonfinite entries, malformed storage and nontridiagonal lines', () => {
  const { matrix, rhs, nx, nt } = makeSystem();
  for (const controls of [{ nx: 1, nt }, { nx, nt: 1 }, { nx, nt, omega: 2 }, { nx, nt, omega: NaN },
    { nx, nt, maxSweeps: -1 }, { nx, nt, tolerance: 0 }, { nx, nt, tolerance: Infinity }])
    assert.throws(() => solveStreamwiseScalarLines(matrix, rhs, controls), /Invalid/);
  const nonfinite = structuredClone(matrix); nonfinite.values[0] = NaN;
  assert.throws(() => solveStreamwiseScalarLines(nonfinite, rhs, { nx, nt }), /finite/);
  const invalidRhs = rhs.slice(); invalidRhs[0] = Infinity;
  assert.throws(() => solveStreamwiseScalarLines(matrix, invalidRhs, { nx, nt }), /right-hand/);
  const malformed = structuredClone(matrix); malformed.colIndex[0] = -1;
  assert.throws(() => solveStreamwiseScalarLines(malformed, rhs, { nx, nt }), /ordered columns/);
  assert.throws(() => solveStreamwiseScalarLines({ n: matrix.n, rowPtr: matrix.rowPtr }, rhs, { nx, nt }), /sparse storage/);
  const distant = sparseMatrix([[0, 1, 2], [0, 1, 2], [0, 1, 2]]);
  for (let i = 0; i < 3; i++) sparseAdd(distant, i, i, 2);
  sparseAdd(distant, 0, 2, .1);
  assert.throws(() => solveStreamwiseScalarLines(distant, [1, 1, 1], { nx: 4, nt: 2 }), /tridiagonal/);
});
