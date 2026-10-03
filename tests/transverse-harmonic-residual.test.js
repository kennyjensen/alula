import test from 'node:test';
import assert from 'node:assert/strict';
import { createTransverseHarmonicResidual } from '../src/geometry/transverse-harmonic-residual.js';
import { quadLaplaceMatrix } from '../src/geometry/harmonic-grid-audit.js';
import { sparseDense, sparseProduct } from '../src/numerics/sparse.js';

const maximum = a => Math.max(...a.map(Math.abs));
const close = (a, b, tolerance = 2e-12) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);
const normalize = d => { const length = Math.hypot(d.x, d.y); return { x: d.x / length, y: d.y / length }; };
const fixture = (affine = false) => {
  const x = [0, .03, .25, .7, 1], eta = [0, .13, .48, 1], nx = x.length - 1, nt = eta.length - 1;
  const nodes = x.map((u, i) => eta.map((v, j) => {
    const px = u + .15 * v + (affine ? 0 : .04 * Math.sin(Math.PI * u) * Math.sin(Math.PI * v));
    return { x: px, y: .2 * px + v + (affine ? 0 : .035 * Math.sin(2 * u + .3 * v)) };
  }));
  const directions = nodes.map((row, i) => row.map((_, j) => ({ x: .15 * Math.sin(i + j), y: .9 + .1 * Math.cos(i - j) })));
  return { nodes, directions, massFlows: eta.slice(1).map((v, j) => v - eta[j]), eta, nx, nt };
};
// Assemble separately with the existing independent physical-space element
// stiffness, including fixed-boundary contributions in the residual.
const stiffness = ({ nodes, eta, nx, nt }) => {
  const n = (nx - 1) * (nt - 1), residual = new Float64Array(n), matrix = new Float64Array(n * n);
  const id = (i, j) => !i || i === nx || !j || j === nt ? -1 : (i - 1) * (nt - 1) + j - 1;
  for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) {
    const ij = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]], k = quadLaplaceMatrix(ij.map(([a, b]) => nodes[a][b]));
    for (let a = 0; a < 4; a++) {
      const row = id(...ij[a]); if (row < 0) continue;
      for (let b = 0; b < 4; b++) {
        residual[row] += k[4 * a + b] * (eta[ij[b][1]] - eta[ij[a][1]]);
        const col = id(...ij[b]); if (col >= 0) matrix[row * n + col] += k[4 * a + b];
      }
    }
  }
  return { residual, matrix };
};

test('inverse harmonic residual matches separate physical-space Q1 stiffness assembly', () => {
  const input = fixture(), before = structuredClone(input), system = createTransverseHarmonicResidual(input);
  const actual = system.evaluate(input.nodes), expected = stiffness(input);
  assert.equal(system.n, (input.nx - 1) * (input.nt - 1));
  actual.residual.forEach((v, k) => close(v, expected.residual[k], 2e-14));
  const withoutJacobian = system.evaluate(input.nodes, { linearize: false });
  assert.deepEqual(withoutJacobian.residual, actual.residual); assert.equal(withoutJacobian.matrix, undefined);
  for (let row = 0; row < system.n; row++) {
    let sum = 0;
    for (let p = actual.matrix.rowPtr[row]; p < actual.matrix.rowPtr[row + 1]; p++) sum += Math.abs(actual.matrix.values[p]);
    close(actual.rowScale[row], sum);
  }
  assert.deepEqual(input, before);
});

test('affine physical flow gives zero residual and exact -K diag(grad eta dot direction) Jacobian', () => {
  const input = fixture(true), system = createTransverseHarmonicResidual(input), actual = system.evaluate(input.nodes);
  const expected = stiffness(input), dense = sparseDense(actual.matrix), projection = [];
  for (let i = 1; i < input.nx; i++) for (let j = 1; j < input.nt; j++) {
    const d = normalize(input.directions[i][j]); projection.push(-.2 * d.x + d.y);
  }
  assert.ok(maximum(actual.residual) < 3e-14);
  for (let row = 0; row < system.n; row++) for (let col = 0; col < system.n; col++)
    close(dense[row * system.n + col], -expected.matrix[row * system.n + col] * projection[col], 3e-12);
});

test('generic analytic Jacobian agrees with directional differences at second order', () => {
  const input = fixture(), system = createTransverseHarmonicResidual(input), actual = system.evaluate(input.nodes);
  const delta = Float64Array.from({ length: system.n }, (_, k) => .15 * Math.sin(k + 1));
  const derivative = sparseProduct(actual.matrix, delta), errors = [];
  const moved = step => input.nodes.map((row, i) => row.map((p, j) => {
    if (!i || i === input.nx || !j || j === input.nt) return { ...p };
    const d = normalize(input.directions[i][j]), value = step * delta[(i - 1) * (input.nt - 1) + j - 1];
    return { x: p.x + value * d.x, y: p.y + value * d.y };
  }));
  for (const h of [.02, .01, .005]) {
    const plus = system.evaluate(moved(h), { linearize: false }).residual;
    const minus = system.evaluate(moved(-h), { linearize: false }).residual;
    errors.push(maximum(derivative.map((v, k) => v - (plus[k] - minus[k]) / (2 * h))));
  }
  assert.ok(errors[1] < .26 * errors[0], JSON.stringify(errors));
  assert.ok(errors[2] < .26 * errors[1], JSON.stringify(errors));
  assert.ok(errors[2] < 2e-5, JSON.stringify(errors));
});

test('inverse harmonic residual normalizes guides, preserves mass scaling and rejects invalid geometry', () => {
  const input = fixture(), a = createTransverseHarmonicResidual(input).evaluate(input.nodes);
  const scaled = { ...input, massFlows: input.massFlows.map(m => 7 * m),
    directions: input.directions.map(row => row.map(d => ({ x: 3 * d.x, y: 3 * d.y }))) };
  const b = createTransverseHarmonicResidual(scaled).evaluate(input.nodes);
  a.residual.forEach((v, k) => close(v, b.residual[k], 1e-14));
  a.matrix.values.forEach((v, k) => close(v, b.matrix.values[k], 2e-12));
  assert.throws(() => createTransverseHarmonicResidual({ ...input, massFlows: [1, 0, 1] }), /Invalid/);
  const zero = structuredClone(input.directions); zero[1][1] = { x: 0, y: 0 };
  assert.throws(() => createTransverseHarmonicResidual({ ...input, directions: zero }), /nonzero/);
  const system = createTransverseHarmonicResidual(input), invalid = structuredClone(input.nodes); invalid[1][1].x = NaN;
  assert.throws(() => system.evaluate(invalid), /state/);
  const folded = input.nodes.map(row => row.map(p => ({ x: -p.x, y: p.y })));
  assert.throws(() => system.evaluate(folded), /quadrilateral/);
});
