import test from 'node:test';
import assert from 'node:assert/strict';
import { createTransverseHarmonicResidual } from '../src/geometry/transverse-harmonic-residual.js';
import { createCoonsGridGeometry } from '../src/geometry/tests/coons-grid-geometry.js';
import { solveCurvedHarmonicReference } from '../src/geometry/tests/curved-harmonic-reference.js';
import { conformalPolynomialFixture } from './fixtures/curved-harmonic-fields.js';
import { sparseProduct } from '../src/numerics/sparse.js';

const maximum = a => Math.max(0, ...a.map(Math.abs));
const normalized = p => { const d = Math.hypot(p.x, p.y); return { x: p.x / d, y: p.y / d }; };
const combine = (a, b, sign = 1) => Object.fromEntries(['point', 'ds', 'dt'].map(key => [key,
  { x: a[key].x + sign * b[key].x, y: a[key].y + sign * b[key].y }]));
const bilinear = (nodes, i, j, s, t) => {
  const p = [nodes[i][j], nodes[i + 1][j], nodes[i + 1][j + 1], nodes[i][j + 1]];
  const shape = [(1 - s) * (1 - t), s * (1 - t), s * t, (1 - s) * t];
  const ds = [t - 1, 1 - t, t, -t], dt = [s - 1, -s, s, 1 - s];
  const out = { point: { ...p[0] }, ds: { x: 0, y: 0 }, dt: { x: 0, y: 0 } };
  for (let k = 1; k < 4; k++) for (const key of ['x', 'y']) {
    const d = p[k][key] - p[0][key];
    out.point[key] += shape[k] * d; out.ds[key] += ds[k] * d; out.dt[key] += dt[k] * d;
  }
  return out;
};
const moved = (input, delta, step) => input.nodes.map((row, i) => row.map((p, j) => {
  if (!i || !j || i === input.nx || j === input.nt) return { ...p };
  const d = normalized(input.directions[i][j]), a = step * delta[(i - 1) * (input.nt - 1) + j - 1];
  return { x: p.x + a * d.x, y: p.y + a * d.y };
}));
function fixture({ coons = false, perturb = true } = {}) {
  const data = conformalPolynomialFixture({ nx: 4, nt: 3 });
  const original = structuredClone(data.nodes), analytic = data.geometry;
  const geometryCorrection = coons ? createCoonsGridGeometry(data).correction
    : (i, j, s, t) => combine(analytic.at(i, j, s, t), bilinear(original, i, j, s, t), -1);
  const directions = original.map((row, i) => row.map((_, j) => {
    const angle = .25 * Math.sin(i + .7 * j); return { x: Math.sin(angle), y: Math.cos(angle) };
  }));
  const input = { ...data, directions, geometryCorrection };
  if (perturb) input.nodes = moved(input, Float64Array.from({ length: 6 }, (_, k) => .008 * Math.sin(k + 1)), 1);
  return input;
}
const geometryFor = input => ({ at: (i, j, s, t) => combine(bilinear(input.nodes, i, j, s, t), input.geometryCorrection(i, j, s, t)) });

// Independently integrate scalar physical Laplace stiffness on the fixed
// curved geometry. Same specified 3x3 quadrature, separate dense assembly;
// this checks the discrete equations, not quadrature convergence.
function fieldStiffness(input, geometry) {
  const { nx, nt } = input, n = (nx - 1) * (nt - 1), matrix = new Float64Array(n * n);
  const id = (i, j) => !i || !j || i === nx || j === nt ? -1 : (i - 1) * (nt - 1) + j - 1;
  const q = [.5 - Math.sqrt(3 / 5) / 2, .5, .5 + Math.sqrt(3 / 5) / 2], w = [5 / 18, 4 / 9, 5 / 18];
  for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) {
    const ids = [id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)];
    for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) {
      const s = q[a], t = q[b], map = geometry.at(i, j, s, t);
      const ds = [t - 1, 1 - t, t, -t], dt = [s - 1, -s, s, 1 - s];
      const J = map.ds.x * map.dt.y - map.ds.y * map.dt.x;
      assert.ok(J > 0);
      const grad = ds.map((d, k) => ({ x: (map.dt.y * d - map.ds.y * dt[k]) / J,
        y: (map.ds.x * dt[k] - map.dt.x * d) / J }));
      for (let r = 0; r < 4; r++) for (let c = 0; c < 4; c++) if (ids[r] >= 0 && ids[c] >= 0)
        matrix[ids[r] * n + ids[c]] += w[a] * w[b] * J * (grad[r].x * grad[c].x + grad[r].y * grad[c].y);
    }
  }
  return matrix;
}

test('fixed all-cell conformal corrections give roundoff harmonic residual at the exact map', () => {
  const input = fixture({ perturb: false }), system = createTransverseHarmonicResidual(input);
  const actual = system.evaluate(input.nodes);
  assert.ok(maximum(actual.residual) < 1e-13, `Residual ${maximum(actual.residual)}`);
  assert.ok(actual.matrix.values.every(Number.isFinite));
  assert.deepEqual(system.evaluate(input.nodes, { linearize: false }).residual, actual.residual);
});

for (const coons of [false, true]) test(`${coons ? 'outer-edge Coons' : 'all-cell conformal'} curved Jacobian matches material finite differences`, t => {
  const input = fixture({ coons }), system = createTransverseHarmonicResidual(input), actual = system.evaluate(input.nodes);
  assert.ok(maximum(actual.residual) > 1e-4, 'Differentiate at a nonstationary perturbed state.');
  const direction = Float64Array.from({ length: system.n }, (_, k) => .15 * Math.cos(k + .4));
  const derivative = sparseProduct(actual.matrix, direction), errors = [];
  for (const h of [.02, .01, .005]) {
    const plus = system.evaluate(moved(input, direction, h), { linearize: false }).residual;
    const minus = system.evaluate(moved(input, direction, -h), { linearize: false }).residual;
    errors.push(maximum(derivative.map((v, k) => v - (plus[k] - minus[k]) / (2 * h))));
  }
  assert.ok(errors[1] < .27 * errors[0], JSON.stringify(errors));
  assert.ok(errors[2] < .27 * errors[1], JSON.stringify(errors));
  assert.ok(errors[2] < 2e-5, JSON.stringify(errors)); t.diagnostic(JSON.stringify({ errors }));
});

test('separate curved field correction solves K deltaEta = -rawResidual on the perturbed map', () => {
  const input = fixture(), geometry = geometryFor(input), system = createTransverseHarmonicResidual(input);
  const residual = system.evaluate(input.nodes).residual, matrix = fieldStiffness(input, geometry);
  const reference = solveCurvedHarmonicReference({ ...input, geometry }, { refinement: 1 });
  const delta = [];
  for (let i = 1; i < input.nx; i++) for (let j = 1; j < input.nt; j++) delta.push(reference.values[i][j] - j / input.nt);
  assert.ok(maximum(delta) > 1e-4);
  const defect = residual.map((v, row) => v + delta.reduce((sum, value, col) => sum + matrix[row * system.n + col] * value, 0));
  assert.ok(maximum(defect) < 2e-13, `Independent field defect ${maximum(defect)}`);
});

test('curved residual preserves rigid motion, length scaling and total mass scaling', () => {
  const input = fixture(), a = createTransverseHarmonicResidual(input).evaluate(input.nodes);
  const angle = .71, c = Math.cos(angle), s = Math.sin(angle), scale = 2.3;
  const rotate = p => ({ x: c * p.x - s * p.y, y: s * p.x + c * p.y });
  const vector = p => { const q = rotate(p); return { x: scale * q.x, y: scale * q.y }; };
  const mapped = { ...input, nodes: input.nodes.map(row => row.map(p => { const q = vector(p); return { x: q.x + 3, y: q.y - 4 }; })),
    directions: input.directions.map(row => row.map(rotate)), massFlows: input.massFlows.map(m => 7 * m),
    geometryCorrection: (...args) => Object.fromEntries(Object.entries(input.geometryCorrection(...args)).map(([key, value]) => [key, vector(value)])),
  };
  const b = createTransverseHarmonicResidual(mapped).evaluate(mapped.nodes);
  a.residual.forEach((v, k) => assert.ok(Math.abs(v - b.residual[k]) < 2e-13));
  a.matrix.values.forEach((v, k) => assert.ok(Math.abs(v - scale * b.matrix.values[k]) < 2e-12));
});

test('curved samples are frozen at construction and invalid corrections cannot shift observation nodes', () => {
  const input = fixture(), original = input.geometryCorrection; let multiplier = 1;
  const mutable = (...args) => Object.fromEntries(Object.entries(original(...args)).map(([key, value]) =>
    [key, { x: multiplier * value.x, y: multiplier * value.y }]));
  const system = createTransverseHarmonicResidual({ ...input, geometryCorrection: mutable }), before = system.evaluate(input.nodes);
  multiplier = 17;
  const after = system.evaluate(input.nodes);
  assert.deepEqual(after.residual, before.residual); assert.deepEqual(after.matrix.values, before.matrix.values);
  assert.throws(() => createTransverseHarmonicResidual({ ...input, geometryCorrection: {} }), /Invalid/);
  assert.throws(() => createTransverseHarmonicResidual({ ...input, geometryCorrection: (...args) => {
    const value = original(...args); value.ds.x = NaN; return value;
  } }), /Nonfinite/);
  assert.throws(() => createTransverseHarmonicResidual({ ...input, geometryCorrection: (...args) => {
    const value = original(...args); value.point.x += 1e-3; return value;
  } }), /observation node/);
});
