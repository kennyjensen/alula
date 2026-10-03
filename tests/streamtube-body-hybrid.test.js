import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { sparseDense, sparseProduct } from '../src/numerics/sparse.js';

const parameters = elements => ({ ...intrinsicBodyFixture({ elements, bodySegments: 4, tubes: 2, mach: .5, alpha: .25 }),
  streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-3 }, upwind: { mucon: 1, mcrit: .5, boundary: { kind: 'unfiltered-first-two' } } });
const seeded = system => system.initial.map((_, k) => (k < system.layout.densityCount ? 1e-3 : 1e-5) * Math.sin(k + .3));
const fd4 = (system, x, d, h) => {
  const values = [-2, -1, 1, 2].map(m => system.residual(x.map((v, k) => v + m * h * d[k])));
  return values[0].map((v, k) => (v - 8 * values[1][k] + 8 * values[2][k] - values[3][k]) / (12 * h));
};
const error = (a, b) => Math.abs(a - b) / Math.max(1, Math.abs(a), Math.abs(b));

test('hybrid body differentiates the full residual in the transition region', t => {
  const system = createStreamtubeBodySystem(parameters(1)), x = seeded(system), n = system.layout.n;
  const flow = system.evaluate(x), matrix = system.jacobian(x);
  assert.ok(flow.diagnostics.hybrid.blendedCells > 0);
  assert.ok(flow.diagnostics.hybrid.entropyCells > 0);
  assert.deepEqual(sparseDense(system.jacobian(x, { sparse: true })), matrix);
  let maximumError = 0;
  for (let k = 0; k < n; k++) {
    const d = new Float64Array(n); d[k] = 1;
    const numerical = fd4(system, x, d, 1e-6);
    for (let row = 0; row < n; row++) maximumError = Math.max(maximumError, error(matrix[row * n + k], numerical[row]));
  }
  assert.ok(maximumError < 2e-7, String(maximumError));
  t.diagnostic(JSON.stringify({ unknowns: n, comparisons: n * n, maximumError, hybrid: flow.diagnostics.hybrid }));
});

test('two-element hybrid capture, stagnation and moving-chart derivatives remain complete', t => {
  const system = createStreamtubeBodySystem(parameters(2)), n = system.layout.n;
  let x = seeded(system), maximumError = 0, comparisons = 0;
  for (let chart = 0; chart < 2; chart++) {
    const flow = system.evaluate(x), matrix = system.jacobian(x, { sparse: true });
    assert.ok(flow.diagnostics.hybrid.blendedCells > 0);
    const globals = Object.values(system.layout.globals).flat().filter(c => c !== null);
    const directions = globals.map(c => Float64Array.from({ length: n }, (_, k) => k === c ? 1 : 0));
    directions.push(Float64Array.from({ length: n }, (_, k) => Math.cos(k + .7)));
    for (const d of directions) {
      const exact = sparseProduct(matrix, d), numerical = fd4(system, x, d, 5e-7);
      for (let k = 0; k < n; k++) { maximumError = Math.max(maximumError, error(exact[k], numerical[k])); comparisons++; }
    }
    x = system.rebase(x);
  }
  assert.ok(maximumError < 2e-7, String(maximumError));
  t.diagnostic(JSON.stringify({ unknowns: n, comparisons, maximumError }));
});

test('zero speed bias recovers every isentropic body row and exposes momentum departures', () => {
  const base = intrinsicBodyFixture({ bodySegments: 4, tubes: 2, mach: .2 });
  const a = createStreamtubeBodySystem({ ...base, streamwiseMode: 'isentropic' });
  const b = createStreamtubeBodySystem({ ...base, streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 },
    upwind: { mucon: 0, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } });
  const x = seeded(a), old = a.evaluate(x), current = b.evaluate(x);
  assert.deepEqual(current.residual, old.residual);
  assert.deepEqual(current.nodes, old.nodes); assert.deepEqual(current.sections, old.sections);
  assert.equal(current.diagnostics.hybrid.momentumCells, 0);
  assert.equal(current.diagnostics.hybrid.blendedCells, 0);
  assert.ok(current.diagnostics.hybrid.maxMomentumDeparture > 0);
  assert.equal(current.diagnostics.hybrid.entropyCells, (b.layout.nx - 1) * b.layout.nt);
  const ja = a.jacobian(x), jb = b.jacobian(x);
  assert.ok(ja.every((v, k) => error(v, jb[k]) < 3e-12));
});

test('hybrid body controls are explicit, copied and restricted to the selected formulation', () => {
  const input = parameters(1), system = createStreamtubeBodySystem(input), before = system.residual(system.initial);
  input.hybrid.epsilonP = 1; assert.deepEqual(system.residual(system.initial), before);
  for (const epsilonP of [0, -1, NaN, Infinity]) assert.throws(() => createStreamtubeBodySystem({ ...parameters(1), hybrid: { epsilonP } }), /epsilonP/);
  assert.throws(() => createStreamtubeBodySystem({ ...parameters(1), upwind: undefined }), /explicit upwinding/);
  assert.throws(() => createStreamtubeBodySystem({ ...parameters(1), streamwiseMode: 'momentum' }), /Hybrid controls/);
});

test('hybrid transition derivatives include wall and wake displacement', t => {
  const input = parameters(1), nx = input.outerLower.length - 1;
  input.displacement = { surfaces: input.bodies.map(b => Object.fromEntries(['upper', 'lower'].map(side => [side,
    Array.from({ length: b.trailingIndex - b.leadingIndex + 1 }, (_, i) => .0001 * (1 + .02 * i))]))),
  wakes: input.bodies.map(b => Array.from({ length: nx - b.trailingIndex }, (_, i) => .00024 * (1 + .01 * i))) };
  const system = createStreamtubeBodySystem(input), x = seeded(system), original = system.evaluate(x);
  assert.ok(original.diagnostics.hybrid.blendedCells > 0);
  const block = system.jacobian(x, { includeDisplacement: true });
  let maximumError = 0, comparisons = 0;
  for (const [col, p] of block.parameters.entries()) {
    const samples = [], h = 1e-7;
    try {
      for (const multiplier of [-2, -1, 1, 2]) {
        const d = structuredClone(input.displacement);
        if (p.kind === 'wake') d.wakes[p.body][p.index] += multiplier * h;
        else for (const side of p.side === 'both' ? ['upper', 'lower'] : [p.side]) d.surfaces[p.body][side][p.index] += multiplier * h;
        system.setDisplacement(d); samples.push(system.residual(x));
      }
    } finally { system.setDisplacement(input.displacement); }
    for (let row = 0; row < system.layout.n; row++) {
      const numerical = (samples[0][row] - 8 * samples[1][row] + 8 * samples[2][row] - samples[3][row]) / (12 * h);
      maximumError = Math.max(maximumError, error(block.displacement[row].get(col) ?? 0, numerical)); comparisons++;
    }
  }
  assert.ok(maximumError < 3e-7, String(maximumError));
  assert.deepEqual(system.residual(x), original.residual);
  t.diagnostic(JSON.stringify({ columns: block.parameters.length, comparisons, maximumError }));
});
