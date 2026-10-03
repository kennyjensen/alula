import test from 'node:test';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { initializeCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-initializer.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

function fixture(wakeGeometry = 'centerline') {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), wakeGeometry };
  const seeded = createCoupledStreamtubeBody(input), x = seeded.initial.slice(0, seeded.ne);
  return { input, options: { reynolds: 1e6, ncrit: 9,
    initialEuler: { x, ...seeded.euler.decode(x) },
    initialBL: seeded.initial.slice(seeded.ne).map((v, k) => k % 4 === 1 || k % 4 === 2 ? 128 * v : v) } };
}

for (const wakeGeometry of ['centerline', 'independent-banks']) test(`${wakeGeometry}: coupled initialization backtracks only the thickness guess and keeps the original four-surface/two-wake equations unsolved`, t => {
  const { input, options } = fixture(wakeGeometry), before = structuredClone({ input, options }), attempts = [];
  const prepared = initializeCoupledStreamtubeBody(input, options, { onAttempt: h => attempts.push(h) });
  const { system, initialization, mesh } = prepared, factor = initialization.thicknessFactor;
  assert.ok(factor < 1 && factor > 0); assert.equal(mesh.quality.valid, true);
  assert.equal(mesh.initialization.flowSolved, false); assert.equal(initialization.equationsChanged, false);
  assert.equal(attempts.length, initialization.history.length);
  assert.ok(attempts.slice(0, -1).every(h => !h.accepted)); assert.equal(attempts.at(-1).accepted, true);
  assert.deepEqual({ input, options }, before);
  assert.equal(system.conditions.reynolds, options.reynolds); assert.equal(system.bl.scale, .001);
  for (let k = 0; k < options.initialBL.length; k++) assert.equal(system.initial[system.ne + k],
    options.initialBL[k] * (k % 4 === 1 || k % 4 === 2 ? factor : 1));
  // Geometry changes do not replace density/global unknowns with an
  // unrelated gas inversion. This makes zero-displacement backtracking continuous.
  for (let k = 0; k < system.euler.layout.densityCount; k++) assert.equal(system.initial[k], options.initialEuler.x[k]);
  assert.equal(system.bl.surfaces.length, 4); assert.equal(system.bl.wakes.length, 2);
  if (wakeGeometry === 'independent-banks') assert.ok(system.evaluate(system.initial).outer.diagnostics.residualByFamily.wakeGap < 1e-12);
  // Admissibility is the initializer's contract, not convergence. The large
  // thickness seed is deliberately outside the proven Newton attraction
  // region. Do not certify it as a root just because every cell is positive.
  const r = solveCoupledStreamtubeBody(system, { maxIterations: 0, tolerance: 1e-10 });
  assert.equal(r.converged, false); assert.equal(r.mesh.quality.valid, true);
  assert.equal(r.mesh.initialization.flowSolved, false);
  assert.ok(r.residual.every(Number.isFinite));
  assert.deepEqual(r.residual, system.residual(system.initial));
  t.diagnostic(JSON.stringify({ factor, attempts: attempts.length, families: r.families }));
});

test('failed coupled initialization preserves inputs and propagates caller cancellation', () => {
  const { input, options } = fixture(), before = structuredClone(options);
  assert.throws(() => initializeCoupledStreamtubeBody(input, options, { maximumBacktracks: 0 }), error => {
    assert.equal(error.initialization.history.length, 1); assert.equal(error.initialization.history[0].accepted, false);
    return /Coupled initialization failed/.test(error.message);
  });
  assert.deepEqual(options, before);
  const stop = new Error('caller stopped');
  assert.throws(() => initializeCoupledStreamtubeBody(input, options, { onAttempt: () => { throw stop; } }), error => error === stop);
  assert.throws(() => initializeCoupledStreamtubeBody(input, {}), /Supply an Euler state/);
  assert.throws(() => initializeCoupledStreamtubeBody(input, options, { maximumBacktracks: -1 }), /backtracking limit/);
  for (const initialThicknessFactor of [0, -1, 2, NaN])
    assert.throws(() => initializeCoupledStreamtubeBody(input, options, { initialThicknessFactor }), /thickness factor/);
});

test('an explicit thinner guess scales thickness alone and keeps the same target equations', () => {
  const { input, options } = fixture('independent-banks');
  // Undo this fixture's intentionally thick seed, then prescribe a quarter
  // starting profile. This tests the retry primitive without a flow sweep.
  options.initialBL = options.initialBL.map((v, k) => k % 4 === 1 || k % 4 === 2 ? v / 128 : v);
  const before = structuredClone({ input, options });
  const a = initializeCoupledStreamtubeBody(input, options, { initialThicknessFactor: .25, maximumBacktracks: 0 });
  assert.equal(a.initialization.thicknessFactor, .25);
  assert.equal(a.initialization.history.length, 1);
  assert.equal(a.system.conditions.reynolds, options.reynolds);
  for (let k = 0; k < options.initialBL.length; k++) assert.equal(a.system.initial[a.system.ne + k],
    options.initialBL[k] * (k % 4 === 1 || k % 4 === 2 ? .25 : 1));
  const v = a.system.evaluate(a.system.initial);
  const restored = createCoupledStreamtubeBody(input, { ...options, initialEuler: {
    x: a.system.initial.slice(0, a.system.ne), nodes: v.outer.nodes, undisplacedNodes: v.outer.undisplacedNodes },
    initialBL: a.system.initial.slice(a.system.ne) });
  assert.deepEqual(restored.residual(restored.initial), v.residual);
  assert.deepEqual({ input, options }, before);
});
