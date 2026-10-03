import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { initializeStreamtubeBodyFromGrid } from '../src/euler/tests/streamtube-body-restart.js';
import { initializePanelStreamtubeBody } from '../src/euler/streamtube-body-initializer.js';
import { initializeStreamtubeDensities } from '../src/euler/streamtube-initial-state.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { sparseDense, sparseProduct } from '../src/numerics/sparse.js';
import { directBodyConservation } from './oracles/streamtube-body.js';
import { stagnationPressureError } from './oracles/streamtube-entropy.js';

const fixture = (elements = 2, rest = {}) => ({ ...intrinsicBodyFixture({ elements, ...rest }), flowModel: 'incompressible', mach: 0 });
const max = a => Math.max(...Array.from(a, Math.abs));
const solved = new Map();
const bodyCase = (elements, controls = {}) => {
  const key = JSON.stringify({ elements, controls });
  if (!solved.has(key)) {
    const input = fixture(elements, controls), system = createStreamtubeBodySystem(input);
    const result = solveStreamtubeBody(system, { tolerance: 1e-11, maxIterations: 15 });
    assert.equal(result.converged, true, result.reason); solved.set(key, { input, system, result });
  }
  return solved.get(key);
};

test('incompressible body layout and every analytic Jacobian column retain all cut, stagnation, capture and farfield constraints', () => {
  const input = fixture(2, { bodySegments: 4, tubes: 2, alpha: .25 });
  input.bodies.forEach(body => { const row = body.surfaceFractions;
    body.surfaceFractions = Object.fromEntries(['upper', 'lower'].map((side, k) => [side,
      row.map((f, i) => i === 0 || i === row.length - 1 ? f : f + (k ? -.01 : .01) * Math.sin(Math.PI * f))])); });
  const system = createStreamtubeBodySystem(input), { layout } = system;
  assert.equal(layout.densityCount, 0); assert.equal(layout.rows.length, layout.n);
  assert.equal(layout.rowCounts.streamwise, undefined); assert.equal(layout.rowCounts.inletDensity, undefined);
  assert.equal(layout.rowCounts.leadingKutta, 2); assert.equal(layout.rowCounts.trailingKutta, 2);
  assert.equal(layout.globals.capture.filter(c => c !== null).length, 1);
  let state = system.initial.map((_, i) => 1e-5 * Math.sin(i + .3));
  for (const [key, value] of [['circulation', .07], ['source', .015], ['doubletX', .01], ['doubletY', -.02]]) state[layout.globals[key]] = value;
  for (let chart = 0; chart < 2; chart++) {
    const dense = system.jacobian(state), sparse = system.jacobian(state, { sparse: true });
    assert.deepEqual(sparseDense(sparse), dense);
    let worst = 0; const families = new Set();
    for (let col = 0; col < layout.n; col++) {
      const h = 1e-7, plus = state.slice(), minus = state.slice(); plus[col] += h; minus[col] -= h;
      const a = system.residual(plus), b = system.residual(minus);
      for (let row = 0; row < layout.n; row++) {
        const v = dense[row * layout.n + col];
        worst = Math.max(worst, Math.abs(v - (a[row] - b[row]) / (2 * h)) / Math.max(1, Math.abs(v)));
        if (v !== 0) families.add(layout.rows[row].kind);
      }
    }
    assert.ok(worst < 4e-7, `Incompressible all-column error ${worst}`);
    assert.equal(families.size, Object.keys(layout.rowCounts).length);
    const direction = state.map((_, i) => .1 * Math.sin(i + .4)), product = sparseProduct(sparse, direction), h = 1e-7;
    const a = system.residual(state.map((v, i) => v + h * direction[i])), b = system.residual(state.map((v, i) => v - h * direction[i]));
    product.forEach((v, i) => assert.ok(Math.abs(v - (a[i] - b[i]) / (2 * h)) < 1e-6 * Math.max(1, Math.abs(v))));
    state = system.rebase(state);
  }
});

test('small one/two-body incompressible solves preserve mass and Bernoulli and solve free captured mass, stagnation and Kutta', () => {
  for (const elements of [1, 2]) {
    const { input, system, result: r } = bodyCase(elements), initial = system.decode(system.initial);
    assert.equal(r.flowModel, 'incompressible'); assert.equal(system.conditions.h0, null);
    assert.ok(r.linearDiagnostics.maxRelativeResidual <= 1e-10);
    const c = directBodyConservation(r, input.bodies, system.conditions);
    for (const k of [0, 3]) assert.ok(Math.abs(c.balance[k]) < 2e-9);
    assert.ok(max(c.cutTraction) < 2e-9);
    for (const block of c.blocks) {
      assert.ok(max(block.internalCancellation) < 2e-9);
      for (const k of [0, 3]) assert.ok(Math.abs(block.maxLocal[k]) < 2e-9);
    }
    for (const row of r.sections) for (const group of row) for (const s of group) {
      assert.equal(s.rho, 1); assert.ok(Math.abs(s.p + .5 * s.q ** 2) < 1e-13);
    }
    for (let b = 0; b < elements; b++) for (let i = 0; i <= system.layout.nx; i++) if (!system.layout.active(b, i))
      assert.deepEqual(r.nodes[b][i].at(-1), r.nodes[b + 1][i][0]);
    if (elements === 1) assert.ok(Math.abs(r.diagnosticForces[0].cl) < 1e-10);
    else {
      assert.ok(Math.abs(r.captured[2] - initial.captured[2]) > .005);
      assert.equal(r.captured[1], initial.captured[1]);
      assert.ok(r.stagnation.every((s, b) => Math.abs(s - system.initialStagnation[b]) > 1e-4));
    }
    assert.match(r.forceStatus, /Unvalidated/);
  }
});

test('incompressible analytic KLU and independently differenced dense Newton solve the same small two-body system', () => {
  const input = fixture(2, { bodySegments: 4, tubes: 2 });
  const a = solveStreamtubeBody(createStreamtubeBodySystem(input), { tolerance: 1e-11 });
  const b = solveStreamtubeBody(createStreamtubeBodySystem(input), { tolerance: 1e-11, jacobianBackend: 'finite-difference' });
  assert.equal(a.converged, true, a.reason); assert.equal(b.converged, true, b.reason);
  assert.ok(max(a.x.map((v, i) => v - b.x[i])) < 2e-9);
  a.nodes.forEach((g, k) => g.forEach((row, i) => row.forEach((p, j) => assert.ok(Math.hypot(p.x - b.nodes[k][i][j].x, p.y - b.nodes[k][i][j].y) < 2e-9))));
});

test('relaxed incompressible geometry transfers to and converges at the actual target Mach without changing prescribed mass', () => {
  // This test exercises state transfer on a panel-seeded grid that converges
  // reliably; the 16-segment lifting grid is a separate solver stress case.
  const prepared = initializePanelStreamtubeBody(fixture(2, { alpha: 2, bodySegments: 8, surfaceSpacing: 'cosine', tubes: 5, tubeGrowth: 3 }),
    { crosslinePlacement: 'potential', outerCrosslineSpread: .15, normalSpacing: 'stagnation' });
  const { input } = prepared, result = solveStreamtubeBody(prepared.system, { initial: prepared.initial, tolerance: 3e-11 });
  assert.equal(result.converged, true, result.reason);
  const before = structuredClone(result);
  const target = { ...input, flowModel: 'compressible', mach: .2, streamwiseMode: 'isentropic' };
  const grid = initializeStreamtubeBodyFromGrid(target, result), initial = grid.system.evaluate(grid.initial);
  assert.deepEqual(result, before); assert.deepEqual(initial.nodes, result.nodes);
  assert.deepEqual(initial.captured, result.captured);
  assert.equal(grid.system.conditions.mach, .2); assert.ok(initial.sections.flat(2).some(s => Math.abs(s.rho - 1) > .001));
  const r = solveStreamtubeBody(grid.system, { initial: grid.initial, tolerance: 1e-11 });
  assert.equal(r.converged, true, r.reason); assert.equal(r.flowModel, 'compressible');
  assert.ok(r.linearDiagnostics.maxRelativeResidual <= 1e-10);
  const pressure = stagnationPressureError(r.sections.map(row => row.flat()), { gamma: grid.system.conditions.gamma,
    referencePressure: grid.system.conditions.pInf, freestreamMach: .2 });
  assert.ok(pressure.maxRelativeError < 2e-10);
  const c = directBodyConservation(r, target.bodies, grid.system.conditions);
  for (const k of [0, 3]) assert.ok(Math.abs(c.balance[k]) < 2e-9);
  assert.ok(max(c.cutTraction) < 2e-9);
  const bad = structuredClone(result); bad.captured[1] += .01;
  assert.throws(() => initializeStreamtubeBodyFromGrid(target, bad), /prescribed/);
  assert.deepEqual(result, before);
  const inc = createStreamtubeBodySystem(input);
  assert.throws(() => initializeStreamtubeDensities(inc, inc.initial), /not applicable/);
  assert.throws(() => createStreamtubeBodySystem({ ...input, mach: .1 }), /subcritical/);
  assert.throws(() => createStreamtubeBodySystem({ ...input, streamwiseMode: 'momentum' }), /streamwise/);
});
