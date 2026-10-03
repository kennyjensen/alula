import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { initializePanelStreamtubeBody } from '../src/euler/streamtube-body-initializer.js';
import { sparseDense, sparseProduct } from '../src/numerics/sparse.js';
import { directBodyConservation } from './oracles/streamtube-body.js';

test('complete body Jacobian matches every independent residual column, including all global constraints', () => {
  for (const streamwiseMode of ['momentum', 'isentropic']) for (const [elements, asymmetric] of [[1, false], [2, false], [2, true]]) {
    const input = { ...intrinsicBodyFixture({ elements, alpha: .25 }), streamwiseMode };
    if (asymmetric) input.bodies.forEach(body => {
      const row = body.surfaceFractions;
      body.surfaceFractions = Object.fromEntries(['upper', 'lower'].map((side, k) => [side,
        row.map((f, i) => i === 0 || i === row.length - 1 ? f : f + (k === 0 ? .01 : -.01) * Math.sin(Math.PI * f))]));
    });
    const system = createStreamtubeBodySystem(input), { layout } = system, n = layout.n;
    let state = system.initial.map((_, i) => (i < layout.densityCount ? 1e-3 : 1e-5) * Math.sin(i + .3));
    for (const [key, value] of [['circulation', .07], ['source', .015], ['doubletX', .01], ['doubletY', -.02]]) state[layout.globals[key]] = value;
    for (let chart = 0; chart < 2; chart++) {
      assert.equal(system.admissible(state), true, String([...system.rejectedTrials]));
      const dense = system.jacobian(state), sparse = system.jacobian(state, { sparse: true });
      assert.deepEqual(sparseDense(sparse), dense);
      let worst = { error: 0 }; const families = {};
      for (let col = 0; col < n; col++) {
        const h = 1e-7, plus = state.slice(), minus = state.slice(); plus[col] += h; minus[col] -= h;
        const a = system.residual(plus), b = system.residual(minus);
        for (let row = 0; row < n; row++) {
          const exact = dense[row * n + col], difference = (a[row] - b[row]) / (2 * h);
          const error = Math.abs(exact - difference) / Math.max(1, Math.abs(exact), Math.abs(difference));
          if (error > worst.error) worst = { error, row, col, family: layout.rows[row].kind, exact, difference };
          families[layout.rows[row].kind] = Math.max(families[layout.rows[row].kind] ?? 0, Math.abs(exact));
        }
      }
      assert.ok(worst.error < 4e-7, JSON.stringify({ elements, asymmetric, chart, ...worst }));
      for (const family of Object.keys(layout.rowCounts)) assert.ok(families[family] > 0, family);
      const direction = state.map((_, i) => Math.cos(i + .2)), product = sparseProduct(sparse, direction), h = 1e-7;
      const a = system.residual(state.map((v, i) => v + h * direction[i])), b = system.residual(state.map((v, i) => v - h * direction[i]));
      for (let i = 0; i < n; i++) assert.ok(Math.abs(product[i] - (a[i] - b[i]) / (2 * h)) < 1e-6 * Math.max(1, Math.abs(product[i])), `Jv row ${i}`);
      system.jacobian(state.map(v => v * .99), { sparse: true }); assert.deepEqual(sparseDense(sparse), dense);
      state = system.rebase(state);
    }
  }
});

test('analytic dense and KLU body solves agree with the independently differenced flow/grid solution', () => {
  for (const elements of [1, 2]) {
    const input = intrinsicBodyFixture({ elements, alpha: elements === 1 ? 2 : 0 });
    const reference = solveStreamtubeBody(createStreamtubeBodySystem(input), { jacobianBackend: 'finite-difference', tolerance: 1e-11 });
    assert.equal(reference.converged, true, reference.reason);
    for (const linearBackend of ['dense', 'klu']) {
      const system = createStreamtubeBodySystem(input);
      const r = solveStreamtubeBody(system, { jacobianBackend: 'analytic', linearBackend, tolerance: 1e-11 });
      assert.equal(r.converged, true, r.reason);
      assert.equal(r.linearBackend, linearBackend); assert.equal(r.jacobianBackend, 'analytic');
      for (let i = 0; i < r.x.length; i++) assert.ok(Math.abs(r.x[i] - reference.x[i]) < 2e-9, `State ${i}`);
      for (let g = 0; g < r.nodes.length; g++) for (let i = 0; i < r.nodes[g].length; i++) for (let j = 0; j < r.nodes[g][i].length; j++) {
        const p = r.nodes[g][i][j], q = reference.nodes[g][i][j]; assert.ok(Math.hypot(p.x - q.x, p.y - q.y) < 2e-9);
      }
      const conservation = directBodyConservation(r, input.bodies, system.conditions);
      assert.ok(Math.max(...conservation.balance.map(Math.abs)) < 2e-9);
      if (linearBackend === 'klu') { assert.ok(r.linearDiagnostics.solves > 0); assert.ok(r.linearDiagnostics.maxRelativeResidual <= 1e-10); }
    }
  }
  const system = createStreamtubeBodySystem(intrinsicBodyFixture());
  assert.throws(() => solveStreamtubeBody(system, { jacobianBackend: 'finite-difference', linearBackend: 'klu' }), /backends/);
  assert.throws(() => solveStreamtubeBody(system, { jacobianBackend: 'unknown' }), /backends/);
});

test('analytic body assembly works on the admissible refined grid whose default finite-difference step fails', () => {
  const r = initializePanelStreamtubeBody(intrinsicBodyFixture({ elements: 2, alpha: 2, bodySegments: 32, surfaceSpacing: 'cosine', tubes: 7, tubeGrowth: 3 }));
  const matrix = r.system.jacobian(r.initial, { sparse: true });
  assert.equal(matrix.n, r.system.layout.n); assert.ok(matrix.values.every(Number.isFinite));
  assert.ok(matrix.values.length < 35 * matrix.n, `${matrix.values.length} entries for ${matrix.n} unknowns`);
  const direction = r.initial.map((_, i) => Math.sin(i + .4)), h = 1e-8;
  const a = r.system.residual(r.initial.map((v, i) => v + h * direction[i])), b = r.system.residual(r.initial.map((v, i) => v - h * direction[i]));
  const product = sparseProduct(matrix, direction);
  let worst = 0;
  for (let i = 0; i < matrix.n; i++) worst = Math.max(worst, Math.abs(product[i] - (a[i] - b[i]) / (2 * h)) / Math.max(1, Math.abs(product[i])));
  assert.ok(worst < 2e-5, `Refined-grid Jv error ${worst}`);
});
