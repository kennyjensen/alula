import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeChannel, solveStreamtubeChannel } from '../src/euler/tests/streamtube-channel.js';
import { sparseDense, sparseProduct } from '../src/numerics/sparse.js';
import { nozzleChannel, vortexChannel, vortexErrors, directChannelConservation } from './oracles/streamtube.js';

const max = a => Math.max(...Array.from(a, Math.abs));
const close = (a, b, tolerance = 5e-8) => assert.ok(Math.abs(a - b) < tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);

test('complete analytic channel Jacobian matches independent whole-residual columns and sparse products', () => {
  for (const streamwiseMode of ['momentum', 'isentropic']) for (const input of [nozzleChannel(7, 1), nozzleChannel(8, 3), vortexChannel(8, 3)]) {
    const system = createStreamtubeChannel({ ...input, streamwiseMode }), { n, densityCount: nd } = system;
    const state = system.initial.map((_, i) => .001 * Math.sin(.7 * (i + 1)));
    const dense = system.jacobian(state), sparse = system.jacobian(state, { sparse: true });
    assert.deepEqual(sparseDense(sparse), dense);
    let flowFromGrid = 0, gridFromFlow = 0;
    for (let col = 0; col < n; col++) {
      const h = 1e-6, plus = state.slice(), minus = state.slice(); plus[col] += h; minus[col] -= h;
      const rp = system.residual(plus), rm = system.residual(minus);
      for (let row = 0; row < n; row++) {
        const derivative = dense[row * n + col];
        close(derivative, (rp[row] - rm[row]) / (2 * h));
        if (row < nd && col >= nd) flowFromGrid = Math.max(flowFromGrid, Math.abs(derivative));
        if (row >= nd && col < nd) gridFromFlow = Math.max(gridFromFlow, Math.abs(derivative));
      }
    }
    if (system.nt > 1) { assert.ok(flowFromGrid > .01); assert.ok(gridFromFlow > .01); }
    const direction = state.map((_, i) => Math.cos(i + .2)), product = sparseProduct(sparse, direction), h = 2e-6;
    const rp = system.residual(state.map((v, i) => v + h * direction[i]));
    const rm = system.residual(state.map((v, i) => v - h * direction[i]));
    product.forEach((v, i) => close(v, (rp[i] - rm[i]) / (2 * h)));
    // Reusing the compiled stencil must not overwrite an earlier matrix.
    system.jacobian(state.map(v => v * .9), { sparse: true });
    assert.deepEqual(sparseDense(sparse), dense);
  }
});

test('analytic dense and KLU channel solves recover the independent finite-difference solution', () => {
  for (const input of [nozzleChannel(12), vortexChannel(12, 4)]) {
    const system = createStreamtubeChannel(input), initial = system.initial.map((_, i) => .001 * Math.sin(i + 1));
    const reference = solveStreamtubeChannel(system, { initial, tolerance: 1e-12, jacobianBackend: 'finite-difference' });
    assert.equal(reference.converged, true, reference.reason);
    for (const linearBackend of ['dense', 'klu']) {
      const r = solveStreamtubeChannel(system, { initial, tolerance: 1e-12, jacobianBackend: 'analytic', linearBackend });
      assert.equal(r.converged, true, r.reason);
      assert.ok(max(r.residual) < 1e-12);
      r.x.forEach((v, i) => close(v, reference.x[i], 1e-9));
      for (const key of ['total', 'external', 'internalCancellation'])
        assert.ok(max(directChannelConservation(r)[key]) < 2e-10);
    }
  }
  const system = createStreamtubeChannel(nozzleChannel(8));
  assert.throws(() => solveStreamtubeChannel(system, { jacobianBackend: 'invalid' }), /backends/);
  assert.throws(() => solveStreamtubeChannel(system, { jacobianBackend: 'finite-difference', linearBackend: 'klu' }), /backends/);
});

test('sparse intrinsic Euler refinement preserves exact-vortex convergence beyond the dense reference grids', () => {
  const errors = [];
  for (const [nx, nt] of [[16, 4], [32, 8], [64, 16]]) {
    const input = vortexChannel(nx, nt), system = createStreamtubeChannel(input);
    const r = solveStreamtubeChannel(system, { tolerance: 1e-11, jacobianBackend: 'analytic', linearBackend: 'klu' });
    assert.equal(r.converged, true, r.reason);
    errors.push(vortexErrors(r, input));
    const matrix = system.jacobian(r.x, { sparse: true });
    assert.ok(matrix.values.length < 15 * system.n, `${matrix.values.length} entries for ${system.n} unknowns`);
    assert.ok(max(directChannelConservation(r).external) < 2e-9);
  }
  for (let i = 1; i < errors.length; i++) for (const key of ['pressureMax', 'positionMax'])
    assert.ok(errors[i][key] < .3 * errors[i - 1][key], JSON.stringify(errors));
  assert.ok(errors.at(-1).pressureMax < .00004, JSON.stringify(errors));
  assert.ok(errors.at(-1).positionMax < 2e-6, JSON.stringify(errors));
});
