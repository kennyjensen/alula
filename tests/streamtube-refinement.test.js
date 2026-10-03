import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { refineStreamtubeBody } from '../src/euler/streamtube-refinement.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directBodyConservation } from './oracles/streamtube-body.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';
import { sparseProduct } from '../src/numerics/sparse.js';

test('nested multielement Euler refinement conserves every parent mass and retains exact material boundaries', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), stagnationMotion: 'walls-only' };
  const source = createStreamtubeBodySystem(input);
  const parent = solveStreamtubeBody(source, { maxIterations: 8 }), before = structuredClone(input);
  assert.equal(parent.converged, true);
  const counts = [[1, 2], [2, 3], [2, 1]];
  for (const normalInterpolation of ['linear', 'streamfunction-quadratic']) for (const streamwiseFactor of [1, 2]) {
    const r = refineStreamtubeBody(input, source, { initial: parent.x, streamwiseFactor, normalSubdivisions: counts, normalInterpolation });
    const value = r.system.evaluate(r.initial), boundary = r.system.decode(r.initial).nodes;
    assert.deepEqual(value.captured, parent.captured);
    assert.equal(r.system.conditions.lengthScale, source.conditions.lengthScale);
    for (let g = 0; g < counts.length; g++) {
      let child = 0;
      for (let j = 0; j < counts[g].length; j++) {
        const mass = value.allocation.groups[g].slice(child, child + counts[g][j]).reduce((s, v) => s + v.massFlow, 0);
        assert.ok(Math.abs(mass - parent.allocation.groups[g][j].massFlow) < 1e-14);
        for (let i = 0; i <= source.layout.nx; i++) {
          const a = value.nodes[g][i * streamwiseFactor][child], b = parent.nodes[g][i][j];
          assert.ok(Math.hypot(a.x - b.x, a.y - b.y) < 2e-14);
        }
        child += counts[g][j];
      }
    }
    for (let b = 0; b < r.system.layout.elements; b++) {
      const body = r.system.layout.bodies[b];
      for (let i = 0; i <= r.system.layout.nx; i++) {
        if (i >= body.leadingIndex && i <= body.trailingIndex) {
          for (const side of ['upper', 'lower']) {
            const p = r.system.curves[b].branch(side, r.system.fractions[b][side][i - body.leadingIndex], value.stagnation[b]).point;
            const actual = side === 'upper' ? value.nodes[b + 1][i][0] : value.nodes[b][i].at(-1);
            assert.ok(Math.hypot(actual.x - p.x, actual.y - p.y) < 2e-14);
          }
        } else assert.deepEqual(boundary[b][i].at(-1), boundary[b + 1][i][0]);
      }
    }
    assert.equal(directStreamtubeVolumeGeometry(value.nodes).valid, true);
    assert.ok(value.diagnostics.maxMach < 1);
  }
  assert.deepEqual(input, before);
  assert.throws(() => refineStreamtubeBody(input, source, { maxNodes: 1 }), /budget/);
  assert.throws(() => refineStreamtubeBody(input, source, { normalSubdivisions: [[0]] }), /subdivisions/);
  assert.throws(() => refineStreamtubeBody(input, source, { normalInterpolation: 'unknown' }), /interpolation/);
  for (const normalInterpolation of ['linear', 'streamfunction-quadratic']) {
    const preview = refineStreamtubeBody(input, { ...source, evaluate: () => assert.fail('Geometry preview must not evaluate the gas state') },
      { initial: parent.x, streamwiseFactor: 1, normalSubdivisions: counts, initializeFlow: false, normalInterpolation });
    const initialized = refineStreamtubeBody(input, source, { initial: parent.x, streamwiseFactor: 1, normalSubdivisions: counts, normalInterpolation });
    assert.deepEqual(preview.initialEuler.nodes, initialized.initialEuler.nodes);
    assert.equal(preview.diagnostics.flowInitialized, false);
    assert.equal(preview.diagnostics.residual, undefined);
    assert.equal(preview.diagnostics.maxMach, undefined);
  }
});

for (const controls of [
  { normalInterpolation: 'linear' },
  ...[2, 4].map(n => ({ normalInterpolation: 'streamfunction-quadratic', streamwiseFactor: 1,
    normalSubdivisions: [[1, n], [n, n], [n, 1]] })),
]) test(`${JSON.stringify(controls)}: controlled two-element Euler refinement converges with independent conservation checks`, () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), stagnationMotion: 'walls-only', normalStencil: 'body-stations' };
  const source = createStreamtubeBodySystem(input), parent = solveStreamtubeBody(source, { maxIterations: 8 });
  assert.equal(parent.converged, true);
  const refined = refineStreamtubeBody(input, source, { initial: parent.x, ...controls });
  const h = 2e-6, direction = refined.initial.map((v, i) => .001 * Math.sin(.43 * i + .2));
  const exact = sparseProduct(refined.system.jacobian(refined.initial, { sparse: true }), direction);
  const plus = refined.system.evaluate(refined.initial.map((v, i) => v + h * direction[i])).residual;
  const minus = refined.system.evaluate(refined.initial.map((v, i) => v - h * direction[i])).residual;
  assert.ok(Math.max(...exact.map((v, i) => {
    const fd = (plus[i] - minus[i]) / (2 * h); return Math.abs(fd - v) / Math.max(1, Math.abs(fd), Math.abs(v));
  })) < 2e-6);
  const flow = solveStreamtubeBody(refined.system, { initial: refined.initial, maxIterations: 12 });
  assert.equal(flow.converged, true, flow.reason);
  const conservation = directBodyConservation(flow, refined.input.bodies, refined.system.conditions);
  assert.ok(Math.max(...conservation.balance.map(Math.abs), ...conservation.cutTraction.map(Math.abs)) < 2e-9);
  assert.equal(directStreamtubeVolumeGeometry(flow.nodes).valid, true);
});
