import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directBodyConservation } from './oracles/streamtube-body.js';
import { initializeStreamtubeBodyFromGrid } from '../src/euler/tests/streamtube-body-restart.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { stagnationPressureError } from './oracles/streamtube-entropy.js';

const max = a => Math.max(...Array.from(a, Math.abs));

test('two-body trust-region solves retain all physical constraints in incompressible and conservative Euler modes', () => {
  for (const flowModel of ['incompressible', 'compressible']) {
    const input = { ...intrinsicBodyFixture({ elements: 2 }), flowModel, mach: flowModel === 'incompressible' ? 0 : .2 };
    const system = createStreamtubeBodySystem(input);
    const r = solveStreamtubeBody(system, { stepMethod: 'dogleg', initialTrustRadius: .01, maxIterations: 50, tolerance: 1e-11 });
    assert.equal(r.converged, true, `${flowModel}: ${r.reason}, R=${r.diagnostics.residual}`);
    assert.equal(r.stepMethod, 'dogleg'); assert.ok(r.diagnostics.residual < 1e-11);
    assert.ok(r.history.slice(1).some(h => h.stepKind !== 'newton'));
    assert.ok(r.history.slice(1).every(h => h.actualReduction > 0 && h.reductionRatio > 1e-4));
    assert.ok(r.linearDiagnostics.maxRelativeResidual < 1e-10);
    const c = directBodyConservation(r, input.bodies, system.conditions);
    for (const k of flowModel === 'incompressible' ? [0, 3] : [0, 1, 2, 3]) assert.ok(Math.abs(c.balance[k]) < 2e-9);
    assert.ok(max(c.cutTraction) < 2e-9);
    for (const block of c.blocks) {
      assert.ok(max(block.internalCancellation) < 2e-9);
      for (const k of [0, 3]) assert.ok(Math.abs(block.maxLocal[k]) < 2e-9);
    }
    assert.equal(r.surfaces.length, 4); assert.match(r.forceStatus, /Unvalidated/);
  }
});

test('analytic dense and KLU trust-region steps agree on the same small body chart', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), flowModel: 'incompressible', mach: 0 };
  const controls = { stepMethod: 'dogleg', initialTrustRadius: .01, maxIterations: 50, tolerance: 1e-11 };
  const a = solveStreamtubeBody(createStreamtubeBodySystem(input), controls);
  const b = solveStreamtubeBody(createStreamtubeBodySystem(input), { ...controls, linearBackend: 'dense' });
  assert.equal(a.converged, true, a.reason); assert.equal(b.converged, true, b.reason);
  a.nodes.forEach((g, k) => g.forEach((row, i) => row.forEach((p, j) => assert.ok(Math.hypot(p.x - b.nodes[k][i][j].x, p.y - b.nodes[k][i][j].y) < 2e-9))));
  assert.ok(max(a.x.map((v, i) => v - b.x[i])) < 2e-9);
});

test('trust-region controls and iteration limits cannot certify an unconverged body flow', () => {
  const input = intrinsicBodyFixture({ elements: 1 }), system = createStreamtubeBodySystem(input);
  assert.throws(() => solveStreamtubeBody(system, { stepMethod: 'unknown' }), /step method/);
  assert.throws(() => solveStreamtubeBody(system, { stepMethod: 'dogleg', jacobianBackend: 'finite-difference' }), /step method/);
  assert.throws(() => solveStreamtubeBody(system, { stepMethod: 'dogleg', initialTrustRadius: 0 }), /trust radius/);
  assert.throws(() => solveStreamtubeBody(system, { projectedSteps: true }), /Projected Euler/);
  assert.throws(() => solveStreamtubeBody(system, { stepMethod: 'dogleg', projectedSteps: 'yes' }), /Projected Euler/);
  assert.throws(() => solveStreamtubeBody(system, { stepMethod: 'dogleg', projectedSteps: false, secondOrderSteps: true }), /Second-order Euler/);
  const r = solveStreamtubeBody(system, { stepMethod: 'dogleg', maxIterations: 0 });
  assert.equal(r.converged, false); assert.equal(r.reason, 'iteration limit'); assert.ok(r.diagnostics.residual > 1e-10);
});

test('inactive Euler corner constraints leave ordinary trust-region steps unchanged', () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  const controls = { stepMethod: 'dogleg', initialTrustRadius: .01, maxIterations: 12, tolerance: 1e-10 };
  const plain = solveStreamtubeBody(createStreamtubeBodySystem(input), { ...controls, projectedSteps: false });
  const projected = solveStreamtubeBody(createStreamtubeBodySystem(input), controls);
  assert.deepEqual(projected.history, plain.history); assert.deepEqual(projected.x, plain.x);
  assert.equal(projected.converged, true); assert.equal(projected.projectedSteps, true);
});

test('trust-region relaxation recovers the failed cold lifting start and transfers it to actual Mach 0.2 Euler', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 8, tubes: 5, alpha: 2 }), flowModel: 'incompressible', mach: 0 };
  const system = createStreamtubeBodySystem(input), start = structuredClone(input);
  const relaxed = solveStreamtubeBody(system, { stepMethod: 'dogleg', tolerance: 1e-11, maxIterations: 25 });
  assert.equal(system.layout.n, 405); assert.equal(relaxed.converged, true, relaxed.reason);
  assert.deepEqual(input, start); assert.ok(relaxed.history.some(h => h.stepKind === 'gradient'));
  const before = structuredClone(relaxed);
  for (const streamwiseMode of ['momentum', 'isentropic']) {
    const target = { ...input, flowModel: 'compressible', mach: .2, streamwiseMode };
    const prepared = initializeStreamtubeBodyFromGrid(target, relaxed);
    assert.deepEqual(prepared.nodes, relaxed.nodes); assert.deepEqual(prepared.system.decode(prepared.initial).captured, relaxed.captured);
    const euler = solveStreamtubeBody(prepared.system, { initial: prepared.initial, tolerance: 1e-11, maxIterations: 20 });
    assert.deepEqual(relaxed, before); assert.equal(euler.converged, true, euler.reason); assert.equal(prepared.system.conditions.mach, .2);
    for (const [r, s, bodies] of [[relaxed, system, input.bodies], [euler, prepared.system, target.bodies]]) {
      assert.ok(r.diagnostics.residual < 1e-11); assert.ok(r.linearDiagnostics.maxRelativeResidual < 1e-10);
      assert.equal(streamtubeMeshSnapshot({ system: s, nodes: r.nodes }).quality.valid, true);
      const c = directBodyConservation(r, bodies, s.conditions);
      for (const k of r.streamwiseMode === 'momentum' ? [0, 1, 2, 3] : [0, 3]) assert.ok(Math.abs(c.balance[k]) < 2e-9);
      assert.ok(max(c.cutTraction) < 2e-9); assert.equal(r.surfaces.length, 4); assert.match(r.forceStatus, /Unvalidated/);
    }
    if (streamwiseMode === 'isentropic') {
      const p = stagnationPressureError(euler.sections.map(row => row.flat()), { gamma: prepared.system.conditions.gamma,
        referencePressure: prepared.system.conditions.pInf, freestreamMach: .2 });
      assert.ok(p.maxRelativeError < 2e-10);
    }
  }
});
