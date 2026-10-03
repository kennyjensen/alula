import test from 'node:test';
import assert from 'node:assert/strict';
import { solveStreamtubeIses } from '../src/euler/streamtube-ises-update.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directBodyConservation } from './oracles/streamtube-body.js';
import { streamtubeEquationOrder } from '../src/euler/streamtube-body-layout.js';
import { solveSparseDirect, solveSparseDirectAligned } from '../src/numerics/klu.js';
import { sparseProduct } from '../src/numerics/sparse.js';

test('physical equation alignment reduces Euler LU fill without changing the Newton variables', () => {
  const system = createStreamtubeBodySystem(intrinsicBodyFixture({ elements: 2, bodySegments: 16, tubes: 4 }));
  const matrix = system.jacobian(system.initial, { sparse: true });
  const exact = Float64Array.from({ length: matrix.n }, (_, i) => Math.sin(.37*i));
  const rhs = sparseProduct(matrix, exact), before = structuredClone(matrix);
  const original = solveSparseDirect(matrix, rhs);
  const aligned = solveSparseDirectAligned(matrix, rhs, streamtubeEquationOrder(system.layout));
  assert.equal(aligned.equationOrdering, 'aligned');
  assert.ok(aligned.factorNonzeros < .8*original.factorNonzeros);
  assert.ok(aligned.relativeResidual <= 1e-10);
  for (let i=0; i<exact.length; i++) assert.ok(Math.abs(aligned.x[i]-exact[i]) < 1e-8);
  assert.deepEqual(matrix, before);
});

test('complete research ISES update sequence closes a conservative two-element root and publishes its grids', () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), before = structuredClone(input), meshes = [];
  const r = solveStreamtubeIses(input, { tolerance: 1e-11, onMesh: m => meshes.push({ iteration: m.iteration, nodes: m.nodes }) });
  assert.deepEqual(input, before); assert.equal(r.converged, true, r.reason); assert.equal(r.surfaces.length, 4);
  assert.equal(r.initialRedistribution.accepted, true); assert.equal(r.initialRedistribution.passages.length, 3);
  assert.ok(r.initialRedistribution.passages.every(p => p.pairs === 5));
  assert.ok(r.history.slice(1).every(h => h.stepKind === 'density-newton' && h.maintenance.inlet));
  assert.deepEqual(meshes.map(m => m.iteration.iteration), r.history.map(h => h.iteration));
  assert.deepEqual(meshes.at(-1).nodes, r.nodes);
  assert.ok(r.linearDiagnostics.maxRelativeResidual <= 1e-10);
  const system = createStreamtubeBodySystem(r.solverInput), balance = directBodyConservation(r, input.bodies, system.conditions);
  for (const v of [...balance.balance, ...balance.cutTraction]) assert.ok(Math.abs(v) < 2e-9);
  const resumed = solveStreamtubeIses(input, { initialEuler: { x: system.initial, nodes: system.decode(system.initial).nodes }, tolerance: 1e-11 });
  assert.equal(resumed.converged, true); assert.ok(Math.abs(resumed.diagnostics.residual - r.diagnostics.residual) < 1e-11);
});

test('research update retains explicit limits and rejects incompatible motion controls', () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  const r = solveStreamtubeIses(input, { maxIterations: 0 });
  assert.equal(r.converged, false); assert.equal(r.reason, 'iteration limit'); assert.equal(r.history.length, 1);
  assert.equal(r.initialRedistribution.accepted, true); assert.equal(r.linearDiagnostics.solves, 0);
  assert.throws(() => solveStreamtubeIses(input, { maxIterations: -1 }), /controls/);
  assert.throws(() => solveStreamtubeIses({ ...input, stagnationMotion: 'unknown' }), /stagnation motion/);
});

test('interpolated stagnation motion converges and replays a two-element Euler checkpoint', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), stagnationMotion: 'interpolated' };
  const before = structuredClone(input);
  const controls = { stepAcceptance: 'armijo', iterationGeometry: 'ises-sampled', tolerance: 1e-10 };
  const r = solveStreamtubeIses(input, { ...controls, maxIterations: 30, retainCheckpoint: true });
  assert.equal(r.converged, true, r.reason);
  assert.equal(r.finalQuality.valid, true);
  assert.equal(r.solverInput.stagnationMotion, 'interpolated');
  assert.deepEqual(input, before);
  const resumed = solveStreamtubeIses(undefined, { ...controls, resume: r.checkpoint, maxIterations: 0 });
  assert.equal(resumed.converged, true);
  assert.deepEqual(resumed.residual, r.residual);
  assert.deepEqual(resumed.nodes, r.nodes);
});
