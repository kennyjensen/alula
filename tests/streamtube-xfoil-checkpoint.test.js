import test from 'node:test';
import assert from 'node:assert/strict';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { coupledAssemblyRestart } from '../src/euler/streamtube-coupled-assembly.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';

const serial = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));

test('a real two-element automatic XFOIL root retains its policy and exact complete restart', t => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  const inputBefore = serial(input), observed = [];
  const controls = { reynolds: 1e6, ncrit: 9, transitionMode: 'automatic',
    edgeMatching: 'section-velocity', stepAcceptance: 'admissible',
    blUpdate: 'xfoil', maxIterations: 12, tolerance: 1e-10 };
  const result = solveCoupledStreamtubeIses(input, { ...controls,
    onCheckpoint: checkpoint => observed.push(serial(checkpoint)) });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.mesh.quality.valid, true);
  assert.ok(Math.max(...result.residual.map(Math.abs)) < controls.tolerance);
  assert.ok(result.history.length > 1 && result.history.length <= controls.maxIterations + 1);
  assert.equal(result.linearDiagnostics.solves, result.history.length - 1);
  assert.equal(result.boundaryLayer.surfaces.length, 4);
  assert.equal(result.boundaryLayer.wakes.length, 2);
  assert.equal(result.blUpdate, 'xfoil');
  assert.equal(result.checkpoint.continuation.blUpdate, 'xfoil');
  assert.equal(observed.length, result.history.length);
  assert.ok(observed.every(c => c.continuation.blUpdate === 'xfoil'));
  const checkpoint = serial(result.checkpoint), before = serial(checkpoint);
  assert.deepEqual(observed.at(-1), checkpoint);

  // Omitted policy inherits the checkpoint, avoiding accidental Giles
  // updates after serialization. Replaying a root performs no linear solve
  // and skips the already-completed initial grid redistribution.
  const { blUpdate, ...resumeControls } = controls;
  const zero = solveCoupledStreamtubeIses(undefined, { ...resumeControls, resume: checkpoint, maxIterations: 0 });
  assert.equal(zero.converged, true); assert.equal(zero.blUpdate, 'xfoil');
  assert.equal(zero.linearDiagnostics.solves, 0);
  assert.equal(zero.initialRedistribution.resumed, true);
  assert.deepEqual(zero.initialRedistribution.passages, []);
  for (const key of ['x', 'residual', 'families']) assert.deepEqual(serial(zero[key]), serial(result[key]), key);
  assert.deepEqual(serial(zero.flow.nodes), serial(result.flow.nodes));
  assert.deepEqual(serial(zero.flow.undisplacedNodes), serial(result.flow.undisplacedNodes));
  assert.deepEqual(serial(zero.checkpoint), checkpoint);
  assert.deepEqual(serial(zero.boundaryLayer.transitionState), serial(result.boundaryLayer.transitionState));
  assert.throws(() => solveCoupledStreamtubeIses(undefined,
    { ...resumeControls, resume: checkpoint, maxIterations: 0, blUpdate: 'giles' }), /BL update controls do not match/);
  assert.deepEqual(serial(checkpoint), before); assert.deepEqual(serial(input), inputBefore);
  // An unfinished initial redistribution has no checkpoint. Recover the
  // complete current layout without borrowing any previous mesh dimensions.
  const fallback = coupledAssemblyRestart({ ...result, checkpoint: undefined });
  const rebuilt = createCoupledStreamtubeBody(fallback.input, { ...fallback.options,
    initialEuler: fallback.initialEuler, initialBL: fallback.initialBL });
  assert.deepEqual(rebuilt.initial, result.x);
  assert.deepEqual(rebuilt.evaluate(rebuilt.initial).residual, result.residual);
  assert.equal(fallback.initialEuler.x.length, rebuilt.ne);
  assert.equal(fallback.initialBL.length, 4 * rebuilt.bl.stations.length);
  assert.equal(coupledAssemblyRestart(result), result.checkpoint.restart);
  t.diagnostic(JSON.stringify({ unknowns: result.x.length, iterations: result.history.length - 1,
    families: result.families, checkpointPolicy: checkpoint.continuation.blUpdate, zeroReplayLinearSolves: zero.linearDiagnostics.solves }));
});
