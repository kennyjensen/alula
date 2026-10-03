import test from 'node:test';
import assert from 'node:assert/strict';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { finishSubsonicPressureCoupling } from '../src/euler/streamtube-coupled-startup.js';

test('subsonic pressure handoff retains conditions and geometry and accounts for both solves', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    wakeGeometry: 'independent-banks', wakeOutlet: 'banks' };
  const source = solveCoupledStreamtubeIses(input, { edgeMatching: 'section-velocity',
    maxIterations: 8, tolerance: 1e-10, stepAcceptance: 'admissible', transitionMode: 'automatic',
    reynolds: 1e6, ncrit: 9 });
  assert.equal(source.converged, true, source.reason);
  const before = structuredClone(source), frames = [], checkpoints = [];
  const result = finishSubsonicPressureCoupling(source, { maxIterations: 20, tolerance: 1e-10,
    onIteration: h => frames.push(h), onCheckpoint: (cp, details) => checkpoints.push({ cp, details }) });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.checkpoint.restart.options.edgeMatching, 'pressure');
  assert.equal(result.checkpoint.restart.input.mach, input.mach);
  assert.equal(result.checkpoint.restart.options.reynolds, 1e6);
  assert.equal(result.checkpoint.restart.options.ncrit, 9);
  assert.deepEqual(result.solverInput.bodies, source.solverInput.bodies);
  assert.deepEqual(result.mesh.cells, source.mesh.cells);
  assert.equal(result.mesh.quality.valid, true);
  assert(result.history.at(-1).iteration <= 20);
  assert.deepEqual(result.history.map(h => h.iteration), Array.from({ length: result.history.length }, (_, i) => i));
  assert.equal(result.couplingStartup.startupIterations, source.history.length - 1);
  assert.equal(frames[0].iteration, source.history.length - 1);
  assert.equal(checkpoints.at(-1).details.history.at(-1).iteration, result.history.at(-1).iteration);
  assert(result.linearDiagnostics.solves >= source.linearDiagnostics.solves);
  assert.equal(result.checkpoint.continuation.linearOrdering, source.checkpoint.continuation.linearOrdering);
  const aligned = finishSubsonicPressureCoupling(source, { maxIterations: 20, linearOrdering: 'aligned-auto' });
  assert.equal(aligned.converged, true, aligned.reason);
  assert.equal(aligned.checkpoint.continuation.linearOrdering, 'aligned-auto');
  assert.equal(aligned.mesh.quality.valid, true);
  assert(aligned.history.at(-1).residual <= 1e-10);
  assert.deepEqual(source, before);
  assert.equal(finishSubsonicPressureCoupling(source, { maxIterations: source.history.length - 1 }), source);
  const sonic = { ...source, flow: { ...source.flow, diagnostics: { ...source.flow.diagnostics, maxMach: 1 } } };
  assert.equal(finishSubsonicPressureCoupling(sonic, { maxIterations: 20 }), sonic);
});
