// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { createHash } from 'node:crypto';
import { solveStreamtubeIses, continueEulerDissipationFromCheckpoint } from '../src/euler/streamtube-ises-update.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';

const serialize = v => JSON.parse(JSON.stringify(v, (_, x) => ArrayBuffer.isView(x) ? Array.from(x) : x));
const archivedPath = 'docs/ises-resume/streamtube-ises-update.before.js.txt';
const archived = fs.readFileSync(archivedPath, 'utf8');
assert.equal(createHash('sha256').update(archived).digest('hex'), '8772724c17bee48825767351e021d988e9d4748e66faf0fcf1ded4e81c9ea72e');
const original = await import('data:text/javascript;base64,' + Buffer.from(archived.replace(/from '([^']+)'/g,
  (_, p) => `from '${pathToFileURL(path.resolve('src/euler', p)).href}'`)).toString('base64'));
const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
const controls = { tolerance: 1e-11, stepAcceptance: 'admissible', maxIterations: 12 };
let complete, saved;

test('aligned Euler ordering retains the archived cold root and convex acceptance', () => {
  const old = original.solveStreamtubeIses(input, controls);
  complete = solveStreamtubeIses(input, controls);
  assert.equal(complete.converged, true, complete.reason);
  assert.equal('checkpoint' in complete, false);
  const { gridAcceptance, ...unchanged } = complete;
  assert.equal(gridAcceptance, 'convex');
  // Elimination order changes floating-point rounding and factor diagnostics,
  // but must preserve the physical root and accepted iteration sequence.
  const compare = (a, b) => {
    if (typeof b === 'number') {
      assert.ok(Math.abs(a-b) <= 1e-10*Math.max(1,Math.abs(b)), `${a} versus ${b}`);
    } else if (b && typeof b === 'object') {
      assert.deepEqual(Object.keys(a), Object.keys(b));
      for (const key of Object.keys(b)) compare(a[key], b[key]);
    } else assert.equal(a,b);
  };
  const { linearDiagnostics: currentLinear, ...current } = serialize(unchanged);
  const { linearDiagnostics: oldLinear, ...previous } = serialize(old);
  compare(current, previous);
  assert.equal(currentLinear.solves, oldLinear.solves);
  assert.ok(currentLinear.maxRelativeResidual <= 1e-10);
  assert.ok(currentLinear.maxFactorNonzeros > 0);
});

test('detached atomic checkpoints resume the exact maintained state and interrupted trajectory', () => {
  const sourceBefore = serialize(input), stop = new Error('cancel after accepted checkpoint');
  assert.throws(() => solveStreamtubeIses(input, { ...controls, onCheckpoint: (checkpoint, details) => {
    saved = serialize(checkpoint);
    const iteration = details.history.at(-1).iteration;
    // Mutations must not modify subsequent values, history, controls or banks.
    checkpoint.input.mach = 99; checkpoint.initialEuler.x.fill(99);
    checkpoint.initialEuler.nodes[0][0][0].x = 99;
    checkpoint.residual.fill(99); checkpoint.continuation.fractions[0][0] = 99;
    checkpoint.continuation.lastRedistributedStagnation[0] = 99;
    details.history.length = 0; details.initialRedistribution.passages.length = 0;
    if (iteration === 2) throw stop;
  } }), error => error === stop);
  const before = serialize(saved);
  assert.deepEqual(serialize(input), sourceBefore);
  const zero = solveStreamtubeIses(undefined, { ...controls, resume: saved, maxIterations: 0 });
  assert.equal(zero.initialRedistribution.resumed, true);
  assert.deepEqual(zero.initialRedistribution.passages, []);
  assert.equal(zero.linearDiagnostics.solves, 0);
  assert.deepEqual(zero.checkpoint, saved);
  assert.deepEqual(serialize(zero.nodes), saved.initialEuler.nodes);
  assert.deepEqual(serialize(zero.residual), saved.residual);
  const rest = solveStreamtubeIses(undefined, { ...controls, resume: saved, maxIterations: 10 });
  assert.equal(rest.converged, true, rest.reason);
  assert.deepEqual(serialize(saved), before);
  for (const key of ['x', 'nodes', 'residual', 'diagnostics', 'surfaces', 'diagnosticForces', 'finalQuality'])
    assert.deepEqual(serialize(rest[key]), serialize(complete[key]), key);
  const withoutIndex = h => { const { iteration, ...rest } = h; return rest; };
  assert.deepEqual(serialize(rest.history.slice(1).map(withoutIndex)), serialize(complete.history.slice(3).map(withoutIndex)));
  assert.ok(rest.history.slice(1).every(h => h.maintenance.inlet));
  assert.deepEqual(rest.checkpoint.continuation.fractions, saved.continuation.fractions);
});

test('resume rejects mismatched controls, corrupted residual, nodes and maintenance history', () => {
  assert.throws(() => solveStreamtubeIses(input, { ...controls, resume: saved }), /no separate initial state/);
  assert.throws(() => solveStreamtubeIses(undefined, { ...controls, resume: saved, iterationGeometry: 'ises-sampled' }), /controls do not match/);
  const broken = mutate => { const copy = structuredClone(saved); mutate(copy); return copy; };
  assert.throws(() => solveStreamtubeIses(undefined, { ...controls, resume: broken(c => c.residual[0] += 1) }), /residual does not replay exactly/);
  assert.throws(() => solveStreamtubeIses(undefined, { ...controls, resume: broken(c => c.initialEuler.nodes[0][1][1].y += 1e-4) }), /does not replay exactly/);
  assert.throws(() => solveStreamtubeIses(undefined, { ...controls, resume: broken(c => c.continuation.fractions[0][1] = 0) }), /inlet fractions/);
  assert.throws(() => solveStreamtubeIses(undefined, { ...controls, resume: broken(c => c.continuation.lastRedistributedStagnation[0] = -1) }), /redistribution history/);
  const retained = solveStreamtubeIses(undefined, { ...controls, resume: saved, maxIterations: 0, retainCheckpoint: true });
  assert.deepEqual(retained.checkpoint, saved);
});

test('Euler dissipation startup restores the requested equations and preserves its source checkpoint', t => {
  const source = JSON.parse(fs.readFileSync(new URL(
    '../docs/solver-reliability/rae-inviscid-audit/fixtures/rae32-dissipation-restored-root.json', import.meta.url)));
  const before = structuredClone(source), published = [];
  const corrupt = structuredClone(source); corrupt.residual[0] += 1;
  assert.throws(() => continueEulerDissipationFromCheckpoint(corrupt), /residual does not replay exactly/);
  const limited = continueEulerDissipationFromCheckpoint(source, { maxIterations: 1, tolerance: 1e-10 });
  assert.equal(limited.converged, false, 'A temporary-law or budget-limited state is not a requested-equation root.');
  assert.ok(limited.history.length <= 2);
  const result = continueEulerDissipationFromCheckpoint(source, { maxIterations: 20, tolerance: 1e-10,
    onCheckpoint: checkpoint => { published.push(structuredClone(checkpoint)); checkpoint.residual.fill(99); },
  });
  assert.equal(result.converged, true, result.reason);
  assert.deepEqual(source, before);
  assert.deepEqual(result.checkpoint.input.upwind, source.input.upwind);
  assert.equal(result.checkpoint.input.mach, .74);
  assert.equal(result.checkpoint.input.alpha, 2.68);
  assert.ok(result.history.some(h => h.mucon < 0));
  assert.ok(result.history.some(h => h.mucon > 0));
  assert.ok(result.history.every(h => h.mcrit >= .75 && h.mcrit <= .99 && h.nextMcrit >= .75 && h.nextMcrit <= .99));
  assert.ok(result.history.length <= 40);
  assert.deepEqual(published.at(-1), result.checkpoint);
  const replay = solveStreamtubeIses(undefined, { resume: result.checkpoint,
    ...result.checkpoint.continuation, maxIterations: 0, tolerance: 1e-10 });
  assert.equal(replay.converged, true, replay.reason);
  assert.equal(replay.finalQuality.valid, true);
  assert.ok(replay.diagnostics.residual <= 1e-10);
  t.diagnostic(JSON.stringify({ iterations: result.history.length, residual: replay.diagnostics.residual }));
});

test('adaptive MCRIT checkpoints preserve the requested law and the interrupted Newton trajectory', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    streamwiseMode: 'momentum', upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
  const before = structuredClone(input), controls = { maxIterations: 12, tolerance: 1e-11,
    adaptiveMcrit: true, stepAcceptance: 'armijo', retainCheckpoint: true };
  const complete = solveStreamtubeIses(input, controls);
  const partial = solveStreamtubeIses(input, { ...controls, maxIterations: 1 });
  const checkpoint = structuredClone(partial.checkpoint);
  assert.ok(checkpoint.input.upwind.mcrit < .99);
  assert.equal(checkpoint.continuation.targetMcrit, .99);
  const zero = solveStreamtubeIses(undefined, { ...controls, ...checkpoint.continuation, resume: checkpoint, maxIterations: 0 });
  assert.deepEqual(zero.checkpoint, checkpoint);
  const rest = solveStreamtubeIses(undefined, { ...controls, ...checkpoint.continuation, resume: checkpoint });
  assert.equal(rest.converged, true, rest.reason);
  assert.deepEqual(rest.checkpoint, complete.checkpoint);
  assert.deepEqual(rest.history.slice(1).map(({ iteration, ...h }) => h),
    complete.history.slice(2).map(({ iteration, ...h }) => h));
  assert.deepEqual(partial.checkpoint, checkpoint);
  assert.deepEqual(input, before);
  assert.throws(() => solveStreamtubeIses(undefined, { ...controls, ...checkpoint.continuation,
    resume: checkpoint, adaptiveMcrit: false }), /controls do not match/);

  // A root under a temporary law must not be called a requested-law root,
  // even when the small subsonic dissipation difference is below tolerance.
  const temporary = structuredClone(complete.checkpoint);
  temporary.input.upwind.mcrit = .8;
  const system = createStreamtubeBodySystem(temporary.input);
  const state = system.adoptGeometry(Float64Array.from(temporary.initialEuler.x), temporary.initialEuler.nodes);
  temporary.residual = Array.from(system.residual(state));
  const limited = solveStreamtubeIses(undefined, { ...controls, ...temporary.continuation, resume: temporary, maxIterations: 0 });
  assert.equal(limited.converged, false);
  assert.match(limited.reason, /requested MCRIT/);
  const restored = solveStreamtubeIses(undefined, { ...controls, ...temporary.continuation, resume: temporary, maxIterations: 1 });
  assert.equal(restored.converged, true, restored.reason);
  assert.equal(restored.checkpoint.input.upwind.mcrit, .99);
});

for (const broadShockMucon of [2, 4]) test(`broad shock startup ${broadShockMucon} checkpoints retain their law and cannot certify it as the target`, () => {
  const input = { ...intrinsicBodyFixture({ bodySegments: 4, tubes: 2 }),
    streamwiseMode: 'momentum', upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
  const before = structuredClone(input), controls = { maxIterations: 1, tolerance: 1e-11,
    adaptiveMcrit: true, firstOrderStartup: true, broadShockStartup: true, broadShockMucon, stopOnGridStagnation: true,
    gridCorrectionBacktracking: 'before-damping', stepAcceptance: 'armijo', retainCheckpoint: true };
  const partial = solveStreamtubeIses(input, controls), checkpoint = partial.checkpoint;
  assert.equal(partial.converged, false);
  assert.equal(checkpoint.input.upwind.mucon, -broadShockMucon);
  assert.equal(checkpoint.input.upwind.mcrit, .75);
  assert.equal(checkpoint.continuation.targetMucon, 1);
  assert.equal(checkpoint.continuation.targetMcrit, .99);
  const zero = solveStreamtubeIses(undefined, { ...checkpoint.continuation, resume: checkpoint,
    maxIterations: 0, tolerance: 1e-11, retainCheckpoint: true });
  assert.equal(zero.converged, false);
  assert.deepEqual(zero.checkpoint, checkpoint);
  for (const mismatch of [{ broadShockMucon: broadShockMucon === 2 ? 4 : 2 }, { stopOnGridStagnation: false }, { broadShockStartup: false }, { gridCorrectionBacktracking: 'after-search' }])
    assert.throws(() => solveStreamtubeIses(undefined, { ...checkpoint.continuation,
      resume: checkpoint, ...mismatch, maxIterations: 0 }), /controls do not match/);
  assert.deepEqual(input, before);
});

test('first-order startup resumes its phase and must restore the requested second-order law', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    streamwiseMode: 'momentum', upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
  const before = structuredClone(input), controls = { maxIterations: 16, tolerance: 1e-11,
    adaptiveMcrit: true, firstOrderStartup: true, stepAcceptance: 'armijo', retainCheckpoint: true };
  const complete = solveStreamtubeIses(input, controls);
  assert.equal(complete.converged, true, complete.reason);
  assert.deepEqual(complete.solverInput.upwind, input.upwind);
  assert.equal(complete.history[0].dissipation.mucon, -1);
  assert.ok(complete.history.some(h => h.dissipation?.restoredSecondOrder));
  assert.ok(complete.linearDiagnostics.solves <= controls.maxIterations);

  const partial = solveStreamtubeIses(input, { ...controls, maxIterations: 1 });
  const checkpoint = structuredClone(partial.checkpoint);
  assert.equal(partial.converged, false);
  assert.equal(checkpoint.input.upwind.mucon, -1);
  assert.equal(checkpoint.continuation.targetMucon, 1);
  const replay = (saved, maxIterations) => solveStreamtubeIses(undefined,
    { ...controls, ...saved.continuation, resume: saved, maxIterations });
  assert.deepEqual(replay(checkpoint, 0).checkpoint, checkpoint);
  const rest = replay(checkpoint, controls.maxIterations - 1);
  assert.equal(rest.converged, true, rest.reason);
  assert.deepEqual(rest.checkpoint, complete.checkpoint);
  assert.deepEqual(rest.history.slice(1).map(({ iteration, ...h }) => h),
    complete.history.slice(2).map(({ iteration, ...h }) => h));
  assert.deepEqual(replay(rest.checkpoint, 0).checkpoint, complete.checkpoint,
    'Resuming the requested-law root must not restart first-order startup.');
  assert.deepEqual(input, before);
  assert.deepEqual(partial.checkpoint, checkpoint);
  for (const mismatch of [{ firstOrderStartup: false }, { targetMucon: 2 }])
    assert.throws(() => solveStreamtubeIses(undefined, { ...controls, ...checkpoint.continuation,
      resume: checkpoint, ...mismatch }), /controls do not match/);

  // A small residual alone cannot certify completion of a temporary phase.
  const temporary = structuredClone(complete.checkpoint);
  temporary.input.upwind.mucon = -1;
  const system = createStreamtubeBodySystem(temporary.input);
  const state = system.adoptGeometry(Float64Array.from(temporary.initialEuler.x), temporary.initialEuler.nodes);
  temporary.residual = Array.from(system.residual(state));
  const limited = replay(temporary, 0);
  assert.equal(limited.converged, false);
  assert.equal(limited.residualConverged, false);
  assert.match(limited.reason, /second-order dissipation/);
  const snapshots = [], events = [], checkpoints = [];
  const restored = solveStreamtubeIses(undefined, { ...controls, ...temporary.continuation,
    resume: temporary, maxIterations: 1, onIteration: h => events.push(structuredClone(h)),
    onMesh: frame => snapshots.push(structuredClone(frame.iteration)),
    onCheckpoint: checkpoint => checkpoints.push(checkpoint) });
  assert.equal(restored.converged, true, restored.reason);
  assert.equal(restored.checkpoint.input.upwind.mucon, 1);
  assert.equal(events.at(-1).dissipation.mucon, 1);
  assert.deepEqual(events.at(-1), snapshots.at(-1));
  assert.deepEqual(checkpoints.at(-1), restored.checkpoint);
  const cancel = new Error('cancel after restoring the dissipation law');
  assert.throws(() => solveStreamtubeIses(undefined, { ...controls, ...temporary.continuation,
    resume: temporary, maxIterations: 1, onCheckpoint: checkpoint => {
      if (checkpoint.input.upwind.mucon > 0) throw cancel;
    } }), error => error === cancel);
});
