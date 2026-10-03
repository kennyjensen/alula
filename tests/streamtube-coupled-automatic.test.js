// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { createHash } from 'node:crypto';
import { pathToFileURL } from 'node:url';
import { solveCoupledStreamtubeAutomatic } from '../src/euler/tests/streamtube-coupled-automatic.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const tolerance = 1e-10;
const serial = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
let cached, sourceError;
function source() {
  if (sourceError) throw sourceError;
  if (!cached) try {
    const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2, mach: .2 }),
      streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 },
      upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
    const result = solveCoupledStreamtubeIses(input, { edgeMatching: 'section-velocity',
      blThermodynamics: 'historical-common-isentrope', transitionMode: 'automatic',
      maxIterations: 12, tolerance, stepAcceptance: 'admissible', iterationRecovery: false });
    assert.equal(result.converged, true, result.reason);
    assert.equal(result.mesh.quality.valid, true);
    assert.ok(result.x.length < 500);
    cached = { result, checkpoint: serial(result.checkpoint) };
  } catch (error) { sourceError = error; throw error; }
  return { result: cached.result, initialCheckpoint: structuredClone(cached.checkpoint) };
}

test('same-Mach automatic coupled entry exactly retains a complete accepted state with no stage or LU', () => {
  const { result: original, initialCheckpoint } = source(), before = structuredClone(initialCheckpoint);
  const result = solveCoupledStreamtubeAutomatic(.2, { initialCheckpoint, stageMaxIterations: 0 });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.linearDiagnostics.solves, 0);
  assert.equal(result.continuation.attempts.length, 0);
  assert.equal(result.continuation.reachedTarget, true);
  assert.deepEqual(result.checkpoint, initialCheckpoint);
  assert.deepEqual(result.residual, original.residual);
  assert.deepEqual(serial(result.x), serial(original.x));
  assert.deepEqual(result.boundaryLayer.transitionState, original.boundaryLayer.transitionState);
  assert.deepEqual(result.flow.nodes, original.flow.nodes);
  assert.deepEqual(initialCheckpoint, before);
  result.checkpoint.restart.initialBL[0] = -5;
  assert.deepEqual(initialCheckpoint, before);
});

test('automatic coupled continuation reaches a new Mach with all four BLs and two wakes and detached live output', t => {
  const { initialCheckpoint } = source(), before = structuredClone(initialCheckpoint);
  const seen = { stages: [], iterations: [], meshes: 0, checkpoints: 0, accepted: 0 };
  const result = solveCoupledStreamtubeAutomatic(.21, { initialCheckpoint, stageMaxIterations: 8, maxStages: 4,
    includeFlowState: true,
    onStage: info => { seen.stages.push(info.mach); info.mach = 99; },
    onIteration: entry => { seen.iterations.push(entry.iteration); entry.residual = 99; },
    onMesh: (snapshot, state) => {
      seen.meshes++; assert.equal(snapshot.mach, .21); assert.ok(snapshot.mesh.vertices.length);
      assert.equal(snapshot.system, undefined);
      snapshot.mesh.vertices[0].x = 99;
      assert.ok(state.flow.sections[0][0][0].p > 0);
      assert.equal(state.iteration.iteration, snapshot.iteration.iteration);
      state.flow.nodes[0][0][0].x = -99;
      state.coupledFamilies.euler = 99;
    },
    onCheckpoint: (checkpoint, details) => {
      seen.checkpoints++; if (details.kind === 'accepted') seen.accepted++;
      assert.equal(checkpoint.restart.input.mach, .21);
      checkpoint.restart.initialBL[0] = -99;
      checkpoint.continuation.fractions[0][1] = -99;
    },
  });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.conditions.mach, .21);
  assert.equal(result.continuation.currentMach, .21);
  assert.equal(result.continuation.reachedTarget, true);
  assert.equal(result.mesh.quality.valid, true);
  assert.ok(Math.max(...result.residual.map(Math.abs)) < tolerance);
  assert.equal(result.boundaryLayer.surfaces.length, 4);
  assert.equal(result.boundaryLayer.wakes.length, 2);
  assert.equal(result.checkpoint.restart.options.transitionMode, 'automatic');
  assert.equal(result.checkpoint.restart.options.blThermodynamics, 'historical-common-isentrope');
  assert.equal(result.initialRedistribution.resumed, true);
  assert.deepEqual(result.initialRedistribution.passages, []);
  assert.deepEqual(initialCheckpoint, before);
  assert.deepEqual(seen.stages, [.21]);
  assert.ok(seen.iterations.length > 1);
  assert.equal(seen.meshes, seen.iterations.length);
  assert.equal(seen.checkpoints, seen.iterations.length + 1);
  assert.equal(seen.accepted, 1);
  assert.equal(result.physicalAcceptance, false);
  assert.equal(result.fullSolverComplete, false);
  t.diagnostic(JSON.stringify({ unknowns: result.x.length, updates: result.history.length - 1,
    residual: Math.max(...result.residual.map(Math.abs)), stages: result.continuation.attempts.length }));
});

test('omitting detached physical-state callbacks retains the entire archived automatic result exactly', async () => {
  const { initialCheckpoint } = source();
  const [entry] = JSON.parse(fs.readFileSync('docs/coupled-automatic/before-flow-state.json')).files;
  const bytes = fs.readFileSync(entry.archive);
  assert.equal(createHash('sha256').update(bytes).digest('hex'), entry.sha256);
  const code = bytes.toString().replace(/from '([^']+)'/g, (_, specifier) => {
    const original = path.resolve(path.dirname(entry.path), specifier);
    const resolved = fs.existsSync(original) ? original : path.resolve('src/euler/tests', specifier);
    return `from '${pathToFileURL(resolved).href}'`;
  });
  const old = await import(`data:text/javascript;base64,${Buffer.from(code).toString('base64')}`);
  const options = { initialCheckpoint, stageMaxIterations: 8, maxStages: 4 };
  const before = old.solveCoupledStreamtubeAutomatic(.21, options);
  const after = solveCoupledStreamtubeAutomatic(.21, { ...options,
    onMesh: (snapshot, state) => { assert.equal(state, undefined); } });
  // Wall clocks are observational and differ between independent runs.
  const clocks = new Set(['jacobianMilliseconds', 'matchingMilliseconds', 'stationSolveMilliseconds', 'autoSolveMilliseconds', 'totalMilliseconds']);
  const stripClocks = value => ArrayBuffer.isView(value) ? value
    : Array.isArray(value) ? value.map(stripClocks)
    : value && typeof value === 'object' ? Object.fromEntries(Object.entries(value).filter(([k]) => !clocks.has(k)).map(([k,v]) => [k,stripClocks(v)])) : value;
  assert.deepEqual(stripClocks(after), stripClocks(before));
});

test('failed increments halve from the last accepted checkpoint and cannot relabel its Mach as the target', () => {
  const { initialCheckpoint } = source(), before = structuredClone(initialCheckpoint);
  const result = solveCoupledStreamtubeAutomatic(.21, { initialCheckpoint, stageMaxIterations: 0,
    maxSubdivisions: 2, maxStages: 8 });
  assert.equal(result.converged, false);
  assert.equal(result.continuation.reachedTarget, false);
  assert.equal(result.status, 'research-coupled-target-not-reached');
  assert.equal(result.stateConverged, true);
  assert.equal(result.conditions.mach, .2);
  assert.equal(result.continuation.currentMach, .2);
  assert.equal(result.continuation.targetMach, .21);
  assert.equal(result.continuation.attempts.length, 3);
  for (const [i, attempt] of result.continuation.attempts.entries()) {
    assert.ok(Math.abs(attempt.mach - (.2 + .01 / 2 ** i)) < 1e-15);
    assert.equal(attempt.lastAcceptedMach, .2);
    assert.equal(attempt.converged, false);
    assert.equal(attempt.globalLinearSolves, 0);
    assert.equal(attempt.initialRedistributionSkipped, true);
  }
  assert.match(result.reason, /exhausted 2 subdivisions/);
  assert.deepEqual(result.checkpoint, before);
  assert.deepEqual(initialCheckpoint, before);
});

test('stage cap is independent of subdivision budget and reverse Mach continuation retains its sign', () => {
  const { initialCheckpoint } = source();
  const result = solveCoupledStreamtubeAutomatic(.19, { initialCheckpoint, stageMaxIterations: 0,
    maxSubdivisions: 5, maxStages: 2 });
  assert.equal(result.converged, false);
  assert.equal(result.continuation.attempts.length, 2);
  assert.equal(result.continuation.attempts[0].mach, .19);
  assert.equal(result.continuation.attempts[1].mach, .195);
  assert.match(result.reason, /2-stage limit/);
  assert.deepEqual(result.checkpoint, initialCheckpoint);
});

test('observer cancellation propagates without starting recovery, including after a committed Newton step', () => {
  const { initialCheckpoint } = source(), before = structuredClone(initialCheckpoint);
  for (const hook of ['onStage', 'onIteration', 'onMesh', 'onCheckpoint']) {
    const cancelled = new Error(`cancel from ${hook}`); let stages = 0;
    const options = { initialCheckpoint, stageMaxIterations: 8, onStage: () => { stages++; } };
    options[hook] = (value, details) => {
      if (hook === 'onIteration' && value.iteration === 0) return;
      if (hook === 'onCheckpoint' && details.history.at(-1).iteration === 0) return;
      if (hook === 'onStage') stages++;
      throw cancelled;
    };
    assert.throws(() => solveCoupledStreamtubeAutomatic(.21, options), error => error === cancelled);
    assert.equal(stages, 1);
    assert.deepEqual(initialCheckpoint, before);
  }
  let caught = false, stages = 0;
  try {
    solveCoupledStreamtubeAutomatic(.21, { initialCheckpoint, onStage: () => { stages++; },
      onIteration: () => { throw null; } });
  } catch (error) { caught = true; assert.equal(error, null); }
  assert.equal(caught, true);
  assert.equal(stages, 1);
});

test('invalid controls, stale source and changed equations fail before any target attempt', () => {
  const { initialCheckpoint } = source(); let stages = 0;
  for (const bad of [{ targetMach: 1 }, { stageMaxIterations: -1 }, { maxStages: 0 }, { maxStages: 257 },
    { maxSubdivisions: 21 }, { maxBacktracks: 21 }, { tolerance: NaN }, { onMesh: true }, { includeFlowState: 'yes' }]) {
    const { targetMach = .21, ...options } = bad;
    assert.throws(() => solveCoupledStreamtubeAutomatic(targetMach, { initialCheckpoint, ...options }), /controls/);
  }
  for (const mutate of [c => { c.families.euler = .001; },
    c => { c.restart.initialEuler.x[0] += .001; },
    c => { c.restart.input.streamwiseMode = 'isentropic'; },
    c => { c.continuation.fractions[0][1] = -1; }]) {
    const bad = structuredClone(initialCheckpoint); mutate(bad);
    assert.throws(() => solveCoupledStreamtubeAutomatic(.21, { initialCheckpoint: bad, onStage: () => { stages++; } }));
  }
  assert.equal(stages, 0);
});

test('bounded Mach steps reach the target from successive converged roots in either direction', () => {
  for (const target of [.21, .19]) {
    const { initialCheckpoint } = source();
    const result = solveCoupledStreamtubeAutomatic(target, { initialCheckpoint, maxMachStep: .004,
      stageMaxIterations: 8, maxStages: 6 });
    assert.equal(result.converged, true, result.reason);
    let previous = .2;
    for (const attempt of result.continuation.attempts) {
      assert.ok(Math.abs(attempt.mach - previous) <= .004 + 1e-15);
      assert.equal(attempt.converged, true);
      previous = attempt.mach;
    }
    assert.equal(previous, target);
    assert.ok(result.continuation.attempts.length >= 3);
  }
});
