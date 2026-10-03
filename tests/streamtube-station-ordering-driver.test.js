// SPDX-License-Identifier: GPL-2.0-or-later
// One supplied two-body update: no panel solve, mesh startup or Mach sweep.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';

const plain = x => JSON.parse(JSON.stringify(x, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const fixture = JSON.parse(fs.readFileSync(new URL('../docs/coupled-current-profile-preparation/two-element-six-update/initial.json', import.meta.url))).checkpoint;
const originalURL = new URL('../src/euler/streamtube-coupled-ises.js', import.meta.url);
const archivedText = fs.readFileSync(new URL('../docs/linear-ordering-integration/before/streamtube-coupled-ises.js.txt', import.meta.url), 'utf8');
const archived = await import(`data:text/javascript;base64,${Buffer.from(archivedText.replace(/from\s+(['"])(\.[^'"]+)\1/g,
  (_, quote, path) => `from ${quote}${new URL(path, originalURL).href}${quote}`)).toString('base64')}`);
const run = (checkpoint, options = {}, solver = solveCoupledStreamtubeIses) => {
  const { iterationGeometry, stepAcceptance, stagnationLimiter } = checkpoint.continuation;
  return solver(undefined, { resume: checkpoint, iterationGeometry, stepAcceptance, stagnationLimiter,
    tolerance: 1e-10, maxIterations: 1, ...options });
};
let defaultResult;

test('default driver retains the complete archived result and checkpoint for a real two-body update', () => {
  const cp = plain(fixture), before = plain(cp);
  const old = run(cp, {}, archived.solveCoupledStreamtubeIses);
  defaultResult = run(cp);
  assert.equal(defaultResult.linearDiagnostics.solves, 1);
  assert.equal(defaultResult.history.length, 2, defaultResult.reason);
  assert.equal(defaultResult.mesh.quality.valid, true);
  assert.equal(defaultResult.gridAcceptance, 'convex');
  assert.deepEqual(plain(defaultResult), { ...plain(old), gridAcceptance: 'convex' });
  assert.equal(Object.hasOwn(defaultResult, 'linearOrdering'), false);
  assert.equal(Object.hasOwn(defaultResult.checkpoint.continuation, 'linearOrdering'), false);
  assert.deepEqual(cp, before);
});

test('station ordering executes coupled Newton and retains original accuracy and accepted-step decisions', t => {
  const cp = plain(fixture), published = [];
  cp.continuation.linearOrdering = 'station';
  cp.continuation.pivotTolerance = .001;
  const before = plain(cp), r = run(cp, { onCheckpoint: c => published.push(plain(c)) });
  assert.equal(r.x.length, 349);
  assert.equal(r.history.length, 2, r.reason);
  assert.equal(r.linearDiagnostics.solves, 1);
  assert.equal(r.mesh.quality.valid, true);
  assert.equal(r.linearOrdering, 'station');
  const d = r.linearDiagnostics.iterations[0];
  assert.equal(d.ordering, 'given');
  assert.equal(d.pivotTolerance, .001);
  assert.equal(d.attempts.length, 1);
  assert.equal(d.attempts[0].btf, false);
  assert.equal(d.stationOrdering.matchedNonzeroDiagonals, r.x.length);
  assert.equal(r.checkpoint.continuation.preferredOrdering, cp.continuation.preferredOrdering);
  assert(r.linearDiagnostics.maxRelativeResidual <= 1e-10);
  assert(published.length === 2 && published.every(c => c.continuation.linearOrdering === 'station'));
  assert.deepEqual(r.boundaryLayer.transitionState, defaultResult.boundaryLayer.transitionState);
  assert.equal(r.history[1].backtracks, defaultResult.history[1].backtracks);
  assert(Math.abs(r.history[1].step - defaultResult.history[1].step) < 1e-12);
  const maxDifference = Math.max(...r.x.map((v, i) => Math.abs(v - defaultResult.x[i])));
  assert(maxDifference < 1e-10, `Same Newton update differs by ${maxDifference}`);
  const replay = run(plain(r.checkpoint), { maxIterations: 0 });
  assert.equal(replay.linearDiagnostics.solves, 0);
  assert.deepEqual(plain(replay.checkpoint), plain(r.checkpoint));
  for (const key of ['x', 'residual', 'families']) assert.deepEqual(plain(replay[key]), plain(r[key]), key);
  assert.deepEqual(cp, before);
  t.diagnostic(JSON.stringify({ unknowns: r.x.length, maxDifference, linearRelativeResidual: r.linearDiagnostics.maxRelativeResidual,
    factorNonzeros: d.attempts[0].factorNonzeros, families: r.families }));
});

test('a true fresh supplied-state349 update selects guarded station ordering and agrees with explicit automatic ordering', t => {
  const saved = plain(fixture), before = plain(saved), f = saved.restart, c = saved.continuation;
  const options = { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL,
    iterationGeometry: c.iterationGeometry, stepAcceptance: c.stepAcceptance, stagnationLimiter: c.stagnationLimiter,
    ...(c.blUpdate === undefined ? {} : { blUpdate: c.blUpdate }),
    ...(c.projectionGeometry === undefined ? {} : { projectionGeometry: c.projectionGeometry }),
    maxIterations: 1, tolerance: 1e-10 };
  const first = [], fresh = solveCoupledStreamtubeIses(f.input, { ...options, onCheckpoint: cp => first.push(plain(cp)) });
  const automatic = solveCoupledStreamtubeIses(f.input, { ...options, linearOrdering: 'auto' });
  assert.equal(fresh.x.length, 349); assert.equal(fresh.history.length, 2, fresh.reason);
  assert.equal(automatic.history.length, 2, automatic.reason);
  assert.notEqual(fresh.initialRedistribution.resumed, true);
  assert.equal(fresh.initialRedistribution.accepted, true);
  assert.equal(fresh.linearOrdering, 'station-auto');
  assert.ok(first.every(cp => cp.continuation.linearOrdering === 'station-auto'));
  assert.equal(fresh.linearDiagnostics.solves, 1); assert.equal(automatic.linearDiagnostics.solves, 1);
  assert.equal(Object.hasOwn(automatic.checkpoint.continuation, 'linearOrdering'), false);
  assert.deepEqual(fresh.solverInput, automatic.solverInput);
  assert.deepEqual(fresh.coupledOptions, automatic.coupledOptions);
  assert.deepEqual(fresh.boundaryLayer.transitionState, automatic.boundaryLayer.transitionState);
  assert.equal(fresh.history[1].backtracks, automatic.history[1].backtracks);
  assert.ok(Math.abs(fresh.history[1].step - automatic.history[1].step) < 1e-12);
  const maximumStateDifference = Math.max(...fresh.x.map((v, i) => Math.abs(v - automatic.x[i])));
  assert.ok(maximumStateDifference < 1e-10, `Fresh Newton states differ by ${maximumStateDifference}`);
  assert.ok(fresh.linearDiagnostics.maxRelativeResidual <= 1e-10);
  assert.equal(fresh.mesh.quality.valid, true); assert.equal(automatic.mesh.quality.valid, true);
  assert.deepEqual(saved, before);
  t.diagnostic(JSON.stringify({ unknowns: fresh.x.length, maximumStateDifference,
    families: fresh.families, stationPolicy: fresh.linearDiagnostics.iterations[0].stationPolicy }));
});

test('ordering restart controls reject mismatches and malformed values before any checkpoint publication', () => {
  const station = plain(fixture); station.continuation.linearOrdering = 'station';
  let reports = 0;
  const onCheckpoint = () => reports++;
  for (const [cp, linearOrdering] of [[fixture, 'station'], [station, 'auto']])
    assert.throws(() => run(cp, { linearOrdering, maxIterations: 0, onCheckpoint }), /ordering controls do not match/);
  for (const linearOrdering of ['given', '', 1, false]) {
    const bad = plain(fixture); bad.continuation.linearOrdering = linearOrdering;
    assert.throws(() => run(bad, { maxIterations: 0, onCheckpoint }), /Invalid coupled ISES iteration controls/);
  }
  assert.equal(reports, 0);
});
