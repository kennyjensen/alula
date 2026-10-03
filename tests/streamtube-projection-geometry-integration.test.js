// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { transferStreamtubeGeometry } from '../src/euler/streamtube-geometry.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';

const copy = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const retained = JSON.parse(fs.readFileSync(new URL(
  '../docs/coupled-current-profile-preparation/two-element-six-update/initial.json', import.meta.url))).checkpoint;
const driverURL = new URL('../src/euler/streamtube-coupled-ises.js', import.meta.url);
const archiveURL = new URL('../docs/coupled-projection-grid-response/runtime-draft/streamtube-coupled-ises.before.js.txt', import.meta.url);
const archivedSource = fs.readFileSync(archiveURL, 'utf8').replace(/from (['"])(\.[^'"]+)\1/g,
  (_, quote, path) => `from ${quote}${new URL(path, driverURL).href}${quote}`);
const previous = await import(`data:text/javascript;base64,${Buffer.from(archivedSource).toString('base64')}`);
const controls = cp => ({ iterationGeometry: cp.continuation.iterationGeometry,
  stepAcceptance: cp.continuation.stepAcceptance, stagnationLimiter: cp.continuation.stagnationLimiter,
  tolerance: 1e-10 });
const withoutResponsePolicy = checkpoint => {
  const cp = copy(checkpoint); delete cp.continuation.projectionGeometry; return cp;
};
function roundoffDifference(a, b, path = '') {
  if (typeof a === 'number' && typeof b === 'number') {
    const relative = Math.abs(a - b) / Math.max(1, Math.abs(a), Math.abs(b));
    assert.ok(relative < 1e-12, `${path}: ${a} versus ${b} (${relative})`);
    return relative;
  }
  if (a && b && typeof a === 'object' && typeof b === 'object') {
    assert.deepEqual(Object.keys(a), Object.keys(b), path);
    return Math.max(0, ...Object.keys(a).map(k => roundoffDifference(a[k], b[k], `${path}.${k}`)));
  }
  assert.equal(a, b, path); return 0;
}

// Geometry-only chart conversion. This does not march a BL or solve a flow;
// bank coordinates are added while preserving the existing physical fields.
function independentFixture() {
  const cp = copy(retained), f = cp.restart;
  const old = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  const input = { ...f.input, wakeGeometry: 'independent-banks' };
  const target = createStreamtubeBodySystem({ ...input, displacement: old.bl.thicknesses(old.initial.subarray(old.ne)) });
  const packed = transferStreamtubeGeometry(old.euler, old.initial.subarray(0, old.ne), target);
  const decoded = target.decode(packed);
  const system = createCoupledStreamtubeBody(input, { ...f.options, initialBL: f.initialBL,
    initialEuler: { x: packed, nodes: decoded.nodes, undisplacedNodes: decoded.undisplacedNodes } });
  const value = system.evaluate(system.initial);
  const oldValue = old.euler.decode(old.initial.subarray(0, old.ne));
  assert.equal(system.n, 362);
  assert.deepEqual(Array.from(system.initial.slice(0, target.layout.densityCount)),
    Array.from(old.initial.slice(0, old.euler.layout.densityCount)));
  assert.deepEqual(decoded.allocation.groups.map(g => g.map(t => t.massFlow)),
    oldValue.allocation.groups.map(g => g.map(t => t.massFlow)));
  const maximumChange = Math.max(...oldValue.nodes.flatMap((g, gi) => g.flatMap((r, i) => r.map((p, j) =>
    Math.hypot(p.x - value.outer.nodes[gi][i][j].x, p.y - value.outer.nodes[gi][i][j].y)))));
  assert.ok(maximumChange < 1e-14);
  cp.families = value.families;
  cp.restart = { ...f, input, initialEuler: { x: system.initial.slice(0, system.ne), nodes: value.outer.nodes,
    undisplacedNodes: value.outer.undisplacedNodes } };
  assert.deepEqual(cp.restart.initialBL, retained.restart.initialBL);
  return copy(cp);
}

test('the default driver retains the complete legacy trajectory without new checkpoint fields', () => {
  const cp = copy(retained), before = copy(cp), options = { ...controls(cp), resume: cp, maxIterations: 2 };
  const old = previous.solveCoupledStreamtubeIses(undefined, options);
  const current = solveCoupledStreamtubeIses(undefined, options);
  assert.equal(current.history.length, 3, current.reason);
  for (const key of ['checkpoint', 'x', 'residual', 'families', 'history', 'flow', 'boundaryLayer', 'lastRejectedStep'])
    assert.deepEqual(copy(current[key]), copy(old[key]), key);
  assert.equal(current.projectionGeometry, undefined);
  assert.equal(current.checkpoint.continuation.projectionGeometry, undefined);
  assert.deepEqual(cp, before);
});

test('projection geometry controls reject incompatible wake, update and restart policies', () => {
  const cp = copy(retained);
  assert.throws(() => solveCoupledStreamtubeIses(undefined, { ...controls(cp), resume: cp,
    maxIterations: 0, projectionGeometry: 'boundary-increment' }), /projection geometry controls do not match/);
  cp.continuation.projectionGeometry = 'boundary-increment';
  assert.throws(() => solveCoupledStreamtubeIses(undefined, { ...controls(cp), resume: cp,
    maxIterations: 0 }), /requires independent wake banks/);
  cp.continuation.blUpdate = 'giles';
  assert.throws(() => solveCoupledStreamtubeIses(undefined, { ...controls(cp), resume: cp,
    maxIterations: 0 }), /Invalid coupled ISES iteration controls/);
});

test('a real multielement independent-bank update preserves the response policy through serialized restart', t => {
  const cp = independentFixture();
  cp.continuation.projectionGeometry = 'boundary-increment';
  const before = copy(cp), options = { ...controls(cp), maxIterations: 2 };
  const combined = solveCoupledStreamtubeIses(undefined, { ...options, resume: cp });
  assert.equal(combined.history.length, 3, combined.reason);
  const first = solveCoupledStreamtubeIses(undefined, { ...options, resume: cp, maxIterations: 1 });
  assert.equal(first.history.length, 2, first.reason);
  const checkpoint = copy(first.checkpoint), checkpointBefore = copy(checkpoint);
  assert.equal(checkpoint.continuation.projectionGeometry, 'boundary-increment');
  const zero = solveCoupledStreamtubeIses(undefined, { ...options, resume: checkpoint, maxIterations: 0 });
  assert.deepEqual(copy(zero.checkpoint), checkpoint);
  for (const key of ['x', 'residual', 'families', 'boundaryLayer'])
    assert.deepEqual(copy(zero[key]), copy(first[key]), `zero-update ${key}`);
  const { geometryChart: zeroChart, ...zeroFlow } = copy(zero.flow);
  const { geometryChart: firstChart, ...firstFlow } = copy(first.flow);
  assert.deepEqual(zeroFlow, firstFlow, 'zero-update physical flow');
  const maximumChartDifference = roundoffDifference(zeroChart, firstChart, 'zero-update chart');
  const second = solveCoupledStreamtubeIses(undefined, { ...options, resume: checkpoint, maxIterations: 1 });
  assert.equal(second.history.length, 2, second.reason);
  // Independent-bank chart reconstruction can change Jacobian directions at
  // roundoff while preserving the exact packed state and residual at restart.
  // Bind BOTH trajectories to the previous driver: this must not be a new
  // response-policy difference hidden by a numerical comparison tolerance.
  assert.ok(combined.history.slice(1).every(h => !h.maintenance.projectionGeometry.active));
  const oldCombined = previous.solveCoupledStreamtubeIses(undefined, { ...options, resume: withoutResponsePolicy(cp) });
  const oldSecond = previous.solveCoupledStreamtubeIses(undefined, { ...options,
    resume: withoutResponsePolicy(checkpoint), maxIterations: 1 });
  assert.deepEqual(withoutResponsePolicy(combined.checkpoint), copy(oldCombined.checkpoint));
  assert.deepEqual(withoutResponsePolicy(second.checkpoint), copy(oldSecond.checkpoint));
  let maximumRestartDifference = 0;
  for (const key of ['checkpoint', 'x', 'residual', 'families', 'flow', 'boundaryLayer'])
    maximumRestartDifference = Math.max(maximumRestartDifference,
      roundoffDifference(copy(second[key]), copy(combined[key]), key));
  const historyEnd = h => { const { iteration, ...rest } = h.at(-1); return rest; };
  roundoffDifference(copy(historyEnd(second.history)), copy(historyEnd(combined.history)), 'last history');
  assert.equal(second.boundaryLayer.surfaces.length, 4); assert.equal(second.boundaryLayer.wakes.length, 2);
  assert.equal(second.projectionGeometry, 'boundary-increment');
  assert.throws(() => solveCoupledStreamtubeIses(undefined, { ...controls(cp), resume: checkpoint,
    maxIterations: 0, projectionGeometry: 'fixed' }), /projection geometry controls do not match/);
  assert.deepEqual(cp, before); assert.deepEqual(checkpoint, checkpointBefore);
  t.diagnostic(JSON.stringify({ unknowns: combined.x.length, acceptedUpdates: combined.history.length - 1,
    maximumRestartDifference, maximumChartDifference, zeroUpdatePhysicalReplayExact: true, legacyTrajectoriesExact: true,
    families: combined.families, responses: combined.history.slice(1).map(h => h.maintenance.projectionGeometry),
    scope: 'Driver/restart integration, not convergence or physical accuracy validation.' }));
});

test('rejection after an active trial response restores the accepted grid, BL and transition state', async t => {
  const cp = independentFixture(); cp.continuation.projectionGeometry = 'boundary-increment';
  const before = copy(cp), options = { ...controls(cp), resume: cp, maxIterations: 1, maxBacktracks: 0 };
  const initial = solveCoupledStreamtubeIses(undefined, { ...options, maxIterations: 0 });
  const originalImport = "import { respondToCoupledProjectionGeometry } from './streamtube-coupled-projection-geometry.js';";
  const source = fs.readFileSync(driverURL, 'utf8');
  assert.equal(source.split(originalImport).length, 2);
  const injected = source.replace(originalImport,
    `import { respondToCoupledProjectionGeometry as actualResponse } from './streamtube-coupled-projection-geometry.js';
     let checkTrial;
     export function installTrialCheck(check) { checkTrial = check; }
     const respondToCoupledProjectionGeometry = (...args) => checkTrial(actualResponse, ...args);`)
    .replace(/from (['"])(\.[^'"]+)\1/g, (_, quote, path) => `from ${quote}${new URL(path, driverURL).href}${quote}`);
  const driver = await import(`data:text/javascript;base64,${Buffer.from(injected).toString('base64')}`);
  let captured, response, calls = 0;
  driver.installTrialCheck((actual, system, x, projection, decoded) => {
    calls++; captured = system;
    // Controlled failure injection, not a claimed Newton trajectory: a raw
    // thickness 1% below the projected thickness exercises the real geometry
    // response before rejecting the trial. The accepted state is untouched.
    const station = system.bl.stations.find(s => s.kind === 'surface' && s.side === 'upper'
      && s.i > system.euler.layout.bodies[s.body].leadingIndex + 1);
    const deltaStar = x[system.ne + 4 * station.id + 2];
    response = actual(system, x, { ...projection, displacementChanges: [{ ...station,
      beforeDeltaStar: .99 * deltaStar, deltaStar }] }, decoded);
    assert.ok(response.diagnostics.active);
    throw new Error('Controlled rejection after active projection response');
  });
  const result = driver.solveCoupledStreamtubeIses(undefined, options);
  assert.equal(calls, 1); assert.equal(result.linearDiagnostics.solves, 1);
  assert.equal(result.lastRejectedStep.stage, 'projection grid response');
  assert.match(result.reason, /Controlled rejection after active projection response/);
  assert.equal(result.history.length, 1);
  for (const key of ['x', 'residual', 'families', 'flow', 'boundaryLayer', 'checkpoint'])
    assert.deepEqual(copy(result[key]), copy(initial[key]), `retained ${key}`);
  assert.deepEqual(copy(captured.bl.snapshotActive()), cp.restart.options.transitionState);
  assert.deepEqual(copy(captured.euler.decode(captured.initial.slice(0, captured.ne)).nodes), copy(initial.flow.nodes));
  assert.deepEqual(copy(captured.evaluate(captured.initial).residual), copy(initial.residual));
  assert.deepEqual(cp, before);
  t.diagnostic(JSON.stringify({ scope: 'Controlled failure after a manufactured active response, not a flow result.',
    response: response.diagnostics, acceptedStateRestoredExactly: true }));
});
