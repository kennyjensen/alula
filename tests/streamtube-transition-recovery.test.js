import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { transitionRecoveryPlan } from '../src/euler/streamtube-transition-recovery.js';

const fixture = JSON.parse(fs.readFileSync(new URL('./fixtures/nlr-transition-cycle.json', import.meta.url)));
// Only the detector's scalar inputs are reproduced. This mock cannot be used
// as a numerical restart and never creates a mesh, kernel or flow system.
function result() {
  const data = structuredClone(fixture), stations = [];
  for (const station of data.currentStations) stations[station.id] = { i: station.i };
  return { converged: false, reason: 'iteration limit', conditions: data.conditions,
    mesh: { quality: { valid: true } }, history: data.history,
    checkpoint: { version: 1, restart: { input: { ...data.physicalInput, bodies: data.bodies,
      outerLower: Array(data.nx + 1).fill(null), weights: data.tubeCounts.map(n => Array(n).fill(1)) } } },
    boundaryLayer: { stations, surfaces: data.surfaces } };
}
const refinedIntervals = plan => plan.streamwiseSubdivisions.flatMap((count, i) => count === 4 ? [i] : []);
const natural = (from, to, extra = {}) => ({ body: 1, side: 'lower', from, to, kind: 'natural', ...extra });
const alternating = () => [natural(20, 19), natural(19, 20), natural(20, 19), natural(19, 20)];
function setEvents(r, events) {
  r.history = events.map((change, k) => ({ iteration: 2 * k + 1, changes: [change] }));
}

test('archived sparse accepted trajectory selects the observed cycling surface and exact tested interval window', () => {
  const r = result(), before = structuredClone(r), plan = transitionRecoveryPlan(r);
  assert.equal(plan.reason, 'repeated-adjacent-natural-transition');
  assert.equal(plan.surfaces.length, 1);
  assert.deepEqual([plan.surfaces[0].body, plan.surfaces[0].side], [1, 'lower']);
  assert.deepEqual(plan.surfaces[0].events.map(e => e.iteration), [14, 16, 18, 20]);
  assert.equal(plan.surfaces[0].transitionStation, 53);
  assert.deepEqual(refinedIntervals(plan), [49, 50, 51, 52, 53]);
  assert.equal(plan.parentNx, 130);
  assert.equal(plan.refinedNx, 145);
  assert.equal(plan.normalFactor, 1);
  assert.equal(plan.nodeCount, 146 * (11 + 14 + 11));
  assert.deepEqual(r, before, 'Planning must leave all geometry and operating inputs unchanged.');
});

test('four connected accepted alternations are required; progress on other surfaces does not manufacture a cycle', () => {
  for (const events of [alternating().slice(0, 3),
    [natural(18, 19), natural(19, 20), natural(20, 21), natural(21, 20)],
    [natural(20, 19), natural(19, 20), natural(22, 21), natural(19, 20)]]) {
    const r = result(); setEvents(r, events);
    assert.equal(transitionRecoveryPlan(r), null);
  }
  const r = result(); setEvents(r, alternating());
  r.history.splice(2, 0, { iteration: 4, changes: [{ body: 0, side: 'upper', from: 3, to: 4, kind: 'natural' }] });
  assert.equal(transitionRecoveryPlan(r).surfaces.length, 1);
});

test('forced, nonadjacent and leading/terminal events reset the accepted natural cycle', () => {
  for (const changed of [natural(20, 19, { kind: 'forced' }),
    natural(20, 19, { kind: 'trailing-edge' }), natural(20, 19, { kind: 'leading' }), natural(22, 19)]) {
    const r = result(), events = alternating(); events[2] = changed; setEvents(r, events);
    assert.equal(transitionRecoveryPlan(r), null);
  }
});

test('auxiliary reconciliation and rejected proposals are not counted as natural changes', () => {
  const r = result(); setEvents(r, alternating());
  r.history.splice(2, 0, { iteration: 4, changes: [natural(20, 20,
    { kind: 'amplification-reconciliation', auxiliaryOnly: true })],
    rejections: [{ changes: [natural(10, 11)] }] });
  assert.equal(transitionRecoveryPlan(r).surfaces.length, 1);
  setEvents(r, [natural(20, 20, { kind: 'amplification-reconciliation', auxiliaryOnly: true })]);
  r.history[0].rejections = alternating().map(change => ({ changes: [change] }));
  assert.equal(transitionRecoveryPlan(r), null);
});

test('converged, invalid-grid, fixed-transition, incomplete and stale histories are ineligible', () => {
  for (const change of [r => { r.converged = true; }, r => { r.mesh.quality.valid = false; },
    r => { r.conditions.transitionMode = 'fixed-trip'; }, r => { delete r.checkpoint; },
    r => { delete r.checkpoint.restart; }, r => { r.history.push({ iteration: 25, changes: [] }); }]) {
    const r = result(); change(r); assert.equal(transitionRecoveryPlan(r), null);
  }
});

test('current checkpoint active interval must be the endpoint of the accepted cycle', () => {
  const r = result(); r.boundaryLayer.surfaces[3].transition = 19;
  r.boundaryLayer.stations[r.boundaryLayer.surfaces[3].ids[19]] = { i: 52 };
  assert.equal(transitionRecoveryPlan(r), null);
});

test('a prior local-transition recovery cannot schedule another attempt', () => {
  const r = result(); r.automaticRefinement = { kind: 'transition-local' };
  assert.equal(transitionRecoveryPlan(r), null);
});

test('multiple cycling surfaces share interval subdivisions instead of multiplying them', () => {
  const r = result(); setEvents(r, alternating());
  r.boundaryLayer.surfaces[2].transition = 20;
  r.boundaryLayer.stations[r.boundaryLayer.surfaces[2].ids[20]] = { i: 53 };
  for (const h of r.history) h.changes.push({ ...h.changes[0], side: 'upper' });
  const plan = transitionRecoveryPlan(r);
  assert.deepEqual(new Set(plan.surfaces.map(s => `${s.body}/${s.side}`)), new Set(['1/upper', '1/lower']));
  assert.deepEqual(refinedIntervals(plan), [49, 50, 51, 52, 53]);
  assert.equal(plan.refinedNx - plan.parentNx, 15);
  assert.equal(plan.normalFactor, 1);
});

test('exact node budget is honored and invalid budgets are rejected', () => {
  const r = result(), exact = transitionRecoveryPlan(r).nodeCount;
  assert.equal(transitionRecoveryPlan(r, { maxNodes: exact }).nodeCount, exact);
  assert.equal(transitionRecoveryPlan(r, { maxNodes: exact - 1 }), null);
  for (const maxNodes of [0, -1, 2.5, 50001, Infinity, NaN])
    assert.throws(() => transitionRecoveryPlan(r, { maxNodes }), /node budget/);
});

test('a surface-edge window never splits inlet or wake intervals', () => {
  const r = result(), body = r.checkpoint.restart.input.bodies[1], branch = r.boundaryLayer.surfaces[3];
  for (const i of [body.leadingIndex + 1, body.trailingIndex]) {
    r.boundaryLayer.stations[branch.ids[branch.transition]].i = i;
    const plan = transitionRecoveryPlan(r), intervals = refinedIntervals(plan);
    assert.ok(intervals.length > 0 && intervals.length <= 5);
    assert.ok(intervals.every(k => k >= body.leadingIndex && k < body.trailingIndex));
  }
});
