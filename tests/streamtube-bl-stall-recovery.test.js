import test from 'node:test';
import assert from 'node:assert/strict';
import { retainBestCoupledCheckpoint, transitionRecoveryPlan } from '../src/euler/streamtube-transition-recovery.js';

// Scalar planning fixture only: these placeholder vectors are not a flow seed.
function fixture() {
  const input = { mach: .2, alpha: 4, hybrid: { ismom: 4 }, bodies: [{ leadingIndex: 5, trailingIndex: 25 }],
    outerLower: Array(31).fill(null), weights: [Array(7).fill(1), Array(7).fill(1)] };
  const options = { reynolds: 1e6, ncrit: 9, transitionMode: 'automatic', tripFractions: [[1, 1]], transitionState: [10, 10] };
  const checkpoint = { version: 1, families: { euler: .0001, boundaryLayer: .02, edgeMatching: .001 },
    restart: { input, options, initialEuler: { x: [1], nodes: [], undisplacedNodes: [] }, initialBL: [1, 2, 3, 4] },
    continuation: { linearOrdering: 'station-auto' } };
  const surfaces = ['upper', 'lower'].map((side, b) => ({ body: 0, side, transition: 10,
    ids: Array.from({ length: 20 }, (_, k) => 20 * b + k) }));
  const stations = surfaces.flatMap(s => s.ids.map((id, k) => ({ id, body: 0, side: s.side, i: 6 + k })));
  const history = Array.from({ length: 21 }, (_, iteration) => ({ iteration, step: iteration ? .02 : 0,
    viscousLimiter: { kind: 'xfoil-bl-update', station: 31, variable: 'delta-star', normalizedIncrement: -25, bound: -.5 } }));
  const result = { converged: false, reason: 'iteration limit', conditions: { transitionMode: 'automatic' },
    checkpoint: structuredClone(checkpoint), history, mesh: { quality: { valid: true } },
    boundaryLayer: { stations, surfaces,
      transitions: surfaces.map(s => ({ body: s.body, side: s.side, kind: 'natural', forced: false })) } };
  result.checkpoint.families = { euler: .01, boundaryLayer: .2, edgeMatching: .02 };
  return { result, bestState: { iteration: 10, checkpoint } };
}

test('BL-limited stagnation selects local natural-transition support without phase reversals', () => {
  const { result, bestState } = fixture(), before = structuredClone({ result, bestState });
  const plan = transitionRecoveryPlan(result, { bestState });
  assert.equal(plan.reason, 'boundary-layer-transition-stall');
  assert.equal(plan.source, 'best-accepted-checkpoint'); assert.equal(plan.bestIteration, 10);
  assert.equal(plan.selection.acceptanceChanged, false);
  assert.deepEqual(plan.surfaces.map(s => [s.body, s.side]), [[0, 'lower']]);
  assert.deepEqual(plan.surfaces[0].stations, [24, 25, 26, 27, 28, 29, 30, 31]);
  assert.deepEqual(plan.surfaces[0].intervals, [9, 10, 11, 12, 13, 14, 15, 16, 17]);
  assert.equal(plan.refinedNx, 39); assert.equal(plan.normalFactor, 1); assert.equal(plan.nodeCount, 40 * 16);
  assert.deepEqual({ result, bestState }, before);
});

test('best-state retention requires a complete accepted checkpoint and strictly better complete residual', () => {
  const { bestState } = fixture(), cp = bestState.checkpoint;
  assert.equal(retainBestCoupledCheckpoint(undefined, cp, { iteration: 10 }), undefined);
  const best = retainBestCoupledCheckpoint(undefined, cp, { iteration: 10, initialRedistribution: { accepted: true } });
  assert.deepEqual(best, bestState); assert.notEqual(best.checkpoint, cp);
  assert.equal(retainBestCoupledCheckpoint(best, cp, { iteration: 11, initialRedistribution: { accepted: true } }), best);
  for (const mutate of [c => { c.version = 2; }, c => { c.restart.initialBL[0] = NaN; },
    c => { c.restart.initialEuler.x = []; }, c => { delete c.continuation; },
    c => { c.restart.options.transitionState = [1.5]; }, c => { c.families.euler = -1; },
    c => { c.restart.options.transitionMode = 'fixed-trip'; }]) {
    const bad = structuredClone(cp); mutate(bad);
    assert.equal(retainBestCoupledCheckpoint(undefined, bad, { iteration: 10, initialRedistribution: { accepted: true } }), undefined);
  }
  cp.restart.initialBL[0] = 500; assert.equal(best.checkpoint.restart.initialBL[0], 1);
});

test('improvement, insufficient history, non-BL residuals, forced transition and distant limiters do not refine', () => {
  for (const change of [
    (r, b) => { b.iteration = 13; }, (r, b) => { r.checkpoint.families.boundaryLayer = .03; },
    (r, b) => { b.checkpoint.families.euler = .015; }, r => { r.checkpoint.families.euler = .15; },
    r => { r.converged = true; }, r => { r.mesh.quality.valid = false; }, r => { r.conditions.transitionMode = 'fixed-trip'; },
    r => { r.automaticRefinement = {}; }, r => { r.reason = 'linear solve failed'; },
    r => { r.boundaryLayer.transitions[1].forced = true; }, r => { delete r.boundaryLayer.transitions; },
    r => { r.history.slice(-8).forEach(h => { h.viscousLimiter.station = 39; }); },
    r => { r.history.slice(-5).forEach(h => { h.viscousLimiter = { kind: 'full-step' }; }); },
    r => { r.history.at(-1).step = NaN; }, r => { r.history.at(-1).step = -1; },
    r => { r.history.at(-1).iteration -= 1; }, r => { r.history.slice(-8).forEach(h => { h.viscousLimiter.bound = NaN; }); },
  ]) {
    const { result, bestState } = fixture(); change(result, bestState);
    assert.equal(transitionRecoveryPlan(result, { bestState }), null, change.toString());
  }
});

test('source physics, selected ISMOM, grid and kernel policy must be the same as the retained trajectory', () => {
  for (const change of [f => { f.input.mach = .3; }, f => { f.input.hybrid.ismom = 3; },
    f => { f.options.ncrit = 4; }, f => { f.options.reynolds = 2e6; },
    f => { f.options.hkFloorLinearization = 'native'; }, f => { f.options.tripFractions = [[.05, .05]]; },
    f => { f.input.bodies[0].leadingIndex++; }]) {
    const { result, bestState } = fixture(); change(result.checkpoint.restart);
    assert.equal(transitionRecoveryPlan(result, { bestState }), null);
  }
});

test('node budget and tolerance bound recovery work without loosening physical acceptance', () => {
  const { result, bestState } = fixture(), plan = transitionRecoveryPlan(result, { bestState });
  assert.equal(transitionRecoveryPlan(result, { bestState, maxNodes: plan.nodeCount }).nodeCount, plan.nodeCount);
  assert.equal(transitionRecoveryPlan(result, { bestState, maxNodes: plan.nodeCount - 1 }), null);
  assert.equal(transitionRecoveryPlan(result, { bestState, tolerance: .01 }), null);
  for (const tolerance of [0, -1, NaN, Infinity]) assert.throws(() => transitionRecoveryPlan(result, { bestState, tolerance }), /tolerance/);
});
