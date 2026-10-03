import test from 'node:test';
import fs from 'node:fs';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import assert from 'node:assert/strict';
import { coupledStartupExtensionPlan } from '../src/euler/streamtube-coupled-startup.js';
const nodes = [[[{ x: 0, y: 0 }, { x: 0, y: 1 }], [{ x: 1, y: 0 }, { x: 1, y: 1 }]]];
const requestedConditions = { mach: .2, reynolds: 1e6, ncrit: 9, transitionMode: 'automatic', tripFractions: [[1, 1]], ismom: 4 };
const controls = { startupAttempt: 1, maxIterations: 40, remainingIterations: 20, requestedConditions, tolerance: 1e-10 };
function result() {
  const history = Array.from({ length: 41 }, (_, iteration) => ({ iteration, euler: 10 / (iteration + 1),
    boundaryLayer: 30 / (iteration + 1), edgeMatching: 20 / (iteration + 1), step: 1, backtracks: 0, rejections: [],
    maintenance: { triggeredBodies: [], passages: [], dekinkRepairs: [] } }));
  const families = Object.fromEntries(['euler', 'boundaryLayer', 'edgeMatching'].map(k => [k, history.at(-1)[k]]));
  const checkpoint = { version: 1, families: { ...families }, restart: {
    input: { mach: .2, bodies: [{ leadingIndex: 2, trailingIndex: 8 }], hybrid: { ismom: 4 }, wakeGeometry: 'independent-banks' },
    options: { reynolds: 1e6, ncrit: 9, transitionMode: 'automatic', tripFractions: [[1, 1]], transitionState: [2, 3] },
    initialEuler: { x: Float64Array.of(1, 2), nodes: structuredClone(nodes), undisplacedNodes: structuredClone(nodes) }, initialBL: Float64Array.of(3, 4, 5, 6) },
    continuation: { fractions: [[0, .5, 1]], lastRedistributedStagnation: [.2], iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible',
      stagnationLimiter: 'listing', linearOrdering: 'station-auto', stationFallback: true, preferredOrdering: 'amd', pivotTolerance: .001 } };
  return { converged: false, reason: 'iteration limit', initialRedistribution: { accepted: true }, families,
    mesh: { quality: { valid: true, invalidCells: [] } }, conditions: { ...requestedConditions }, history, checkpoint,
    x: Float64Array.of(1, 2, 3, 4, 5, 6), residual: Float64Array.of(families.euler, -families.boundaryLayer, families.edgeMatching, 0, 0, 0),
    flow: { nodes: structuredClone(nodes), undisplacedNodes: structuredClone(nodes) },
    linearDiagnostics: { solves: 40, iterations: Array.from({ length: 40 }, (_, i) => ({ iteration: i + 1 })) }, lastRejectedStep: null };
}
test('only work budget changes; complete accepted checkpoint, phase/chart/maintenance and full original history are preserved', () => {
  const r = result(), before = structuredClone(r), plan = coupledStartupExtensionPlan(r, controls);
  assert.ok(plan); assert.equal(plan.additionalIterations, 20); assert.equal(plan.maximumTotalIterations, 60);
  assert.equal(plan.equationsChanged, false); assert.equal(plan.tolerance, 1e-10);
  assert.deepEqual(plan.resume, r.checkpoint); assert.deepEqual(plan.original.history, r.history);
  assert.equal(plan.original.linearDiagnostics.solves, 40); assert.equal(plan.accounting.continuationIterationOffset, 40);
  assert.equal(plan.accounting.continuationInitialFrameIsAnUpdate, false); assert.deepEqual(r, before);
  plan.resume.restart.initialBL[0] = 99; plan.original.history[0].iteration = -1;
  assert.deepEqual(r, before, 'The extension cannot overwrite evidence of the original retained endpoint.');
});
test('one bounded chunk honors smaller remaining budgets and forbids repeat extensions or any missing/zero cap', () => {
  for (const remainingIterations of [1, 7, 20, 100]) assert.equal(coupledStartupExtensionPlan(result(), { ...controls, remainingIterations }).additionalIterations, Math.min(20, remainingIterations));
  for (const patch of [{ startupAttempt: 2 }, { maxIterations: 0 }, { maxIterations: 39 }, { remainingIterations: 0 }, { remainingIterations: -1 },
    { remainingIterations: 1.5 }, { remainingIterations: Infinity }, { extensionCount: 1 }, { extensionCount: null }, { tolerance: 0 },
    { coarseInitialization: {} }]) assert.equal(coupledStartupExtensionPlan(result(), { ...controls, ...patch }), null);
  assert.equal(coupledStartupExtensionPlan(result()), null);
});
test('prescribed trips receive the same progress budget without changing their physical locations', () => {
  for (const tripFractions of [[[.03, .07]], [[.4, 1]], [[1, .2]]]) {
    const r = result();
    r.checkpoint.restart.options.tripFractions = structuredClone(tripFractions);
    const request = { ...requestedConditions, tripFractions };
    r.conditions = { ...request };
    const before = structuredClone(r);
    const plan = coupledStartupExtensionPlan(r, { ...controls, requestedConditions: request });
    assert.ok(plan);
    assert.equal(plan.additionalIterations, 20);
    assert.deepEqual(plan.resume.restart.options.tripFractions, tripFractions);
    assert.deepEqual(r, before);
    // Trips alone cannot buy more work for a stalled state.
    r.history.at(-1).step = 1e-12;
    assert.equal(coupledStartupExtensionPlan(r, { ...controls, requestedConditions: request }), null);
  }
  for (const trip of [0, -.1, 1.1, NaN, Infinity]) {
    const r = result(), tripFractions = [[trip, .07]];
    r.checkpoint.restart.options.tripFractions = tripFractions;
    assert.equal(coupledStartupExtensionPlan(r, { ...controls,
      requestedConditions: { ...requestedConditions, tripFractions } }), null);
  }
});
test('every family must be strictly positive, finite and decreasing in both final transitions', () => {
  for (const name of ['euler', 'boundaryLayer', 'edgeMatching']) for (const row of [38, 39]) {
    const r = result(); r.history[row][name] = r.history[row + 1][name];
    assert.equal(coupledStartupExtensionPlan(r, controls), null);
  }
  for (const value of [0, -1, NaN, Infinity]) for (const name of ['euler', 'boundaryLayer', 'edgeMatching']) {
    const r = result(); r.history[38][name] = value; assert.equal(coupledStartupExtensionPlan(r, controls), null);
  }
});
test('a productive-looking scalar trend does not override step/event/geometry/domain or incomplete-state guards', () => {
  const mutations = [r => r.converged = true, r => r.reason = 'admissibility', r => r.lastRejectedStep = {},
    r => r.initialRedistribution.accepted = false, r => r.mesh.quality.valid = false, r => r.mesh.quality.invalidCells.push(1),
    r => r.history.at(-1).step = .999, r => r.history.at(-1).backtracks = 1, r => r.history.at(-1).rejections.push({}),
    r => r.history.at(-1).activeChange = true, r => r.history.at(-1).changes = [{}],
    r => r.history.at(-1).maintenance.geometryRedistribution = true, r => r.history.at(-1).maintenance.passages.push({}),
    r => r.history.at(-1).maintenance.triggeredBodies.push(0), r => r.history.at(-1).maintenance.dekinkRepairs.push({}),
    r => r.history.at(-1).eventProfile = {}, r => delete r.checkpoint, r => r.checkpoint.restart.options.transitionState = [],
    r => r.checkpoint.continuation.lastRedistributedStagnation = [], r => r.checkpoint.continuation.fractions = [[NaN]],
    r => r.checkpoint.families.euler *= .5, r => r.x[0] += 1, r => r.flow.nodes[0][0][0].x += 1,
    r => r.flow.undisplacedNodes[0][0][0].y += 1, r => r.residual[0] = Infinity, r => r.history[20].iteration = 21,
    r => r.history.at(-1).euler *= .5, r => r.history.at(-1).maintenance = null];
  for (const mutate of mutations) { const r = result(); mutate(r); assert.equal(coupledStartupExtensionPlan(r, controls), null); }
});
test('different operating points, auxiliary startup/refinement routes and experimental profile policies are excluded', () => {
  for (const field of ['automaticRefinement', 'gridSequence', 'refinement', 'recoveryAttempt', 'ncritContinuation', 'machContinuation', 'startupExtension']) {
    const r = result(); r[field] = {}; assert.equal(coupledStartupExtensionPlan(r, controls), null);
  }
  for (const patch of [{ mach: .3 }, { ncrit: 4 }, { reynolds: 2e6 }, { ismom: 3 }, { transitionMode: 'fixed-trip' }, { tripFractions: [[.1, 1]] }])
    assert.equal(coupledStartupExtensionPlan(result(), { ...controls, requestedConditions: { ...requestedConditions, ...patch } }), null);
  const experimental = result(); experimental.checkpoint.continuation.eventProfile = 'xfoil-mrchdu';
  assert.equal(coupledStartupExtensionPlan(experimental, controls), null);
});
test('a final full step qualifies after earlier phase movement, while earlier-history spikes do not rewrite the endpoint evidence', () => {
  const r = result(); r.history[38].activeChange = true; r.history[38].changes = [{ from: 22, to: 21 }];
  r.history[39].activeChange = true; r.history[39].changes = [{ from: 21, to: 20 }];
  r.history[37].maintenance.geometryRedistribution = true; r.history[20].boundaryLayer = 100;
  assert.ok(coupledStartupExtensionPlan(r, controls));
});


test('an explicit small per-stage cap also bounds the extra stage', () => {
  const r = result(); r.history = r.history.slice(0, 6);
  const fields = ['euler', 'boundaryLayer', 'edgeMatching'];
  r.families = Object.fromEntries(fields.map(k => [k, r.history.at(-1)[k]]));
  r.checkpoint.families = { ...r.families };
  r.residual = Float64Array.of(r.families.euler, r.families.boundaryLayer, r.families.edgeMatching, 0, 0, 0);
  const plan = coupledStartupExtensionPlan(r, { ...controls, maxIterations: 5, remainingIterations: 20 });
  assert.ok(plan); assert.equal(plan.additionalIterations, 5); assert.equal(plan.maximumTotalIterations, 10);
});


test('the captured browser endpoint qualifies after an exact zero-update replay and public history reconstruction', () => {
  const base = new URL('../docs/solver-reliability/rae-default-extension/browser-first-start/', import.meta.url);
  const saved = JSON.parse(fs.readFileSync(new URL('first-start-checkpoint-arrays.json', base))), cp = saved.checkpoint ?? saved;
  const before = structuredClone(cp), c = cp.continuation;
  const history = fs.readFileSync(new URL('events.jsonl', base), 'utf8').trim().split('\n').map(JSON.parse)
    .filter(e => e.type === 'iteration' && e.iteration.stage === 'coupled').map(e => e.iteration);
  assert.deepEqual(history.map(h => h.iteration), Array.from({ length: 41 }, (_, i) => i));
  const replay = solveCoupledStreamtubeIses(undefined, { resume: cp, maxIterations: 0, tolerance: 1e-10,
    iterationGeometry: c.iterationGeometry, stepAcceptance: c.stepAcceptance, stagnationLimiter: c.stagnationLimiter });
  assert.equal(replay.history.length, 1); assert.equal(replay.linearDiagnostics.solves, 0);
  assert.deepEqual(replay.families, cp.families); assert.equal(replay.mesh.quality.valid, true);
  replay.history = history;
  // The capture records only the second20-update chunk's factor diagnostics.
  // Do not invent the missing original40-update linear work log.
  delete replay.linearDiagnostics;
  const plan = coupledStartupExtensionPlan(replay, { startupAttempt: 1, maxIterations: 40, remainingIterations: 20, tolerance: 1e-10,
    requestedConditions: { mach: .2, reynolds: cp.restart.options.reynolds, ncrit: 9, transitionMode: 'automatic', tripFractions: [[1, 1]], ismom: 4 } });
  assert.ok(plan); assert.equal(plan.additionalIterations, 20); assert.equal(plan.maximumTotalIterations, 60);
  assert.equal(plan.accounting.continuationIterationOffset, 40);
  assert.deepEqual(Array.from(plan.resume.restart.initialEuler.x), cp.restart.initialEuler.x);
  assert.deepEqual(Array.from(plan.resume.restart.initialBL), cp.restart.initialBL);
  assert.deepEqual(plan.resume.restart.initialEuler.nodes, cp.restart.initialEuler.nodes);
  assert.deepEqual(plan.resume.restart.initialEuler.undisplacedNodes, cp.restart.initialEuler.undisplacedNodes);
  assert.deepEqual(plan.resume.restart.options, cp.restart.options);
  assert.deepEqual(plan.resume.continuation, cp.continuation); assert.deepEqual(cp, before);
});


test('the saved first NLR attempt fails the full-step and all-family progress gates', () => {
  const saved = JSON.parse(fs.readFileSync(new URL('../docs/solver-reliability/gui-defaults-browser-recovery-v1/nlr7301/result.json', import.meta.url)));
  const attempt = saved.initialization.attempts[0], tail = attempt.history.slice(-3);
  assert.equal(attempt.iterations, 40); assert.equal(attempt.converged, false);
  assert.equal(tail.at(-1).step, .0016546361906289201);
  assert.ok(tail[0].euler < tail[1].euler && tail[1].euler < tail[2].euler);
  assert.ok(tail[0].edgeMatching < tail[1].edgeMatching && tail[1].edgeMatching < tail[2].edgeMatching);
  // The original receipt keeps this attempt's history but not its complete
  // checkpoint. A controlled packet isolates the exact archived history;
  // it does not reconstruct or certify the discarded flow state.
  const r = result(); r.history = structuredClone(attempt.history); r.families = { ...attempt.families };
  r.checkpoint.families = { ...r.families };
  r.residual = Float64Array.of(r.families.euler, r.families.boundaryLayer, r.families.edgeMatching, 0, 0, 0);
  assert.equal(coupledStartupExtensionPlan(r, controls), null);
});

test('the saved successful NACA browser root remains ineligible for any extension', () => {
  const saved = JSON.parse(fs.readFileSync(new URL('../docs/solver-reliability/gui-defaults-browser-qualified-v2/single/result.json', import.meta.url)));
  assert.equal(saved.converged, true); assert.equal(saved.initialization.attempts.length, 1);
  assert.equal(saved.history.at(-1).iteration, 29);
  const plan = coupledStartupExtensionPlan(saved, { ...controls, requestedConditions: {
    ...requestedConditions, reynolds: saved.checkpoint.restart.options.reynolds } });
  assert.equal(plan, null);
});

test('full Armijo steps earn a bounded extension when a shock redistributes family maxima', () => {
  const make = () => {
    const r = result();
    r.checkpoint.restart.options.edgeMatching = 'section-velocity';
    // Last BL maximum rises, but the complete squared residual falls strongly.
    r.history[39].boundaryLayer = .6;
    const norm = r.residual.reduce((sum, v) => sum + v * v, 0);
    r.history[39].residualDecrease = { method: 'armijo-squared-residual', step: 1,
      beforeSquaredNorm: 4 * norm, afterSquaredNorm: 2 * norm, allowedSquaredNorm: 3.9996 * norm };
    r.history[40].residualDecrease = { method: 'armijo-squared-residual', step: 1,
      beforeSquaredNorm: 2 * norm, afterSquaredNorm: norm, allowedSquaredNorm: 1.9998 * norm };
    return r;
  };
  const r = make(), before = structuredClone(r), plan = coupledStartupExtensionPlan(r, controls);
  assert.ok(plan); assert.equal(plan.additionalIterations, 20);
  assert.match(plan.selection, /two full Armijo steps/);
  assert.deepEqual(r, before);
  for (const change of [
    r => { r.history[39].step = .5; },
    r => { r.history[39].activeChange = true; },
    r => { r.history[39].backtracks = 1; },
    r => { r.history[39].residualDecrease.afterSquaredNorm *= 1.01; },
    r => { r.history[40].residualDecrease.afterSquaredNorm *= .99; },
    r => { r.history[40].residualDecrease.allowedSquaredNorm = Infinity; },
    r => { r.history[40].residualDecrease.method = 'unverified'; },
    r => { r.checkpoint.restart.options.edgeMatching = 'pressure'; },
    r => { r.history[39].residualDecrease.beforeSquaredNorm = 1.1 * r.history[40].residualDecrease.afterSquaredNorm; },
  ]) {
    const rejected = make(); change(rejected);
    assert.equal(coupledStartupExtensionPlan(rejected, controls), null);
  }
});

function eventfulProgress() {
  const r = result();
  r.checkpoint.continuation.stepAcceptance = 'event-armijo';
  for (const h of r.history.slice(-20)) {
    h.step = .01; h.activeChange = true; h.changes = [{ from: 2, to: 3 }];
    h.residualDecrease = { method: 'armijo-squared-residual', step: .01,
      beforeSquaredNorm: 1e6, afterSquaredNorm: 1e-6, allowedSquaredNorm: 9e5 };
  }
  for (const i of [23, 39]) {
    const h = r.history[i]; delete h.activeChange; h.changes = [];
    h.residualDecrease = { method: 'armijo-squared-residual', step: .01,
      beforeSquaredNorm: i * 100, afterSquaredNorm: i * 70, allowedSquaredNorm: i * 99 };
  }
  return r;
}

test('bounded work can continue through transitions using only paired unchanged-equation descent', () => {
  const r = eventfulProgress(), before = structuredClone(r), plan = coupledStartupExtensionPlan(r, controls);
  assert.ok(plan);
  assert.match(plan.selection, /measured descent/);
  assert.equal(plan.progress.descendingSteps, 2);
  assert.ok(Math.abs(plan.progress.pairedMeritRatio - .49) < 1e-14);
  assert.equal(plan.progress.comparesAcrossEvents, false);
  assert.deepEqual(plan.resume, r.checkpoint);
  assert.deepEqual(r, before);
  // Apparent decreases across N/shear equation changes buy no extra work.
  r.history[23].activeChange = true;
  assert.equal(coupledStartupExtensionPlan(r, controls), null);
});

test('measured progress excludes tiny rays, failed trials, incomplete certificates and changing dissipation', () => {
  for (const mutation of [
    r => r.history[23].accepted = false,
    r => r.history[23].step = 1e-12,
    r => r.history[23].residualDecrease.step = .1,
    r => r.history[23].residualDecrease.afterSquaredNorm = NaN,
    r => r.history[23].residualDecrease.allowedSquaredNorm = Infinity,
    r => r.history[23].residualDecrease.method = 'unverified',
    r => r.checkpoint.continuation.dissipationEnhancement = true,
    r => { for (const i of [23, 39]) r.history[i].residualDecrease.afterSquaredNorm *= 1.4; },
  ]) {
    const r = eventfulProgress(); mutation(r);
    assert.equal(coupledStartupExtensionPlan(r, controls), null);
  }
});

test('each new chunk earns its own finite reserve and preserves a consistent cumulative history', () => {
  const r = eventfulProgress();
  for (let count = 0; count < 4; count++) {
    const plan = coupledStartupExtensionPlan(r, { ...controls, remainingIterations: 80,
      extensionCount: count });
    assert.ok(plan);
    assert.equal(plan.originalIterations, 40 + 20 * count);
    assert.equal(plan.initialIterations, 40);
    assert.equal(plan.additionalIterations, 20);
    assert.equal(plan.maximumTotalIterations, 60 + 20 * count);
    const nextRows = structuredClone(r.history.slice(-20));
    nextRows.forEach(h => { h.iteration += 20; }); r.history.push(...nextRows);
    r.startupExtension = { initialIterations: 40, extensionCount: count + 1, iterations: 20 * (count + 1) };
  }
  assert.equal(r.history.length, 121);
  assert.equal(coupledStartupExtensionPlan(r, { ...controls, remainingIterations: 80, extensionCount: 4 }), null);
  const bad = eventfulProgress(); bad.startupExtension = { initialIterations: 40, extensionCount: 1, iterations: 20 };
  assert.equal(coupledStartupExtensionPlan(bad, { ...controls, extensionCount: 1 }), null);
});

test('a finer surface gets a finite reserve and must keep earning progress through its final partial chunk', () => {
  const r = eventfulProgress();
  r.checkpoint.restart.input.bodies[0] = { leadingIndex: 21, trailingIndex: 147 };
  let count = 0, plan;
  while ((plan = coupledStartupExtensionPlan(r, { ...controls, extensionCount: count }))) {
    assert(plan.additionalIterations > 0 && plan.additionalIterations <= 20);
    assert.equal(plan.maximumExtensionIterations, 252);
    const offset = r.history.length - 1;
    const nextRows = structuredClone(r.history.slice(-20).slice(0, plan.additionalIterations));
    nextRows.forEach((h, i) => { h.iteration = offset + i + 1; });
    r.history.push(...nextRows);
    r.startupExtension = { initialIterations: 40, extensionCount: ++count, iterations: r.history.length - 41 };
    if (count === 4) {
      assert.equal(r.history.length, 121);
      assert.ok(coupledStartupExtensionPlan(r, { ...controls, extensionCount: count }),
        'A productive fine-grid state must survive the former coarse-grid cap.');
      const stalled = structuredClone(r);
      stalled.history.slice(-20).forEach(h => { delete h.residualDecrease; h.step = 1e-12; });
      assert.equal(coupledStartupExtensionPlan(stalled, { ...controls, extensionCount: count }), null,
        'More geometric capacity must not buy work for a stalled solve.');
    }
  }
  assert.equal(r.history.length, 293);
  assert.equal(count, 13);
  assert.equal(r.startupExtension.iterations, 252);
});

test('the coarse NACA at eight degrees finishes from its productive transition-limited budget endpoint', async () => {
  const { gunzipSync } = await import('node:zlib');
  const saved = JSON.parse(gunzipSync(fs.readFileSync(new URL('./fixtures/naca8-progress-extension.json.gz', import.meta.url))));
  const original = structuredClone(saved), cp = saved.checkpoint, options = cp.restart.options;
  let r = solveCoupledStreamtubeIses(undefined, { ...cp.continuation, resume: cp, maxIterations: 0, tolerance: 1e-10 });
  r.history = saved.history;
  const conditions = { mach: cp.restart.input.mach, reynolds: options.reynolds, ncrit: options.ncrit,
    transitionMode: options.transitionMode, tripFractions: options.tripFractions, ismom: cp.restart.input.hybrid.ismom };
  let reserve = 80, count = 0;
  while (!r.converged && reserve > 0) {
    const plan = coupledStartupExtensionPlan(r, { startupAttempt: 1, maxIterations: 40,
      remainingIterations: reserve, extensionCount: count, requestedConditions: conditions });
    assert.ok(plan, `No measured progress after ${r.history.length - 1} updates`);
    const before = r, offset = r.history.length - 1;
    r = solveCoupledStreamtubeIses(undefined, { ...plan.resume.continuation, resume: plan.resume,
      maxIterations: plan.additionalIterations, tolerance: 1e-10 });
    r.history = [...before.history, ...r.history.slice(1).map(h => ({ ...h, iteration: h.iteration + offset }))];
    r.startupExtension = { initialIterations: 40, extensionCount: ++count, iterations: r.history.length - 41 };
    reserve -= plan.additionalIterations;
  }
  assert.equal(r.converged, true, r.reason);
  assert.ok(Math.max(...Object.values(r.families)) <= 1e-10);
  assert.equal(r.mesh.quality.valid, true);
  assert.ok(r.history.length <= 121);
  assert.equal(r.conditions.mach, .2);
  assert.equal(r.conditions.reynolds, options.reynolds); // Kernel Reynolds uses the solver's reference length.
  assert.equal(r.conditions.ncrit, 9);
  assert.deepEqual(r.solverInput.upwind, cp.restart.input.upwind);
  assert.deepEqual(saved, original);
});
