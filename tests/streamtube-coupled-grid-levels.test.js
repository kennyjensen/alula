import { coupledConvergenceSatisfied } from '../src/euler/streamtube-coupled-convergence.js';
// SPDX-License-Identifier: GPL-2.0-or-later
// Pure dimension/mass plans are real. Controller tests replace only numerical
// constructors/preparation with declared stubs; no global PDE solve is run.
import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import * as shearPolicy from '../src/euler/streamtube-coupled-shear-policy.js';
import { gridReplayDeparture } from '../src/euler/streamtube-coupled-grid-levels.js';
import { planCoupledGridLevels } from '../src/euler/streamtube-coupled-grid-levels.js';
import { planCoupledGridSequence, COUPLED_GRID_SEQUENCE_MAX_NODES } from '../src/euler/streamtube-coupled-grid-sequence.js';

const clone = structuredClone;
const points = (nx, y) => Array.from({ length: nx + 1 }, (_, i) => ({ x: i / nx, y }));
const input = (tubes = [3, 3], nx = 10, finite = false) => ({
  mach: .185, alpha: 6, streamwiseMode: 'hybrid', hybrid: { ismom: 4, epsilonP: 1e-5 },
  upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } },
  wakeGeometry: 'independent-banks', wakeOutlet: 'banks',
  bodies: tubes.slice(1).map((_, element) => ({ element, leadingIndex: 2, trailingIndex: 7,
    points: [{ x: 1, y: 0 }, { x: 0, y: 0 }, { x: 1, y: finite ? -.001 : 0 }],
    ...(finite ? { trailingEdge: { kind: 'finite-base', upperIndex: 0, lowerIndex: 2 } } : {}) })),
  weights: tubes.map(n => Array.from({ length: n }, (_, j) => j + 1)),
  outerLower: points(nx, -1), outerUpper: points(nx, 1), cutPaths: [points(nx, 0)], gridSpacing: { supplied: true },
});
const options = input => ({ reynolds: 2.51e6, ncrit: 4, edgeMatching: 'section-velocity',
  transitionMode: 'automatic', blThermodynamics: 'historical-common-isentrope',
  tripFractions: input.bodies.map(() => [1, 1]), transitionState: input.bodies.flatMap(() => [1, 1]) });
const args = (requestedGridIntervals = 128, source = input(), target = input([5, 5], 70)) => ({
  sourceInput: source, input: target, sourceGridIntervals: 16, requestedGridIntervals,
});

test('full-path plans retain moderate factors, mass partitions and exact node counts before numerical work', () => {
  const a = args(), before = clone(a), plan = planCoupledGridLevels(a);
  assert.deepEqual(plan.levels.map(p => p.nominalGridIntervals), [32, 64, 128]);
  assert.deepEqual(plan.levels.map(p => p.streamwiseFactor), [2, 2, 2]);
  assert.deepEqual(plan.levels.map(p => p.child.nx), [20, 40, 80]);
  assert.deepEqual(plan.levels.map(p => p.nodeCount), [252, 492, 972]);
  assert.equal(plan.globalSolvesBeforeFinalSeed, 2); assert.equal(plan.finalSeedOnly, true);
  assert.deepEqual(plan.levels.map(p => p.final), [false, false, true]);
  let weights = clone(a.sourceInput.weights);
  for (const level of plan.levels) {
    const child = weights.map((row, g) => row.flatMap((w, j) =>
      Array(level.normalSubdivisions[g][j]).fill(w / level.normalSubdivisions[g][j])));
    weights.forEach((row, g) => assert.ok(Math.abs(row.reduce((a, b) => a + b) - child[g].reduce((a, b) => a + b)) < 1e-13));
    weights = child;
  }
  assert.deepEqual(a, before);
  const large = planCoupledGridLevels(a, { maxStreamwiseFactor: 4 });
  assert.deepEqual(large.levels.map(p => p.nominalGridIntervals), [64, 128]);
  assert.deepEqual(large.levels.map(p => p.streamwiseFactor), [4, 2]);
});

test('the existing 32-level plan is byte-for-byte the same single nested subdivision', () => {
  for (const finite of [false, true]) {
    const a = args(32, input([10, 13, 10], 130, finite), input([12, 15, 12], 231, finite));
    const one = planCoupledGridSequence(a), many = planCoupledGridLevels(a);
    const { index, final, targetWeights, ...same } = many.levels[0];
    assert.deepEqual(same, one); assert.deepEqual(targetWeights, a.input.weights);
    assert.equal(index, 0); assert.equal(final, true); assert.equal(many.globalSolvesBeforeFinalSeed, 0);
    assert.equal(one.child.nx * one.child.tubes.reduce((a, b) => a + b), 10140);
  }
});

test('normal resolution can stage independently while never merging a captured parent tube', () => {
  const a = args(16, input([2, 3], 10), input([20, 2], 20)), before = clone(a);
  const plan = planCoupledGridLevels(a);
  assert.deepEqual(plan.levels.map(p => p.child.tubes), [[8, 3], [20, 3]]);
  assert.deepEqual(plan.levels.map(p => p.streamwiseFactor), [1, 1]);
  assert.deepEqual(plan.levels.map(p => p.nominalGridIntervals), [16, 16]);
  for (const level of plan.levels) level.targetWeights.forEach((row, g) =>
    assert.equal(row.reduce((a, b) => a + b), a.input.weights[g].reduce((a, b) => a + b)));
  assert.deepEqual(a, before);
});

test('whole-path oversized output, malformed controls and no-op requests fail during pure planning', () => {
  assert.throws(() => planCoupledGridLevels(args(), { maxNodes: 500 }), /node budget/);
  for (const bound of [0, 1, 5, 2.5, NaN]) assert.throws(() => planCoupledGridLevels(args(), { maxStreamwiseFactor: bound }), /subdivision bound/);
  for (const bound of [0, 50001, Infinity]) assert.throws(() => planCoupledGridLevels(args(), { maxNodes: bound }), /node budget/);
  for (const counts of [[16, 8], [0, 128], [16, 129]])
    assert.throws(() => planCoupledGridLevels({ ...args(), sourceGridIntervals: counts[0], requestedGridIntervals: counts[1] }), /ordered source/);
  assert.throws(() => planCoupledGridLevels(args(16, input(), input())), /does not refine/);
  assert.throws(() => planCoupledGridLevels({ ...args(), input: input([0, 2]) }), /complete matching/);
  const rounded = planCoupledGridLevels(args(47));
  assert.equal(rounded.finalNominalGridIntervals, 64);
  assert.equal(rounded.requestedGridIntervals, 47);
});

const families = { euler: 1e-12, boundaryLayer: 2e-12, edgeMatching: 3e-12 };
const residual = () => Float64Array.from([1e-12, 2e-12, 3e-12, 0, 0]);
const nodes = nx => [[[{ x: 0, y: 0 }, { x: 0, y: 1 }], [{ x: nx, y: 0 }, { x: nx, y: 1 }]]];
function root(data, controls = options(data)) {
  const nx = data.outerLower.length - 1, initialEuler = { x: [nx], nodes: nodes(nx), undisplacedNodes: nodes(nx) };
  const initialBL = [1, .001, .002, 1], x = [...initialEuler.x, ...initialBL];
  return { converged: true, reason: 'residual', x, residual: Array.from(residual()), families: clone(families),
    history: [{ iteration: 0 }, { iteration: 1 }], boundaryLayer: { transitionState: clone(controls.transitionState) },
    mesh: { quality: { valid: true }, initialization: {} },
    checkpoint: { version: 1, families: clone(families), restart: { input: clone(data), options: clone(controls), initialEuler, initialBL } } };
}
let serial = 0;
async function harness(config = {}) {
  const calls = { constructs: 0, evaluations: 0, preparations: [], solves: [] }, key = `__gridLevelsDraft${serial++}`;
  const create = (data, controls) => {
    calls.constructs++;
    const initial = Float64Array.from([...controls.initialEuler.x, ...controls.initialBL]);
    return { initial, conditions: { ncrit: controls.ncrit, reynolds: controls.reynolds },
      euler: { layout: { nx: data.outerLower.length - 1 }, conditions: { mach: data.mach, hybrid: clone(data.hybrid) } },
      bl: { transitionMode: controls.transitionMode, trips: clone(controls.tripFractions), snapshotActive: () => clone(controls.transitionState) },
      admissible: () => !config.inadmissible,
      evaluate() { calls.evaluations++; return { families: clone(families), residual: residual(),
        outer: { nodes: controls.initialEuler.nodes.map(group => group.map(row => row.map(p =>
          ({ ...p, x: p.x + (config.nodeDeparture ?? 0) })))), undisplacedNodes: controls.initialEuler.undisplacedNodes } }; } };
  };
  globalThis[key] = { ...shearPolicy, coupledConvergenceSatisfied, gridReplayDeparture, planCoupledGridSequence, COUPLED_GRID_SEQUENCE_MAX_NODES,
    createCoupledStreamtubeBody: create, streamtubeMeshSnapshot: () => ({ quality: { valid: !config.invalidMesh } }),
    prepareCoupledGridSequence(a, controls) {
      const source = a.sourceResult.checkpoint.restart.input;
      calls.preparations.push({ sourceNx: source.outerLower.length - 1, target: a.requestedGridIntervals,
        sourceGridIntervals: a.sourceGridIntervals, options: clone(a.options), controls: clone(controls) });
      if (calls.preparations.length === config.rejectPreparation) throw new Error('Controlled child grid rejection.');
      const plan = planCoupledGridSequence({ ...a, sourceInput: source }, controls);
      const data = { ...clone(source), outerLower: points(plan.child.nx, -1), outerUpper: points(plan.child.nx, 1),
        cutPaths: [points(plan.child.nx, 0)], weights: source.weights.map((row, g) => row.flatMap((w, j) =>
          Array(plan.normalSubdivisions[g][j]).fill(w / plan.normalSubdivisions[g][j]))) };
      const value = root(data, { ...clone(a.options), transitionState: a.sourceResult.checkpoint.restart.options.transitionState });
      const f = value.checkpoint.restart, system = create(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
      return { input: data, options: clone(f.options),
        iterationControls: { shearCoordinate: shearPolicy.coupledCheckpointShearCoordinate(a.sourceResult.checkpoint) }, system, initialEuler: f.initialEuler, initialBL: f.initialBL,
        mesh: value.mesh, value: system.evaluate(system.initial), transfer: { ...plan, initialGuessOnly: true }, initialization: { thicknessFactor: 1 } };
    } };
  const draft = fs.readFileSync(new URL('../src/euler/streamtube-coupled-grid-levels.js', import.meta.url), 'utf8')
    .replace(/^import \{([^}]+)\} from '[^']+';/gm, (_, names) => `const {${names}} = globalThis[${JSON.stringify(key)}];`);
  const module = await import('data:text/javascript;base64,' + Buffer.from(draft).toString('base64'));
  delete globalThis[key];
  const a = args(), sourceResult = root(a.sourceInput), request = { sourceResult, input: a.input,
    options: clone(sourceResult.checkpoint.restart.options), sourceGridIntervals: 16, requestedGridIntervals: 128 };
  const solve = (prepared, level) => {
    calls.solves.push(clone(level));
    const candidate = root(prepared.input, prepared.options);
    if (prepared.iterationControls.shearCoordinate === 'logarithmic')
      candidate.checkpoint.continuation = { blUpdate: 'xfoil', shearCoordinate: 'logarithmic' };
    config.changeResult?.(candidate, level);
    return candidate;
  };
  return { ...module, calls, request, solve, create };
}

test('intermediate-only execution returns the final seed and verified64 root without solving128', async () => {
  const h = await harness(), before = clone(h.request), events = [];
  const result = h.prepareCoupledGridLevels(h.request, { solveIntermediateLevel: h.solve, onLevel: e => events.push(e) });
  assert.deepEqual(h.calls.solves.map(s => s.nominalGridIntervals), [32, 64]);
  assert.deepEqual(h.calls.preparations.map(s => s.target), [32, 64, 128]);
  assert.equal(result.sourceGridIntervals, 64); assert.equal(result.sourceResult.checkpoint.restart.input.outerLower.length - 1, 40);
  assert.equal(result.prepared.transfer.nominalGridIntervals, 128); assert.equal(result.prepared.input.outerLower.length - 1, 80);
  assert.equal(result.diagnostics.finalSeedPrepared, true); assert.equal(result.diagnostics.reachedTarget, false);
  assert.equal(result.diagnostics.stateConverged, true); assert.equal(result.diagnostics.actualGridIntervals, 64);
  assert.deepEqual(result.diagnostics.levels.map(l => l.converged), [true, true, false]);
  assert.equal(result.sourceSystem.initial[0], 40); assert.deepEqual(h.request, before);
  assert.deepEqual(events.map(e => e.stage), ['preparing', 'prepared', 'converged', 'preparing', 'prepared', 'converged', 'preparing', 'prepared']);
});

test('32-level execution reuses the supplied source system and calls the unchanged preparer once', async () => {
  const h = await harness(), request = { ...h.request, requestedGridIntervals: 32 };
  const f = request.sourceResult.checkpoint.restart;
  request.sourceSystem = h.create(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  const result = h.prepareCoupledGridLevels(request);
  assert.equal(h.calls.preparations.length, 1); assert.equal(h.calls.solves.length, 0);
  assert.equal(result.sourceResult, request.sourceResult); assert.equal(result.sourceSystem, request.sourceSystem);
  assert.equal(result.sourceGridIntervals, 16); assert.equal(result.prepared.transfer.nominalGridIntervals, 32);
  assert.equal(result.diagnostics.finalSeedPrepared, true); assert.equal(result.diagnostics.reachedTarget, false);
});

test('an oversized final plan fails before source construction or any intermediate solve', async () => {
  const h = await harness();
  assert.throws(() => h.prepareCoupledGridLevels(h.request, { maxNodes: 500, solveIntermediateLevel: h.solve }), /node budget/);
  assert.deepEqual(h.calls, { constructs: 0, evaluations: 0, preparations: [], solves: [] });
  assert.throws(() => h.prepareCoupledGridLevels(h.request), /intermediate coupled solve callback/);
  assert.equal(h.calls.constructs, 0);
});

test('mapping and nonconverged intermediate failures retain the last verified source and its resolution', async () => {
  for (const config of [{ rejectPreparation: 2 }, { rejectPreparation: 3 },
    { changeResult: (result, level) => { if (level.index === 1) { result.converged = false; result.reason = 'iteration limit'; } } }]) {
    const h = await harness(config), result = h.prepareCoupledGridLevels(h.request, { solveIntermediateLevel: h.solve });
    assert.equal(result.prepared, null); assert.equal(result.diagnostics.finalSeedPrepared, false);
    assert.equal(result.diagnostics.reachedTarget, false); assert.equal(result.sourceResult.converged, true);
    const retained = config.rejectPreparation === 3 ? 64 : 32;
    assert.equal(result.sourceGridIntervals, retained); assert.equal(result.diagnostics.actualGridIntervals, retained);
    assert.equal(result.sourceSystem.initial[0], retained === 64 ? 40 : 20);
    assert.equal(result.diagnostics.failure.level, config.rejectPreparation === 3 ? 2 : 1);
  }
});

test('changed physical controls, stale residual/phase and changed topology cannot become retained roots', async () => {
  const changes = [
    r => { r.checkpoint.restart.input.mach = .2; }, r => { r.checkpoint.restart.input.hybrid.ismom = 3; },
    r => { r.checkpoint.restart.input.bodies[0].points[1].x = .01; }, r => { r.checkpoint.restart.options.ncrit = 5; },
    r => { r.checkpoint.restart.options.reynolds = 2.7e6; }, r => { r.checkpoint.restart.options.tripFractions[0][0] = .5; },
    r => { r.checkpoint.restart.options.hkFloorLinearization = 'native'; }, r => { r.residual[0] = 0; },
    r => { r.boundaryLayer.transitionState[0] = 2; }, r => { r.x[0] += 1; },
    r => { r.checkpoint.restart.input.weights[0][0] *= 2; }, r => { r.checkpoint.restart.input.outerLower.push({ x: 1, y: -1 }); },
  ];
  for (const changeResult of changes) {
    const h = await harness({ changeResult }), result = h.prepareCoupledGridLevels(h.request, { solveIntermediateLevel: h.solve });
    assert.equal(result.prepared, null); assert.equal(result.sourceResult, h.request.sourceResult);
    assert.equal(result.sourceGridIntervals, 16); assert.equal(result.diagnostics.failure.stage, 'intermediate-solve');
    assert.equal(h.calls.preparations.length, 1);
  }
});

test('observer and intermediate-solver cancellation propagates the original thrown value', async () => {
  for (const where of ['preparing', 'prepared', 'converged', 'solve']) {
    const h = await harness(), stop = Object.freeze({ stop: where });
    assert.throws(() => h.prepareCoupledGridLevels(h.request, {
      onLevel: e => { if (e.stage === where) throw stop; },
      solveIntermediateLevel: (p, level) => { if (where === 'solve') throw stop; return h.solve(p, level); },
    }), value => value === stop);
    assert.ok(h.calls.preparations.length <= 1);
  }
});

test('the explicit native predictor reaches every level while default interpolation preserves the old preparation controls', async () => {
  for (const blPredictor of ['interpolate', 'xfoil-mrchdu']) {
    const h = await harness();
    h.prepareCoupledGridLevels(h.request, { solveIntermediateLevel: h.solve, blPredictor });
    assert.equal(h.calls.preparations.length, 3);
    for (const call of h.calls.preparations) assert.deepEqual(call.controls,
      { tolerance: 1e-10, maxNodes: 50000, ...(blPredictor === 'xfoil-mrchdu' ? { blPredictor } : {}) });
  }
  const h = await harness();
  for (const blPredictor of ['unknown', null, false])
    assert.throws(() => h.prepareCoupledGridLevels(h.request, { blPredictor, solveIntermediateLevel: h.solve }), /BL predictor/);
  assert.deepEqual(h.calls, { constructs: 0, evaluations: 0, preparations: [], solves: [] });
});

test('a returned numerical failure preserves its code and diagnostics alongside the last true root', async () => {
  const h = await harness(), diagnostics = { cell: { i: 17, group: 1, tube: 2 }, pressure: -1 };
  const result = h.prepareCoupledGridLevels(h.request, { solveIntermediateLevel: () => ({ result: null,
    failure: { reason: 'pressure rejected', code: 'streamtube-interface-pressure', diagnostics } }) });
  assert.equal(result.prepared, null); assert.equal(result.sourceResult, h.request.sourceResult);
  assert.equal(result.diagnostics.failure.code, 'streamtube-interface-pressure');
  assert.deepEqual(result.diagnostics.failure.diagnostics, diagnostics);
  diagnostics.cell.i = 200;
  assert.equal(result.diagnostics.failure.diagnostics.cell.i, 17);
});


test('every intermediate grid retains inherited shear metadata and the final seed carries it separately from physics', async () => {
  for (const coordinate of [undefined, 'linear', 'logarithmic']) {
    const h = await harness();
    if (coordinate !== undefined) h.request.sourceResult.checkpoint.continuation = { blUpdate: 'xfoil', shearCoordinate: coordinate };
    const before = clone(h.request);
    const output = h.prepareCoupledGridLevels(h.request, { solveIntermediateLevel: h.solve });
    assert.equal(output.diagnostics.finalSeedPrepared, true);
    assert.deepEqual(h.calls.solves.map(s => s.nominalGridIntervals), [32, 64]);
    assert.equal(shearPolicy.coupledCheckpointShearCoordinate(output.sourceResult.checkpoint), coordinate ?? 'linear');
    assert.equal(output.prepared.iterationControls.shearCoordinate, coordinate ?? 'linear');
    assert.equal(Object.hasOwn(output.prepared.options, 'shearCoordinate'), false);
    assert.deepEqual(h.request, before);
  }
});

test('an intermediate root cannot silently change or invalidate its inherited shear coordinate', async () => {
  for (const [sourceCoordinate, changedCoordinate] of [['linear', 'logarithmic'], ['logarithmic', 'linear'],
    ['logarithmic', undefined], ['linear', null], ['linear', 'wrong']]) {
    const h = await harness({ changeResult: result => {
      result.checkpoint.continuation = { blUpdate: 'xfoil', ...(changedCoordinate === undefined ? {} : { shearCoordinate: changedCoordinate }) };
    } });
    h.request.sourceResult.checkpoint.continuation = { blUpdate: 'xfoil', shearCoordinate: sourceCoordinate };
    const before = clone(h.request.sourceResult);
    const output = h.prepareCoupledGridLevels(h.request, { solveIntermediateLevel: h.solve });
    assert.equal(output.prepared, null);
    assert.equal(output.sourceResult, h.request.sourceResult);
    assert.equal(output.sourceGridIntervals, 16);
    assert.equal(output.diagnostics.failure.stage, 'intermediate-solve');
    assert.match(output.diagnostics.failure.reason ?? output.diagnostics.failure.message, /shear-coordinate/);
    assert.equal(h.calls.preparations.length, 1);
    assert.deepEqual(h.request.sourceResult, before);
  }
});

test('grid-level replay accepts coordinate roundoff but rejects changed geometry and exact residual mismatches', async () => {
  const h = await harness({ nodeDeparture: Number.EPSILON });
  const result = h.prepareCoupledGridLevels(h.request, { solveIntermediateLevel: h.solve });
  assert.ok(result.prepared);
  const moved = await harness({ nodeDeparture: 1e-8 });
  assert.throws(() => moved.prepareCoupledGridLevels(moved.request, { solveIntermediateLevel: moved.solve }), e =>
    e.code === 'GRID_LEVEL_REPLAY' && e.diagnostics.nodes.equivalent === false && e.diagnostics.residual === true);
  const stale = await harness({ changeResult: r => { r.residual[0] += Number.EPSILON; } });
  const rejected = stale.prepareCoupledGridLevels(stale.request, { solveIntermediateLevel: stale.solve });
  assert.equal(rejected.diagnostics.failure.code, 'GRID_LEVEL_REPLAY');
  assert.equal(rejected.diagnostics.failure.diagnostics.residual, false);
});

test('coordinate comparison rejects missing/nonfinite nodes and topology changes', () => {
  const a = nodes(10);
  assert.deepEqual(gridReplayDeparture(a, structuredClone(a)), { equivalent: true, maximum: 0 });
  for (const b of [undefined, [], [[[]]], [[[{ x: NaN, y: 0 }]]]])
    assert.equal(gridReplayDeparture(a, b).equivalent, false);
});
