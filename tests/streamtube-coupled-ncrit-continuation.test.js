import { observableFlow, observableBL } from '../src/euler/streamtube-flow-preview.js';
// SPDX-License-Identifier: GPL-2.0-or-later
// Controller tests use deterministic solver seams. Physical/numerical restart
// behavior is exercised independently by streamtube-coupled-ncrit-restart.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';

const original = fs.readFileSync(new URL('../src/euler/streamtube-coupled-ncrit-continuation.js', import.meta.url), 'utf8');
let nextHarness = 0;
const nodes = ncrit => [[[{ x: 0, y: ncrit }, { x: 0, y: ncrit + 1 }],
  [{ x: 1, y: ncrit }, { x: 1, y: ncrit + 1 }]]];
function root(ncrit, { converged = true, policy = 'admissible' } = {}) {
  const families = { euler: 1e-12, boundaryLayer: converged ? 2e-12 : .2, edgeMatching: 3e-13 };
  const checkpoint = { version: 1, families, restart: { input: { mach: .185, alpha: 6, bodies: [{ element: 0 }] },
    options: { ncrit, hkFloorLinearization: 'native' }, initialEuler: { x: [ncrit], nodes: nodes(ncrit), undisplacedNodes: nodes(ncrit) },
    initialBL: [.03, .001, .002, 1] }, continuation: { iterationGeometry: 'ises-sampled', stepAcceptance: policy, stagnationLimiter: 'prose' } };
  return { converged, checkpoint, families, reason: converged ? 'residual' : 'iteration limit',
    history: [{ iteration: 0 }, { iteration: 1, ...families }],
    mesh: { initialization: { gridSmoothing: { enabled: true } }, quality: { valid: true }, nodes: nodes(ncrit) },
    conditions: { ncrit }, rootIdentity: `N${ncrit}` };
}
async function harness({ warnings = () => [], solve = () => true } = {}) {
  const preparations = [], solves = [], systems = [], key = `__ncritContinuationTest${nextHarness++}`;
  const create = (input, options) => {
    const ncrit = options.ncrit;
    const system = { initial: [ncrit], euler: { layout: { bodies: input.bodies } },
      bl: { scale: .001, wakes: [], stations: [{ id: 0, kind: 'surface', body: 0, side: 'upper', i: 1 }], surfaces: [{ body: 0, side: 'upper', ids: [0] }] },
      evaluate: () => ({ outer: { nodes: options.initialEuler.nodes, undisplacedNodes: options.initialEuler.undisplacedNodes,
        cells: [[[{ interfacePressure: { lower: 1, upper: 1 } }]]] } }) };
    systems.push(ncrit); return system;
  };
  globalThis[key] = { observableFlow, observableBL,
    prepareCoupledNcritRestart: (selected, source, options) => {
      preparations.push({ selected, source: structuredClone(source), options });
      const checkpoint = structuredClone(root(selected).checkpoint);
      checkpoint.continuation = structuredClone(source.continuation);
      return { checkpoint, diagnostics: { warnings: warnings(selected, source), targetConverged: false } };
    },
    createCoupledStreamtubeBody: create,
    streamtubeMeshSnapshot: state => {
      assert.ok(state.nodes, 'The real mesh snapshot requires top-level physical nodes.');
      assert.deepEqual(state.nodes, state.flow.nodes, 'Published mesh and pressure use the same physical grid.');
      return { nodes: structuredClone(state.nodes), initialization: {}, quality: { valid: true } };
    },
    solveCoupledStreamtubeIses: (_, options) => {
      const { resume } = options, ncrit = resume.restart.options.ncrit;
      solves.push({ ncrit, resume: structuredClone(resume), maxIterations: options.maxIterations, tolerance: options.tolerance });
      for (const name of ['iterationGeometry', 'stepAcceptance', 'stagnationLimiter'])
        assert.equal(options[name], resume.continuation[name], `Exact saved ${name} is mandatory for a real resume.`);
      const converged = solve(ncrit, resume), value = root(ncrit, { converged, policy: options.stepAcceptance });
      value.checkpoint.continuation = structuredClone(resume.continuation);
      const system = create(value.checkpoint.restart.input, value.checkpoint.restart.options);
      const iteration = { iteration: 1, ...value.families };
      options.onIteration?.(iteration);
      options.onCheckpoint?.(value.checkpoint, { history: value.history });
      options.onMesh?.({ system: system.euler, nodes: nodes(ncrit), flow: { nodes: nodes(ncrit), undisplacedNodes: nodes(ncrit),
        cells: [[[{ interfacePressure: { lower: 1, upper: 1 } }]]] }, iteration });
      return value;
    },
  };
  const text = original.replace(/^import \{([^}]+)\} from '[^']+';/gm,
    (_, exports) => `const {${exports}} = globalThis.${key};`);
  const module = await import('data:text/javascript;base64,' + Buffer.from(text).toString('base64'));
  delete globalThis[key];
  return { ...module, preparations, solves, systems };
}

test('coarse native warnings halve once, then return the last true root for the requested-grid handoff', async () => {
  const h = await harness({ warnings: selected => selected === 8 ? ['native local warning'] : [] });
  const source = root(4), before = structuredClone(source), flows = [], stages = [], meshes = [], iterations = [];
  const p = h.continueCoupledNcrit(source, { targetNcrit: 9, phase: 'coarse', maxIterations: 40,
    onFlow: x => flows.push(x), onStage: x => stages.push(x), onMesh: x => meshes.push(x), onIteration: x => iterations.push(x) });
  assert.deepEqual(source, before);
  assert.deepEqual(h.preparations.map(x => [x.source.restart.options.ncrit, x.selected]), [[4, 5], [5, 6], [6, 7], [7, 8], [7, 7.5], [7.5, 8]]);
  assert.deepEqual(h.solves.map(x => x.ncrit), [5, 6, 7, 7.5]);
  assert.ok(h.preparations.every(x => x.options.blPredictor === 'xfoil-mrchdu' && x.options.requestedNcrit === 9));
  assert.ok(h.solves.every(x => x.maxIterations === 40));
  assert.equal(p.result.converged, true); assert.equal(p.result.checkpoint.restart.options.ncrit, 7.5);
  assert.equal(p.result.rootIdentity, 'N7.5'); assert.equal(p.diagnostics.reachedTarget, false);
  assert.equal(p.diagnostics.actualNcrit, 7.5); assert.equal(p.diagnostics.targetNcrit, 9);
  assert.equal(p.diagnostics.stateConverged, true); assert.equal(p.diagnostics.attempts.filter(a => a.preparation.warnings.length).length, 2);
  const retained = flows.filter(x => x.ncritContinuation.stateConverged);
  assert.deepEqual(retained.map(x => x.actualNcrit), [7, 7.5]);
  for (const packet of flows) {
    assert.equal(packet.actualNcrit, packet.checkpoint.restart.options.ncrit);
    assert.equal(packet.targetNcrit, 9); assert.equal(packet.ncritContinuation.reachedTarget, false);
    assert.deepEqual(packet.flow.nodes, packet.checkpoint.restart.initialEuler.nodes);
  }
  assert.ok(iterations.every(x => x.stage === 'coupled-coarse-initialization' && x.targetNcrit === 9));
  assert.ok(stages.at(-1).ncritContinuation.stateConverged);
  assert.ok(meshes.every(x => x.initialization.gridSmoothing.enabled));
});

test('fine preserve retries keep converged intermediate sources and return honest labels when the target remains unsolved', async () => {
  const h = await harness({ solve: selected => selected !== 9 }), source = root(8, { policy: 'armijo' });
  source.gridSequence = { method: 'nested', actualCells: 10140 };
  const flows = [], p = h.continueCoupledNcrit(source, { targetNcrit: 9, phase: 'fine', onFlow: x => flows.push(x) });
  assert.deepEqual(h.preparations.map(x => [x.source.restart.options.ncrit, x.selected]),
    [[8, 9], [8, 8.5], [8.5, 9], [8.5, 8.75], [8.75, 9], [8.75, 8.875], [8.875, 9]]);
  assert.ok(h.preparations.every(x => x.options.blPredictor === 'preserve'));
  assert.equal(p.result.checkpoint.restart.options.ncrit, 8.875); assert.equal(p.result.converged, true);
  assert.equal(p.diagnostics.reachedTarget, false); assert.equal(p.diagnostics.stateConverged, true);
  assert.equal(p.result.checkpoint.continuation.stepAcceptance, 'armijo');
  assert.deepEqual(p.result.mesh.initialization.gridRefinement, source.gridSequence);
  assert.deepEqual(flows.filter(x => x.ncritContinuation.stateConverged).map(x => x.actualNcrit), [8, 8.5, 8.75, 8.875]);
  assert.equal(flows.at(-1).actualNcrit, 8.875); assert.equal(flows.at(-1).ncritContinuation.reachedTarget, false);
  assert.ok(flows.some(x => x.actualNcrit === 9 && !x.ncritContinuation.stateConverged));
});

test('integer and final-target clipping halve the actual failed increment without repeating identical seeds', async () => {
  const h = await harness({ solve: () => false });
  const p = h.continueCoupledNcrit(root(7.5), { targetNcrit: 9, phase: 'fine' });
  assert.deepEqual(h.preparations.map(x => x.selected), [8, 7.75, 7.625]);
  assert.equal(p.result.checkpoint.restart.options.ncrit, 7.5);
  const clipped = await harness({ solve: () => false });
  clipped.continueCoupledNcrit(root(8), { targetNcrit: 8.05, phase: 'fine' });
  assert.deepEqual(clipped.preparations.map(x => x.selected), [8.05]);
});

test('successful fine continuation certifies only the final requested criterion', async () => {
  const h = await harness(), source = root(7.5), before = structuredClone(source);
  const p = h.continueCoupledNcrit(source, { targetNcrit: 9, phase: 'fine' });
  assert.deepEqual(h.solves.map(x => x.ncrit), [8, 9]); assert.deepEqual(source, before);
  assert.equal(p.result.checkpoint.restart.options.ncrit, 9); assert.equal(p.result.conditions.ncrit, 9);
  assert.equal(p.result.converged, true); assert.equal(p.diagnostics.reachedTarget, true);
  assert.equal(p.diagnostics.stateConverged, true);
});

test('a bounded half-unit policy advances only complete roots and preserves a clipped target', async () => {
  const h = await harness(), source = root(5), before = structuredClone(source);
  const continued = h.continueCoupledNcrit(source, { targetNcrit: 7.2, phase: 'fine', maximumStep: .5 });
  assert.deepEqual(h.preparations.map(x => [x.source.restart.options.ncrit, x.selected]),
    [[5, 5.5], [5.5, 6], [6, 6.5], [6.5, 7], [7, 7.2]]);
  assert.deepEqual(source, before);
  assert.equal(continued.diagnostics.maximumStep, .5);
  assert.equal(continued.diagnostics.reachedTarget, true);
  assert.equal(continued.result.checkpoint.restart.options.ncrit, 7.2);
  const failing = await harness({ solve: n => n !== 5.5 });
  const retained = failing.continueCoupledNcrit(root(5), { targetNcrit: 5.5, phase: 'fine', maximumStep: .5 });
  assert.deepEqual(failing.preparations.map(x => [x.source.restart.options.ncrit, x.selected]),
    [[5, 5.5], [5, 5.25], [5.25, 5.5], [5.25, 5.375], [5.375, 5.5]]);
  assert.equal(retained.result.checkpoint.restart.options.ncrit, 5.375);
  assert.equal(retained.diagnostics.reachedTarget, false);
});

test('local refinement provenance follows successful, rejected and restored Ncrit frames', async () => {
  for (const location of ['automaticRefinement', 'mesh']) {
    const h = await harness({ solve: n => n !== 6 }), source = root(5), meshes = [];
    const provenance = { kind: 'transition-local', parentNx: 128, refinedNx: 170, nodeCount: 6156 };
    if (location === 'automaticRefinement') source.automaticRefinement = provenance;
    else source.mesh.initialization.gridRefinement = provenance;
    const continued = h.continueCoupledNcrit(source, { targetNcrit: 6, phase: 'fine', maximumStep: .5,
      onMesh: mesh => meshes.push(mesh) });
    assert.ok(meshes.length > 2);
    assert.ok(meshes.every(mesh => JSON.stringify(mesh.initialization.gridRefinement) === JSON.stringify(provenance)));
    assert.deepEqual(continued.result.mesh.initialization.gridRefinement, provenance);
    assert.equal(continued.diagnostics.actualNcrit, 5.875);
  }
});

test('invalid maximum Ncrit steps are rejected before preparing a candidate', async () => {
  const h = await harness();
  for (const [phase, maximumStep] of [['fine', null], ['fine', NaN], ['fine', Infinity],
    ['fine', 0], ['fine', .1], ['fine', 1.01], ['coarse', .25]])
    assert.throws(() => h.continueCoupledNcrit(root(5), { targetNcrit: 9, phase, maximumStep }),
      /Invalid coupled Ncrit continuation/);
  assert.deepEqual(h.preparations, []); assert.deepEqual(h.solves, []);
});

test('Ncrit continuation carries saved linear policy and fallback state through every resumed root', async () => {
  for (const policy of [undefined, 'station', 'station-auto']) {
    const h = await harness(), source = root(7.5);
    if (policy !== undefined) source.checkpoint.continuation.linearOrdering = policy;
    if (policy === 'station-auto') source.checkpoint.continuation.stationFallback = true;
    const before = structuredClone(source), checkpoints = [];
    const prepared = h.continueCoupledNcrit(source, { targetNcrit: 9, phase: 'fine',
      onIterationCheckpoint: value => checkpoints.push(value) });
    assert.deepEqual(h.solves.map(s => s.ncrit), [8, 9]);
    for (const cp of [prepared.result.checkpoint, ...checkpoints, ...h.solves.map(s => s.resume)]) {
      assert.deepEqual(cp.continuation, before.checkpoint.continuation);
      assert.equal(Object.hasOwn(cp.restart.options, 'linearOrdering'), false);
      assert.equal(Object.hasOwn(cp.restart.input, 'linearOrdering'), false);
    }
    assert.deepEqual(source, before);
  }
});

test('observer cancellation propagates unchanged without launching a retry', async () => {
  for (const name of ['onStage', 'onIteration', 'onIterationCheckpoint', 'onMesh', 'onFlow']) {
    const h = await harness(), cancelled = new Error(`Cancelled at ${name}`); cancelled.name = 'AbortError';
    const source = root(7), before = structuredClone(source);
    assert.throws(() => h.continueCoupledNcrit(source, { targetNcrit: 9, phase: 'fine', [name]: () => { throw cancelled; } }),
      error => error === cancelled);
    assert.deepEqual(source, before); assert.ok(h.preparations.length <= 1); assert.ok(h.solves.length <= 1);
  }
});

test('invalid control/source packets reject before any preparation or global solver call', async () => {
  const h = await harness();
  for (const [source, options] of [[undefined, {}], [root(8, { converged: false }), { targetNcrit: 9, phase: 'fine' }],
    [root(8), { targetNcrit: Infinity, phase: 'fine' }], [root(8), { targetNcrit: 7, phase: 'coarse' }],
    [root(8), { targetNcrit: 9, phase: 'fine', maxIterations: -1 }], [root(8), { targetNcrit: 9, phase: 'fine', tolerance: 0 }]])
    assert.throws(() => h.continueCoupledNcrit(source, options), /Invalid coupled Ncrit continuation/);
  assert.deepEqual(h.preparations, []); assert.deepEqual(h.solves, []);
});
