// SPDX-License-Identifier: GPL-2.0-or-later
// Synthetic metadata tests of the public gate and controlled adapter seams.
// They do not manufacture a numerical root or execute an evaluator/Newton/LU.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { coupledAssemblyConditions } from '../src/euler/streamtube-coupled-assembly.js';
import { coupledMachPlan } from '../src/euler/tests/streamtube-coupled-mach-assembly.js';

const runtime = 'src/euler/tests/streamtube-coupled-mach-assembly.js';
const archive = 'docs/surface-pchip-refinement/mach-continuation/before/streamtube-coupled-mach-assembly.js';
const dataURL = code => 'data:text/javascript;base64,' + Buffer.from(code).toString('base64');
let serial = 0;
async function adapter(file, hooks) {
  const key = `__coupledRegionTest${serial++}`;
  if (hooks) globalThis[key] = hooks;
  const imports = hooks ? {
    './streamtube-coupled-flow-restart.js': dataURL(`export const initializeCoupledStreamtubeFromFlow=globalThis['${key}'].source;`),
    './streamtube-coupled-automatic.js': dataURL(`export const solveCoupledStreamtubeAutomatic=globalThis['${key}'].automatic;`),
    './streamtube-coupled-ises.js': dataURL('export const solveCoupledStreamtubeIses=()=>{throw Error("Unexpected ISES call");};'),
    './streamtube-coupled-hybrid-certification.js': dataURL('export const certifyCoupledStreamtubeHybrid=()=>{throw Error("Unexpected certification call");};'),
  } : {};
  try {
    const text = fs.readFileSync(file, 'utf8').replace(/from '([^']+)'/g, (_, specifier) => {
      // The archived adapter predates consolidation of equation controls;
      // current adapters also moved their test-only neighbors one level down.
      // Preserve the archived implementation, resolving only those imports.
      const name = path.basename(specifier);
      const alias = name === 'streamtube-equation-controls.js'
        ? path.resolve('src/euler/streamtube-equation-selection.js') : undefined;
      const base = file === runtime ? 'src/euler/tests' : 'src/euler';
      const original = path.resolve(base, specifier);
      const resolved = alias ?? (fs.existsSync(original) ? original : path.resolve('src/euler/tests', name));
      return `from '${imports[specifier] ?? imports[`./${name}`] ?? pathToFileURL(resolved).href}'`;
    });
    return await import(dataURL(text));
  } finally { if (hooks) delete globalThis[key]; }
}
const old = await adapter(archive);

function parentFixture(elements = 1, regional = true, ismom = 3) {
  const input = intrinsicBodyFixture({ elements, bodySegments: 4, tubes: 2, mach: .2 });
  input.bodies = input.bodies.map((body, element) => ({ ...body, element }));
  const sourceCase = { elements: input.bodies.map(b => ({ points: b.points })), flowModel: 'streamtube-grid',
    quadBoundaryLayers: true, mach: .2, alpha: 0, referenceChord: 2, reynolds: 2e6,
    transitionMode: 'automatic', ncrit: 9, materialTrips: input.bodies.map(() => [1, 1]),
    ...(ismom === undefined ? {} : { eulerIsmom: ismom }) };
  const settings = coupledAssemblyConditions(sourceCase, input.bodies, 1);
  const topology = { nx: input.outerLower.length - 1, tubes: input.weights.map(row => row.length),
    bodies: input.bodies.map(({ leadingIndex, trailingIndex }) => ({ leadingIndex, trailingIndex })) };
  let cut = 0;
  const region = { version: 1, parentCheckpointSha256: 'ab'.repeat(32), topology,
    regions: input.bodies.map((b, body) => {
      cut += topology.tubes[body];
      return { body, throughRow: b.leadingIndex + 2, lowerTube: cut - 1, upperTube: cut + 1 };
    }) };
  input.streamwiseMode = 'hybrid';
  input.hybrid = { epsilonP: 1e-5, ...(ismom === undefined ? {} : { ismom }), ...(regional ? { entropyRegions: region } : {}) };
  input.upwind = { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } };
  const families = { euler: 1e-12, boundaryLayer: 2e-12, edgeMatching: 3e-12 };
  const checkpoint = { version: 1, families, restart: { input,
    options: { ...settings.options, transitionMode: 'automatic', transitionState: Array(2 * elements).fill(1),
      blThermodynamics: 'historical-common-isentrope' },
    // Placeholder finite arrays are deliberate: only metadata/seam contracts
    // are exercised. The real restart adapter performs full replay validation.
    initialEuler: { x: [0], nodes: [[[ { x: 0, y: 0 } ]]], undisplacedNodes: [[[ { x: 0, y: 0 } ]]] }, initialBL: [0, 1, 2, 1] },
    continuation: { fractions: input.bodies.map(b => Array.from({ length: b.leadingIndex + 1 }, (_, i) => i / b.leadingIndex)),
      lastRedistributedStagnation: input.bodies.map(() => .5), preferredOrdering: 'amd', pivotTolerance: .001,
      linearOrdering: 'station', iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible',
      stagnationLimiter: 'listing', blUpdate: 'xfoil' } };
  return { model: 'research-streamtube-euler-bl', converged: true, mesh: { quality: { valid: true } },
    sourceCase, checkpoint, families, ...settings.normalization,
    initialization: { euler: { gridSmoothing: { enabled: true } } },
    solverSettings: { tolerance: 1e-10, streamwiseMode: 'hybrid', hybrid: structuredClone(input.hybrid),
      upwind: structuredClone(input.upwind), linearOrdering: 'station', pivotTolerance: .001 } };
}

test('exact explicit regions admit one/two-body warm plans without mutating or aliasing caller data', () => {
  for (const elements of [1, 2]) {
    const parent = parentFixture(elements), before = structuredClone(parent);
    const expected = structuredClone(parent.checkpoint.restart.input.hybrid.entropyRegions);
    const plan = coupledMachPlan({ ...parent.sourceCase, mach: .24 }, parent, { entropyRegions: expected });
    assert.equal(plan.route, 'warm-hybrid'); assert.equal(plan.sourceMach, .2); assert.equal(plan.targetMach, .24);
    assert.deepEqual(plan.entropyRegions, expected);
    plan.entropyRegions.regions[0].throughRow++;
    assert.deepEqual(parent, before);
    assert.deepEqual(expected, before.checkpoint.restart.input.hybrid.entropyRegions);
  }
});

test('regions require explicit warm ISMOM3 and exact requested topology, bounds and provenance', () => {
  const parent = parentFixture(2), request = { ...parent.sourceCase, mach: .24 };
  const region = parent.checkpoint.restart.input.hybrid.entropyRegions;
  assert.throws(() => coupledMachPlan(request, parent), /unchanged explicit historical hybrid/);
  assert.throws(() => coupledMachPlan(request, undefined, { entropyRegions: region }), /warm parent/);
  assert.throws(() => coupledMachPlan({ ...request, eulerIsmom: 4 }, parent, { entropyRegions: region }), /explicit ISMOM3/);
  const auto = { ...request }; delete auto.eulerIsmom;
  assert.throws(() => coupledMachPlan(auto, parent, { entropyRegions: region }), /explicit ISMOM3/);
  for (const bad of [null, {}, { ...region, version: 2 }, { ...region, parentCheckpointSha256: 'not-a-hash' }])
    assert.throws(() => coupledMachPlan(request, parent, { entropyRegions: bad }), /entropy regions/);
  for (const mutate of [r => r.topology.nx++, r => {
    r.topology.tubes[0]++; r.regions.forEach(b => { b.lowerTube++; b.upperTube++; });
  },
    r => r.topology.bodies[0].leadingIndex++, r => r.topology.bodies[0].trailingIndex++]) {
    const bad = structuredClone(region); mutate(bad);
    assert.throws(() => coupledMachPlan(request, parent, { entropyRegions: bad }), /topology does not match/);
  }
  for (const mutate of [r => r.regions[0].throughRow++, r => r.regions[0].lowerTube--,
    r => { r.parentCheckpointSha256 = 'cd'.repeat(32); }]) {
    const different = structuredClone(region); mutate(different);
    assert.throws(() => coupledMachPlan(request, parent, { entropyRegions: different }), /unchanged explicit historical hybrid/);
  }
  const alteredParent = structuredClone(parent); alteredParent.checkpoint.restart.input.hybrid.entropyRegions.regions[0].throughRow++;
  assert.throws(() => coupledMachPlan(request, alteredParent, { entropyRegions: region }), /unchanged explicit historical hybrid/);
  const noRegion = parentFixture(2, false);
  assert.throws(() => coupledMachPlan(request, noRegion, { entropyRegions: region }), /unchanged explicit historical hybrid/);
  const wrongMode = structuredClone(parent); wrongMode.checkpoint.restart.input.hybrid.ismom = 4;
  assert.throws(() => coupledMachPlan(request, wrongMode, { entropyRegions: region }), /ISMOM3 checkpoint/);
});

test('omitted regions retain archived default plans, schema and regional-parent rejection', () => {
  for (const mode of [undefined, 1, 2, 3, 4]) {
    const parent = parentFixture(2, false, mode);
    if (mode === undefined) {
      delete parent.sourceCase.eulerIsmom; delete parent.checkpoint.restart.input.hybrid.ismom;
      delete parent.solverSettings.hybrid.ismom;
    }
    for (const mach of [.2, .24, .4]) {
      const request = { ...parent.sourceCase, mach };
      assert.deepEqual(coupledMachPlan(request, parent), old.coupledMachPlan(request, parent));
      assert.deepEqual(coupledMachPlan(request), old.coupledMachPlan(request));
      assert.equal(Object.hasOwn(coupledMachPlan(request, parent), 'entropyRegions'), false);
    }
  }
  const parent = parentFixture(), request = { ...parent.sourceCase, mach: .24 };
  let previous; try { old.coupledMachPlan(request, parent); } catch (error) { previous = error.message; }
  assert(previous); assert.throws(() => coupledMachPlan(request, parent), error => error.message === previous);
});

function controlledHooks(parent) {
  const calls = { source: [], automatic: [] };
  const hooks = {
    source(mach, cp, options) {
      calls.source.push({ mach, checkpoint: structuredClone(cp), options: structuredClone(options) });
      return { system: { euler: { conditions: { lengthScale: parent.solverLength } }, bl: { transitionMode: 'automatic' } } };
    },
    automatic(targetMach, options) {
      const cp = structuredClone(options.initialCheckpoint);
      calls.automatic.push({ targetMach, checkpoint: cp, stageMaxIterations: options.stageMaxIterations,
        maxWakeRecoveries: options.maxWakeRecoveries, blPredictor: options.blPredictor });
      options.onStage?.({ mach: targetMach });
      const event = { mach: cp.restart.input.mach, targetMach, stage: 'coupled-wake-grid', wakeRecovery: 1 };
      options.onIteration?.({ ...event, iteration: 0 });
      options.onCheckpoint?.(structuredClone(cp), { ...event, kind: 'accepted', reachedTarget: false });
      assert.deepEqual(options.initialCheckpoint, cp, 'Callback mutation must not change the forwarded checkpoint.');
      return { checkpoint: cp, solverInput: cp.restart.input, families: cp.families,
        conditions: { mach: cp.restart.input.mach }, mesh: { quality: { valid: true } },
        converged: false, stateConverged: true, continuation: { boundaryLayerReinitializedOnWarmRestart: false } };
    },
  };
  return { hooks, calls };
}

test('controlled public continuation preserves exact region/station/pivot through source, callbacks and automatic boundary', async () => {
  const parent = parentFixture(2);
  parent.solverSettings.linearOrdering = 'auto'; parent.solverSettings.pivotTolerance = 1;
  const before = structuredClone(parent), request = { ...parent.sourceCase, mach: .24 };
  const expected = structuredClone(parent.checkpoint.restart.input.hybrid.entropyRegions), { hooks, calls } = controlledHooks(parent);
  const current = await adapter(runtime, hooks), events = [];
  const result = current.solveCoupledStreamtubeMach(request, { parentResult: parent, entropyRegions: expected,
    maxStages: 1, maxWakeRecoveries: 0,
    onPrepared: packet => {
      assert.deepEqual(packet.checkpoint, parent.checkpoint);
      assert.deepEqual(packet.parentResult.solverSettings.hybrid, parent.checkpoint.restart.input.hybrid);
      assert.equal(packet.parentResult.solverSettings.linearOrdering, 'station');
      assert.equal(packet.parentResult.solverSettings.pivotTolerance, .001);
      assert.equal(current.coupledMachPlan(request, packet.parentResult, { entropyRegions: expected }).route, 'warm-hybrid');
      packet.checkpoint.restart.input.hybrid.entropyRegions.regions[0].throughRow++;
      packet.parentResult.solverSettings.hybrid.entropyRegions.regions[0].throughRow++;
    },
    onStage: event => events.push(event), onIteration: event => events.push(event),
    onCheckpoint: (cp, event) => {
      assert.deepEqual(cp.restart.input.hybrid.entropyRegions, expected);
      assert.equal(cp.continuation.linearOrdering, 'station'); assert.equal(cp.continuation.pivotTolerance, .001);
      events.push(event); cp.restart.input.hybrid.entropyRegions.regions[0].throughRow++; cp.continuation.pivotTolerance = 1;
    },
  });
  assert.equal(calls.source.length, 1); assert.equal(calls.automatic.length, 1);
  assert.deepEqual(calls.source[0].checkpoint, before.checkpoint);
  assert.deepEqual(calls.automatic[0].checkpoint, before.checkpoint);
  assert.equal(calls.automatic[0].blPredictor, 'xfoil-mrchdu'); assert.equal(calls.automatic[0].maxWakeRecoveries, 0);
  assert.deepEqual(result.checkpoint, before.checkpoint); assert.deepEqual(result.solverSettings.hybrid, before.checkpoint.restart.input.hybrid);
  assert.equal(result.solverSettings.linearOrdering, 'station'); assert.equal(result.solverSettings.pivotTolerance, .001);
  assert.equal(result.converged, false); assert.equal(result.stateConverged, true);
  assert.equal(result.actualMach, .2); assert.equal(result.targetMach, .24);
  assert(events.some(event => event.stage === 'coupled-wake-grid' && event.actualMach === .2));
  assert.deepEqual(parent, before); assert.deepEqual(expected, before.checkpoint.restart.input.hybrid.entropyRegions);
  // A controlled returned checkpoint can record a different pivot after an
  // automatic-ordering solve. Output labels must follow that returned CP,
  // not the source metadata prepared before continuation.
  const autoParent = structuredClone(parent); delete autoParent.checkpoint.continuation.linearOrdering;
  autoParent.solverSettings.linearOrdering = 'station';
  const next = controlledHooks(autoParent), ordinary = next.hooks.automatic;
  next.hooks.automatic = (mach, options) => {
    const value = ordinary(mach, options); value.checkpoint.continuation.pivotTolerance = 1; return value;
  };
  const nextAdapter = await adapter(runtime, next.hooks);
  const changed = nextAdapter.solveCoupledStreamtubeMach(request, { parentResult: autoParent, entropyRegions: expected,
    onPrepared: packet => {
      assert.equal(packet.parentResult.solverSettings.linearOrdering, 'auto');
      assert.equal(packet.parentResult.solverSettings.pivotTolerance, .001);
    } });
  assert.equal(changed.solverSettings.linearOrdering, 'auto'); assert.equal(changed.solverSettings.pivotTolerance, 1);
});

test('public ISMOM4 Mach continuation preserves private station history without adding physical or GUI controls', async () => {
  for (const policy of [undefined, 'station', 'station-auto']) {
    const parent = parentFixture(2, false, 4);
    delete parent.checkpoint.continuation.linearOrdering;
    if (policy !== undefined) parent.checkpoint.continuation.linearOrdering = policy;
    if (policy === 'station-auto') parent.checkpoint.continuation.stationFallback = true;
    const before = structuredClone(parent), request = { ...parent.sourceCase, mach: .24 };
    const { hooks, calls } = controlledHooks(parent), current = await adapter(runtime, hooks);
    const checkpoints = [], result = current.solveCoupledStreamtubeMach(request, { parentResult: parent,
      maxStages: 1, maxWakeRecoveries: 0, onCheckpoint: cp => checkpoints.push(cp) });
    for (const cp of [result.checkpoint, ...checkpoints, ...calls.source.map(s => s.checkpoint),
      ...calls.automatic.map(s => s.checkpoint)]) {
      assert.deepEqual(cp.continuation, before.checkpoint.continuation);
      assert.deepEqual(cp.restart.input.hybrid, before.checkpoint.restart.input.hybrid);
      assert.deepEqual(cp.restart.input.upwind, before.checkpoint.restart.input.upwind);
      assert.equal(cp.restart.input.hybrid.ismom, 4);
      assert.equal(Object.hasOwn(cp.restart.input, 'linearOrdering'), false);
      assert.equal(Object.hasOwn(cp.restart.options, 'linearOrdering'), false);
    }
    assert.deepEqual(parent, before); assert.equal(request.eulerIsmom, 4);
    assert.equal(Object.hasOwn(request, 'linearOrdering'), false);
  }
});

test('controlled omitted-option execution preserves archived output with current policy metadata', async () => {
  const parent = parentFixture(2, false), request = { ...parent.sourceCase, mach: .24 };
  const a = controlledHooks(parent), b = controlledHooks(parent);
  const current = await adapter(runtime, a.hooks), previous = await adapter(archive, b.hooks);
  const controls = { parentResult: parent, maxStages: 1, maxWakeRecoveries: 0 };
  const expected = previous.solveCoupledStreamtubeMach(request, controls);
  // These explicit policy labels were added after the archived adapter. No
  // numerical state or continuation call may differ in this controlled case.
  Object.assign(expected.solverSettings, { shearCoordinate: 'linear', hkFloorLinearization: 'exact',
    maxMachStep: .05, maxFirstOrderRecoveries: 2 });
  assert.deepEqual(current.solveCoupledStreamtubeMach(request, controls), expected);
  assert.deepEqual(a.calls, b.calls);
});

test('corrupt source replay is propagated before automatic continuation or recovery', async () => {
  const parent = parentFixture(), expected = parent.checkpoint.restart.input.hybrid.entropyRegions;
  const failure = new Error('Controlled source replay rejection'), { hooks, calls } = controlledHooks(parent);
  hooks.source = () => { throw failure; };
  const current = await adapter(runtime, hooks);
  assert.throws(() => current.solveCoupledStreamtubeMach({ ...parent.sourceCase, mach: .24 }, {
    parentResult: parent, entropyRegions: expected }), error => error === failure);
  assert.equal(calls.automatic.length, 0);
});
