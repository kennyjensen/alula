// SPDX-License-Identifier: GPL-2.0-or-later
// Adapter contracts only: explicit numerical/mesh stubs, no global flow solve.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import * as shearPolicy from '../src/euler/streamtube-coupled-shear-policy.js';
import { observableFlow, observableBL } from '../src/euler/streamtube-flow-preview.js';

const clone = structuredClone;
let serial = 0;
async function harness(config = {}) {
  const calls = { solves: [], meshes: [] }, key = `__gridLevelSolve${serial++}`;
  const prepared = { input: { mach: .185, alpha: 6, streamwiseMode: 'hybrid',
    hybrid: { ismom: 4, epsilonP: 1e-5 }, upwind: { mcrit: .99, mucon: 1 } },
    options: { reynolds: 2449255.3983699363, ncrit: 7.5, transitionMode: 'automatic',
      tripFractions: [[1, 1], [1, 1]], transitionState: [7, 9, 11, 13], hkFloorLinearization: 'native',
      blThermodynamics: 'historical-common-isentrope', edgeMatching: 'section-velocity' },
    initialEuler: { x: [2, 3], nodes: [[[{ x: 0, y: 1 }]]], undisplacedNodes: [[[{ x: 0, y: 0 }]]] },
    initialBL: [.1, .001, .002, 1], transfer: { nominalGridIntervals: 32, child: { nx: 260, tubes: [12, 15, 12] } },
    system: { bl: { scale: .001, wakes: [], stations: [{ id: 31, kind: 'surface', body: 1, side: 'upper', i: 120 }],
      surfaces: [{ body: 1, side: 'upper', ids: [31] }] } } };
  const level = { index: 0, nominalGridIntervals: 32, requestedGridIntervals: 32, targetGridIntervals: 128, targetNcrit: 9 };
  const normalization = { referenceReynolds: 2.7e6, kernelReynolds: prepared.options.reynolds,
    referenceChord: 1, solverLength: .9071316290259023 };
  const stubs = { ...shearPolicy, observableFlow, observableBL,
    streamtubeMeshSnapshot(state) {
      calls.meshes.push(clone(state.nodes));
      return { nodes: clone(state.nodes), quality: { valid: true }, initialization: {}, iteration: clone(state.iteration) };
    },
    solveCoupledStreamtubeIses(input, options) {
      calls.solves.push({ input: clone(input), options: Object.fromEntries(Object.entries(options).filter(([, v]) => typeof v !== 'function')) });
      if (Object.hasOwn(config, 'error')) throw config.error;
      let checkpoint, state;
      for (const iteration of [0, 1]) {
        const nodes = [[[{ x: iteration, y: 10 * iteration }]]], undisplacedNodes = [[[{ x: iteration, y: -1 }]]];
        checkpoint = { version: 1, marker: iteration, restart: { input: clone(input), options: clone(prepared.options),
          initialEuler: { x: [iteration], nodes: clone(nodes), undisplacedNodes: clone(undisplacedNodes) }, initialBL: clone(prepared.initialBL) } };
        const h = { iteration, residual: iteration ? 1e-12 : 1, marker: iteration };
        options.onIteration?.(h);
        options.onCheckpoint?.(checkpoint, { kind: iteration ? 'accepted' : 'initial', marker: iteration });
        state = { system: { layout: { bodies: [{ element: 1 }, { element: 0 }] } }, nodes,
          flow: { nodes, undisplacedNodes, cells: [[[{ interfacePressure: { lower: iteration + 1, upper: iteration + 2 } }]]] }, iteration: h };
        options.onMesh?.(state);
      }
      return { converged: config.converged ?? true, reason: config.converged === false ? 'iteration limit' : 'residual',
        checkpoint, mesh: stubs.streamtubeMeshSnapshot(state) };
    },
  };
  globalThis[key] = stubs;
  const source = fs.readFileSync(new URL('../src/euler/streamtube-coupled-grid-levels.js', import.meta.url), 'utf8');
  const isolated = `const { coupledFreshShearCoordinate, observableFlow, observableBL, streamtubeMeshSnapshot, solveCoupledStreamtubeIses } = globalThis[${JSON.stringify(key)}];\n`
    + source.slice(source.indexOf('export function solveCoupledGridLevel('));
  const module = await import('data:text/javascript;base64,' + Buffer.from(isolated).toString('base64'));
  delete globalThis[key];
  return { ...module, calls, prepared, level, controls: { maxIterations: 40, tolerance: 1e-10, normalization, startupAttempt: 1 } };
}

test('each published intermediate pressure frame pairs its own checkpoint, physical mesh and refined BL station topology', async () => {
  const h = await harness(), before = clone(h.prepared), stages = [], iterations = [], checkpoints = [], meshes = [], flows = [];
  const result = h.solveCoupledGridLevel(h.prepared, h.level, { ...h.controls,
    onStage: e => stages.push(e), onIteration: e => iterations.push(e),
    onIterationCheckpoint: (cp, details) => checkpoints.push({ cp, details }),
    onMesh: (mesh, phase, stage) => meshes.push({ mesh, phase, stage }), onFlow: frame => flows.push(frame) });
  assert.equal(result.converged, true); assert.equal(h.calls.solves.length, 1);
  assert.deepEqual(meshes.map(m => m.phase), ['initial', 'solving']);
  assert.equal(flows.length, 2); assert.equal(checkpoints.length, 2);
  for (let i = 0; i < 2; i++) {
    assert.equal(flows[i].checkpoint.marker, i); assert.equal(flows[i].iteration.marker, i);
    assert.deepEqual(flows[i].flow.nodes, checkpoints[i].cp.restart.initialEuler.nodes);
    assert.deepEqual(flows[i].flow.undisplacedNodes, checkpoints[i].cp.restart.initialEuler.undisplacedNodes);
    assert.deepEqual(flows[i].flow.nodes, meshes[i].mesh.nodes);
    assert.deepEqual(flows[i].bl, h.prepared.system.bl);
    assert.deepEqual(flows[i].bodies, [{ element: 1 }, { element: 0 }]);
    assert.deepEqual(flows[i].normalization, h.controls.normalization);
    assert.equal(flows[i].flow.cells[0][0][0].interfacePressure.lower, i + 1);
    assert.equal(flows[i].gridLevel, 32); assert.equal(flows[i].requestedGridIntervals, 128);
    assert.equal(flows[i].actualNcrit, 7.5); assert.equal(flows[i].targetNcrit, 9);
    assert.equal(flows[i].stage, 'coupled-grid-refinement');
  }
  for (const value of [...stages, ...iterations, ...checkpoints.map(c => c.details)]) {
    assert.equal(value.gridLevel, 32); assert.equal(value.requestedGridIntervals, 128);
    assert.equal(value.startupAttempt, 1); assert.equal(value.stage, 'coupled-grid-refinement');
    assert.equal(value.actualNcrit, 7.5); assert.equal(value.targetNcrit, 9);
  }
  assert.deepEqual(result.mesh.initialization.gridRefinement, h.prepared.transfer);
  assert.deepEqual(h.prepared, before);
  flows[0].checkpoint.restart.options.ncrit = 99;
  assert.equal(checkpoints[0].cp.restart.options.ncrit, 7.5);
  assert.equal(result.checkpoint.restart.options.ncrit, 7.5);
});

test('ordinary ISES receives all prepared physical, transition and native-Hk controls unchanged', async () => {
  const h = await harness();
  h.solveCoupledGridLevel(h.prepared, h.level, h.controls);
  const call = h.calls.solves[0];
  assert.deepEqual(call.input, h.prepared.input);
  for (const [key, value] of Object.entries(h.prepared.options)) assert.deepEqual(call.options[key], value);
  assert.deepEqual(call.options.initialEuler, h.prepared.initialEuler); assert.deepEqual(call.options.initialBL, h.prepared.initialBL);
  assert.equal(call.options.maxIterations, 40); assert.equal(call.options.tolerance, 1e-10);
  assert.equal(call.options.shearCoordinate, 'linear');
  assert.equal(call.options.blUpdate, 'xfoil'); assert.equal(call.options.projectionGeometry, 'boundary-increment');
  assert.equal(call.options.iterationGeometry, 'ises-sampled'); assert.equal(call.options.stepAcceptance, 'admissible');
});

test('ordinary numerical exceptions are explicit failures and iteration-limit results remain unconverged', async () => {
  const error = Object.assign(new Error('Strict gas/grid update rejected.'), { code: 'controlled-domain', diagnostics: { cell: { i: 20 } } });
  const h = await harness({ error }), result = h.solveCoupledGridLevel(h.prepared, h.level, h.controls);
  assert.equal(result.result, null); assert.equal(result.failure.reason, error.message); assert.equal(result.failure.code, error.code);
  assert.deepEqual(result.failure.diagnostics, error.diagnostics);
  const limited = await harness({ converged: false }), stopped = limited.solveCoupledGridLevel(limited.prepared, limited.level, limited.controls);
  assert.equal(stopped.converged, false); assert.equal(stopped.reason, 'iteration limit'); assert.ok(stopped.checkpoint);
});

test('every observer cancellation preserves the thrown identity including undefined', async () => {
  for (const observer of ['onStage', 'onIteration', 'onIterationCheckpoint', 'onMesh', 'onFlow'])
    for (const stop of [undefined, null, 'user stop', Object.freeze({ cancelled: observer }), new Error('Cancel')]) {
      const h = await harness();
      assert.throws(() => h.solveCoupledGridLevel(h.prepared, h.level,
        { ...h.controls, [observer]: () => { throw stop; } }), value => Object.is(value, stop));
      assert.ok(h.calls.solves.length <= 1);
    }
});


test('each grid-level solve inherits explicit coordinates while native Hk alone still selects linear', async () => {
  for (const coordinate of [undefined, 'linear', 'logarithmic']) {
    const h = await harness();
    if (coordinate !== undefined) h.prepared.iterationControls = { shearCoordinate: coordinate };
    const before = clone(h.prepared);
    h.solveCoupledGridLevel(h.prepared, h.level, h.controls);
    assert.equal(h.calls.solves.length, 1);
    assert.equal(h.calls.solves[0].options.shearCoordinate, coordinate ?? 'linear');
    assert.equal(h.calls.solves[0].options.hkFloorLinearization, 'native');
    assert.equal(Object.hasOwn(h.prepared.options, 'shearCoordinate'), false);
    assert.deepEqual(h.prepared, before);
  }
});

test('invalid prepared coordinates are reported before invoking ISES', async () => {
  for (const iterationControls of [null, false, [], { shearCoordinate: null }, { shearCoordinate: 'wrong' }]) {
    const h = await harness(); h.prepared.iterationControls = iterationControls;
    const output = h.solveCoupledGridLevel(h.prepared, h.level, h.controls);
    assert.equal(output.result, null);
    assert.match(output.failure.reason, /Invalid/);
    assert.equal(h.calls.solves.length, 0);
  }
});
