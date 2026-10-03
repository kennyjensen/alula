import { observableFlow, observableBL } from '../src/euler/streamtube-flow-preview.js';
import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
const source = fs.readFileSync(new URL('../src/euler/streamtube-transition-recovery.js', import.meta.url), 'utf8');
let serial = 0;

async function harness({ native = false, initializationValid = true, initializationFailure } = {}) {
  const families = native ? { euler: 1e-12, boundaryLayer: 2e-12, edgeMatching: 3e-12 }
    : { euler: .001, boundaryLayer: .1, edgeMatching: .001 };
  const oldInput = { hybrid: { ismom: 4 }, mach: .2 }, options = {
    reynolds: 1e6, ncrit: 9, transitionMode: 'automatic', hkFloorLinearization: 'native', transitionState: [1], tripFractions: [[1, 1]] };
  const nodes = [[[{ x: 0, y: 0 }, { x: 0, y: 1 }], [{ x: 1, y: 0 }, { x: 1, y: 1 }]]];
  const cp = { version: 1, families, restart: { input: oldInput, options, initialEuler: { x: [1], nodes }, initialBL: [1, 2, 3, 4] },
    continuation: { projectionGeometry: 'boundary-increment', linearOrdering: 'station-auto',
      ...(native ? { iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', stagnationLimiter: 'listing', blUpdate: 'xfoil' } : {}) } };
  const result = { converged: native, checkpoint: cp, families, mesh: { initialization: { gridSmoothing: { enabled: true } } },
    boundaryLayer: { transitions: [{ body: 0, side: 'upper', kind: 'natural', forced: false }] } };
  const plan = { reason: 'boundary-layer-transition-stall', surfaces: [{ body: 0, side: 'upper' }],
    streamwiseSubdivisions: [2], normalFactor: 1 };
  const stations = [{ id: 44, kind: 'surface', body: 0, side: 'upper', i: 2 }], surfaces = [{ body: 0, side: 'upper', ids: [44] }];
  const captured = [], profiles = [], refinedInput = { ...oldInput, refined: true }, refinedMesh = {
    quality: { valid: initializationValid }, cells: [[0, 1, 2, 3], [1, 2, 3, 4]], initialization: {} };
  const stubs = { observableFlow, observableBL,
    createCoupledStreamtubeBody: () => ({ n: 5, initial: [0], euler: { layout: { nx: 1 } }, evaluate: () => ({ families }) }),
    refineCoupledStreamtubeBody: () => ({ input: refinedInput, options, initialEuler: { x: [2], nodes }, initialBL: [2, 3, 4, 5],
      system: { n: 9, euler: { layout: { nx: 2 } }, bl: { scale: .001, wakes: [], stations, surfaces } }, diagnostics: { quality: { valid: true } } }),
    prepareRefinedCoupledGridProfile(prepared, controls) {
      profiles.push({ prepared, controls });
      return { ...prepared, initialBL: [7, 8, 9, 10], value: { families },
        transfer: { profilePreparation: { warnings: [], physicalEulerInvariantsPreserved: true } } };
    },
    streamtubeMeshSnapshot: () => structuredClone(refinedMesh),
    solveCoupledStreamtubeIses(input, controls) {
      captured.push({ input, controls });
      if (initializationFailure && controls.maxIterations === 0) return initializationFailure;
      const checkpoint = controls.resume ?? { ...cp, restart: { ...cp.restart, input,
        initialEuler: controls.initialEuler, initialBL: controls.initialBL, options: { ...options, transitionState: [2] } } };
      const iteration = { iteration: 1, residual: 1e-12 };
      controls.onIteration?.(iteration); controls.onCheckpoint?.(checkpoint, { history: [iteration] });
      controls.onMesh?.({ system: { layout: { bodies: [{ leadingIndex: 1, trailingIndex: 3 }] } }, iteration,
        flow: { nodes, undisplacedNodes: nodes, cells: [[[{ interfacePressure: { lower: 2, upper: 3 } }]]] } });
      return { converged: true, checkpoint, mesh: structuredClone(refinedMesh), initialRedistribution: { accepted: initializationValid } };
    },
  };
  const key = `__recoveryObserver${++serial}`; globalThis[key] = stubs;
  const text = source.replace(/import \{([^}]+)\} from '[^']+';/g, (_, names) => `const {${names}} = globalThis[${JSON.stringify(key)}];`);
  const module = await import('data:text/javascript;base64,' + Buffer.from(`${text}\n//# sourceURL=recovery-observers-${serial}.js`).toString('base64'));
  return { ...module, result, plan, captured, profiles, stations, surfaces, release: () => { delete globalThis[key]; } };
}

test('recovery publishes paired refined station topology, checkpoint, pressure and mesh', async () => {
  const h = await harness(), events = [], normalization = { referenceChord: 1 };
  try {
    const result = h.recoverCoupledTransition(h.result, { plan: h.plan, maxIterations: 12, tolerance: 1e-10,
      normalization, startupAttempt: 1, onIterationCheckpoint: (cp, details) => events.push({ type: 'checkpoint', cp, details }),
      onMesh: mesh => events.push({ type: 'mesh', mesh }), onFlow: frame => events.push({ type: 'flow', frame }) });
    assert.deepEqual(events.map(e => e.type), ['checkpoint', 'mesh', 'flow']);
    const [checkpoint, mesh, flow] = events;
    assert.deepEqual(flow.frame.checkpoint, checkpoint.cp);
    assert.deepEqual(flow.frame.bl, { scale: .001, wakes: [], stations: h.stations, surfaces: h.surfaces });
    assert.equal(flow.frame.checkpoint.restart.input.refined, true);
    assert.equal(mesh.mesh.cells.length, 2);
    assert.deepEqual(flow.frame.automaticRefinement, result.automaticRefinement);
    assert.deepEqual(mesh.mesh.initialization.gridRefinement, result.automaticRefinement);
    assert.deepEqual(mesh.mesh.initialization.gridSmoothing, { enabled: true });
    assert.equal(flow.frame.actualNcrit, 9); assert.equal(flow.frame.mach, .2);
    assert.deepEqual(flow.frame.normalization, normalization);
  } finally { h.release(); }
});

test('recovery retains explicit physical and native-kernel controls and caller iteration cap', async () => {
  const h = await harness();
  try {
    h.recoverCoupledTransition(h.result, { plan: h.plan, maxIterations: 3, tolerance: 2e-11 });
    const { input, controls } = h.captured[0];
    assert.equal(input.hybrid.ismom, 4); assert.equal(input.mach, .2);
    for (const key of ['reynolds', 'ncrit', 'transitionMode', 'hkFloorLinearization', 'tripFractions'])
      assert.deepEqual(controls[key], h.result.checkpoint.restart.options[key]);
    assert.equal(controls.linearOrdering, 'station-auto'); assert.equal(controls.projectionGeometry, 'boundary-increment');
    assert.equal(controls.maxIterations, 3); assert.equal(controls.tolerance, 2e-11);
    delete h.result.checkpoint.continuation.linearOrdering;
    h.recoverCoupledTransition(h.result, { plan: h.plan, maxIterations: 1 });
    assert.equal(h.captured[1].controls.linearOrdering, 'auto');
  } finally { h.release(); }
});

test('all recovery observers propagate cancellation values, including undefined', async () => {
  for (const callback of ['onStage', 'onIteration', 'onIterationCheckpoint', 'onMesh', 'onFlow']) {
    for (const value of [new Error('cancel'), undefined]) {
      const h = await harness(); let caught = false;
      try {
        try { h.recoverCoupledTransition(h.result, { plan: h.plan, [callback]: () => { throw value; } }); }
        catch (error) { caught = true; assert.equal(error, value); }
        assert.equal(caught, true, callback);
      } finally { h.release(); }
    }
  }
});

test('a forced transition at the retained best state is not refined as a natural bubble', async () => {
  const h = await harness();
  try {
    h.result.boundaryLayer.transitions[0].forced = true;
    assert.throws(() => h.recoverCoupledTransition(h.result, { plan: h.plan }), /natural transition/);
    assert.equal(h.captured.length, 0);
  } finally { h.release(); }
});

test('native refinement establishes its chart once, prepares physical profiles, then resumes the same complete state', async () => {
  const h = await harness({ native: true });
  try {
    const result = h.recoverCoupledTransition(h.result, { plan: h.plan, maxIterations: 12, blPredictor: 'xfoil-mrchdu' });
    assert.equal(h.captured.length, 2); assert.equal(h.captured[0].controls.maxIterations, 0);
    assert.equal(h.captured[1].input, undefined); assert.equal(h.captured[1].controls.maxIterations, 12);
    assert.deepEqual(h.profiles[0].controls, { wakeWidthIncrement: true, reinitializeAmplification: true });
    assert.deepEqual(h.captured[1].controls.resume.restart.initialBL, [7, 8, 9, 10]);
    assert.deepEqual(h.captured[1].controls.resume.continuation, h.result.checkpoint.continuation);
    for (const name of ['iterationGeometry', 'stepAcceptance', 'stagnationLimiter', 'blUpdate', 'linearOrdering', 'projectionGeometry'])
      assert.equal(h.captured[1].controls[name], h.result.checkpoint.continuation[name]);
    assert.equal(result.automaticRefinement.parentNx, 1); assert.equal(result.automaticRefinement.refinedNx, 2);
    assert.deepEqual(result.automaticRefinement.profilePreparation.warnings, []);
  } finally { h.release(); }
});

test('native startup refinement refuses a nonroot or invalid new chart before any Newton solve', async () => {
  const invalid = await harness({ native: true, initializationValid: false });
  try {
    assert.throws(() => invalid.recoverCoupledTransition(invalid.result, { plan: invalid.plan, blPredictor: 'xfoil-mrchdu' }), /initialize its new grid/);
    assert.equal(invalid.captured.length, 1); assert.equal(invalid.captured[0].controls.maxIterations, 0);
    assert.equal(invalid.profiles.length, 0);
  } finally { invalid.release(); }
  const nonroot = await harness();
  try {
    assert.throws(() => nonroot.recoverCoupledTransition(nonroot.result, { plan: nonroot.plan, blPredictor: 'xfoil-mrchdu' }), /converged parent/);
    assert.equal(nonroot.captured.length, 0); assert.equal(nonroot.profiles.length, 0);
  } finally { nonroot.release(); }
});

test('a checkpointless native initializer preserves its original numerical rejection and stops preparation', async () => {
  const failure = { converged: false, reason: 'Coupled ISES initial redistribution rejected: invalid corner',
    families: { euler: .2, boundaryLayer: .3, edgeMatching: .4 }, lastRejectedStep: null,
    initialRedistribution: { accepted: false, rejection: 'invalid corner',
      admissibility: { kind: 'grid-convexity', firstCell: { group: 1, interval: 17, tube: 0 } } },
    mesh: { quality: { valid: true, invalidCells: [] } } };
  const before = structuredClone(failure), h = await harness({ native: true, initializationFailure: failure });
  try {
    assert.throws(() => h.recoverCoupledTransition(h.result, { plan: h.plan, blPredictor: 'xfoil-mrchdu' }), error => {
      assert.equal(error.code, 'coupled-refinement-initialization');
      assert.equal(error.diagnostics.reason, failure.reason);
      assert.deepEqual(error.diagnostics.initialRedistribution, failure.initialRedistribution);
      assert.deepEqual(error.diagnostics.families, failure.families);
      assert.deepEqual(error.diagnostics.meshQuality, failure.mesh.quality);
      assert.equal(error.diagnostics.checkpointAvailable, false);
      error.diagnostics.initialRedistribution.admissibility.firstCell.interval = 99;
      return true;
    });
    assert.deepEqual(failure, before);
    assert.equal(h.captured.length, 1); assert.equal(h.profiles.length, 0);
  } finally { h.release(); }
});
