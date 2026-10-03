import * as shearPolicy from '../src/euler/streamtube-coupled-shear-policy.js';
import { coupledCheckpointHkPolicy } from '../src/euler/streamtube-coupled-assembly.js';
// SPDX-License-Identifier: GPL-2.0-or-later
// Public controller with numerical dependencies replaced by explicit stubs.
// No physical systems, gas, BL, Jacobians, LU or Newton calls are performed.
// Driver rejection fixtures use the real cheap convex-grid predicates.
import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import { checkpointDataEqual } from '../src/euler/streamtube-nested-checkpoint.js';
import { streamtubeEquationControls } from '../src/euler/streamtube-equation-selection.js';
import { colderStreamtubeStartupMach } from '../src/euler/tests/streamtube-cold-startup.js';
import { assertConvexStreamtubeGrid } from '../src/geometry/streamtube-convex-step.js';
import { requireConvexGridUpdate, requireConvexPublishedGrid } from '../src/euler/streamtube-grid-update.js';
import { quadCoupledNcrit, quadCoupledNcritResult } from '../src/ui/quad-coupled-ncrit.js';

const fixtureLayout = (grid, bodies) => {
  let column = 0;
  return { bodies, densityCount: 1, elements: grid.length - 1, nx: grid[0].length - 1,
    tubes: grid.map(group => group[0].length - 1),
    nodes: grid.map(group => group.map(row => row.map(() => ({ kind: 'free', column: column++ })))) };
};

const copy = structuredClone, tolerance = 1e-10;
const caseData = () => ({ flowModel: 'streamtube-grid', quadBoundaryLayers: true, mach: .74, alpha: 2.68,
  elements: [{ name: 'controlled body', points: [{ x: 1, y: 0 }, { x: 0, y: 0 }, { x: 1, y: 0 }] }],
  referenceChord: 1, reynolds: 2.7e6, ncrit: 4, transitionMode: 'automatic', materialTrips: [[1, 1]],
  gridIntervals: 128, gridTubes: 11, gridEllipticSmoothing: true });
const families = { euler: 0, boundaryLayer: 0, edgeMatching: 0 };
const normal = { referenceChord: 1, referenceReynolds: 2.7e6, solverLength: 1, kernelReynolds: 2.7e6 };
const settings = () => ({ options: { reynolds: 2.7e6, ncrit: 4, edgeMatching: 'section-velocity',
  transitionMode: 'automatic', tripFractions: [[1, 1]] }, normalization: copy(normal), materialTrips: [[1, 1]], elementOrder: [0] });
const nodes = [[[{ x: 0, y: 0 }, { x: 1, y: 1 }]]];
const mesh = () => ({ quality: { valid: true }, vertices: [{ x: 0, y: 0 }, { x: 1, y: 1 }],
  initialization: { gridSmoothing: { converged: true } }, tag: 'one-original-mesh' });
function acceptedSource(mach) {
  const input = { mach, alpha: 2.68, bodies: [{ element: 0 }], flowModel: 'compressible', streamwiseMode: 'isentropic',
    wakeGeometry: 'independent-banks', wakeOutlet: 'banks' };
  const checkpoint = { version: 1, restart: { input, options: { ...settings().options, transitionState: [1, 1] },
    initialEuler: { x: [0], nodes: copy(nodes), undisplacedNodes: copy(nodes) }, initialBL: [1, 1, 2, 1] },
    continuation: { blUpdate: 'xfoil', iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', stagnationLimiter: 'listing' }, families: copy(families) };
  return { model: 'research-streamtube-euler-bl', converged: true, stateConverged: true, mach,
    mesh: mesh(), checkpoint, families: copy(families), ...copy(normal), solverInput: copy(input),
    sourceCase: { ...caseData(), mach }, initialization: { euler: { gridSmoothing: { converged: true } } },
    solverSettings: { tolerance, streamwiseMode: 'isentropic' } };
}
function capacity(stage = 'gas-initialization') {
  return Object.assign(new Error('original typed startup capacity failure'), { code: 'streamtube-sonic-capacity', stage,
    diagnostics: { stage: 'gas-initialization', capacityRatio: 1.104,
      stagnationDensityFallback: { code: 'streamtube-interface-pressure', reason: 'nonpositive interface pressure' } } });
}
let sequence = 0;
async function harness(config = {}) {
  const calls = { cold: [], retarget: [], reconstructed: [], automatic: [], certified: [], warmReplay: 0 };
  const originalPacket = { identity: 'original-geometry', systemFunction() {} };
  const request = caseData(), originalError = config.error ?? capacity();
  let certifications = 0;
  const stubs = { targetMachGridPlan: () => null, ...shearPolicy, coupledCheckpointHkPolicy, checkpointDataEqual, streamtubeEquationControls, colderStreamtubeStartupMach,
    coupledAssemblyConditions: settings,
    solveCoupledStreamtubeAssembly(input, options) {
      calls.cold.push({ caseData: copy(input), prepared: options.preparedEuler });
      options.onStage?.({ stage: 'euler', startupAttempt: 0 });
      options.onMesh?.(mesh(), 'initial', 'euler');
      options.onEulerPrepared?.(options.preparedEuler ?? originalPacket);
      if (calls.cold.length <= (config.failures ?? 2)) throw originalError;
      const result = acceptedSource(input.mach);
      config.mutateColdResult?.(result, input);
      return result;
    },
    retargetPreparedStreamtubeAssembly(packet, nextCase) {
      calls.retarget.push({ packet, caseData: copy(nextCase) });
      if (config.retargetError) throw config.retargetError;
      return { identity: `retarget-${nextCase.mach}`, systemFunction() {} };
    },
    certifyCoupledStreamtubeHybrid(checkpoint) {
      calls.certified.push(checkpoint.restart.input.mach);
      if (config.warmRecovery && certifications++ === 0)
        throw Object.assign(new Error('controlled incompatible source'), { code: 'coupled-hybrid-certification-speed-bias', diagnostics: { bias: .1 } });
      const next = copy(checkpoint);
      Object.assign(next.restart.input, { streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 },
        upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } });
      next.restart.options.blThermodynamics = 'historical-common-isentrope';
      return { checkpoint: next, diagnostics: { sourceMach: next.restart.input.mach } };
    },
    solveCoupledStreamtubeIses(_input, options) {
      calls.warmReplay++;
      return { converged: true, checkpoint: copy(options.resume), mesh: mesh(), conditions: { mach: options.resume.restart.input.mach, referenceChord: 1 },
        solverInput: copy(options.resume.restart.input), families: copy(options.resume.families), residual: [0], flow: { nodes: copy(nodes) } };
    },
    initializeCoupledStreamtubeFromFlow(mach, checkpoint) {
      calls.reconstructed.push(mach);
      assert.equal(mach, checkpoint.restart.input.mach, 'Source reconstruction must use its actual Mach.');
      return { system: { euler: { conditions: { lengthScale: 1 } }, bl: { transitionMode: 'automatic' } } };
    },
    solveCoupledStreamtubeAutomatic(targetMach, options) {
      calls.automatic.push({ targetMach, checkpoint: copy(options.initialCheckpoint) });
      // End the simulated target attempt retaining the certified colder root.
      const checkpoint = copy(options.initialCheckpoint);
      return { converged: false, stateConverged: true, reason: 'controlled target stop', checkpoint,
        conditions: { mach: checkpoint.restart.input.mach }, solverInput: checkpoint.restart.input,
        mesh: mesh(), families: copy(families), continuation: { boundaryLayerReinitializedOnWarmRestart: false } };
    } };
  const key = `__coldMachRecovery${++sequence}`; globalThis[key] = stubs;
  const source = fs.readFileSync(new URL('../src/euler/tests/streamtube-coupled-mach-assembly.js', import.meta.url), 'utf8')
    .replace(/import \{([^}]+)\} from '[^']+';/g, (_, names) => `const {${names.replace(/\bas\b/g, ':')}} = globalThis[${JSON.stringify(key)}];`);
  const module = await import('data:text/javascript;base64,' + Buffer.from(source + `\n//# sourceURL=${key}`).toString('base64'));
  return { module: { ...module, solveCoupledStreamtubeMach: (input, options) =>
    module.solveCoupledStreamtubeMach(input, { alphaContinuation: false, targetInviscidStartup: false, ...options }) },
    calls, originalPacket, originalError, request, release: () => delete globalThis[key] };
}

test('repeated colder startup reuses the original packet and carries actual source Mach into continuation', async () => {
  const h = await harness(), before = copy(h.request), stages = [], accepted = [], prepared = [];
  try {
    const result = h.module.solveCoupledStreamtubeMach(h.request, { onStage: e => stages.push(e),
      onPrepared: e => prepared.push(e), onCheckpoint: (cp, event) => accepted.push({ cp, event }) });
    assert.deepEqual(h.calls.cold.map(c => c.caseData.mach), [.2, .1, .05]);
    assert.equal(h.calls.cold[0].prepared, undefined);
    assert.equal(h.calls.retarget.length, 2);
    for (const call of h.calls.retarget) assert.equal(call.packet, h.originalPacket);
    for (const call of h.calls.cold) assert.deepEqual({ ...call.caseData, mach: .74 }, before);
    assert.deepEqual(h.calls.reconstructed, [.05]); assert.deepEqual(h.calls.certified, [.05]);
    assert.equal(h.calls.automatic[0].targetMach, .74);
    assert.equal(h.calls.automatic[0].checkpoint.restart.input.mach, .05);
    assert.equal(prepared[0].actualMach, .05); assert.equal(prepared[0].targetMach, .74);
    assert.equal(accepted[0].event.actualMach, .05); assert.equal(accepted[0].event.targetMach, .74);
    assert.deepEqual(stages.filter(s => s.stage === 'euler-cold-startup').map(s => [s.previousMach, s.actualMach, s.targetMach]), [[.2, .1, .74], [.1, .05, .74]]);
    assert.equal(result.actualMach, .05); assert.equal(result.targetMach, .74);
    assert.equal(result.machContinuation.sourceMach, .05); assert.equal(result.requestedCase.mach, .74);
    assert.equal(result.converged, false); assert.equal(result.stateConverged, true);
    assert.deepEqual(result.initialization.coldMachStartup.attempts.map(a => a.mach), [.2, .1]);
    assert.deepEqual(h.request, before);
  } finally { h.release(); }
});

test('capacity failures from later phases, ordinary errors and low target routes never retry colder', async () => {
  for (const config of [{ error: capacity('boundary-layer-initialization') }, { error: capacity('transition-refinement') },
    { error: new Error('untyped failure') }, { target: .2 }]) {
    const h = await harness(config);
    try {
      const result = h.module.solveCoupledStreamtubeMach({ ...h.request, mach: config.target ?? .74 });
      assert.equal(h.calls.cold.length, 1); assert.equal(h.calls.retarget.length, 0);
      assert.equal(h.calls.automatic.length, 0); assert.equal(result.stateConverged, false);
      assert.equal(result.reason, h.originalError.message); assert.equal(result.actualMach, .2);
    } finally { h.release(); }
  }
});

test('secondary static-enthalpy rejection reuses the mesh and preserves explicit ISMOM4 through the cold floor', async () => {
  const error = capacity();
  error.diagnostics.capacityRatio = 3.881;
  error.diagnostics.stagnationDensityFallback = { attempted: true, admissible: false,
    code: 'streamtube-static-enthalpy', reason: 'Body transport chain group=1, tube=1: Nonpositive streamtube static enthalpy.',
    diagnostics: { staticEnthalpy: -1, group: 1, tube: 1 } };
  error.diagnostics.pressureDomainDensityFallback = { attempted: true, admissible: false,
    code: 'streamtube-interface-pressure', reason: 'Controlled complete pressure-domain seed rejection.' };
  for (const targetMach of [.74, .2]) {
    const h = await harness({ error, failures: 99 });
    const request = { ...h.request, mach: targetMach, eulerIsmom: 4, gridTubes: 7, reynolds: 2.51e6 }, before = copy(request);
    try {
      const result = h.module.solveCoupledStreamtubeMach(request);
      const expected = targetMach === .74 ? [.2, .1, .05, .025] : [.2];
      assert.deepEqual(h.calls.cold.map(c => c.caseData.mach), expected);
      for (const call of h.calls.cold) assert.deepEqual({ ...call.caseData, mach: targetMach }, before);
      for (const call of h.calls.retarget) {
        assert.equal(call.packet, h.originalPacket);
        assert.deepEqual({ ...call.caseData, mach: targetMach }, before);
        assert.deepEqual(streamtubeEquationControls(call.caseData.eulerIsmom), streamtubeEquationControls(4));
      }
      assert.equal(h.calls.automatic.length, 0); assert.equal(h.calls.certified.length, 0);
      assert.equal(result.converged, false); assert.equal(result.stateConverged, false); assert.equal(result.availableFlow, false);
      assert.equal(result.actualMach, expected.at(-1)); assert.equal(result.targetMach, targetMach);
      assert.equal(result.requestedCase.eulerIsmom, 4); assert.equal(result.sourceCase.eulerIsmom, 4);
      assert.equal(result.mesh.tag, 'one-original-mesh'); assert.deepEqual(result.failure.diagnostics, error.diagnostics);
      assert.deepEqual(request, before);
    } finally { h.release(); }
  }
});

test('a later unknown or geometric density fallback prevents colder retry despite an earlier typed enthalpy rejection', async () => {
  for (const fallback of [{ reason: 'unclassified failure' },
    { code: 'streamtube-grid-domain', reason: 'Crossed quadrilateral polygon.' },
    { code: 'streamtube-pressure-domain-transport', reason: 'Biased transport invalidates the density proof.' }]) {
    const error = capacity();
    error.diagnostics.stagnationDensityFallback = { code: 'streamtube-static-enthalpy' };
    error.diagnostics.pressureDomainDensityFallback = fallback;
    const h = await harness({ error });
    try {
      const result = h.module.solveCoupledStreamtubeMach({ ...h.request, eulerIsmom: 4 });
      assert.equal(h.calls.cold.length, 1); assert.equal(h.calls.retarget.length, 0);
      assert.equal(result.actualMach, .2); assert.equal(result.converged, false);
      assert.deepEqual(result.failure.diagnostics.pressureDomainDensityFallback, fallback);
    } finally { h.release(); }
  }
});

test('a converged intermediate Ncrit root is retained at its actual condition without Mach continuation', async () => {
  const h = await harness({ failures: 0, mutateColdResult(result) {
    result.checkpoint.restart.options.ncrit = 8;
    result.conditions = { mach: result.mach, ncrit: 8 };
    result.actualNcrit = 8; result.targetNcrit = 9;
    result.ncritContinuation = { actualNcrit: 8, targetNcrit: 9, reachedTarget: false };
    // Deliberately leave converged:true: the adapter must not promote it.
  } });
  try {
    const request = { ...h.request, mach: .185, ncrit: 9 };
    const result = h.module.solveCoupledStreamtubeMach(request);
    assert.equal(result.converged, false); assert.equal(result.stateConverged, true);
    assert.equal(result.actualNcrit, 8); assert.equal(result.targetNcrit, 9);
    assert.equal(result.sourceCase.ncrit, 8); assert.equal(result.requestedCase.ncrit, 9);
    assert.equal(result.checkpoint.restart.options.ncrit, 8);
    assert.equal(result.machContinuation.reachedTarget, true, 'Mach was reached; Ncrit was not.');
    assert.equal(result.ncritContinuation.reachedTarget, false);
    assert.equal(h.calls.cold.length, 1); assert.equal(h.calls.automatic.length, 0);
    assert.equal(h.calls.certified.length, 0); assert.equal(h.calls.warmReplay, 0);
  } finally { h.release(); }
});

test('observer cancellation retains its identity at initial mesh, retry stage and accepted-source preparation', async () => {
  for (const kind of ['stage', 'mesh', 'retry', 'prepared']) {
    const h = await harness(), stop = kind === 'stage' ? 'primitive cancellation' : Object.freeze({ cancelled: kind });
    try {
      const options = kind === 'mesh' ? { onMesh: () => { throw stop; } }
        : kind === 'prepared' ? { onPrepared: () => { throw stop; } }
          : { onStage: e => { if (kind === 'stage' || e.stage === 'euler-cold-startup') throw stop; } };
      assert.throws(() => h.module.solveCoupledStreamtubeMach(h.request, options), e => e === stop);
      assert.equal(h.calls.automatic.length, 0);
      assert.equal(h.calls.cold.length, kind === 'prepared' ? 3 : 1);
    } finally { h.release(); }
  }
});

test('failed retarget preserves the rejected startup mesh and both causes in a structured failure', async () => {
  const h = await harness({ retargetError: new Error('controlled retarget binding rejection') });
  try {
    const result = h.module.solveCoupledStreamtubeMach(h.request);
    assert.equal(result.converged, false); assert.equal(result.stateConverged, false); assert.equal(result.availableFlow, false);
    assert.equal(result.actualMach, .2); assert.equal(result.targetMach, .74); assert.equal(result.requestedCase.mach, .74);
    assert.equal(result.mesh.tag, 'one-original-mesh');
    const text = JSON.stringify(result);
    assert.ok(text.includes(h.originalError.message)); assert.ok(text.includes('controlled retarget binding rejection'));
    assert.deepEqual(result.initialization.coldMachStartup.attempts.map(a => a.mach), [.2]);
    assert.equal(h.calls.cold.length, 1); assert.equal(h.calls.automatic.length, 0);
  } finally { h.release(); }
});

test('certification recovery records the actual colder baseline and rejects simultaneous warm/cold inputs', async () => {
  const h = await harness({ warmRecovery: true }), parent = acceptedSource(.25), before = copy(parent);
  try {
    const result = h.module.solveCoupledStreamtubeMach(h.request, { parentResult: parent });
    assert.equal(h.calls.warmReplay, 1); assert.deepEqual(h.calls.certified, [.25, .05]);
    assert.equal(result.machContinuation.sourceMach, .05);
    assert.equal(result.machContinuation.coldRecovery.baselineMach, .05);
    assert.equal(result.targetMach, .74); assert.deepEqual(parent, before);
    const count = h.calls.cold.length;
    assert.throws(() => h.module.solveCoupledStreamtubeMach(h.request, { parentResult: parent, preparedEuler: h.originalPacket }), /cannot also supply/);
    assert.equal(h.calls.cold.length, count);
  } finally { h.release(); }
});

// Run the actual Worker message handler against a supplied controller result.
// Dynamic numerical imports are stubs; no native Worker is created.
async function workerError(result, request, thrown) {
  const messages = [], self = { postMessage: data => messages.push(copy(data)) };
  const key = `__coldMachErrorWorker${++sequence}`;
  globalThis[key] = { self, quadCoupledNcrit, quadCoupledNcritResult, modules: {
    '../euler/streamtube-coupled-mach-assembly.js': { solveCoupledStreamtubeMach() { if (thrown) throw thrown; return result; } },
    // The Worker now starts the requested direct coupled solve through the
    // assembly entry point. Keep this Worker-only error harness aligned with
    // that public route; both stubs deliberately return the supplied result.
    '../euler/streamtube-coupled-assembly.js': { solveCoupledStreamtubeAssembly() { if (thrown) throw thrown; return result; } },
    '../ui/quad-coupled-result.js': { quadCoupledResultForDisplay: value => value }, '../ui/quad-coupled-transonic-coefficients.js': {},
    '../ui/quad-mesh-progress.js': { createQuadMeshProgress: () => ({ stage() {}, mesh: value => value }) },
  } };
  const source = fs.readFileSync(new URL('../src/worker/solver.js', import.meta.url), 'utf8')
    .replace(/import \{([^}]+)\} from '[^']+';/g, (_, names) => `const {${names.replace(/\bas\b/g, ':')}} = globalThis[${JSON.stringify(key)}];`)
    .replace(/await import\('([^']+)'\)/g, (_, path) => `globalThis[${JSON.stringify(key)}].modules[${JSON.stringify(path)}]`);
  try {
    await import('data:text/javascript;base64,' + Buffer.from(`const self = globalThis[${JSON.stringify(key)}].self;\n${source}`).toString('base64'));
    await self.onmessage({ data: { id: 41, task: 'continue-coupled', caseData: request } });
    assert.equal(messages.length, 1);
    assert.equal(messages[0].type, 'error');
    return messages[0];
  } finally { delete globalThis[key]; }
}

test('terminal gas diagnostics and compact cause survive controller and Worker without a retry or at the cold floor', async () => {
  for (const target of [.2, .74]) {
    const cause = Object.assign(new Error('nonpositive interface pressure'), {
      code: 'streamtube-interface-pressure', stage: 'gas-initialization',
      diagnostics: { interfacePressure: -.002, cell: { i: 128, group: 1, tube: 1 } },
      unrelatedLargeData: { shouldNotBeCopied: true },
    });
    const error = Object.assign(capacity(), { cause });
    error.diagnostics.terminalMarker = 'terminal gas diagnostic';
    const h = await harness({ error, failures: 99 }), request = { ...h.request, mach: target };
    try {
      const result = h.module.solveCoupledStreamtubeMach(request);
      assert.deepEqual(h.calls.cold.map(c => c.caseData.mach), target === .2 ? [.2] : [.2, .1, .05, .025]);
      assert.deepEqual(result.failure.diagnostics, error.diagnostics);
      assert.notEqual(result.failure.diagnostics, error.diagnostics);
      assert.deepEqual(result.failure.cause, { message: cause.message, code: cause.code, stage: cause.stage, diagnostics: cause.diagnostics });
      const message = await workerError(result, request);
      assert.equal(message.message, error.message);
      assert.equal(message.code, error.code);
      assert.equal(message.stage, 'gas-initialization');
      assert.equal(message.actualMach, target === .2 ? .2 : .025);
      assert.equal(message.targetMach, target);
      assert.equal(message.diagnostics.terminalMarker, 'terminal gas diagnostic');
      assert.equal(message.diagnostics.capacityRatio, 1.104);
      assert.deepEqual(message.cause, result.failure.cause);
      assert.equal(message.cause.diagnostics.cell.i, 128);
      assert.equal(Object.hasOwn(message.cause, 'unrelatedLargeData'), false);
      if (target === .74) assert.deepEqual(message.diagnostics.coldMachStartup.attempts.map(a => a.mach), [.2, .1, .05]);
      error.diagnostics.terminalMarker = 'changed after return';
      cause.diagnostics.interfacePressure = 0;
      assert.equal(result.failure.diagnostics.terminalMarker, 'terminal gas diagnostic');
      assert.equal(result.failure.cause.diagnostics.interfacePressure, -.002);
    } finally { h.release(); }
  }
});

test('both Worker precursor failure paths read the actual diagnostic stop reason and rejected pressure location', async () => {
  const precursor = { mach: .1, diagnostics: { reason: 'Euler grid report', solverStopReason: 'update rejected' },
    flow: { reason: 'fallback flow reason', lastRejectedStep: { stage: 'admissibility',
      message: 'Nonpositive interface pressure at cell i=128, group=1, tube=1.' } }, mesh: { quality: { valid: true } } };
  const error = Object.assign(new Error('Coupled Euler precursor did not converge'), {
    code: 'coupled-euler-precursor', stage: 'euler', initialization: { stage: 'euler', result: precursor },
  });
  const returned = await workerError({ reason: error.message, actualMach: .1, targetMach: .74,
    failure: { code: error.code, stage: error.stage }, initialization: error.initialization }, caseData());
  assert.equal(returned.precursor.solverStopReason, 'update rejected');
  assert.deepEqual(returned.precursor.lastRejectedStep, precursor.flow.lastRejectedStep);
  const caught = await workerError(undefined, caseData(), error);
  assert.equal(caught.actualMach, .1);
  assert.equal(caught.targetMach, .74);
  assert.equal(caught.precursor.solverStopReason, 'update rejected');
  assert.deepEqual(caught.precursor.lastRejectedStep, precursor.flow.lastRejectedStep);
});

test('Euler ISES keeps typed pressure diagnostics on a rejected update without changing accepted state', async () => {
  const error = Object.assign(new Error('Nonpositive interface pressure.'), {
    code: 'streamtube-interface-pressure', diagnostics: { interfacePressure: -.002, cell: { i: 128, group: 1, tube: 1 } },
  });
  const grid = Array.from({ length: 2 }, (_, g) => Array.from({ length: 2 }, (_, i) => [{ x: i, y: g }, { x: i, y: g + 1 }]));
  const input = { bodies: [{ leadingIndex: 0 }] }, state = Float64Array.of(0);
  let systems = 0, proposals = 0;
  const value = () => ({ nodes: copy(grid), stagnation: [.5], residual: Float64Array.of(1), diagnostics: { residual: 1 } });
  const stubs = {
    assertConvexStreamtubeGrid, requireConvexGridUpdate, requireConvexPublishedGrid,
    createStreamtubeBodySystem() {
      const index = systems++;
      return { initial: state.slice(), layout: fixtureLayout(grid, input.bodies),
        adoptGeometry: x => x.slice(), decode: value, jacobian: () => ({}),
        evaluate() { if (index >= 2) throw error; return value(); } };
    },
    solveStreamtubeBody: (_system, options) => ({ ...value(), x: options.initial.slice(), converged: false }),
    streamtubeMeshSnapshot: () => ({ quality: { valid: true } }),
    captureStreamtubeInletFractions: () => [[0, 1]],
    redistributeStreamtubeTangentially: (nodes, options) => ({ nodes, referenceBank: options.referenceBank,
      fixedBanks: options.fixedBanks, solution: { pairs: 5, relativeResidual: 0 }, maxDisplacement: 0 }),
    adjustStreamtubeInlets: nodes => ({ nodes, maxDisplacement: 0, maxArcErrorBefore: 0, maxArcErrorAfter: 0 }),
    dekinkStreamtubeInteriors: nodes => ({ nodes, repairs: [] }),
    solveSparseDirect: () => ({ x: Float64Array.of(.1), relativeResidual: 0, refinements: 0 }),
    proposeDensityNewton: () => { proposals++; return { x: Float64Array.of(.1), step: 1 }; },
  };
  const key = `__isesErrorPreservation${++sequence}`; globalThis[key] = stubs;
  const source = fs.readFileSync(new URL('../src/euler/streamtube-ises-update.js', import.meta.url), 'utf8')
    .replace(/import \{([^}]+)\} from '[^']+';/g, (_, names) => `const {${names.replace(/\bas\b/g, ':')}} = globalThis[${JSON.stringify(key)}];`);
  try {
    const { solveStreamtubeIses } = await import('data:text/javascript;base64,' + Buffer.from(source).toString('base64'));
    const result = solveStreamtubeIses(input, { maxIterations: 1, maxBacktracks: 0, stepAcceptance: 'admissible' });
    const failed = result.lastRejectedStep;
    assert.equal(failed.code, error.code);
    assert.deepEqual(failed.diagnostics, error.diagnostics);
    assert.equal(result.converged, false);
    assert.deepEqual(result.x, state);
    assert.equal(result.history.length, 1);
    assert.equal(proposals, 1);
    assert.equal(failed.message, error.message);
    error.diagnostics.interfacePressure = 1;
    assert.equal(failed.diagnostics.interfacePressure, -.002);
  } finally { delete globalThis[key]; }
});

test('initial SMOVE pressure failures retain the original state and reach both Worker precursor paths before any Newton update', async () => {
  for (const stepAcceptance of ['listing', 'admissible']) {
    const failure = Object.assign(new Error('Body cell i=124, group=0, tube=2: Nonpositive interface pressure.'), {
      code: 'streamtube-interface-pressure', diagnostics: { interfacePressure: { lower: -.2, upper: Infinity },
        pressureSum: Infinity, pressureDifference: NaN, normalInertia: -4, cell: { i: 124, group: 0, tube: 2 } },
    });
    const expected = copy(failure.diagnostics), state = Float64Array.of(.137);
    const grid = Array.from({ length: 2 }, (_, g) => Array.from({ length: 2 }, (_, i) => [{ x: i, y: g }, { x: i, y: g + 1 }]));
    const input = { bodies: [{ leadingIndex: 0 }] }, before = copy(input);
    let systems = 0, evaluations = 0;
    const value = () => ({ nodes: copy(grid), stagnation: [.5], residual: Float64Array.of(1), diagnostics: { residual: 1 } });
    const stubs = {
      assertConvexStreamtubeGrid, requireConvexGridUpdate, requireConvexPublishedGrid,
      createStreamtubeBodySystem() {
        const index = systems++;
        return { initial: state.slice(), layout: fixtureLayout(grid, input.bodies),
          adoptGeometry: x => x.slice(), decode: value,
          jacobian: () => assert.fail('Initialization failure must not form a Jacobian.'),
          evaluate() { evaluations++; if (index > 0) throw failure; return value(); } };
      },
      solveStreamtubeBody: (_system, options) => {
        assert.equal(options.maxIterations, 0);
        return { ...value(), x: options.initial.slice(), converged: false };
      },
      streamtubeMeshSnapshot: () => ({ quality: { valid: true } }),
      captureStreamtubeInletFractions: () => [[0, 1]],
      redistributeStreamtubeTangentially: (nodes, options) => ({ nodes, referenceBank: options.referenceBank,
        fixedBanks: options.fixedBanks, solution: { pairs: 5, relativeResidual: 0 }, maxDisplacement: 0 }),
      solveSparseDirect: () => assert.fail('Initialization failure must not solve a linear system.'),
      proposeDensityNewton: () => assert.fail('Initialization failure must not propose a Newton update.'),
    };
    const key = `__isesInitialErrorPreservation${++sequence}`; globalThis[key] = stubs;
    const source = fs.readFileSync(new URL('../src/euler/streamtube-ises-update.js', import.meta.url), 'utf8')
      .replace(/import \{([^}]+)\} from '[^']+';/g, (_, names) => `const {${names.replace(/\bas\b/g, ':')}} = globalThis[${JSON.stringify(key)}];`);
    try {
      const { solveStreamtubeIses } = await import('data:text/javascript;base64,' + Buffer.from(source).toString('base64'));
      const result = solveStreamtubeIses(input, { maxIterations: 1, maxBacktracks: 2, stepAcceptance,
        onCheckpoint: () => assert.fail('Rejected initialization cannot publish a checkpoint.') });
      assert.equal(result.initialRedistribution.accepted, false);
      assert.equal(result.initialRedistribution.rejection, failure.message);
      assert.equal(result.reason, `ISES initial redistribution rejected: ${failure.message}`);
      assert.deepEqual(result.lastRejectedStep, { stage: 'initial redistribution', message: failure.message,
        code: failure.code, diagnostics: expected });
      assert.equal(result.converged, false); assert.deepEqual(result.x, state); assert.deepEqual(result.nodes, grid);
      assert.equal(result.history.length, 1); assert.equal(result.linearDiagnostics.solves, 0);
      assert.equal(systems, stepAcceptance === 'listing' ? 2 : 4); assert.equal(evaluations, systems);
      assert.deepEqual(input, before);
      const precursor = { mach: .2, diagnostics: { reason: result.reason, solverStopReason: result.reason },
        flow: result, mesh: { quality: { valid: true } } };
      const error = Object.assign(new Error(`Coupled Euler precursor did not converge: ${result.reason}`), {
        code: 'coupled-euler-precursor', stage: 'euler', initialization: { stage: 'euler', result: precursor },
      });
      const returned = await workerError({ reason: error.message, actualMach: .2, targetMach: .74,
        failure: { code: error.code, stage: error.stage }, initialization: error.initialization }, caseData());
      const caught = await workerError(undefined, caseData(), error);
      for (const packet of [returned.precursor, caught.precursor]) {
        assert.equal(packet.solverStopReason, result.reason);
        assert.deepEqual(packet.lastRejectedStep, result.lastRejectedStep);
      }
      failure.diagnostics.cell.i = 999;
      assert.deepEqual(result.lastRejectedStep.diagnostics, expected);
    } finally { delete globalThis[key]; }
  }
});
