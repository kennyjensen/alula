// SPDX-License-Identifier: GPL-2.0-or-later
// Actual ISES controller with synthetic numerical boundaries; no flow/Jacobian/LU.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createIterationProgress, physicalIterationVector, residualMerit } from '../src/euler/streamtube-iteration-progress.js';

let serial = 0;
async function harness({ terminalFailure = false } = {}) {
  const nodes = Array.from({ length: 2 }, (_, g) => Array.from({ length: 3 }, (_, i) => [{ x: i, y: g }, { x: i, y: g + 1 }]));
  const outer = { nodes, undisplacedNodes: nodes, stagnation: [0], diagnostics: { maxMach: .2 } };
  const families = { euler: .1, boundaryLayer: 0, edgeMatching: 0 }, calls = { constructs: 0, helper: [], direct: [] };
  let terminalFailures = terminalFailure ? 1 : 0;
  const linear = () => ({ x: new Float64Array(5), relativeResidual: 0, refinements: 0,
    ordering: 'colamd', pivotTolerance: 1, attempts: [{ ordering: 'colamd', pivotTolerance: 1 }] });
  const bindings = {
    createIterationProgress, physicalIterationVector, residualMerit,
    createCoupledStreamtubeBody() {
      const ordinal = ++calls.constructs;
      return { ordinal, n: 5, ne: 1, initial: new Float64Array(5),
        euler: { layout: { densityCount: 1, bodies: [{ leadingIndex: 1 }], rows: [{ kind: 'streamwise' }] },
          conditions: { pressureScale: 1 }, decode: () => outer,
          setDisplacement() {}, adoptGeometry: x => x.slice() },
        bl: { stations: [], thicknesses: () => ({}), snapshotActive: () => [] },
        jacobian: () => ({ ordinal }), admissible: () => true,
        admissibleValue(x, controls) { return this.admissible(x, controls) ? this.evaluate(x) : null; },
        evaluate: () => ({ residual: Float64Array.of(.1, 0, 0, 0, 0), families, outer }) };
    },
    coupledStreamtubeTripEvents: () => ({ snapshot: () => [], restore() {}, prepare: () => ({ changed: false }) }),
    coupledStreamtubeResult: (system, x, details) => ({ x, ...details, families, lastSystem: system.ordinal }),
    proposeCoupledDensityNewton: (_, x) => ({ x: x.slice(), step: 1 }),
    proposeCoupledXfoilBLUpdate: (_, x) => ({ x: x.slice(), step: 1 }),
    solveCoupledLinearSystem: (matrix, rhs, options) => {
      calls.helper.push({ ordinal: matrix.ordinal, rhs: Array.from(rhs), ...options });
      if (terminalFailures-- > 0) throw Object.assign(new Error('Controlled station and automatic linear rejection.'), {
        code: 'KLU_RESIDUAL_LIMIT', diagnostics: { stationPolicy: {
          mode: 'station-auto', recommendation: 'auto', attempted: true, fallback: true,
          autoFailure: { code: 'KLU_RESIDUAL_LIMIT' } } } });
      return { ordering: null, linear: linear(), stationPolicy: {
        mode: 'station-auto', recommendation: 'auto', attempted: options.stationEnabled,
        fallback: options.stationEnabled, selected: { ordering: 'colamd', pivotTolerance: 1 } } };
    },
    solveSparseDirect: (_, __, options) => { calls.direct.push(options); return linear(); },
    createStreamtubeStationOrdering: () => ({ Puser: [0, 1, 2, 3, 4], Quser: [0, 1, 2, 3, 4],
      diagnostics: { n: 5, borderPairCount: 0, maximumLocalStationMismatch: 0 } }),
    redistributeStreamtubeTangentially: group => ({ nodes: group, solution: { pairs: 5, relativeResidual: 0 }, maxDisplacement: 0 }),
    captureStreamtubeInletFractions: () => [[0, 1]],
    adjustStreamtubeInlets: nodes => ({ nodes }), dekinkStreamtubeInteriors: nodes => ({ nodes, repairs: [] }),
    assertConvexStreamtubeGrid() {}, prepareConvexWakeGridUpdate: (_, nodes) => ({ nodes }),
    respondToCoupledProjectionGeometry: () => assert.fail('No projection response requested.'),
    requireCoupledResidualDecrease: () => assert.fail('No Armijo policy requested.'),
  };
  const key = `__stationContinuation${serial++}`; globalThis[key] = bindings;
  const source = fs.readFileSync(new URL('../src/euler/streamtube-coupled-ises.js', import.meta.url), 'utf8')
    .replace(/import \{([^}]+)\} from '[^']+';/g, (_, names) => `const {${names.replace(/\s+as\s+/g, ': ')}}=globalThis[${JSON.stringify(key)}];`);
  try {
    const module = await import('data:text/javascript;base64,' + Buffer.from(source).toString('base64'));
    return { ...module, calls, input: { bodies: [{ leadingIndex: 1 }] } };
  } finally { delete globalThis[key]; }
}

test('fresh ordering fallback persists across accepted chart rebases and chunk resume, while a new grid starts fresh', async () => {
  const h = await harness(), checkpoints = [];
  const result = h.solveCoupledStreamtubeIses(h.input, { maxIterations: 2, onCheckpoint: cp => checkpoints.push(structuredClone(cp)) });
  assert.equal(result.linearOrdering, 'station-auto'); assert.equal(result.history.length, 3);
  assert.deepEqual(h.calls.helper.map(c => c.stationEnabled), [true, false]);
  assert.ok(h.calls.helper[1].ordinal > h.calls.helper[0].ordinal, 'The accepted step creates a new numerical chart.');
  assert.equal(h.calls.direct.length, 0);
  assert.equal(Object.hasOwn(checkpoints[0].continuation, 'stationFallback'), false);
  assert.ok(checkpoints.slice(1).every(cp => cp.continuation.stationFallback === true));
  const cp = result.checkpoint, before = structuredClone(cp);
  const resumed = h.solveCoupledStreamtubeIses(undefined, { resume: cp, maxIterations: 1 });
  assert.equal(h.calls.helper.at(-1).stationEnabled, false);
  assert.equal(resumed.checkpoint.continuation.stationFallback, true);
  assert.equal(resumed.checkpoint.continuation.linearOrdering, 'station-auto');
  assert.deepEqual(cp, before);
  h.solveCoupledStreamtubeIses(h.input, { maxIterations: 1 });
  assert.equal(h.calls.helper.at(-1).stationEnabled, true);
});

test('terminal dual-linear rejection retains the fallback marker in the unchanged accepted checkpoint and the next resume skips station', async () => {
  const h = await harness({ terminalFailure: true });
  const stopped = h.solveCoupledStreamtubeIses(h.input, { maxIterations: 1 });
  assert.equal(stopped.history.length, 1, 'No Newton candidate was accepted after linear failure.');
  assert.equal(stopped.linearDiagnostics.solves, 0);
  assert.equal(stopped.lastRejectedStep.stage, 'linear solve');
  assert.equal(stopped.lastRejectedStep.code, 'KLU_RESIDUAL_LIMIT');
  assert.equal(stopped.lastRejectedStep.diagnostics.stationPolicy.recommendation, 'auto');
  assert.equal(stopped.checkpoint.continuation.linearOrdering, 'station-auto');
  assert.equal(stopped.checkpoint.continuation.stationFallback, true);
  assert.deepEqual(Array.from(stopped.x), [0, 0, 0, 0, 0]);
  const before = structuredClone(stopped.checkpoint);
  const resumed = h.solveCoupledStreamtubeIses(undefined, { resume: stopped.checkpoint, maxIterations: 1 });
  assert.deepEqual(h.calls.helper.map(c => c.stationEnabled), [true, false]);
  assert.equal(resumed.history.length, 2); assert.equal(resumed.checkpoint.continuation.stationFallback, true);
  assert.deepEqual(stopped.checkpoint, before);
});

test('legacy omitted and explicit station resumes keep their original dispatch and reject malformed fallback before construction', async () => {
  const h = await harness(), zero = h.solveCoupledStreamtubeIses(h.input, { maxIterations: 0, linearOrdering: 'auto' });
  const old = structuredClone(zero.checkpoint), station = structuredClone(old);
  station.continuation.linearOrdering = 'station';
  const ordinary = h.solveCoupledStreamtubeIses(undefined, { resume: old, maxIterations: 1 });
  assert.equal(Object.hasOwn(ordinary.checkpoint.continuation, 'linearOrdering'), false);
  assert.equal(Object.hasOwn(ordinary.checkpoint.continuation, 'stationFallback'), false);
  assert.equal(h.calls.direct.at(-1).preferredOrdering, old.continuation.preferredOrdering);
  const strict = h.solveCoupledStreamtubeIses(undefined, { resume: station, maxIterations: 1 });
  assert.equal(strict.checkpoint.continuation.linearOrdering, 'station');
  assert.equal(h.calls.direct.at(-1).ordering, 'given'); assert.equal(h.calls.direct.at(-1).pivotFallback, false);
  assert.equal(h.calls.helper.length, 0);
  const constructs = h.calls.constructs;
  for (const [policy, fallback] of [['station-auto', null], ['station-auto', 'true'], ['station-auto', 1],
    ['station', true], ['auto', true], [undefined, true]]) {
    const bad = structuredClone(old);
    if (policy !== undefined) bad.continuation.linearOrdering = policy;
    bad.continuation.stationFallback = fallback;
    assert.throws(() => h.solveCoupledStreamtubeIses(undefined, { resume: bad, maxIterations: 0 }), /Invalid coupled ISES iteration controls/);
  }
  assert.equal(h.calls.constructs, constructs); assert.equal(h.calls.helper.length, 0);
});
