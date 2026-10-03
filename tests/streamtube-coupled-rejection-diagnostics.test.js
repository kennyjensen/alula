// SPDX-License-Identifier: GPL-2.0-or-later
// Execute the driver state machine with controlled flow boundaries and the
// real cheap convex-grid predicates; no flow solve or Jacobian is evaluated.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { assertConvexStreamtubeGrid } from '../src/geometry/streamtube-convex-step.js';
import { prepareConvexWakeGridUpdate, requireConvexGridUpdate } from '../src/euler/streamtube-grid-update.js';
import { createIterationProgress, physicalIterationVector, residualMerit } from '../src/euler/streamtube-iteration-progress.js';

const currentSource = fs.readFileSync(new URL('../src/euler/streamtube-coupled-ises.js', import.meta.url), 'utf8');
const previousSource = fs.readFileSync(new URL('../docs/surface-pchip-refinement/mach-continuation/before/streamtube-coupled-ises.js.boundary-diagnostics.txt', import.meta.url), 'utf8');
let sequence = 0;
async function run(source, failure) {
  const key = `__coupledRejectionMetadata${++sequence}`;
  const grid = Array.from({ length: 2 }, (_, g) => Array.from({ length: 4 }, (_, i) => [{ x: i, y: g }, { x: i, y: g + 1 }]));
  const bodies = [{ leadingIndex: 1, trailingIndex: 2 }], x = Float64Array.of(0, .03, 1, 2, 1);
  const input = { bodies }, families = { euler: 1, boundaryLayer: 0, edgeMatching: 0 };
  const continuation = { fractions: [[0, 1]], lastRedistributedStagnation: [.5], preferredOrdering: 'amd', pivotTolerance: .001,
    iterationGeometry: 'convex', stepAcceptance: 'admissible', stagnationLimiter: 'listing' };
  const resume = { version: 1, families, continuation, restart: { input,
    options: { transitionMode: 'automatic', transitionState: [0] }, initialEuler: { x: x.slice(0, 1), nodes: grid,
      undisplacedNodes: grid }, initialBL: x.slice(1) } };
  let phase = [0], constructors = 0;
  const calls = { matrixStubs: 0, linearStubs: 0, proposalStubs: 0 };
  const decoded = () => ({ nodes: structuredClone(grid), undisplacedNodes: structuredClone(grid), stagnation: [.5],
    diagnostics: { maxMach: .2 } });
  const stubs = {
    assertConvexStreamtubeGrid, prepareConvexWakeGridUpdate, requireConvexGridUpdate,
    createIterationProgress, physicalIterationVector, residualMerit,
    createCoupledStreamtubeBody(_input, options) {
      const index = constructors++;
      return { n: 5, ne: 1, initial: Float64Array.from([...options.initialEuler.x, ...options.initialBL]),
        euler: { layout: { densityCount: 1, bodies }, setDisplacement() {}, decode: decoded, adoptGeometry: state => state.slice() },
        bl: { stations: [{ id: 0 }], thicknesses: state => state.slice(), snapshotActive: () => phase.slice() },
        evaluate: state => ({ residual: Float64Array.of(1, 0, 0, 0, 0), families: { ...families }, outer: decoded(), state: state.slice() }),
        admissible(_state, { onFailure } = {}) {
          if (index === 0 || !failure) return true;
          onFailure?.({ kind: 'evaluation-error', stage: 'evaluation', message: failure.message,
            ...(failure.code === undefined ? {} : { code: failure.code }),
            ...(failure.diagnostics === undefined ? {} : { diagnostics: structuredClone(failure.diagnostics) }) });
          return false;
        },
        admissibleValue(state, controls) { return this.admissible(state, controls) ? this.evaluate(state) : null; },
        jacobian() { calls.matrixStubs++; return {}; },
      };
    },
    coupledStreamtubeTripEvents: () => ({ snapshot: () => phase.slice(), restore: saved => { phase = saved.slice(); },
      prepare() { phase = [1]; return { changed: false }; } }),
    coupledStreamtubeResult: (_system, state, options) => ({ ...options, x: Array.from(state), phase: phase.slice() }),
    solveSparseDirect() { calls.linearStubs++; return { x: Float64Array.of(.1, 0, 0, 0, 0), relativeResidual: 0,
      refinements: 0, ordering: 'amd', pivotTolerance: .001, attempts: [] }; },
    proposeCoupledDensityNewton(_system, state, direction, { maximumStep }) {
      calls.proposalStubs++; return { x: state.map((v, i) => v + maximumStep * direction[i]), step: maximumStep };
    },
    adjustStreamtubeInlets: nodes => ({ nodes, maxDisplacement: 0, maxArcErrorBefore: 0, maxArcErrorAfter: 0 }),
    dekinkStreamtubeInteriors: nodes => ({ nodes, repairs: [] }),
  };
  globalThis[key] = stubs;
  const patched = source.replace(/import \{([^}]+)\} from '[^']+';/g,
    (_, names) => `const {${names.replace(/\s+as\s+/g, ': ')}} = globalThis[${JSON.stringify(key)}];`);
  try {
    const { solveCoupledStreamtubeIses } = await import(`data:text/javascript;base64,${Buffer.from(patched).toString('base64')}`);
    const result = solveCoupledStreamtubeIses(undefined, { resume, maxIterations: 1, maxBacktracks: 1,
      iterationGeometry: 'convex', stepAcceptance: 'admissible', stagnationLimiter: 'listing', iterationRecovery: false });
    if (source === currentSource) for (const step of result.history) {
      assert.equal(step.shearCoordinate, 'linear');
      assert.equal(step.hkFloorLinearization, 'exact');
    }
    return { result, calls };
  } finally { delete globalThis[key]; }
}
const withoutNewMetadata = value => JSON.parse(JSON.stringify(value, (key, entry) =>
  // The archived driver predates these method labels; validate them above
  // while retaining exact comparison of every numerical state and retry.
  ['code', 'diagnostics', 'gridAcceptance', 'shearCoordinate', 'hkFloorLinearization'].includes(key) ? undefined : entry));

test('coupled retry and final rejection preserve actual boundary diagnostics and accepted state', async () => {
  const failure = { message: 'Body inlet/outlet section must remain subsonic: i=3, group=1, tube=0.',
    code: 'streamtube-boundary-subsonic', diagnostics: { i: 3, group: 1, tube: 0, boundary: 'outlet',
      machSquared: 1.21, rho: .4, p: 1.1, q: 2, enthalpy: 4,
      machSquaredUpperBound: 1, subsonicMargin: -.21 } };
  const current = await run(currentSource, failure), old = await run(previousSource, failure);
  assert.deepEqual(withoutNewMetadata(current), withoutNewMetadata(old));
  assert.equal(current.result.gridAcceptance, 'convex');
  const rejected = current.result.lastRejectedStep;
  assert.equal(rejected.code, failure.code); assert.deepEqual(rejected.diagnostics, failure.diagnostics);
  assert.equal(rejected.rejections.length, 2);
  for (const trial of rejected.rejections) {
    assert.equal(trial.code, failure.code); assert.deepEqual(trial.diagnostics, failure.diagnostics);
    assert.deepEqual(trial.admissibility.diagnostics, failure.diagnostics);
  }
  assert.deepEqual(current.result.x, [0, .03, 1, 2, 1]); assert.deepEqual(current.result.phase, [0]);
  assert.equal(current.result.history.length, 1);
  assert.deepEqual(current.calls, { matrixStubs: 1, linearStubs: 1, proposalStubs: 2 });
  failure.diagnostics.machSquared = 99;
  assert.equal(rejected.diagnostics.machSquared, 1.21);
  for (const trial of rejected.rejections) assert.equal(trial.diagnostics.machSquared, 1.21);
});

test('successful and untyped-rejection driver outputs retain archived default behavior and shape', async () => {
  for (const failure of [null, { message: 'controlled untyped domain rejection' }]) {
    const current = await run(currentSource, failure), previous = await run(previousSource, failure);
    assert.equal(current.result.gridAcceptance, 'convex');
    assert.deepEqual(withoutNewMetadata(current), withoutNewMetadata(previous));
  }
});
