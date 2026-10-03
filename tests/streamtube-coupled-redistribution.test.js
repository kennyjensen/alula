// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createCoupledStreamtubeBody, coupledStreamtubeTripEvents, coupledStreamtubeResult } from '../src/euler/streamtube-coupled.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { redistributeStreamtubeTangentially } from '../src/geometry/streamtube-tangential-redistribution.js';
import { prepareConvexWakeGridUpdate } from '../src/euler/streamtube-grid-update.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const plain = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const fixture = () => ({ input: { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
  wakeGeometry: 'independent-banks', wakeOutlet: 'banks' },
options: { edgeMatching: 'section-velocity', transitionMode: 'automatic', reynolds: 1e6, ncrit: 9,
  // Controlled directions below replace the ordinary direct solve. Keep
  // station-ordering and iteration-recovery machinery outside this harness.
  linearOrdering: 'auto', iterationRecovery: false } });

// Retain real geometry, coupled gas/BL evaluation and active transitions.
// Inject candidate-domain failures only AFTER its ordinary domain check.
// Trigger tests supply a controlled direction: no global Jacobian or LU.
async function driver({ rejectedCandidates = 0, failMovement = false, direction, rejectedRawProposals = 0 } = {}) {
  const url = new URL('../src/euler/streamtube-coupled-ises.js', import.meta.url);
  const key = `coupled-redistribution-test-${Math.random()}`;
  const calls = [], systems = [], supplied = [], rejected = [], linearCalls = [], gridBounds = [];
  globalThis[key] = {
    prepareConvexWakeGridUpdate(previous, proposed, layout, allocation, minimumUncorrectedFraction) {
      gridBounds.push(minimumUncorrectedFraction);
      const result = prepareConvexWakeGridUpdate(previous, proposed, layout, allocation, minimumUncorrectedFraction);
      if (rejectedRawProposals-- > 0) throw Object.assign(new Error('Controlled raw coupled grid rejection.'), {
        code: 'streamtube-grid-step', stepFraction: .5,
      });
      return result;
    },
    createCoupledStreamtubeBody(input, options) {
      supplied.push(plain({ input, options }));
      const system = createCoupledStreamtubeBody(input, options), index = systems.length, admissibleValue = system.admissibleValue;
      systems.push(system);
      system.jacobian = x => ({ controlledDirection: true, state: Array.from(x) });
      system.admissibleValue = (x, controls = {}) => {
        const value = admissibleValue(x, controls);
        if (!value) return false;
        if (index > 0 && index <= rejectedCandidates) {
          rejected.push({ index, phase: system.bl.snapshotActive(), state: Array.from(x) });
          controls.onFailure?.({ kind: 'controlled-candidate-domain', stage: 'gas', message: 'Controlled coupled pressure rejection.' });
          return false;
        }
        return value;
      };
      return system;
    },
    coupledStreamtubeTripEvents,
    coupledStreamtubeResult,
    solveSparseDirect(matrix) {
      assert.equal(matrix.controlledDirection, true);
      assert.ok(direction, 'This test must not enter Newton without a controlled direction.');
      linearCalls.push(matrix.state);
      return { x: direction.slice(), relativeResidual: 0, refinements: 0, ordering: 'amd', pivotTolerance: .001, attempts: [] };
    },
    redistributeStreamtubeTangentially(nodes, options) {
      calls.push({ nodes: structuredClone(nodes), options: structuredClone(options) });
      if (failMovement && options.correctionScale > .5)
        throw new Error('Controlled nonpositive coordinate span.');
      return redistributeStreamtubeTangentially(nodes, options);
    },
  };
  let text = fs.readFileSync(url, 'utf8');
  for (const statement of [
    "import { createCoupledStreamtubeBody, coupledStreamtubeTripEvents, coupledStreamtubeResult } from './streamtube-coupled.js';",
    "import { solveSparseDirect } from '../numerics/klu.js';",
    "import { redistributeStreamtubeTangentially } from '../geometry/streamtube-tangential-redistribution.js';",
    "import { prepareConvexWakeGridUpdate } from './streamtube-grid-update.js';",
  ]) {
    assert.ok(text.includes(statement));
    const names = statement.slice(statement.indexOf('{'), statement.indexOf('}') + 1);
    text = text.replace(statement, `const ${names} = globalThis[${JSON.stringify(key)}];`);
  }
  text = text.replace(/from '(\.[^']+)'/g, (_, relative) => `from '${new URL(relative, url).href}'`);
  try {
    const { solveCoupledStreamtubeIses: solve } = await import(`data:text/javascript;base64,${Buffer.from(text).toString('base64')}`);
    return { solve, calls, systems, supplied, rejected, linearCalls, gridBounds };
  } finally { delete globalThis[key]; }
}

test('initial coupled SMOVE retries the same complete Euler/BL state across every passage', async () => {
  const f = fixture(), before = plain(f), d = await driver({ rejectedCandidates: 2 }), checkpoints = [];
  const result = d.solve(f.input, { ...f.options, maxIterations: 0, stepAcceptance: 'admissible', maxBacktracks: 2,
    onCheckpoint: checkpoint => checkpoints.push(plain(checkpoint)) });
  assert.equal(result.initialRedistribution.accepted, true, result.reason);
  assert.equal(result.initialRedistribution.correctionScale, .25);
  assert.equal(result.initialRedistribution.backtracks, 2);
  assert.deepEqual(result.initialRedistribution.rejections.map(r => [r.correctionScale, r.stage]),
    [[1, 'admissibility'], [.5, 'admissibility']]);
  assert.deepEqual(d.calls.map(c => c.options.correctionScale), [1, 1, 1, .5, .5, .5, .25, .25, .25]);
  for (let i = 3; i < d.calls.length; i++) assert.deepEqual(d.calls[i].nodes, d.calls[i % 3].nodes);
  const original = d.systems[0], value = original.evaluate(original.initial), ne = original.ne;
  assert.deepEqual(plain(result.x.slice(0, original.euler.layout.densityCount)),
    plain(original.initial.slice(0, original.euler.layout.densityCount)));
  assert.deepEqual(plain(result.x.slice(ne)), plain(original.initial.slice(ne)));
  assert.deepEqual(result.boundaryLayer.transitionState, original.bl.snapshotActive());
  assert.deepEqual(result.flow.allocation, value.outer.allocation);
  for (const trial of d.supplied.slice(1)) {
    assert.deepEqual(trial.options.initialEuler.undisplacedNodes, value.outer.undisplacedNodes);
    assert.deepEqual(trial.options.initialBL, plain(original.initial.slice(ne)));
    assert.deepEqual(trial.options.transitionState, original.bl.snapshotActive());
  }
  assert.equal(result.mesh.quality.valid, true);
  assert.equal(d.rejected.length, 2);
  assert.equal(checkpoints.length, 1);
  assert.deepEqual(checkpoints[0].restart.initialEuler.nodes, result.flow.nodes);
  assert.equal(result.linearDiagnostics.solves, 0);
  assert.equal(d.linearCalls.length, 0);
  assert.deepEqual(plain(f), before);
});

test('failed initial coupled redistribution preserves phase, BL fields, density, mass and both geometry descriptions', async () => {
  for (const [policy, attempts] of [['listing', 1], ['admissible', 2]]) {
    const f = fixture(), d = await driver({ rejectedCandidates: 10 });
    const result = d.solve(f.input, { ...f.options, maxIterations: 0, stepAcceptance: policy, maxBacktracks: 1,
      onCheckpoint: () => assert.fail('A rejected initializer cannot publish a restart.') });
    const source = d.systems[0], original = source.evaluate(source.initial);
    assert.equal(result.initialRedistribution.accepted, false);
    assert.match(result.reason, /initial redistribution rejected/);
    assert.equal(d.rejected.length, attempts);
    assert.equal('checkpoint' in result, false);
    assert.deepEqual(plain(result.x), plain(source.initial));
    assert.deepEqual(plain(result.residual), plain(original.residual));
    assert.deepEqual(result.flow.nodes, original.outer.nodes);
    assert.deepEqual(result.flow.undisplacedNodes, original.outer.undisplacedNodes);
    assert.deepEqual(result.flow.allocation, original.outer.allocation);
    assert.deepEqual(result.boundaryLayer.transitionState, source.bl.snapshotActive());
    assert.equal(result.linearDiagnostics.solves, 0);
  }
});

test('initial full coupled SMOVE remains exact for both listing and admissible policies; resume does not repeat it', async () => {
  const f = fixture(), outcomes = [];
  for (const policy of ['listing', 'admissible']) {
    const d = await driver(), result = d.solve(f.input, { ...f.options, maxIterations: 0, stepAcceptance: policy });
    assert.equal(result.initialRedistribution.accepted, true, result.reason);
    assert.deepEqual(d.calls.map(c => c.options.correctionScale), [1, 1, 1]);
    assert.equal('correctionScale' in result.initialRedistribution, false);
    assert.ok(result.initialRedistribution.passages.every(p => !('correctionScale' in p)));
    outcomes.push(result);
  }
  for (const key of ['x', 'flow', 'residual', 'families', 'boundaryLayer', 'initialRedistribution'])
    assert.deepEqual(plain(outcomes[1][key]), plain(outcomes[0][key]), key);
  const d = await driver({ failMovement: true });
  const resumed = d.solve(undefined, { resume: outcomes[1].checkpoint, maxIterations: 0, stepAcceptance: 'admissible' });
  assert.equal(d.calls.length, 0);
  assert.deepEqual(plain(resumed.checkpoint), plain(outcomes[1].checkpoint));
});

test('initial coupled coordinate-span rejection uses the same bounded correction scale', async () => {
  const f = fixture(), d = await driver({ failMovement: true });
  const result = d.solve(f.input, { ...f.options, maxIterations: 0, stepAcceptance: 'admissible', maxBacktracks: 1 });
  assert.equal(result.initialRedistribution.accepted, true, result.reason);
  assert.equal(result.initialRedistribution.correctionScale, .5);
  assert.equal(result.initialRedistribution.rejections[0].stage, 'SMOVE');
  assert.deepEqual(d.calls.map(c => c.options.correctionScale), [1, .5, .5, .5]);
});

let triggerSource;
function triggerFixture() {
  const f = fixture();
  triggerSource ??= solveCoupledStreamtubeIses(f.input, { ...f.options, maxIterations: 0,
    stepAcceptance: 'admissible' }).checkpoint;
  const checkpoint = plain(triggerSource), r = checkpoint.restart;
  const system = createCoupledStreamtubeBody(r.input, { ...r.options, initialEuler: r.initialEuler,
    initialBL: Float64Array.from(r.initialBL) });
  const i = system.euler.layout.bodies[0].leadingIndex + 1, nodes = r.initialEuler.nodes;
  const a = nodes[0][i].at(-1), b = nodes[1][i][0], spacing = .5 * Math.hypot(a.x - b.x, a.y - b.y);
  checkpoint.continuation.lastRedistributedStagnation[0] -= (.5 - 1e-8) * spacing;
  const direction = new Float64Array(system.n), column = system.euler.layout.globals.stagnation[0];
  direction[column] = .001 * spacing / system.euler.curves[0].length;
  return { checkpoint, direction, column };
}

test('triggered coupled SMOVE halves one coordinate correction for every passage with the common Newton step', async () => {
  const { checkpoint, direction, column } = triggerFixture(), before = plain(checkpoint);
  const d = await driver({ direction, rejectedCandidates: 2 }), checkpoints = [];
  const result = d.solve(undefined, { resume: checkpoint, maxIterations: 1, stepAcceptance: 'admissible', maxBacktracks: 2,
    onCheckpoint: cp => checkpoints.push(plain(cp)) });
  assert.equal(result.history.length, 2, result.reason);
  const update = result.history[1];
  assert.equal(update.step, .25);
  assert.equal(update.backtracks, 2);
  assert.deepEqual(update.maintenance.triggeredBodies, [0]);
  assert.deepEqual(d.calls.map(c => c.options.correctionScale), [1, 1, 1, .5, .5, .5, .25, .25, .25]);
  assert.deepEqual(update.maintenance.passages.map(p => p.correctionScale), [.25, .25, .25]);
  assert.equal(result.x[column], checkpoint.restart.initialEuler.x[column] + .25 * direction[column]);
  assert.equal(d.linearCalls.length, 1);
  assert.equal(result.linearDiagnostics.solves, 1);
  assert.equal(checkpoints.length, 2);
  assert.deepEqual(checkpoints[0], before);
  assert.deepEqual(result.checkpoint.continuation.lastRedistributedStagnation, result.flow.stagnation);
  assert.deepEqual(plain(checkpoint), before);
});

test('exhausted triggered coupled retries restore the complete accepted phase, fields and redistribution history', async () => {
  const { checkpoint, direction } = triggerFixture(), before = plain(checkpoint);
  const d = await driver({ direction, rejectedCandidates: 10 }), checkpoints = [];
  const result = d.solve(undefined, { resume: checkpoint, maxIterations: 1, stepAcceptance: 'admissible', maxBacktracks: 1,
    onCheckpoint: cp => checkpoints.push(plain(cp)) });
  assert.equal(result.history.length, 1);
  assert.equal(result.lastRejectedStep.stage, 'admissibility');
  assert.deepEqual(d.calls.map(c => c.options.correctionScale), [1, 1, 1, .5, .5, .5]);
  assert.deepEqual(plain(result.checkpoint), before);
  assert.deepEqual(plain(result.x), [...before.restart.initialEuler.x, ...before.restart.initialBL]);
  assert.deepEqual(result.flow.nodes, before.restart.initialEuler.nodes);
  assert.deepEqual(result.flow.undisplacedNodes, before.restart.initialEuler.undisplacedNodes);
  assert.deepEqual(result.boundaryLayer.transitionState, before.restart.options.transitionState);
  assert.equal(checkpoints.length, 1);
  assert.deepEqual(plain(checkpoint), before);
});

test('raw coupled grid retries retain the full correction on the first actual SMOVE attempt', async () => {
  const { checkpoint, direction } = triggerFixture(), before = plain(checkpoint);
  const d = await driver({ direction, rejectedRawProposals: 1 }), checkpoints = [];
  const result = d.solve(undefined, { resume: checkpoint, maxIterations: 1, stepAcceptance: 'admissible', maxBacktracks: 1,
    onCheckpoint: cp => checkpoints.push(plain(cp)) });
  assert.equal(result.history.length, 2, result.reason);
  const update = result.history[1];
  assert.equal(update.step, .5); assert.equal(update.backtracks, 1);
  assert.equal(update.rejections[0].stage, 'Newton grid step');
  assert.deepEqual(update.maintenance.triggeredBodies, [0]);
  assert.deepEqual(d.calls.map(c => c.options.correctionScale), [1, 1, 1]);
  assert.ok(update.maintenance.passages.every(p => !('correctionScale' in p)));
  assert.equal(checkpoints.length, 2); assert.equal(result.mesh.quality.valid, true);
  assert.equal(d.linearCalls.length, 1); assert.deepEqual(plain(checkpoint), before);
  const reference = await driver({ direction: direction.map(v => .5 * v) });
  const expected = reference.solve(undefined, { resume: checkpoint, maxIterations: 1, stepAcceptance: 'admissible' });
  for (const key of ['x', 'flow', 'residual', 'boundaryLayer']) assert.deepEqual(plain(result[key]), plain(expected[key]), key);
});

test('a resumed wake recovery retains its coordinate policy across backtracking and checkpoints', async () => {
  const { checkpoint, direction } = triggerFixture();
  checkpoint.continuation.wakeCoordinateRecovery = true;
  const before = plain(checkpoint), d = await driver({ direction, rejectedRawProposals: 1 });
  const result = d.solve(undefined, { resume: checkpoint, maxIterations: 1,
    stepAcceptance: 'admissible', maxBacktracks: 2 });
  assert.equal(result.history.length, 2, result.reason);
  assert.deepEqual(d.gridBounds, [1, 1]);
  assert.equal(result.checkpoint.continuation.wakeCoordinateRecovery, true);
  const resumed = d.solve(undefined, { resume: result.checkpoint, maxIterations: 0, stepAcceptance: 'admissible' });
  assert.deepEqual(plain(resumed.checkpoint), plain(result.checkpoint));
  assert.deepEqual(plain(checkpoint), before);
  assert.throws(() => d.solve(undefined, { resume: { ...checkpoint,
    continuation: { ...checkpoint.continuation, wakeCoordinateRecovery: 'yes' } },
  maxIterations: 0, stepAcceptance: 'admissible' }), /Invalid coupled ISES iteration controls/);
});

test('coupled correction damping counts rejected SMOVE candidates separately from raw proposals', async () => {
  const { checkpoint, direction } = triggerFixture();
  const d = await driver({ direction, rejectedRawProposals: 1, rejectedCandidates: 1 });
  const result = d.solve(undefined, { resume: checkpoint, maxIterations: 1, stepAcceptance: 'admissible', maxBacktracks: 2 });
  assert.equal(result.history.length, 2, result.reason);
  const update = result.history[1];
  assert.equal(update.step, .25);
  assert.deepEqual(update.rejections.map(r => r.stage), ['Newton grid step', 'admissibility']);
  assert.deepEqual(d.calls.map(c => c.options.correctionScale), [1, 1, 1, .5, .5, .5]);
  assert.deepEqual(update.maintenance.passages.map(p => p.correctionScale), [.5, .5, .5]);
  assert.equal(result.mesh.quality.valid, true); assert.equal(d.linearCalls.length, 1);
});
