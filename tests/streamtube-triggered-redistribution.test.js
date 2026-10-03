// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { solveStreamtubeIses } from '../src/euler/streamtube-ises-update.js';
import { redistributeStreamtubeTangentially } from '../src/geometry/streamtube-tangential-redistribution.js';
import { requireConvexGridUpdate, requireConvexPublishedGrid } from '../src/euler/streamtube-grid-update.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

// Control-flow regression: use a real tiny multielement grid, coordinate
// redistribution and gas evaluation, but supply one controlled Newton
// direction instead of assembling a Jacobian or doing a global linear solve.
// The redistribution event is deliberately close, as in the RAE64 failure.
let source;
function fixture(policy = 'admissible') {
  source ??= solveStreamtubeIses(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), {
    maxIterations: 0, retainCheckpoint: true, stepAcceptance: 'admissible',
  }).checkpoint;
  const checkpoint = structuredClone(source), system = createStreamtubeBodySystem(checkpoint.input);
  const i = system.layout.bodies[0].leadingIndex + 1;
  const lower = checkpoint.initialEuler.nodes[0][i].at(-1), upper = checkpoint.initialEuler.nodes[1][i][0];
  const spacing = .5 * Math.hypot(upper.x - lower.x, upper.y - lower.y);
  checkpoint.continuation.stepAcceptance = policy;
  checkpoint.continuation.lastRedistributedStagnation[0] -= (.5 - 1e-8) * spacing;
  const column = system.layout.globals.stagnation[0], direction = new Float64Array(system.layout.n);
  direction[column] = .001 * spacing / system.curves[0].length;
  return { checkpoint, direction, column };
}

async function driver(direction, { rejectedCandidates = 0, failMovement = false, rejectedRawProposals = 0 } = {}) {
  const url = new URL('../src/euler/streamtube-ises-update.js', import.meta.url);
  const key = `ises-triggered-redistribution-test-${Math.random()}`;
  const calls = [], systems = [], evaluations = [], linearCalls = [];
  const seam = {
    requireConvexPublishedGrid,
    requireConvexGridUpdate(previous, proposed) {
      const result = requireConvexGridUpdate(previous, proposed);
      if (rejectedRawProposals-- > 0) throw Object.assign(new Error('Controlled raw Newton grid rejection.'), {
        code: 'streamtube-grid-step', stepFraction: .5,
      });
      return result;
    },
    createStreamtubeBodySystem(input) {
      const system = createStreamtubeBodySystem(input), index = systems.length, evaluate = system.evaluate;
      systems.push(system);
      system.jacobian = x => ({ controlledDirection: true, x: Array.from(x) });
      system.evaluate = x => {
        const value = evaluate(x), rejected = index > 0 && index <= rejectedCandidates;
        evaluations.push({ index, rejected, nodes: structuredClone(value.nodes) });
        if (rejected) throw new Error('Controlled candidate interface pressure rejection.');
        return value;
      };
      return system;
    },
    solveStreamtubeBody,
    solveSparseDirect(matrix, rhs) {
      assert.equal(matrix.controlledDirection, true);
      linearCalls.push({ state: matrix.x, rhs: Array.from(rhs) });
      return { x: direction.slice(), relativeResidual: 0, refinements: 0 };
    },
    redistributeStreamtubeTangentially(nodes, options) {
      calls.push({ nodes: structuredClone(nodes), options: structuredClone(options) });
      if (failMovement && options.correctionScale > .5)
        throw new Error('Controlled nonpositive tangential coordinate span.');
      return redistributeStreamtubeTangentially(nodes, options);
    },
  };
  globalThis[key] = seam;
  let text = fs.readFileSync(url, 'utf8');
  for (const statement of [
    "import { createStreamtubeBodySystem, solveStreamtubeBody } from './streamtube-body.js';",
    "import { solveSparseDirect } from '../numerics/klu.js';",
    "import { redistributeStreamtubeTangentially } from '../geometry/streamtube-tangential-redistribution.js';",
    "import { requireConvexGridUpdate, requireConvexPublishedGrid } from './streamtube-grid-update.js';",
  ]) {
    assert.ok(text.includes(statement));
    const names = statement.slice(statement.indexOf('{'), statement.indexOf('}') + 1);
    text = text.replace(statement, `const ${names} = globalThis[${JSON.stringify(key)}];`);
  }
  text = text.replace(/from '(\.[^']+)'/g, (_, relative) => `from '${new URL(relative, url).href}'`);
  try {
    const { solveStreamtubeIses: solve } = await import(`data:text/javascript;base64,${Buffer.from(text).toString('base64')}`);
    return { solve, calls, systems, evaluations, linearCalls };
  } finally { delete globalThis[key]; }
}

test('triggered SMOVE backtracks its coordinate correction and the Newton step together in every passage', async () => {
  const { checkpoint, direction, column } = fixture(), before = structuredClone(checkpoint);
  const d = await driver(direction, { rejectedCandidates: 2 }), publications = [];
  const result = d.solve(undefined, { resume: checkpoint, maxIterations: 1, stepAcceptance: 'admissible', maxBacktracks: 2,
    onCheckpoint: (cp, details) => publications.push({ cp, details }) });
  assert.equal(result.history.length, 2, result.reason);
  const update = result.history[1];
  assert.equal(update.step, .25);
  assert.equal(update.backtracks, 2);
  assert.deepEqual(update.rejections.map(r => [r.step, r.stage]), [[1, 'admissibility'], [.5, 'admissibility']]);
  assert.deepEqual(update.maintenance.triggeredBodies, [0]);
  assert.deepEqual(d.calls.map(c => c.options.correctionScale), [1, 1, 1, .5, .5, .5, .25, .25, .25]);
  assert.deepEqual(update.maintenance.passages.map(p => p.correctionScale), [.25, .25, .25]);
  assert.equal(result.x[column], checkpoint.initialEuler.x[column] + .25 * direction[column]);
  assert.deepEqual(Array.from(result.x).slice(0, d.systems[0].layout.densityCount),
    checkpoint.initialEuler.x.slice(0, d.systems[0].layout.densityCount));
  assert.equal(d.linearCalls.length, 1);
  assert.equal(result.linearDiagnostics.solves, 1);
  assert.equal(result.finalQuality.valid, true);
  assert.equal(publications.length, 2, 'Rejected movements must not publish checkpoints.');
  assert.deepEqual(publications[0].cp, before);
  assert.deepEqual(publications[1].cp.initialEuler.nodes, result.nodes);
  assert.deepEqual(publications[1].cp.residual, Array.from(result.residual));
  assert.deepEqual(result.checkpoint.continuation.lastRedistributedStagnation, result.stagnation);
  assert.equal(d.evaluations.filter(e => e.rejected).length, 2);
  assert.deepEqual(checkpoint, before);
});

test('successful first triggered movement retains the listed full correction and exact physical result', async () => {
  const outcomes = [];
  for (const policy of ['listing', 'admissible']) {
    const { checkpoint, direction } = fixture(policy), d = await driver(direction);
    const result = d.solve(undefined, { resume: checkpoint, maxIterations: 1, stepAcceptance: policy });
    assert.equal(result.history.length, 2, result.reason);
    assert.equal(result.history[1].step, 1);
    assert.deepEqual(result.history[1].maintenance.triggeredBodies, [0]);
    assert.deepEqual(d.calls.map(c => c.options.correctionScale), [1, 1, 1]);
    assert.ok(result.history[1].maintenance.passages.every(p => !('correctionScale' in p)));
    outcomes.push(result);
  }
  for (const key of ['x', 'nodes', 'residual', 'diagnostics', 'stagnation', 'finalQuality'])
    assert.deepEqual(outcomes[1][key], outcomes[0][key], key);
  assert.deepEqual(outcomes[1].history[1].maintenance, outcomes[0].history[1].maintenance);
});

test('a rejected triggered movement cannot commit its stagnation history, geometry, density or checkpoint', async () => {
  for (const [policy, maxBacktracks, scales] of [['listing', 2, [1, 1, 1]], ['admissible', 1, [1, 1, 1, .5, .5, .5]]]) {
    const { checkpoint, direction } = fixture(policy), before = structuredClone(checkpoint);
    const d = await driver(direction, { rejectedCandidates: 10 }), publications = [];
    const result = d.solve(undefined, { resume: checkpoint, maxIterations: 1, stepAcceptance: policy, maxBacktracks,
      onCheckpoint: cp => publications.push(cp) });
    assert.equal(result.history.length, 1);
    assert.match(result.reason, /rejected during admissibility.*Controlled candidate interface pressure/);
    assert.deepEqual(d.calls.map(c => c.options.correctionScale), scales);
    assert.equal(publications.length, 1);
    assert.deepEqual(result.checkpoint, before);
    assert.deepEqual(result.nodes, before.initialEuler.nodes);
    assert.deepEqual(Array.from(result.x), before.initialEuler.x);
    assert.deepEqual(Array.from(result.residual), before.residual);
    assert.deepEqual(checkpoint, before);
  }
});

test('a failed coordinate movement is retried from the retained state with a reduced coordinate correction', async () => {
  const { checkpoint, direction } = fixture(), d = await driver(direction, { failMovement: true });
  const result = d.solve(undefined, { resume: checkpoint, maxIterations: 1, stepAcceptance: 'admissible', maxBacktracks: 1 });
  assert.equal(result.history.length, 2, result.reason);
  assert.equal(result.history[1].step, .5);
  assert.equal(result.history[1].rejections[0].stage, 'SMOVE');
  assert.deepEqual(d.calls.map(c => c.options.correctionScale), [1, .5, .5, .5]);
  assert.deepEqual(result.history[1].maintenance.passages.map(p => p.correctionScale), [.5, .5, .5]);
  assert.equal(d.linearCalls.length, 1);
});

test('raw Euler grid retries retain the full correction on the first actual SMOVE attempt', async () => {
  const { checkpoint, direction } = fixture(), before = structuredClone(checkpoint);
  const d = await driver(direction, { rejectedRawProposals: 1 }), checkpoints = [];
  const result = d.solve(undefined, { resume: checkpoint, maxIterations: 1, stepAcceptance: 'admissible', maxBacktracks: 1,
    onCheckpoint: cp => checkpoints.push(cp) });
  assert.equal(result.history.length, 2, result.reason);
  const update = result.history[1];
  assert.equal(update.step, .5); assert.equal(update.backtracks, 1);
  assert.equal(update.rejections[0].stage, 'Newton grid step');
  assert.deepEqual(update.maintenance.triggeredBodies, [0]);
  assert.deepEqual(d.calls.map(c => c.options.correctionScale), [1, 1, 1]);
  assert.ok(update.maintenance.passages.every(p => !('correctionScale' in p)));
  assert.equal(checkpoints.length, 2); assert.equal(result.finalQuality.valid, true);
  assert.equal(d.linearCalls.length, 1); assert.deepEqual(checkpoint, before);
  // The retained geometry must match the listed full redistribution of
  // this same smaller Newton proposal, not a pre-damped coordinate field.
  const reference = await driver(direction.map(v => .5 * v));
  const expected = reference.solve(undefined, { resume: checkpoint, maxIterations: 1, stepAcceptance: 'admissible' });
  for (const key of ['x', 'nodes', 'residual']) assert.deepEqual(result[key], expected[key], key);
});

test('Euler correction damping counts rejected SMOVE candidates separately from raw proposals', async () => {
  const { checkpoint, direction } = fixture(), d = await driver(direction, { rejectedRawProposals: 1, rejectedCandidates: 1 });
  const result = d.solve(undefined, { resume: checkpoint, maxIterations: 1, stepAcceptance: 'admissible', maxBacktracks: 2 });
  assert.equal(result.history.length, 2, result.reason);
  const update = result.history[1];
  assert.equal(update.step, .25);
  assert.deepEqual(update.rejections.map(r => r.stage), ['Newton grid step', 'admissibility']);
  assert.deepEqual(d.calls.map(c => c.options.correctionScale), [1, 1, 1, .5, .5, .5]);
  assert.deepEqual(update.maintenance.passages.map(p => p.correctionScale), [.5, .5, .5]);
  assert.equal(result.finalQuality.valid, true); assert.equal(d.linearCalls.length, 1);
});
