// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { createCoupledStreamtubeBody, coupledStreamtubeTripEvents, coupledStreamtubeResult } from '../src/euler/streamtube-coupled.js';
import { solveStreamtubeIses } from '../src/euler/streamtube-ises-update.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { redistributeStreamtubeTangentially } from '../src/geometry/streamtube-tangential-redistribution.js';
import { limitStreamtubeGridStep, interpolateStreamtubeGridNodes } from '../src/geometry/streamtube-convex-step.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const plain = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const clone = structuredClone;
// An independent signed-corner check. The regression must not certify a
// guard with the same production predicate that it is meant to exercise.
function convex(nodes) {
  for (const group of nodes) for (let i = 0; i < group.length - 1; i++) for (let j = 0; j < group[i].length - 1; j++) {
    const p = [group[i][j], group[i + 1][j], group[i + 1][j + 1], group[i][j + 1]];
    for (let k = 0; k < 4; k++) {
      const a = p[k], b = p[(k + 1) % 4], c = p[(k + 2) % 4];
      if (!((b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x) > 0)) return false;
    }
  }
  return true;
}
const input = () => ({ ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
  wakeGeometry: 'independent-banks', wakeOutlet: 'banks' });
const options = { edgeMatching: 'section-velocity', transitionMode: 'automatic', reynolds: 1e6, ncrit: 9,
  iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible' };
const sources = new Map();
function fixture(kind, blUpdate = 'giles') {
  const key = `${kind}:${blUpdate}`;
  if (!sources.has(key)) {
    const value = kind === 'euler'
      ? solveStreamtubeIses(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), {
        ...options, maxIterations: 0, retainCheckpoint: true })
      : solveCoupledStreamtubeIses(input(), { ...options, maxIterations: 0, blUpdate });
    assert(value.checkpoint, value.reason);
    sources.set(key, plain(value.checkpoint));
  }
  return clone(sources.get(key));
}
function restored(kind, checkpoint) {
  if (kind === 'euler') {
    const system = createStreamtubeBodySystem(checkpoint.input);
    const state = system.adoptGeometry(Float64Array.from(checkpoint.initialEuler.x), checkpoint.initialEuler.nodes);
    return { system, euler: system, state, nodes: checkpoint.initialEuler.nodes, ne: system.layout.n };
  }
  const r = checkpoint.restart;
  const system = createCoupledStreamtubeBody(r.input, { ...r.options, initialEuler: r.initialEuler,
    initialBL: Float64Array.from(r.initialBL) });
  return { system, euler: system.euler, state: system.initial, nodes: r.initialEuler.nodes, ne: system.ne };
}

// Real tiny Euler/BL systems and their evaluations are retained. Only the
// global factorization is replaced with a prescribed geometric challenge.
// This isolates update policy without solving an unrelated flow problem.
async function driver(kind, direction, { foldRedistribution = false } = {}) {
  const url = new URL(kind === 'euler' ? '../src/euler/streamtube-ises-update.js' : '../src/euler/streamtube-coupled-ises.js', import.meta.url);
  const key = `convex-driver-${Math.random()}`, systems = [], linearCalls = [], redistributionCalls = [];
  const seam = { solveStreamtubeBody, coupledStreamtubeTripEvents, coupledStreamtubeResult,
    createStreamtubeBodySystem(data) {
      const system = createStreamtubeBodySystem(data); systems.push(system);
      system.jacobian = state => ({ controlledDirection: true, state: Array.from(state) }); return system;
    },
    createCoupledStreamtubeBody(data, controls) {
      const system = createCoupledStreamtubeBody(data, controls); systems.push(system);
      system.jacobian = state => ({ controlledDirection: true, state: Array.from(state) }); return system;
    },
    solveSparseDirect(matrix) {
      assert.equal(matrix.controlledDirection, true); assert(direction);
      linearCalls.push(matrix.state);
      return { x: direction.slice(), relativeResidual: 0, refinements: 0, ordering: 'amd', pivotTolerance: .001, attempts: [] };
    },
    solveCoupledLinearSystem(matrix) {
      return { ordering: null, linear: seam.solveSparseDirect(matrix), stationPolicy: { recommendation: 'station-auto' } };
    },
    solveSparseDirectAligned(matrix) {
      return { ...seam.solveSparseDirect(matrix), factorNonzeros: 0, equationOrdering: 'aligned' };
    },
    redistributeStreamtubeTangentially(nodes, controls) {
      redistributionCalls.push(clone(controls));
      const result = redistributeStreamtubeTangentially(nodes, controls);
      if (foldRedistribution) {
        const moved = clone(result.nodes);
        moved[1][1] = { x: moved[1][0].x, y: moved[1][0].y - 1 };
        assert.equal(convex([moved]), false);
        return { ...result, nodes: moved };
      }
      return result;
    },
  };
  globalThis[key] = seam;
  let source = fs.readFileSync(url, 'utf8');
  const statements = kind === 'euler'
    ? ["import { createStreamtubeBodySystem, solveStreamtubeBody } from './streamtube-body.js';"]
    : ["import { createCoupledStreamtubeBody, coupledStreamtubeTripEvents, coupledStreamtubeResult } from './streamtube-coupled.js';"];
  statements.push(kind === 'euler' ? "import { solveSparseDirectAligned } from '../numerics/klu.js';"
    : "import { solveSparseDirect } from '../numerics/klu.js';",
    "import { redistributeStreamtubeTangentially } from '../geometry/streamtube-tangential-redistribution.js';");
  if (kind === 'coupled') statements.push("import { solveCoupledLinearSystem } from './streamtube-coupled-linear-solve.js';");
  for (const statement of statements) {
    assert(source.includes(statement));
    const names = statement.slice(statement.indexOf('{'), statement.indexOf('}') + 1);
    source = source.replace(statement, `const ${names} = globalThis[${JSON.stringify(key)}];`);
  }
  source = source.replace(/from '(\.[^']+)'/g, (_, relative) => `from '${new URL(relative, url).href}'`);
  try {
    const module = await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);
    return { solve: module[kind === 'euler' ? 'solveStreamtubeIses' : 'solveCoupledStreamtubeIses'], systems, linearCalls, redistributionCalls };
  } finally { delete globalThis[key]; }
}

test('the real RAE64 Newton proposal is shortened before the first concave corner', () => {
  const { rawNewton } = JSON.parse(fs.readFileSync(new URL('./fixtures/rae64-grid-update.json', import.meta.url)));
  const previous = [rawNewton.before], proposed = [rawNewton.proposed];
  assert.equal(convex(previous), true);
  assert.equal(convex(proposed), false);
  const limited = limitStreamtubeGridStep(previous, proposed);
  assert.equal(limited.limited, true);
  assert(limited.step > 0 && limited.step < .0031);
  assert.equal(convex(interpolateStreamtubeGridNodes(previous, proposed, limited.step)), true);
  assert.equal(convex(interpolateStreamtubeGridNodes(previous, proposed, Math.min(1, limited.limiter.boundaryStep * 1.001))), false);
});

for (const [kind, blUpdate] of [['euler', 'giles'], ['coupled', 'giles'], ['coupled', 'xfoil']]) {
  test(`${kind}/${blUpdate}: one Newton direction preserves every published grid and scales density and BL together`, async () => {
    const checkpoint = fixture(kind, blUpdate), r = restored(kind, checkpoint);
    // A nonzero but sub-threshold reference offset makes the maintenance
    // reset observable without triggering Giles's stagnation criterion.
    checkpoint.continuation.lastRedistributedStagnation = checkpoint.continuation.lastRedistributedStagnation.map(s => s - 1e-8);
    const before = clone(checkpoint);
    const direction = new Float64Array(r.state.length);
    const column = r.euler.layout.nodes[0][1][1].column;
    direction[column] = -3;
    direction[0] = .01;
    let thetaColumn;
    if (kind === 'coupled') {
      thetaColumn = r.ne + 4 * r.system.bl.stations[0].id + 1;
      direction[thetaColumn] = .01 * r.state[thetaColumn];
    }
    const raw = r.euler.decode(Float64Array.from(r.state.subarray(0, r.ne), (v, k) => v + direction[k])).nodes;
    assert.equal(convex(raw), false, 'The prescribed full Newton direction must actually fold the physical grid.');
    const d = await driver(kind, direction), publications = [], meshes = [];
    const result = d.solve(undefined, { resume: checkpoint, ...options, maxIterations: 1, maxBacktracks: 12,
      onCheckpoint: cp => publications.push(plain(cp)), onMesh: update => meshes.push(clone(update.nodes)) });
    assert.equal(result.history.length, 2, result.reason);
    const update = result.history[1];
    assert(update.step > 0 && update.step < 1);
    assert(update.rejections.some(row => row.stage === 'Newton grid step' && row.code === 'streamtube-grid-step'));
    // This extension retains Giles's five-pair SMOVE, but additionally calls
    // it when the geometric step bound fires before his stagnation trigger.
    assert.deepEqual(update.maintenance.triggeredBodies, []);
    assert.equal(update.maintenance.geometryRedistribution, true);
    assert.equal(update.maintenance.passages.length, r.euler.layout.tubes.length);
    assert(update.maintenance.passages.every(p => p.pairs === 5));
    assert.equal(d.linearCalls.length, 1);
    assert.equal(result.linearDiagnostics.solves, 1);
    assert.equal(publications.length, 2);
    assert(meshes.length >= 2);
    for (const nodes of meshes) assert.equal(convex(nodes), true);
    for (const cp of publications) assert.equal(convex(kind === 'euler' ? cp.initialEuler.nodes : cp.restart.initialEuler.nodes), true);
    assert.deepEqual(publications.at(-1).continuation.lastRedistributedStagnation,
      Array.from(kind === 'euler' ? result.stagnation : result.flow.stagnation));
    assert.notDeepEqual(publications.at(-1).continuation.lastRedistributedStagnation,
      before.continuation.lastRedistributedStagnation);
    assert.equal(kind === 'euler' ? result.finalQuality.valid : result.mesh.quality.valid, true);
    assert(Math.abs(result.x[0] - r.state[0] - Math.log1p(update.step * direction[0])) < 2e-16);
    if (thetaColumn !== undefined)
      assert.equal(result.x[thetaColumn], r.state[thetaColumn] + update.step * direction[thetaColumn]);
    assert.deepEqual(checkpoint, before, 'Retries must not mutate the supplied accepted checkpoint.');
  });
}

for (const kind of ['euler', 'coupled']) {
  test(`${kind}: a gas-admissible concave initial state cannot publish a checkpoint`, async () => {
    const checkpoint = fixture(kind), r = restored(kind, checkpoint), proposed = clone(r.nodes);
    proposed[0][1][1].y += 1;
    const limit = limitStreamtubeGridStep(r.nodes, proposed);
    assert.equal(limit.limited, true);
    const nodes = interpolateStreamtubeGridNodes(r.nodes, proposed, limit.limiter.boundaryStep * 1.001);
    assert.equal(convex(nodes), false);
    const d = await driver(kind, null);
    if (kind === 'euler') {
      const x = r.system.adoptGeometry(r.state, nodes);
      assert.doesNotThrow(() => r.system.evaluate(x), 'Sampled conservation/gas checks alone allow this actual concave grid.');
      assert.throws(() => d.solve(checkpoint.input, { ...options, maxIterations: 0,
        initialEuler: { x, nodes }, onCheckpoint: () => assert.fail('An invalid initial state cannot publish.') }),
      error => error.code === 'streamtube-grid-nonconvex');
    } else {
      const seed = checkpoint.restart;
      const x = r.euler.adoptGeometry(r.state.subarray(0, r.ne), nodes);
      const initialEuler = { x, ...r.euler.decode(x) };
      const candidate = createCoupledStreamtubeBody(seed.input, { ...seed.options,
        initialEuler, initialBL: Float64Array.from(seed.initialBL) });
      assert.equal(convex(candidate.euler.decode(candidate.initial.subarray(0, candidate.ne)).nodes), false);
      assert.equal(candidate.admissible(candidate.initial, { requireConvex: false }), true);
      assert.throws(() => d.solve(seed.input, { ...seed.options, ...options, maxIterations: 0,
        initialEuler, initialBL: Float64Array.from(seed.initialBL),
        onCheckpoint: () => assert.fail('An invalid initial state cannot publish.') }), /Inadmissible coupled ISES initial state/);
    }
    assert.equal(d.linearCalls.length, 0);
    assert.equal(d.redistributionCalls.length, 0);
  });

  test(`${kind}: a concave initial SMOVE proposal is rejected before creating or publishing its candidate`, async () => {
    const d = await driver(kind, null, { foldRedistribution: true }), checkpoints = [];
    const data = kind === 'euler' ? intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }) : input();
    const result = d.solve(data, { ...options, maxIterations: 0, maxBacktracks: 1, retainCheckpoint: true,
      onCheckpoint: cp => checkpoints.push(cp) });
    assert.equal(result.initialRedistribution.accepted, false);
    assert.match(result.reason, /Nonconvex grid cell/);
    assert.equal(checkpoints.length, 0);
    assert.equal('checkpoint' in result, false);
    assert.equal(d.systems.length, 1, 'A known concave proposal must not reach physical candidate construction.');
    assert.equal(d.linearCalls.length, 0);
    assert.equal(convex(kind === 'euler' ? result.nodes : result.flow.nodes), true);
  });

  test(`${kind}: a safe full update retains the exact listing result`, async () => {
    const original = fixture(kind), r = restored(kind, original), direction = new Float64Array(r.state.length), results = [];
    direction[0] = 1e-4;
    for (const policy of ['listing', 'admissible']) {
      const checkpoint = clone(original); checkpoint.continuation.stepAcceptance = policy;
      const d = await driver(kind, direction);
      const result = d.solve(undefined, { resume: checkpoint, ...options, stepAcceptance: policy, maxIterations: 1 });
      assert.equal(result.history.length, 2, result.reason);
      assert.equal(result.history[1].step, 1);
      assert.equal(result.history[1].maintenance.geometryRedistribution, undefined);
      assert.equal(d.redistributionCalls.length, 0);
      assert.deepEqual(result.checkpoint.continuation.lastRedistributedStagnation,
        checkpoint.continuation.lastRedistributedStagnation);
      results.push(result);
    }
    for (const key of kind === 'euler' ? ['x', 'nodes', 'residual', 'stagnation'] : ['x', 'flow', 'residual', 'boundaryLayer'])
      assert.deepEqual(plain(results[1][key]), plain(results[0][key]), key);
  });
}
