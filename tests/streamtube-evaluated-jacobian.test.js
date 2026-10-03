// SPDX-License-Identifier: GPL-2.0-or-later
// Tiny complete matrices only: no LU, Newton, mesh initialization or flow solve.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { archivedBody, archivedCoupled } from './helpers/archived-streamtube-assembly.js';

const fixture = JSON.parse(fs.readFileSync(new URL('../docs/coupled-current-profile-preparation/two-element-six-update/initial.json', import.meta.url))).checkpoint;
const plainValue = value => Object.fromEntries(Object.entries(value).filter(([, v]) => typeof v !== 'function'));
const snapshot = value => structuredClone(plainValue(value));
const thickness = input => ({ surfaces: input.bodies.map(b => Object.fromEntries(['upper', 'lower'].map(side => [side,
  Array.from({ length: b.trailingIndex - b.leadingIndex + 1 }, (_, i) => .0001 * (1 + .02 * i))]))),
  wakes: input.bodies.map(b => Array.from({ length: input.outerLower.length - 1 - b.trailingIndex }, (_, i) => .00024 * (1 + .01 * i))) });
function blockSnapshot(block, system) {
  return { state: block.state, displacement: block.displacement, parameters: block.parameters,
    boundaries: system.layout.bodies.flatMap((_, body) => ['lower', 'upper'].flatMap(side =>
      Array.from({ length: system.layout.nx }, (_, k) => ({ body, side, i: k + 1,
        pressure: block.boundaryPressure(body, k + 1, side), velocity: block.boundaryEdgeVelocity(body, k + 1, side),
        distanceWeightedVelocity: block.boundaryEdgeVelocity(body, k + 1, side, 'distance-weighted') })))),
    outlet: system.layout.bodies.flatMap((_, body) => ['lower', 'upper'].map(side => block.boundaryOutletTangency(body, side))) };
}
function countEvaluations(system) {
  const evaluate = system.evaluate; let count = 0;
  system.evaluate = (...args) => { count++; return evaluate(...args); };
  return { get count() { return count; }, reset() { count = 0; } };
}
const modes = [
  ['momentum', { streamwiseMode: 'momentum' }],
  ['isentropic', { streamwiseMode: 'isentropic' }],
  ...[1, 2, 3, 4].map(ismom => [`ISMOM${ismom}`, { streamwiseMode: 'hybrid', hybrid: { ismom, epsilonP: 1e-3 },
    upwind: { mucon: 1, mcrit: .1, boundary: { kind: 'unfiltered-first-two' } } }]),
  ['incompressible', { mach: 0, flowModel: 'incompressible' }],
];

test('returned evaluation data cannot contaminate a later atomic Jacobian call', () => {
  const input = intrinsicBodyFixture({ bodySegments: 4, tubes: 2 });
  const current = createStreamtubeBodySystem(input), old = archivedBody.createStreamtubeBodySystem(input), x = current.initial;
  const first = current.evaluateJacobian(x, { sparse: true });
  first.value.residual.fill(NaN); first.value.nodes[0][0][0].x = NaN;
  const second = current.evaluateJacobian(x, { sparse: true });
  assert.deepEqual(snapshot(second.value), snapshot(old.evaluate(x)));
  assert.deepEqual(second.jacobian, old.jacobian(x, { sparse: true }));
});

for (const [label, controls] of modes) test(`${label}: atomic Euler evaluation is identical to archived assembly across mutable state, thickness and chart changes`, () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2, mach: .2, alpha: .25 }), ...controls };
  input.displacement = thickness(input);
  const current = createStreamtubeBodySystem(structuredClone(input)), old = archivedBody.createStreamtubeBodySystem(structuredClone(input));
  const counts = countEvaluations(current), oldCounts = countEvaluations(old);
  let x = current.initial.map((_, i) => (i < current.layout.densityCount ? 1e-3 : 1e-5) * Math.sin(i + .3));
  const check = () => {
    const originalState = x.slice(), chart = current.geometryChart(), displacement = structuredClone(current.displacement);
    counts.reset(); oldCounts.reset();
    const expectedValue = old.evaluate(x), expectedJacobian = old.jacobian(x, { sparse: true, includeDisplacement: true });
    const pair = current.evaluateJacobian(x, { sparse: true, includeDisplacement: true });
    assert.equal(counts.count, 1, 'Atomic operation must evaluate Euler exactly once.');
    assert.equal(oldCounts.count, 2, 'Archived reference includes the duplicate evaluation.');
    assert.deepEqual(snapshot(pair.value), snapshot(expectedValue));
    assert.deepEqual(blockSnapshot(pair.jacobian, current), blockSnapshot(expectedJacobian, old));
    assert.deepEqual(x, originalState); assert.deepEqual(current.geometryChart(), chart); assert.deepEqual(current.displacement, displacement);
    // Existing public methods retain their return shape and exact derivative values.
    assert.deepEqual(current.jacobian(x), old.jacobian(x));
    assert.equal(Object.hasOwn(current.evaluate(x), 'jacobian'), false);
    return snapshot(pair.value);
  };
  const initial = check();
  x[0] += 2e-5; // Same typed-array identity, changed numerical state.
  const stateChanged = check(); assert.notDeepEqual(stateChanged.residual, initial.residual);
  const d = structuredClone(input.displacement); d.surfaces[0].upper[1] *= 1.2; d.wakes[1][0] *= 1.1;
  current.setDisplacement(d); old.setDisplacement(d);
  const thicknessChanged = check(); assert.notDeepEqual(thicknessChanged.nodes, stateChanged.nodes);
  const physicalNodes = current.decode(x).nodes, next = current.adoptGeometry(x, physicalNodes), nextOld = old.adoptGeometry(x, physicalNodes);
  assert.deepEqual(next, nextOld); assert.deepEqual(current.geometryChart(), old.geometryChart()); x = next; check();
  const invalid = x.slice(); invalid[0] = NaN;
  assert.throws(() => current.evaluateJacobian(invalid), /Invalid intrinsic body state/);
});

for (const mode of ['momentum', 'isentropic', 'hybrid']) test(`${mode}: coupled assembly drops only the duplicate Euler evaluation and retains complete BL/edge/phase rows`, () => {
  const f = structuredClone(fixture.restart);
  f.input.streamwiseMode = mode;
  if (mode === 'hybrid') {
    f.input.streamwiseMode = 'hybrid'; f.input.hybrid = { ismom: 4, epsilonP: 1e-5 };
    f.input.upwind = { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } };
    f.options.blThermodynamics = 'historical-common-isentrope';
  }
  const build = create => create(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: Float64Array.from(f.initialBL) });
  const current = build(createCoupledStreamtubeBody), old = build(archivedCoupled.createCoupledStreamtubeBody);
  assert.equal(current.n, 349);
  const counts = countEvaluations(current.euler), oldCounts = countEvaluations(old.euler);
  let x = current.initial.slice(); assert.deepEqual(x, old.initial);
  const check = () => {
    const original = x.slice(), phase = current.bl.snapshotActive(), chart = current.euler.geometryChart();
    counts.reset(); oldCounts.reset();
    const expected = old.jacobian(x), matrix = current.jacobian(x);
    assert.equal(counts.count, 1); assert.equal(oldCounts.count, 2);
    assert.deepEqual(matrix, expected);
    assert.deepEqual(x, original); assert.deepEqual(current.bl.snapshotActive(), phase); assert.deepEqual(current.euler.geometryChart(), chart);
    const a = current.evaluate(x, { jacobian: true }), b = old.evaluate(x, { jacobian: true });
    assert.deepEqual(Object.keys(a), Object.keys(b)); assert.equal(Object.hasOwn(a, 'eulerJacobian'), false);
    assert.deepEqual(a.residual, b.residual); assert.deepEqual(a.layers, b.layers); assert.deepEqual(a.families, b.families);
    assert.deepEqual(snapshot(a.outer), snapshot(b.outer));
    return matrix;
  };
  const first = check();
  x[0] += 1e-6; x[current.ne + 2] *= 1.001;
  assert.notDeepEqual(check().values, first.values);
  const next = current.rebase(x), nextOld = old.rebase(x); assert.deepEqual(next, nextOld); x = next; check();
  const phase = current.bl.snapshotActive(), changed = phase.slice(), changingSurface = phase.findIndex(p => p > 1);
  assert.ok(changingSurface >= 0); changed[changingSurface]--; assert.notDeepEqual(changed, phase);
  current.bl.restoreActive(changed); old.bl.restoreActive(changed);
  // Relabelling alone leaves this saved physical profile outside the new
  // transition interval. Both implementations must reject it, then assemble
  // the restored phase from current data rather than a previous cached block.
  assert.throws(() => old.jacobian(x), /Transition is outside this active interval/);
  assert.throws(() => current.jacobian(x), /Transition is outside this active interval/);
  assert.deepEqual(current.bl.snapshotActive(), changed); assert.deepEqual(old.bl.snapshotActive(), changed);
  current.bl.restoreActive(phase); old.bl.restoreActive(phase); check();
  const invalid = x.slice(); invalid[current.ne + 1] = -1;
  counts.reset(); oldCounts.reset();
  assert.throws(() => current.jacobian(invalid), /Inadmissible coupled streamtube state/);
  assert.throws(() => old.jacobian(invalid), /Inadmissible coupled streamtube state/);
  assert.equal(counts.count, 0); assert.equal(oldCounts.count, 0);
});
