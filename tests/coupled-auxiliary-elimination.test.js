import test from 'node:test';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { prepareCoupledAuxiliaries } from '../scripts/validation/coupled-auxiliary-elimination.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { solveLinear } from '../src/numerics/linear.js';
import { sparseProduct } from '../src/numerics/sparse.js';

const make = extra => createCoupledStreamtubeBody(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
  { reynolds: 1e6, ncrit: 9, transitionMode: 'automatic', edgeMatching: 'section-velocity', ...extra });
const maximum = values => Math.max(0, ...Array.from(values, Math.abs));

test('auxiliary elimination closes four surface and two wake transport systems without moving geometry or thickness', () => {
  const s = make(), phase = s.bl.snapshotActive(), x = s.initial.slice();
  for (const p of s.bl.stations) x[s.ne + 4 * p.id] *= ['laminar', 'similarity'].includes(p.regime) ? .8 : .6;
  const copy = x.slice(), before = s.evaluate(x), r = prepareCoupledAuxiliaries(s, x), after = s.evaluate(r.x);
  assert.deepEqual(x, copy); assert.deepEqual(s.bl.snapshotActive(), phase);
  assert.ok(maximum(s.bl.stations.map(p => after.residual[s.ne + 4 * p.id])) < 2e-12);
  r.x.forEach((v, i) => { if (i < s.ne || (i - s.ne) % 4 !== 0) assert.equal(v, x[i]); });
  assert.deepEqual(after.outer.nodes, before.outer.nodes);
  assert.equal(after.families.euler, before.families.euler); assert.equal(after.families.edgeMatching, before.families.edgeMatching);
  for (const w of s.bl.wakes) {
    const [a, b] = s.bl.surfaces.filter(p => p.body === w.body).map(p => after.layers.states[p.ids.at(-1)]), z = after.layers.states[w.ids[0]];
    assert.ok(Math.abs(z.aux * (a.theta + b.theta) - a.aux * a.theta - b.aux * b.theta) < 1e-14);
  }
  assert.deepEqual(prepareCoupledAuxiliaries(s, r.x).x, r.x);
});

test('eliminated auxiliaries and residuals have the implicit derivatives of the complete coupled Jacobian', t => {
  const s = make(), x = prepareCoupledAuxiliaries(s, s.initial).x, aux = s.bl.stations.map(p => s.ne + 4 * p.id);
  const matrix = s.jacobian(x), dense = s.jacobian(x, { sparse: false });
  const a = Float64Array.from(aux.flatMap(row => aux.map(col => dense[row * s.n + col])));
  const direction = x.map((v, i) => i < s.ne ? .001 * Math.sin(.7 * i + .3)
    : (i - s.ne) % 4 === 0 ? 0 : .02 * Math.max(.01, Math.abs(v)) * Math.sin(.4 * i + .1));
  const product = sparseProduct(matrix, direction), da = solveLinear(a, Float64Array.from(aux, i => -product[i]));
  const tangent = direction.slice(); aux.forEach((i, j) => { tangent[i] = da[j]; });
  const expected = sparseProduct(matrix, tangent), checks = [];
  for (const h of [1e-4, 5e-5]) {
    const plus = prepareCoupledAuxiliaries(s, x.map((v, i) => v + h * direction[i])).x;
    const minus = prepareCoupledAuxiliaries(s, x.map((v, i) => v - h * direction[i])).x;
    const p = s.residual(plus), m = s.residual(minus); let stateError = 0, residualError = 0;
    for (let i = 0; i < s.n; i++) {
      const fdState = (plus[i] - minus[i]) / (2 * h), fdResidual = (p[i] - m[i]) / (2 * h);
      stateError = Math.max(stateError, Math.abs(fdState - tangent[i]) / Math.max(1, Math.abs(fdState), Math.abs(tangent[i])));
      residualError = Math.max(residualError, Math.abs(fdResidual - expected[i]) / Math.max(1, Math.abs(fdResidual), Math.abs(expected[i])));
    }
    checks.push({ h, stateError, residualError });
    assert.ok(stateError < 5e-6, JSON.stringify(checks)); assert.ok(residualError < 5e-6, JSON.stringify(checks));
  }
  t.diagnostic(JSON.stringify({ unknowns: s.n, eliminated: aux.length, checks }));
});

test('leading trips and laminar-to-TE surfaces retain their existing equations', () => {
  // A material trip must lie downstream of the displaced stagnation point.
  // Five percent is before this coarse fixture's first station but remains
  // on the upper branch after its multielement inviscid precursor.
  const s = make({ ncrit: 14, reynolds: 3e5, tripFractions: [[.05, 1], [1, 1]] });
  assert.equal(s.bl.surfaces[0].transition, 0);
  const r = prepareCoupledAuxiliaries(s, s.initial), value = s.evaluate(r.x);
  assert.ok(maximum(s.bl.stations.map(p => value.residual[s.ne + 4 * p.id])) < 2e-12);
  assert.ok(s.bl.surfaces.some(p => p.transition === p.ids.length - 1));
});

test('an unprepared interval change cannot partially modify the caller state or phases', () => {
  const s = make(), phase = s.bl.snapshotActive(), wrong = phase.slice(); wrong[wrong.length - 1]--;
  s.bl.restoreActive(wrong); const before = s.initial.slice();
  assert.throws(() => prepareCoupledAuxiliaries(s, s.initial), /Prepare the transition interval/);
  assert.deepEqual(s.initial, before); assert.deepEqual(s.bl.snapshotActive(), wrong);
});
