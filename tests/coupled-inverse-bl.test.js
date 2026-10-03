import test from 'node:test';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { prepareCoupledInverseBL } from '../scripts/validation/coupled-inverse-bl.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { solveLinear } from '../src/numerics/linear.js';
import { sparseProduct } from '../src/numerics/sparse.js';

const make = extra => createCoupledStreamtubeBody(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
  { reynolds: 1e6, ncrit: 9, transitionMode: 'automatic', edgeMatching: 'section-velocity', ...extra });

test('conditional inverse BL closes all four surfaces, TE matches and two wakes at unchanged displacement', () => {
  const s = make(), before = s.evaluate(s.initial), phase = s.bl.snapshotActive(), copy = s.initial.slice();
  const p = prepareCoupledInverseBL(s, s.initial);
  assert.deepEqual(s.initial, copy); assert.deepEqual(s.bl.snapshotActive(), phase);
  assert.deepEqual(p.x.slice(0, s.ne), copy.slice(0, s.ne));
  for (const row of s.bl.thicknessMap) for (const col of row.keys()) assert.equal(p.x[col], copy[col]);
  s.bl.restoreActive(p.phase); const after = s.evaluate(p.x);
  assert.ok(after.families.boundaryLayer < 3e-11); assert.equal(after.families.euler, before.families.euler);
  assert.deepEqual(after.outer.nodes, before.outer.nodes);
  assert.equal(p.summary.surfaces, 4); assert.equal(p.summary.wakes, 2);
  assert.equal(p.intervals.length, s.bl.stations.length);
  const again = prepareCoupledInverseBL(s, p.x);
  assert.equal(again.summary.updates, 0); assert.deepEqual(again.phase, p.phase);
});

test('complete Jacobian includes the implicit inverse BL response to retained grid/displacement variables', t => {
  const s = make(), prepared = prepareCoupledInverseBL(s, s.initial), x = prepared.x;
  s.bl.restoreActive(prepared.phase);
  const te = new Set(s.bl.wakes.map(w => w.ids[0]));
  const aux = s.bl.stations.flatMap(p => (te.has(p.id) ? [0, 1, 2] : [0, 1, 3]).map(k => s.ne + 4 * p.id + k));
  const rows = s.bl.stations.flatMap(p => [0, 1, 2].map(k => s.ne + 4 * p.id + k)), eliminated = new Set(aux);
  const dense = s.jacobian(x, { sparse: false }), matrix = s.jacobian(x);
  const a = Float64Array.from(rows.flatMap(r => aux.map(c => dense[r * s.n + c]))), checks = [];
  for (const seed of [.3, 1.7]) {
    const direction = x.map((v, i) => eliminated.has(i) ? 0 : i < s.ne ? .0001 * Math.sin(.7 * i + seed)
      : .005 * Math.max(.01, Math.abs(v)) * Math.sin(.4 * i + seed));
    const jb = sparseProduct(matrix, direction), da = solveLinear(a, Float64Array.from(rows, i => -jb[i]));
    const tangent = direction.slice(); aux.forEach((i, j) => { tangent[i] = da[j]; });
    const expected = sparseProduct(matrix, tangent);
    for (const h of [1e-4, 5e-5]) {
      const plus = prepareCoupledInverseBL(s, x.map((v, i) => v + h * direction[i]));
      const minus = prepareCoupledInverseBL(s, x.map((v, i) => v - h * direction[i]));
      assert.deepEqual(plus.phase, prepared.phase); assert.deepEqual(minus.phase, prepared.phase);
      const p = s.residual(plus.x), m = s.residual(minus.x); let stateError = 0, residualError = 0;
      for (let i = 0; i < s.n; i++) {
        const fdState = (plus.x[i] - minus.x[i]) / (2 * h), fdResidual = (p[i] - m[i]) / (2 * h);
        stateError = Math.max(stateError, Math.abs(fdState - tangent[i]) / Math.max(1, Math.abs(fdState), Math.abs(tangent[i])));
        residualError = Math.max(residualError, Math.abs(fdResidual - expected[i]) / Math.max(1, Math.abs(fdResidual), Math.abs(expected[i])));
      }
      checks.push({ seed, h, stateError, residualError });
      assert.ok(stateError < 5e-6, JSON.stringify(checks)); assert.ok(residualError < 5e-6, JSON.stringify(checks));
    }
  }
  t.diagnostic(JSON.stringify({ unknowns: s.n, eliminated: aux.length, checks }));
});

test('conditional inverse BL handles a leading trip and laminar surfaces to the trailing edge', () => {
  const s = make({ ncrit: 14, reynolds: 3e5, tripFractions: [[.05, 1], [1, 1]] });
  const p = prepareCoupledInverseBL(s, s.initial); s.bl.restoreActive(p.phase);
  assert.equal(p.phase[0], 0); assert.ok(s.bl.surfaces.some(p => p.transition === p.ids.length - 1));
  assert.ok(s.evaluate(p.x).families.boundaryLayer < 3e-11);
});

test('failed conditional inverse preparation leaves the supplied full state and active intervals intact', () => {
  const s = make(), x = s.initial.slice(), phase = s.bl.snapshotActive();
  x[s.ne + 1] *= .8; const copy = x.slice();
  assert.throws(() => prepareCoupledInverseBL(s, x, { maxIterations: 0 }), /Conditional inverse BL failed/);
  assert.deepEqual(x, copy); assert.deepEqual(s.bl.snapshotActive(), phase);
});
