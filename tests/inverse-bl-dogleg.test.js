import test from 'node:test';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { prepareCoupledInverseBL } from '../scripts/validation/coupled-inverse-bl.js';
import { createInverseBLDogleg } from '../scripts/validation/inverse-bl-dogleg.js';
import { solveLinear } from '../src/numerics/linear.js';
import { sparseDense, sparseProduct } from '../src/numerics/sparse.js';

test('inverse dogleg model closes both elements BL rows and predicts the complete residual, with unit model agreement', t => {
  const s = createCoupledStreamtubeBody(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    { transitionMode: 'automatic', edgeMatching: 'section-velocity', reynolds: 1e6, ncrit: 9 });
  const prepared = prepareCoupledInverseBL(s, s.initial); s.bl.restoreActive(prepared.phase);
  const matrix = s.jacobian(prepared.x), residual = s.residual(prepared.x);
  const newton = solveLinear(sparseDense(matrix), residual.map(v => -v));
  const trust = createInverseBLDogleg(s, matrix, residual, newton);
  const norms = [trust.model.cauchyNorm / 2, (trust.model.cauchyNorm + trust.model.newtonNorm) / 2, trust.model.newtonNorm * 1.01];
  let maximumBL = 0, maximumMeritError = 0;
  for (const radius of norms) {
    const p = trust.propose(radius), jd = sparseProduct(matrix, p.direction), next = residual.map((v, i) => v + jd[i]);
    for (const station of s.bl.stations) for (const eq of [0, 1, 2])
      maximumBL = Math.max(maximumBL, Math.abs(next[s.ne + 4 * station.id + eq]));
    const a = trust.assess(p.reducedDirection, next, radius);
    assert.equal(a.accepted, true); assert.ok(a.scaledStepNorm <= radius * (1 + 1e-12));
    maximumMeritError = Math.max(maximumMeritError, Math.abs(a.ratio - 1));
    const reject = trust.assess(p.reducedDirection, residual.map(v => v * 2), radius);
    assert.equal(reject.accepted, false); assert.equal(reject.nextRadius, .25 * radius);
    const weak = trust.assess(p.reducedDirection, residual.map((v, i) => v + .01 * jd[i]), radius);
    assert.equal(weak.accepted, true); assert.ok(weak.ratio < .25); assert.equal(weak.nextRadius, .25 * radius);
  }
  assert.ok(maximumBL < 1e-10 && maximumMeritError < 1e-12);
  assert.equal(s.bl.surfaces.length, 4); assert.equal(s.bl.wakes.length, 2);
  t.diagnostic(JSON.stringify({ unknowns: s.n, maximumBL, maximumMeritError }));
});
