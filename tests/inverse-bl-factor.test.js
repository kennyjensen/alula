import test from 'node:test';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { prepareCoupledInverseBL } from '../scripts/validation/coupled-inverse-bl.js';
import { factorInverseBL } from '../scripts/validation/inverse-bl-factor.js';
import { solveLinear } from '../src/numerics/linear.js';
import { sparseProduct } from '../src/numerics/sparse.js';

test('station inverse factor agrees with a dense solve of every BL equation on both elements and wakes', t => {
  const s = createCoupledStreamtubeBody(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    { transitionMode: 'automatic', edgeMatching: 'section-velocity', reynolds: 1e6, ncrit: 9 });
  const p = prepareCoupledInverseBL(s, s.initial); s.bl.restoreActive(p.phase);
  const matrix = s.jacobian(p.x), dense = s.jacobian(p.x, { sparse: false }), factor = factorInverseBL(s, matrix);
  const columns = factor.blocks.flatMap(b => b.columns), rows = factor.blocks.flatMap(b => b.rows);
  const a = Float64Array.from(rows.flatMap(row => columns.map(col => dense[row * s.n + col])));
  let maxError = 0, maxResidual = 0;
  for (const seed of [.3, 1.7]) {
    const retained = p.x.map((v, i) => factor.owner[i] >= 0 ? 0 : .003 * Math.sin(i * .37 + seed));
    const rhs = p.x.map((v, i) => .002 * Math.sin(i * .47 + seed)), product = sparseProduct(matrix, retained);
    const reference = solveLinear(a, Float64Array.from(rows, row => rhs[row] - product[row]));
    const x = factor.lift(retained, rhs), actual = sparseProduct(matrix, x);
    columns.forEach((col, k) => { maxError = Math.max(maxError, Math.abs(x[col] - reference[k])); });
    rows.forEach(row => { maxResidual = Math.max(maxResidual, Math.abs(actual[row] - rhs[row])); });
    factor.retained.forEach(col => { assert.equal(x[col], retained[col]); });
  }
  assert.ok(maxError < 1e-10 && maxResidual < 1e-10);
  assert.equal(columns.length, 93); assert.equal(s.bl.wakes.length, 2);
  t.diagnostic(JSON.stringify({ blocks: factor.blocks.length, maxError, maxResidual }));
});
