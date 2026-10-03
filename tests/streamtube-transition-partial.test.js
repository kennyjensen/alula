import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { boundedResidualPartial } from '../src/numerics/bounded-partial.js';

test('the retained near-node forced trip has finite accurate derivatives inside its active interval', t => {
  const { input } = JSON.parse(fs.readFileSync(new URL('fixtures/near-node-forced-trip.json', import.meta.url)));
  const kernel = createIntegralKernel({ reynolds: 1e6, mach: .2, ncrit: 9, exactJacobian: true });
  const base = kernel.interval(input), s = input.downstream.s, h = Math.cbrt(Number.EPSILON) * s;
  assert.throws(() => kernel.interval({ ...input, downstream: { ...input.downstream, s: s - h } }), /outside this active interval/);
  const f = position => {
    const r = kernel.interval({ ...input, downstream: { ...input.downstream, s: position } });
    assert.equal(r.transition.forced, true); return r.residual;
  };
  const d = boundedResidualPartial(f, s, { step: h, lower: input.tripS, base: base.residual });
  // Independent fourth-order forward formula; all samples lie inside the
  // forced-trip branch. No global solve or use of the production derivative.
  const values = Array.from({ length: 5 }, (_, k) => f(s + k * h));
  const reference = values[0].map((v, i) => (-25 * v + 48 * values[1][i] - 36 * values[2][i] + 16 * values[3][i] - 3 * values[4][i]) / (12 * h));
  const errors = d.map((v, i) => Math.abs(v - reference[i]) / Math.max(1, Math.abs(v), Math.abs(reference[i])));
  assert.ok(errors.every(e => e < 2e-6), JSON.stringify({ d, reference, errors }));
  assert.equal(kernel.interval(input).transition.forced, true);
  t.diagnostic(JSON.stringify({ gap: s - input.tripS, originalStep: h, errors }));
});
