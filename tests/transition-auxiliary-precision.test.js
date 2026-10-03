// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { evaluateTransitionInterval } from '../src/viscous/transition-interval.js';

const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/transition-auxiliary-precision.json', import.meta.url)));
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const relative = (a, b) => Math.abs(a - b) / Math.max(1, Math.abs(a), Math.abs(b));
const close = (a, b, tolerance = 1e-9, label = '') => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && relative(a, b) < tolerance, `${label}: ${a} != ${b}, error=${relative(a, b)}`);
const kernel = () => createIntegralKernel(fixture.parameters);
const pack = r => [...r.residual, r.transition.s];
const column = (r, side) => [...r.partials[side].map(row => row[0]), r.partials.location[side][0]];

// Independent sixth-order central difference. For shear use log coordinates
// so all samples remain positive, including coefficients much smaller than
// the former absolute finite-difference scale of .01.
function sixthOrder(f, h) {
  const samples = [-3, -2, -1, 1, 2, 3].map(offset => f(offset * h));
  const coefficients = [-1, 9, -45, 45, -9, 1];
  return samples[0].map((_, row) => samples.reduce((s, v, j) => s + coefficients[j] * v[row], 0) / (60 * h));
}

test('failed-browser transition fixture binds original Fortran with only root tolerance changed', () => {
  for (const [path, expected] of Object.entries(fixture.provenance.hashes))
    assert.equal(hash(fs.readFileSync(new URL('../' + path, import.meta.url))), expected, path);
  const original = fs.readFileSync(new URL('../third_party/Xfoil/src/xblsys.f', import.meta.url), 'utf8');
  const { from, to, sourceSha256 } = fixture.provenance.patch;
  assert.equal(original.split(from).length, 2);
  assert.equal(hash(original.replace(from, to)), sourceSha256);
  assert.match(fixture.provenance.scope, /only amplification stopping tolerance/);
  assert.equal(fixture.input.downstream.id, 161);
  const before = structuredClone(fixture.input), k = kernel();
  for (const jacobian of [false, true]) {
    const r = evaluateTransitionInterval(k, fixture.input, { jacobian });
    assert.deepEqual(r.residual, fixture.retained.residual);
    assert.deepEqual(r.transition, fixture.retained.transition);
    r.residual.forEach((v, i) => close(v, fixture.expected.residual[i], 1e-14));
  }
  assert.deepEqual(fixture.input, before);
});

test('near-endpoint natural transition auxiliary partials match two independent Fortran steps', t => {
  const r = evaluateTransitionInterval(kernel(), fixture.input, { jacobian: true });
  let maximum = 0;
  for (const reference of fixture.derivatives) {
    const actual = column(r, reference.side);
    actual.forEach((v, row) => {
      maximum = Math.max(maximum, relative(v, reference.derivative[row]));
      close(v, reference.derivative[row], 1e-9, `${reference.side}[${row}], step ${reference.step}`);
    });
  }
  for (let row = 0; row < 3; row++)
    close(r.partials.downstream[row][0], fixture.expected.downstream[row][0], 2e-14, `native Ctau row ${row}`);
  assert.equal(r.partials.location.downstream[0], 0);
  t.diagnostic(JSON.stringify({ maximum, upstream: column(r, 'upstream'), downstream: column(r, 'downstream') }));
});

test('turbulent auxiliary has no effect on natural or forced transition position across six decades', t => {
  let maximum = 0;
  for (const tripS of [Number.MAX_VALUE, .21, fixture.input.downstream.s]) {
    const k = kernel(), original = evaluateTransitionInterval(k, { ...fixture.input, tripS });
    for (const aux of [1e-8, 1e-6, fixture.input.downstream.aux, .01, .05]) {
      const input = { ...fixture.input, tripS, downstream: { ...fixture.input.downstream, aux } };
      const before = structuredClone(input), r = evaluateTransitionInterval(k, input, { jacobian: true });
      assert.equal(r.transition.s, original.transition.s);
      assert.equal(r.transition.amplification, original.transition.amplification);
      assert.equal(r.transition.forced, original.transition.forced);
      assert.equal(r.partials.location.downstream[0], 0);
      for (const h of [1e-3, 5e-4]) {
        const reference = sixthOrder(offset => {
          const sampled = evaluateTransitionInterval(k, { ...input, downstream: { ...input.downstream, aux: aux * Math.exp(offset) } });
          assert.equal(sampled.transition.s, original.transition.s);
          assert.equal(sampled.transition.forced, original.transition.forced);
          return pack(sampled);
        }, h);
        column(r, 'downstream').forEach((v, row) => {
          const error = relative(v * aux, reference[row]); maximum = Math.max(maximum, error);
          close(v * aux, reference[row], 2e-10, `${tripS}, Ctau=${aux}, row=${row}`);
        });
      }
      assert.deepEqual(input, before);
    }
  }
  t.diagnostic(JSON.stringify({ maximum }));
});

test('natural/trip switch keeps the fourth-order amplification derivative on its selected side', () => {
  const k = kernel(), base = evaluateTransitionInterval(k, fixture.input, { jacobian: true });
  for (const direction of [-1, 1]) {
    const input = { ...fixture.input, tripS: base.transition.s + direction * 1e-10 };
    const r = evaluateTransitionInterval(k, input, { jacobian: true });
    assert.equal(r.transition.forced, direction < 0);
    if (direction > 0) column(r, 'upstream').forEach((v, row) => close(v, column(base, 'upstream')[row], 5e-9));
    else {
      // A pinned forced root and its interpolated laminar/turbulent state
      // have no dependence on the upstream amplification inside this branch.
      column(r, 'upstream').forEach(v => close(v, 0, 5e-9));
    }
  }
});
