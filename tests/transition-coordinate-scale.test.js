import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { evaluateTransitionInterval } from '../src/viscous/transition-interval.js';

const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/terminal-trip-coordinate-scale.json', import.meta.url)));
const fields = ['aux', 'theta', 'deltaStar', 'ue', 's'];
const error = (a, b) => Math.abs(a - b) / Math.max(1, Math.abs(a), Math.abs(b));
const close = (a, b, limit = 2e-6, label = '') => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && error(a, b) < limit, `${label}: ${a} != ${b}; error=${error(a, b)}, limit=${limit}`);
const kernel = () => createIntegralKernel(fixture.parameters);
const packed = r => [...r.residual, r.transition.s];
const value = (k, input) => {
  const result = evaluateTransitionInterval(k, input);
  assert.equal(result.transition.forced, true);
  return packed(result);
};

// Independent five-point fourth-order formulas. The reference step is larger
// than the runtime local step and uses no production derivative helper.
function fourthOrder(f, x, h, sign = 0) {
  if (!sign) {
    const p = f(x + h), pp = f(x + 2 * h), m = f(x - h), mm = f(x - 2 * h);
    return p.map((v, i) => (8 * (v - m[i]) - (pp[i] - mm[i])) / (12 * h));
  }
  const b = f(x), samples = [1, 2, 3, 4].map(j => f(x + sign * j * h));
  return b.map((v, i) => (48 * (samples[0][i] - v) - 36 * (samples[1][i] - v)
    + 16 * (samples[2][i] - v) - 3 * (samples[3][i] - v)) / (12 * sign * h));
}

test('terminal-trip fixture preserves the exact retained interval and pre-change residual', () => {
  const bytes = fs.readFileSync(new URL(`../${fixture.provenance.source}`, import.meta.url));
  assert.equal(createHash('sha256').update(bytes).digest('hex'), fixture.provenance.sha256);
  const original = JSON.parse(bytes).fixture;
  assert.deepEqual(fixture.input, original.input); assert.deepEqual(fixture.parameters, original.parameters);
  assert.deepEqual(fixture.expected.residual, original.base.slice(0, 3));
  const actual = evaluateTransitionInterval(kernel(), fixture.input);
  assert.deepEqual(actual.residual, fixture.expected.residual);
  assert.equal(actual.transition.s, fixture.expected.transitionS); assert.equal(actual.transition.forced, fixture.expected.forced);
});

test('terminal-trip individual partials match independent fourth-order differences on the short interval', t => {
  const input = fixture.input, k = kernel(), result = evaluateTransitionInterval(k, input, { jacobian: true });
  const span = input.downstream.s - input.upstream.s; let worst = { error: 0 };
  const compare = (actual, reference, label) => actual.forEach((v, row) => {
    const e = error(v, reference[row]); if (e > worst.error) worst = { error: e, label, row };
    close(v, reference[row], 2e-6, `${label}[${row}]`);
  });
  for (const side of ['upstream', 'downstream']) for (const [column, key] of fields.entries()) {
    const h = 1e-4 * (key === 's' ? span : Math.max(Math.abs(input[side][key]), key === 'aux' ? .01 : 1e-7));
    const f = v => value(k, { ...input, [side]: { ...input[side], [key]: v } });
    const reference = fourthOrder(f, input[side][key], h, side === 'downstream' && key === 's' ? 1 : 0);
    const actual = [...result.partials[side].map(row => row[column]), result.partials.location[side][column]];
    compare(actual, reference, `${side}.${key}`);
  }
  compare([...result.partials.trip, result.partials.location.trip],
    fourthOrder(tripS => value(k, { ...input, tripS }), input.tripS, 1e-4 * span, -1), 'tripS');
  t.diagnostic(JSON.stringify({ span, worst }));
});

test('terminal shear row is invariant when the trip moves with its downstream station', t => {
  const input = fixture.input, k = kernel(), before = structuredClone(input);
  const result = evaluateTransitionInterval(k, input, { jacobian: true }), span = input.downstream.s - input.upstream.s;
  // At zero turbulent length the shear condition is local to the TE state.
  // This metamorphic check does not use a derivative step or fitted reference.
  for (const fraction of [-.1, -.01, .01, .1]) {
    const s = input.downstream.s + fraction * span;
    const moved = evaluateTransitionInterval(k, { ...input, downstream: { ...input.downstream, s }, tripS: s });
    assert.equal(moved.transition.s, s); close(moved.residual[0], result.residual[0], 2e-12, `motion ${fraction}`);
  }
  const sum = result.partials.downstream[0][4] + result.partials.trip[0];
  close(sum, 0, 2e-6, 'downstream.s + tripS shear derivative');
  close(result.partials.upstream[0][4], 0, 2e-6, 'upstream.s shear derivative');
  close(result.partials.location.downstream[4] + result.partials.location.trip, 1, 2e-6, 'moving terminal location');
  assert.deepEqual(input, before); t.diagnostic(JSON.stringify({ terminalShearCoordinateSum: sum }));
});

test('an interior forced trip retains its real location and shear dependence', t => {
  const span = fixture.input.downstream.s - fixture.input.upstream.s;
  const input = { ...fixture.input, tripS: fixture.input.upstream.s + .6 * span }, k = kernel();
  const result = evaluateTransitionInterval(k, input, { jacobian: true });
  assert.equal(result.transition.forced, true); close(result.transition.s, input.tripS, 2e-12);
  const reference = fourthOrder(tripS => value(k, { ...input, tripS }), input.tripS, 1e-4 * span);
  const actual = [...result.partials.trip, result.partials.location.trip];
  actual.forEach((v, row) => close(v, reference[row], 2e-6, `interior trip row${row}`));
  assert.ok(Math.abs(actual[0]) > .1, `Missing physical trip dependence: ${actual[0]}`);
  close(actual[3], 1, 2e-6, 'interior forced-trip location');
  const moved = evaluateTransitionInterval(k, { ...input, tripS: input.tripS + .01 * span });
  assert.ok(Math.abs(moved.residual[0] - result.residual[0]) > 1e-8);
  t.diagnostic(JSON.stringify({ tripDerivative: actual }));
});
