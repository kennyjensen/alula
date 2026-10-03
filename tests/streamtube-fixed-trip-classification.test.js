// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { evaluateFixedTripInterval } from '../src/euler/streamtube-boundary-layers.js';

const capture = JSON.parse(fs.readFileSync('docs/rae2822/bl-transition-audit.json', 'utf8')).comparisons[0];
const fresh = () => createIntegralKernel(capture.kernelParameters);

test('the frozen RAE selected trip corrects its Boolean label with every residual and legacy derivative unchanged', () => {
  const input = structuredClone(capture.input), before = structuredClone(input), kernel = fresh();
  const raw = kernel.interval(input);
  assert.deepEqual(raw, capture.raw);
  assert.equal(raw.transition.s, input.tripS); assert.equal(raw.transition.forced, false);
  const corrected = evaluateFixedTripInterval(kernel, input);
  assert.equal(corrected.transition.forced, true); assert.equal(corrected.transition.nativeForced, false);
  assert.equal(corrected.transition.s, input.tripS);
  assert.deepEqual({ ...corrected, transition: raw.transition }, raw,
    'all three rows, upstream/downstream Jacobians, parameter derivatives and gas properties remain exact');
  assert.deepEqual(input, before);
});

test('a later trip still rejects actual earlier natural transition', () => {
  const input = structuredClone(capture.input), kernel = fresh();
  const natural = kernel.transitionCheck({ ...input, tripS: Number.MAX_VALUE });
  assert.equal(natural.transition, true);
  input.tripS = .5 * (natural.s + input.downstream.s);
  assert.ok(natural.s < input.tripS && input.tripS < input.downstream.s);
  const raw = kernel.interval(input);
  assert.equal(raw.transition.forced, false); assert.ok(raw.transition.s < input.tripS);
  assert.throws(() => evaluateFixedTripInterval(kernel, input), /Natural transition precedes the fixed-trip research interval/);
});

test('an already forced mixed interval preserves its entire original result without reselection', () => {
  // Analytic control: lower upstream amplification while preserving the
  // physical station geometry and Ncrit. This is not a changed RAE solve.
  const input = structuredClone(capture.input), kernel = fresh();
  input.upstream.aux = 0;
  const raw = kernel.interval(input);
  assert.equal(raw.transition.forced, true);
  const interval = kernel.interval; let calls = 0;
  kernel.interval = data => { calls++; return interval(data); };
  kernel.transitionCheck = () => { throw new Error('Already forced intervals must not be reselected.'); };
  assert.deepEqual(evaluateFixedTripInterval(kernel, input), raw);
  assert.equal(calls, 1);
});

test('ordinary laminar blocks preserve their entire original result and bypass transition selection', () => {
  const input = { ...structuredClone(capture.input), regime: 'laminar' }, kernel = fresh();
  input.downstream.aux = input.upstream.aux;
  const raw = kernel.interval(input);
  assert.equal(raw.transition, null);
  kernel.transitionCheck = () => { throw new Error('Laminar blocks must not select transition here.'); };
  assert.deepEqual(evaluateFixedTripInterval(kernel, input), raw);
});
