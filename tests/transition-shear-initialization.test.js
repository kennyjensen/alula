import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { evaluateTransitionInterval } from '../src/viscous/transition-interval.js';
import { initializeMixedTransitionShear } from '../src/viscous/transition-shear-initialization.js';
const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/automatic-transition.json', import.meta.url)));

test('mixed shear initialization closes natural and interior-forced rows from tiny and large stale shear', () => {
  let changedEnergy = false;
  for (const { input } of fixture.cases) for (const seed of [1e-10, .3]) {
    const kernel = createIntegralKernel({ ...input.parameters, exactJacobian: true });
    const before = evaluateTransitionInterval(kernel, input), original = structuredClone(input);
    const expected = { s: before.transition.s, kind: before.transition.forced ? 'forced' : 'natural' };
    const r = initializeMixedTransitionShear(kernel, { ...input, downstream: { ...input.downstream, aux: seed }, expected });
    assert.ok(r.aux > 0); assert.ok(r.diagnostics.evaluations < 100);
    assert.ok(Math.abs(r.diagnostics.residual) <= 1e-14);
    const after = evaluateTransitionInterval(kernel, { ...input, downstream: { ...input.downstream, aux: r.aux } });
    assert.equal(after.transition.s, before.transition.s); assert.equal(after.transition.forced, before.transition.forced);
    assert.ok(Math.abs(after.residual[0]) <= 1e-14);
    assert.equal(after.residual[1], before.residual[1], 'Momentum is independent of shear at fixed physical primitives.');
    changedEnergy ||= after.residual[2] !== before.residual[2];
    assert.deepEqual(input, original);
  }
  assert.equal(changedEnergy, true, 'The test must not incorrectly require mixed energy to be invariant.');
});

test('a scalar preparation rejects a changed transition branch without mutating input', () => {
  const { input } = fixture.cases.find(c => c.input.name === 'natural-M0.2');
  const kernel = createIntegralKernel({ ...input.parameters, exactJacobian: true }), before = structuredClone(input);
  assert.throws(() => initializeMixedTransitionShear(kernel, { ...input, expected: { kind: 'forced', s: input.downstream.s } }), /selected transition/);
  assert.deepEqual(input, before);
});

test('terminal preparation rejects a nonfinite native shear residual', () => {
  const input = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/terminal-transition.json', import.meta.url))).cases[0];
  const native = createIntegralKernel({ ...input.parameters, exactJacobian: true });
  const kernel = { ...native, interval: args => {
    const value = native.interval(args);
    return { ...value, residual: [NaN, ...value.residual.slice(1)] };
  } };
  assert.throws(() => initializeMixedTransitionShear(kernel, { ...input,
    expected: { kind: 'forced', s: input.downstream.s } }), /did not close/);
});
