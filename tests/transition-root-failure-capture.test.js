import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { XFOIL_TRANSITION_ROOT_ACCURACY } from '../src/viscous/xfoil/xblsys.js';

// The original 1e-12 inputs now converge (critical-cycle.test.js). Only these
// diagnostic tests deliberately request 1e-30, below the measured roundoff
// floor, to exercise failure capture without requiring the repaired bug.
const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/transition-root-failures.json', import.meta.url)));
const jsonCopy = value => JSON.parse(JSON.stringify(value));
const diagnosticCases = fixture.cases.map(c => ({ ...c, payload: { ...c.payload,
  parameters: { ...c.payload.parameters, transitionTolerance: 1e-30 } } }));
const first = diagnosticCases[0];
const freezeInput = input => Object.freeze({ ...input,
  upstream: Object.freeze({ ...input.upstream }), downstream: Object.freeze({ ...input.downstream }) });
function capture(kernel, method, input) {
  let caught;
  try { kernel[method](input); } catch (error) { caught = error; }
  assert.ok(caught instanceof Error, 'The captured current-algorithm failure must occur.');
  assert.equal(caught.code, XFOIL_TRANSITION_ROOT_ACCURACY);
  assert.equal(caught.message, first.message);
  return caught.transitionRootFailure;
}

test('below-roundoff accuracy failures replay serializable local diagnostics without caller mutation', () => {
  assert.equal(fixture.cases.length, 13);
  for (const observed of diagnosticCases) {
    const payload = jsonCopy(observed.payload), input = freezeInput(payload.input);
    const snapshot = jsonCopy(input);
    const replay = capture(createIntegralKernel(payload.parameters), payload.method, input);
    assert.deepEqual(jsonCopy(replay), payload, observed.name);
    assert.deepEqual(input, snapshot, observed.name);
    assert.equal(replay.parameters.transitionTolerance, 1e-30);
  }
});

test('failure snapshots copy only whitelisted primitives and remain detached from caller metadata', () => {
  const payload = jsonCopy(first.payload), input = payload.input;
  for (const state of [input.upstream, input.downstream]) {
    state.ctx = { mustNotBeCaptured: true };
    state.alias = () => state;
    state.stationId = 123;
    state.wakeGap = undefined;
  }
  input.upstream.amplification = input.upstream.aux;
  // transitionCheck deliberately replaces the downstream amplification guess;
  // its public replay input still records the caller's optional value exactly.
  input.downstream.amplification = 123;
  const replay = capture(createIntegralKernel(payload.parameters), payload.method, input);
  const saved = jsonCopy(replay);
  assert.deepEqual(Object.keys(replay).sort(), ['version', 'code', 'parameters', 'method', 'input'].sort());
  for (const side of ['upstream', 'downstream']) {
    assert.deepEqual(Object.keys(replay.input[side]).sort(), ['s', 'aux', 'theta', 'deltaStar', 'ue', 'amplification'].sort());
    assert.notEqual(replay.input[side], input[side]);
    input[side].s += 1;
    input[side].aux = 0;
    input[side].ctx.mustNotBeCaptured = false;
  }
  assert.deepEqual(replay, saved);
  assert.deepEqual(jsonCopy(capture(createIntegralKernel(saved.parameters), saved.method, saved.input)), saved);
});

test('the same real failure is captured at the explicit transition-interval call with normalized defaults', () => {
  const payload = jsonCopy(first.payload);
  const input = freezeInput({ upstream: payload.input.upstream, downstream: payload.input.downstream, regime: 'transition' });
  const replay = capture(createIntegralKernel(payload.parameters), 'interval', input);
  assert.deepEqual(jsonCopy(replay), { ...payload, method: 'interval',
    input: { ...payload.input, regime: 'transition', similarityExponent: 1 } });
});

test('an unrelated wake interval does not contaminate a repeated transition-failure capture', () => {
  const payload = jsonCopy(first.payload), kernel = createIntegralKernel(payload.parameters);
  const before = capture(kernel, payload.method, payload.input);
  const a = { s: 2, aux: .04, theta: .002, deltaStar: .004, ue: .9 };
  const wake = kernel.interval({ regime: 'wake', upstream: a, downstream: { ...a, s: 2.1 } });
  assert.ok(wake.residual.every(Number.isFinite));
  const after = capture(kernel, payload.method, payload.input);
  assert.deepEqual(after, before);
  assert.notEqual(after, before);
  assert.notEqual(after.input.upstream, before.input.upstream);
});

test('ordinary physical-domain failures retain their own error and receive no transition-root payload', () => {
  const payload = jsonCopy(first.payload), kernel = createIntegralKernel(payload.parameters);
  const input = { ...payload.input, upstream: { ...payload.input.upstream, theta: -1 } };
  assert.throws(() => kernel.transitionCheck(input), error => {
    assert.equal(error.message, 'Inadmissible integral BL station.');
    assert.equal(error.code, undefined);
    assert.equal(error.transitionRootFailure, undefined);
    return true;
  });
});
