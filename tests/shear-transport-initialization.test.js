import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { initializeTransportShear } from '../src/viscous/shear-transport-initialization.js';

const fixture = JSON.parse(readFileSync(new URL('./fixtures/fortran/event-shear-transport.json', import.meta.url)));
const frozen = value => Object.freeze({ ...value, upstream: Object.freeze({ ...value.upstream }), downstream: Object.freeze({ ...value.downstream }) });

test('native ordinary turbulent/wake scalar roots preserve primitives and satisfy original Fortran rows', () => {
  for (const c of fixture.cases) {
    const kernel = createIntegralKernel(c.parameters), input = frozen(c.input), before = JSON.stringify(input);
    const initial = kernel.interval(input), prepared = initializeTransportShear(kernel, input);
    const after = kernel.interval({ ...input, downstream: { ...input.downstream, aux: prepared.aux } });
    assert.ok(Math.abs(prepared.aux - c.expectedAux) <= 2e-11 * c.expectedAux, c.name);
    assert.ok(Math.abs(after.residual[0]) <= 1e-14, c.name);
    for (let r = 0; r < 3; r++) assert.ok(Math.abs(after.residual[r] - c.expectedResidual[r]) < 2e-13, `${c.name}/row${r}`);
    assert.equal(after.residual[1], initial.residual[1], 'shear preparation preserves the momentum equation');
    assert.equal(JSON.stringify(input), before);
    assert.equal(prepared.diagnostics.equationsChanged, false);
    assert.equal(prepared.diagnostics.initialGuessOnly, true);
    assert.equal(prepared.diagnostics.tolerance, 1e-14);
    assert.equal(initializeTransportShear(kernel, { ...input, downstream: { ...input.downstream, aux: prepared.aux } }).aux,
      prepared.aux, 'already closed seed is retained bit-exact');
  }
});

test('both widely separated positive seeds recover the same unique native root', () => {
  for (const c of fixture.cases) {
    const kernel = createIntegralKernel(c.parameters);
    for (const multiplier of [1e-5, 1e5]) {
      const p = initializeTransportShear(kernel, { ...c.input, downstream: { ...c.input.downstream, aux: c.expectedAux * multiplier } });
      assert.ok(Math.abs(p.aux - c.expectedAux) < 2e-11 * c.expectedAux, `${c.name}/${multiplier}`);
    }
  }
});

test('rejects invalid regimes and nonpositive or nonfinite shear without mutation', () => {
  const c = fixture.cases[0], kernel = createIntegralKernel(c.parameters);
  for (const regime of ['laminar', 'transition', 'similarity', undefined])
    assert.throws(() => initializeTransportShear(kernel, { ...c.input, regime }), /turbulent or wake/);
  for (const side of ['upstream', 'downstream']) for (const aux of [0, -1, Infinity, NaN]) {
    const input = frozen({ ...c.input, [side]: { ...c.input[side], aux } });
    assert.throws(() => initializeTransportShear(kernel, input), /positive finite/);
    assert.ok(Object.is(input[side].aux, aux));
  }
});

test('propagates physical-domain failures and cancellation, and rejects nonmonotone or unresolved native rows', () => {
  const c = fixture.cases[0], kernel = createIntegralKernel(c.parameters);
  const bad = frozen({ ...c.input, downstream: { ...c.input.downstream, ue: 20 } });
  assert.throws(() => initializeTransportShear(kernel, bad), error => error.code === 'BL_EDGE_STATE_DOMAIN');
  assert.throws(() => initializeTransportShear(kernel, { ...c.input, downstream: { ...c.input.downstream, s: c.input.upstream.s } }), /increase downstream/);
  const cancellation = Object.assign(new Error('cancelled'), { name: 'AbortError' });
  assert.throws(() => initializeTransportShear({ interval() { throw cancellation; } }, c.input), error => error === cancellation);
  for (const derivative of [0, 1, NaN]) assert.throws(() => initializeTransportShear({
    interval() { return { residual: [1], downstream: [[derivative]] }; }
  }, c.input), /strictly decreasing/);
  assert.throws(() => initializeTransportShear({ interval() { return { residual: [1], downstream: [[-1]] }; } }, c.input), /Could not bracket/);
});
