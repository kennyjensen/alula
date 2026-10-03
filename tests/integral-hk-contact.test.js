import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { wakeHkContactCase } from './helpers/integral-hk-contact.js';
const fixture = JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/hk-contact.json', import.meta.url)));
const difference = (a, b) => Math.max(...a.flat().map((v, i) => Math.abs(v - b.flat()[i])));

test('one ULP at the wake DSLIM contact preserves residuals and the admissible one-sided derivative', () => {
  const c = wakeHkContactCase(), exact = createIntegralKernel(c.parameters), native = createIntegralKernel({ ...c.parameters, hkFloorLinearization: 'native' });
  assert.deepEqual(c.below, fixture.cases[0].input); assert.deepEqual(c.contact, fixture.cases[1].input);
  const bits = new DataView(new ArrayBuffer(8)); bits.setFloat64(0, c.below.downstream.deltaStar); const belowBits = bits.getBigUint64(0);
  bits.setFloat64(0, c.contact.downstream.deltaStar); assert.equal(bits.getBigUint64(0) - belowBits, 1n);
  const a = exact.interval(c.below), b = exact.interval(c.contact), na = native.interval(c.below), nb = native.interval(c.contact);
  assert(a.properties.rawHk > 1 && a.properties.rawHk < c.minimumHk); assert.equal(b.properties.rawHk, c.minimumHk);
  assert.deepEqual(a.residual, b.residual); assert(difference(a.downstream, b.downstream) < 1e-10);
  for (const [x, y] of [[a, na], [b, nb]]) { assert.deepEqual(x.residual, y.residual); assert.deepEqual(x.properties, y.properties); }
  assert.equal(na.hkFloorLinearizationUsed, true); assert.equal(nb.hkFloorLinearizationUsed, undefined);
  assert(difference(na.downstream, nb.downstream) < 1e-10);
  // At roundoff contact both use the increasing-Hk one-sided derivative.
  assert.deepEqual(a.downstream, na.downstream); assert.deepEqual(b, nb);
  const explicitExact = createIntegralKernel({ ...c.parameters, hkFloorLinearization: 'exact' });
  assert.deepEqual(explicitExact.interval(c.below), a);
});

test('the two contact states retain original Fortran residuals and every native local partial', () => {
  for (const [file, hash] of Object.entries(fixture.provenance.sha256))
    assert.equal(createHash('sha256').update(fs.readFileSync(new URL('../' + file, import.meta.url))).digest('hex'), hash, file);
  const kernel = createIntegralKernel({ ...fixture.parameters, exactJacobian: false });
  for (const c of fixture.cases) {
    const result = kernel.interval(c.input);
    for (const key of Object.keys(c.expected)) c.expected[key].flat().forEach((v, i) => {
      const actual = result[key].flat()[i]; assert(Number.isFinite(v) && Number.isFinite(actual));
      assert(Math.abs(v - actual) <= 3e-10 * Math.max(1, Math.abs(v), Math.abs(actual)), key);
    });
  }
});


test('wake contact derivative agrees with increasing-thickness finite differences; strict clipping remains exact', () => {
  const c = wakeHkContactCase(), kernel = createIntegralKernel(c.parameters);
  for (const input of [c.below, c.contact]) {
    const value = kernel.interval(input), h = input.downstream.theta * 1e-7;
    const samples = [0, 1, 2, 3, 4].map(k => kernel.interval({ ...input,
      downstream: { ...input.downstream, deltaStar: input.downstream.deltaStar + k * h } }).residual);
    for (let r = 0; r < 3; r++) {
      const fd = (-25 * samples[0][r] + 48 * samples[1][r] - 36 * samples[2][r] + 16 * samples[3][r] - 3 * samples[4][r]) / (12 * h);
      assert.ok(Math.abs(fd - value.downstream[r][2]) < 2e-5 * Math.max(1, Math.abs(fd)), `${r}: ${fd} vs ${value.downstream[r][2]}`);
    }
  }
  const clipped = { ...c.below, downstream: { ...c.below.downstream, deltaStar: c.below.downstream.deltaStar - 1e-6 * c.below.downstream.theta } };
  const exact = kernel.interval(clipped), native = createIntegralKernel({...c.parameters,hkFloorLinearization:'native'}).interval(clipped);
  assert.deepEqual(exact.residual,native.residual);
  assert.notDeepEqual(exact.downstream,native.downstream);
});
