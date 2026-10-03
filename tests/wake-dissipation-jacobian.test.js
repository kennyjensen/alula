import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { dilw } from '../src/viscous/xfoil/xblsys.js';
import { createIntegralKernel } from '../src/viscous/integral.js';

const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/wake-jacobian.json', import.meta.url)));
const close = (a, b, tolerance = 3e-10) => assert.ok(Math.abs(a - b) < tolerance * Math.max(1e-8, Math.abs(a), Math.abs(b)), `${a} != ${b}`);

test('native wake dissipation values/approximate partials retain Fortran parity, while exact Hk slopes match differences of the original value', () => {
  for (const [path, hash] of Object.entries(fixture.provenance.sourceHashes))
    assert.equal(createHash('sha256').update(fs.readFileSync(new URL(`../${path}`, import.meta.url))).digest('hex'), hash);
  for (const { hk, rt, expected } of fixture.primitives) {
    const native = dilw(hk, rt), exact = dilw(hk, rt, true);
    for (const key of ['di', 'diHk', 'diRt']) close(native[key], expected[key]);
    assert.equal(exact.di, native.di); assert.equal(exact.diRt, native.diRt);
    for (const difference of expected.differences) close(exact.diHk, difference, 3e-7);
    assert.ok(Math.abs(exact.diHk - native.diHk) > 1e-12);
  }
});

test('retained wake blocks preserve original residuals while corrected state and compressible partials match independent differences', () => {
  let witnessedMismatch = false;
  for (const { parameters, input, expected } of fixture.intervals) {
    const native = createIntegralKernel({ ...parameters, exactJacobian: false }).interval(input);
    const kernel = createIntegralKernel({ ...parameters, exactJacobian: true }), exact = kernel.interval(input);
    for (const key of ['residual', 'upstream', 'downstream']) native[key].flat().forEach((v, i) => close(v, expected[key].flat()[i]));
    assert.deepEqual(exact.residual, native.residual);
    for (const side of ['upstream', 'downstream']) for (const [col, key] of ['aux', 'theta', 'deltaStar', 'ue', 's'].entries()) {
      const value = input[side][key];
      for (const relativeStep of [1e-3, 1e-4, 1e-5]) {
        const h = relativeStep * Math.max(Math.abs(value), 1e-6);
        const samples = [-2, -1, 1, 2].map(m => kernel.interval({ ...input, [side]: { ...input[side], [key]: value + m * h } }).residual);
        for (let row = 0; row < 3; row++) {
          const fd = (samples[0][row] - 8 * samples[1][row] + 8 * samples[2][row] - samples[3][row]) / (12 * h);
          assert.ok(Math.abs(fd - exact[side][row][col]) < 2e-6 * Math.max(1, Math.abs(fd), Math.abs(exact[side][row][col])));
          if (row === 2 && Math.abs(fd - native[side][row][col]) > .1 * Math.max(1, Math.abs(fd))) witnessedMismatch = true;
        }
      }
    }
    for (const [key, value, derivative] of [['mach', parameters.mach ** 2, 'machSquaredDerivative'], ['reynolds', parameters.reynolds, 'reynoldsDerivative']]) {
      const h = 1e-4 * value;
      const f = [-2, -1, 1, 2].map(m => createIntegralKernel({ ...parameters, exactJacobian: true,
        [key]: key === 'mach' ? Math.sqrt(value + m * h) : value + m * h }).interval(input).residual);
      for (let row = 0; row < 3; row++) {
        const fd = (f[0][row] - 8 * f[1][row] + 8 * f[2][row] - f[3][row]) / (12 * h);
        assert.ok(Math.abs((exact[derivative][row] - fd) * value) < 2e-7 * Math.max(1, Math.abs(fd * value)),
          `${key}, row ${row}: ${exact[derivative][row]} != ${fd}`);
      }
    }
  }
  assert.equal(witnessedMismatch, true, 'the retained data must exercise the original inaccurate derivative');
});
