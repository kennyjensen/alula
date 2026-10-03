import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createIntegralKernel } from '../src/viscous/integral.js';

const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/bl-active-limits.json', import.meta.url)));
test('active Hk and upwind limits preserve native values and default partials while exact state derivatives follow the active branch', () => {
  for (const [file, hash] of Object.entries(fixture.provenance.sourceHashes))
    assert.equal(createHash('sha256').update(fs.readFileSync(new URL(`../${file}`, import.meta.url))).digest('hex'), hash);
  for (const { parameters, input, expected } of fixture.cases) {
    const native = createIntegralKernel({ ...parameters, exactJacobian: false }).interval(input);
    for (const key of ['residual', 'upstream', 'downstream']) native[key].flat().forEach((v, i) =>
      assert.ok(Math.abs(v - expected[key].flat()[i]) < 3e-10 * Math.max(1, Math.abs(v))));
    const kernel = createIntegralKernel({ ...parameters, exactJacobian: true }), block = kernel.interval(input);
    assert.deepEqual(block.residual, native.residual);
    if (input.regime === 'laminar') {
      assert.equal(block.properties.hk, 1.05);
      // The former local central stencil actually leaves the physical
      // domain here, despite an admissible base state.
      const s = input.downstream;
      assert.throws(() => kernel.interval({ ...input, downstream: { ...s, theta: s.theta * (1 + Math.cbrt(Number.EPSILON)) } }), /compressible BL edge state/);
    } else {
      const a = kernel.station(input.upstream, 'wake'), b = kernel.station(input.downstream, 'wake');
      assert.ok(Math.log((b.hk - 1) / (a.hk - 1)) ** 2 > 15);
    }
    for (const side of ['upstream', 'downstream']) for (const [col, key] of ['aux', 'theta', 'deltaStar', 'ue', 's'].entries()) {
      for (const relativeStep of [1e-3, 1e-4, 1e-5]) {
        const h = relativeStep * Math.max(Math.abs(input[side][key]), 1e-5);
        const sign = side === 'downstream' && input.regime === 'laminar' && ['theta', 'ue'].includes(key) ? -1 : 1;
        const samples = [0, 1, 2, 3, 4].map(m => kernel.interval({ ...input,
          [side]: { ...input[side], [key]: input[side][key] + m * sign * h } }).residual);
        for (let row = 0; row < 3; row++) {
          const fd = (-25 * samples[0][row] + 48 * samples[1][row] - 36 * samples[2][row] + 16 * samples[3][row] - 3 * samples[4][row]) / (12 * sign * h);
          const a = block[side][row][col];
          assert.ok(Math.abs(fd - a) < 2e-6 * Math.max(1, Math.abs(fd), Math.abs(a)), `${input.regime} ${side}.${key}, row ${row}, h=${h}: ${fd} != ${a}`);
        }
      }
    }
    for (const [key, value, derivative] of [['mach', parameters.mach ** 2, 'machSquaredDerivative'], ['reynolds', parameters.reynolds, 'reynoldsDerivative']]) {
      const h = 1e-4 * value, sign = key === 'mach' ? -1 : 1;
      const f = [0, 1, 2, 3, 4].map(m => createIntegralKernel({ ...parameters, exactJacobian: true,
        [key]: key === 'mach' ? Math.sqrt(value + m * sign * h) : value + m * sign * h }).interval(input).residual);
      for (let row = 0; row < 3; row++) {
        const fd = (-25 * f[0][row] + 48 * f[1][row] - 36 * f[2][row] + 16 * f[3][row] - 3 * f[4][row]) / (12 * sign * h);
        assert.ok(Math.abs((fd - block[derivative][row]) * value) < 2e-6 * Math.max(1, Math.abs(fd * value)), `${key}, row ${row}`);
      }
    }
  }
});
