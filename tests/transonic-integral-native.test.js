import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createIntegralKernel } from '../src/viscous/integral.js';

const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/transonic-integral.json', import.meta.url)));
const error = (a, b) => Math.abs(a - b) / Math.max(1, Math.abs(a), Math.abs(b));

test('RAE64x11 turbulent shock intervals retain complete derivatives along the restricted Newton step', t => {
  const f = JSON.parse(fs.readFileSync(new URL('fixtures/rae64x11-turbulent-curvature-intervals.json', import.meta.url)));
  assert.equal(createHash('sha256').update(fs.readFileSync(new URL(`../${f.provenance.checkpoint}`, import.meta.url))).digest('hex'), f.provenance.sha256);
  const kernel = createIntegralKernel(f.parameters), keys = ['aux', 'theta', 'deltaStar', 'ue', 's'];
  for (const c of f.cases) {
    const original = structuredClone(c.input), base = kernel.interval(c.input);
    base.residual.forEach((r, row) => assert.ok(Math.abs(r - c.expectedResidual[row]) < 1e-12));
    // Scale each independent perturbation to its physical variable. A large
    // unnormalized global direction alone obscures this kernel's local limit.
    for (const side of ['upstream', 'downstream']) for (const [col, key] of keys.entries()) {
      const scale = key === 's' ? c.input.downstream.s - c.input.upstream.s : Math.abs(c.input[side][key]);
      const errors = [1e-4, 1e-5, 1e-6].map(epsilon => {
        const h = scale * epsilon;
        const at = sign => kernel.interval({ ...c.input, [side]: { ...c.input[side], [key]: c.input[side][key] + sign * h } }).residual;
        const plus = at(1), minus = at(-1);
        return Math.max(...plus.map((r, row) => Math.abs((r - minus[row]) / (2 * h) - base[side][row][col])
          * scale / Math.max(1, Math.abs(base[side][row][col]) * scale)));
      });
      assert.ok(Math.min(...errors) < 1e-8, `i=${c.station.i}/${side}/${key}: ${errors}`);
    }
    const derivative = base.residual.map((_, row) => keys.reduce((sum, key, col) => sum
      + base.upstream[row][col] * c.direction.upstream[key]
      + base.downstream[row][col] * c.direction.downstream[key], 0));
    derivative.forEach((d, row) => assert.ok(Math.abs(d + base.residual[row]) < 1e-10));
    const at = h => kernel.interval({ ...c.input, ...Object.fromEntries(['upstream', 'downstream'].map(side =>
      [side, Object.fromEntries(keys.map(key => [key, c.input[side][key] + h * c.direction[side][key]]))])) }).residual;
    const sweep = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8].map(h => {
      const plus = at(h), minus = at(-h);
      return { h, error: Math.max(...derivative.map((d, row) => Math.abs((plus[row] - minus[row]) / (2 * h) - d))) };
    });
    assert.ok(sweep[1].error < .02 * sweep[0].error);
    assert.ok(sweep[2].error < .02 * sweep[1].error);
    assert.ok(Math.min(...sweep.map(s => s.error)) < 1e-8);
    const step = .0029164780714, residual = at(step);
    assert.ok(Math.abs(residual[1]) > Math.abs(base.residual[1]), 'The real limited step has nonlinear momentum growth despite a descending tangent.');
    assert.deepEqual(c.input, original);
    t.diagnostic(JSON.stringify({ i: c.station.i, derivative, sweep }));
  }
});

test('RAE64x11 amplification onset has a consistent Jacobian despite strong Newton-direction curvature', t => {
  const c = JSON.parse(fs.readFileSync(new URL('fixtures/rae64x11-laminar-amplification-interval.json', import.meta.url)));
  assert.equal(createHash('sha256').update(fs.readFileSync(new URL(`../${c.provenance.checkpoint}`, import.meta.url))).digest('hex'), c.provenance.sha256);
  const kernel = createIntegralKernel(c.parameters), before = structuredClone(c.input);
  const base = kernel.interval(c.input), keys = ['aux', 'theta', 'deltaStar', 'ue', 's'];
  const derivative = keys.reduce((sum, key, i) => sum
    + base.upstream[0][i] * c.direction.upstream[key]
    + base.downstream[0][i] * c.direction.downstream[key], 0);
  assert.ok(Math.abs(base.residual[0] - c.expected.amplificationResidual) < 1e-12);
  assert.ok(Math.abs(derivative + base.residual[0]) < 1e-10,
    'This is the actual coupled Newton direction, satisfying J dx = -R.');
  const at = h => {
    const input = { ...c.input };
    for (const side of ['upstream', 'downstream']) {
      input[side] = Object.fromEntries(keys.map(key => [key, c.input[side][key] + h * c.direction[side][key]]));
      assert.ok(kernel.station(input[side], 'laminar').rawHk > 1);
    }
    return kernel.interval(input).residual[0];
  };
  const sweep = [1e-5, 1e-6, 1e-7, 1e-8, 1e-9].map(h => ({ h,
    error: Math.abs((at(h) - at(-h)) / (2 * h) - derivative) }));
  // A random, normalized direction missed this very narrow local scale.
  // Test second-order truncation decay, then agreement before roundoff wins.
  assert.ok(sweep[0].error > Math.abs(derivative));
  assert.ok(sweep[1].error < .02 * sweep[0].error);
  assert.ok(sweep[2].error < .02 * sweep[1].error);
  assert.ok(Math.min(...sweep.map(s => s.error)) < 1e-6);
  assert.deepEqual(c.input, before);
  t.diagnostic(JSON.stringify({ derivative, sweep }));
});

test('transonic prescribed-edge BL residuals and original approximate derivatives match unchanged Fortran', t => {
  for (const [name, hash] of Object.entries(fixture.provenance.sourceHashes))
    assert.equal(createHash('sha256').update(fs.readFileSync(new URL(`../${name}`, import.meta.url))).digest('hex'), hash, name);
  let residualError = 0, approximateError = 0;
  assert.equal(fixture.cases.length, 9);
  for (const c of fixture.cases) {
    const legacy = createIntegralKernel(c.parameters).interval(c.input);
    const complete = createIntegralKernel({ ...c.parameters, exactJacobian: true }).interval(c.input);
    assert.deepEqual(complete.residual, legacy.residual);
    for (const key of ['residual', 'upstream', 'downstream', 'reynoldsDerivative', 'machSquaredDerivative', 'tripDerivative']) {
      const expected = c.expected[key].flat(), actual = legacy[key].flat(); assert.equal(actual.length, expected.length);
      actual.forEach((v, i) => {
        assert.ok(Number.isFinite(v) && Number.isFinite(expected[i]));
        const e = error(v, expected[i]); approximateError = Math.max(approximateError, e);
        if (key === 'residual') residualError = Math.max(residualError, e);
        assert.ok(e < fixture.controls.nativeJacobianParityLimit, `${c.name}/${key}/${i}: ${e}`);
      });
    }
  }
  t.diagnostic(JSON.stringify({ residualError, approximateError }));
});

test('complete transonic BL Jacobian matches independent fourth-order differences of original Fortran residuals', t => {
  const blocks = fixture.cases.map(c => createIntegralKernel({ ...c.parameters, exactJacobian: true }).interval(c.input));
  let maximumError = 0;
  assert.equal(fixture.derivatives.length, 180);
  for (const d of fixture.derivatives) for (let row = 0; row < 3; row++) {
    const native = (d.samples[0][row] - 8 * d.samples[1][row] + 8 * d.samples[2][row] - d.samples[3][row]) / (12 * d.step);
    assert.equal(native, d.expected[row]);
    const actual = blocks[d.caseIndex][d.side][row][d.column], e = error(actual, native);
    assert.ok(Number.isFinite(actual) && Number.isFinite(native)); maximumError = Math.max(maximumError, e);
    assert.ok(e < fixture.controls.derivativeLimit, `${fixture.cases[d.caseIndex].name}/${d.side}/${d.key}/${row}: ${e}`);
  }
  t.diagnostic(JSON.stringify({ scalarDerivativeChecks: 540, maximumError }));
});
