import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { evaluateTransitionInterval } from '../src/viscous/transition-interval.js';

const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/prepared-transition-root-precision.json', import.meta.url)));
const original = JSON.parse(fs.readFileSync(new URL('../docs/current-prepared-transition-shear.json', import.meta.url)));
const hash = s => createHash('sha256').update(s).digest('hex');
const close = (a, b, limit = 3e-10) => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && Math.abs(a - b) <= limit * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const evaluate = (kernel, input, jacobian = false) => input.regime === 'transition'
  ? evaluateTransitionInterval(kernel, input, { jacobian }) : kernel.interval(input);

test('native precision oracle changes only the stopping tolerance and adds an iteration trace', () => {
  for (const [file, expected] of Object.entries(fixture.provenance.sha256))
    assert.equal(hash(fs.readFileSync(new URL('../' + file, import.meta.url))), expected, file);
  for (const v of fixture.provenance.variants) {
    let code = fs.readFileSync(new URL('../third_party/Xfoil/src/xblsys.f', import.meta.url), 'utf8');
    assert.equal(v.patches.length, v.tolerance === 5e-5 ? 1 : 2);
    for (const p of v.patches) { assert.equal(code.split(p.from).length, 2); code = code.replace(p.from, p.to); }
    assert.equal(hash(code), v.sourceSHA256);
  }
  const columns = fixture.originalBranch.columns;
  assert.ok(columns.filter(c => c.relativeStep === 2e-5).every(c => new Set(c.samples.map(v => v.trace.length)).size > 1));
  assert.ok(columns.filter(c => c.relativeStep <= 1e-5).every(c => new Set(c.samples.map(v => v.trace.length)).size === 1));
});

test('resolved transition roots match the precision Fortran blocks on stalled and two-element states', () => {
  for (const c of fixture.cases) {
    const k = createIntegralKernel(c.parameters), value = evaluate(k, c.input);
    value.residual.forEach((v, r) => close(v, c.expected.residual[r]));
    if (value.transition) close(value.transition.s, c.expected.transitionS, 3e-12);
    // An unrelated wake call must not change the next surface calculation.
    const a = { s: 2, aux: .04, theta: .002, deltaStar: .004, ue: .9 };
    k.interval({ regime: 'wake', upstream: a, downstream: { ...a, s: 2.1 } });
    assert.deepEqual(evaluate(k, c.input).residual, value.residual);
  }
});

test('resolved transition Jacobians match two independent native difference steps', t => {
  const values = fixture.cases.map(c => evaluate(createIntegralKernel(c.parameters), c.input, true));
  const keys = ['aux', 'theta', 'deltaStar', 'ue', 's']; let maximum = 0;
  for (const c of fixture.derivatives) {
    const v = values[c.case], matrix = v.partials?.[c.side] ?? v[c.side];
    matrix.forEach((row, r) => {
      const a = row[keys.indexOf(c.key)], b = c.expected[r]; close(a, b, 2e-6);
      maximum = Math.max(maximum, Math.abs(a - b) / Math.max(1, Math.abs(a), Math.abs(b)));
    });
  }
  t.diagnostic(JSON.stringify({ cases: values.length, columns: fixture.derivatives.length, maximum }));
});

test('the default integral kernel retains original native stopping accuracy', () => {
  for (const c of original.cases) {
    const value = evaluate(createIntegralKernel(c.parameters), c.input);
    value.residual.forEach((v, r) => close(v, c.expected.residual[r]));
  }
});
