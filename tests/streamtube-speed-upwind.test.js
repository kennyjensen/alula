// SPDX-License-Identifier: GPL-2.0-or-later
// MSES manual, printed p. 12: local speed upwinding only. These tests do
// not compute a shock, change thermodynamic speeds, or validate an airfoil.
import test from 'node:test';
import assert from 'node:assert/strict';
import { evaluateStreamtubeSpeedUpwind, linearizeStreamtubeSpeedUpwind } from '../src/euler/streamtube-speed-upwind.js';

const evaluate = evaluateStreamtubeSpeedUpwind;
const linearize = linearizeStreamtubeSpeedUpwind;
const close = (actual, expected, tolerance = 3e-13) => assert.ok(
  Math.abs(actual - expected) <= tolerance * Math.max(1, Math.abs(actual), Math.abs(expected)),
  `${actual} != ${expected}`);
const base = { speeds: [1.1, 1.3, 1.15], machSquared: [.88, 1.17], spacing: [.32, .47],
  mucon: 1.25, mcrit: .91, gamma: 1.4 };
const arrays = { speeds: 3, machSquared: 2, spacing: 2 };
const scalars = ['mucon', 'mcrit', 'gamma'];
const shift = (parameters, tangent, step) => {
  const p = structuredClone(parameters);
  for (const [key, n] of Object.entries(arrays)) if (tangent[key])
    p[key] = Array.from({ length: n }, (_, i) => p[key][i] + step * tangent[key][i]);
  for (const key of scalars) if (tangent[key]) p[key] += step * tangent[key];
  return p;
};
const fourthDifference = (p, tangent, h) => {
  const a = evaluate(shift(p, tangent, -2 * h)), b = evaluate(shift(p, tangent, -h));
  const c = evaluate(shift(p, tangent, h)), d = evaluate(shift(p, tangent, 2 * h));
  return Object.fromEntries(['speed', 'coefficient', 'correction'].map(key =>
    [key, (a[key] - 8 * b[key] + 8 * c[key] - d[key]) / (12 * h)]));
};

test('uniform speed is unchanged across first/second order and subsonic, sonic and supersonic activation', () => {
  for (const q of [0, .4, 3]) for (const mucon of [-2, 0, 2])
    for (const machSquared of [[0, 0], [.2, .8], [1, 1], [1.5, 4]]) {
      const p = { ...base, speeds: [q, q, q], machSquared, mucon };
      const value = evaluate(p);
      assert.equal(value.speed, q); close(value.correction, 0);
      if (mucon !== 0) assert.equal(value.secondOrder, mucon > 0);
      const d = linearize(p).apply({ machSquared: [.2, -.1], spacing: [.03, -.04],
        mcrit: .1, gamma: .2, mucon: mucon === 0 ? 0 : .1 });
      close(d.speed, 0); close(d.correction, 0);
    }
});

test('second-order upwinding exactly reproduces affine speed on nonuniform physical arc stations', () => {
  for (const s of [[.1, .4, 1.7], [2, 2.01, 3.1], [0, .7, .72]])
    for (const slope of [-.3, .8]) for (const machSquared of [[.4, 1.8], [1, 1], [2.5, 3]]) {
      const speeds = s.map(x => 2 + slope * x), spacing = [s[1] - s[0], s[2] - s[1]];
      const value = evaluate({ ...base, speeds, spacing, machSquared });
      close(value.speed, speeds[2]); close(value.correction, 0);
    }
});

test('quadratic speed has the exact curvature correction on unequal arc intervals', () => {
  // For q(s)=a+b*s+c*s^2, unequal divided differences give
  // qtilde_i-q_i = -mu*c*h_i*(h_i+h_previous). At mean M²=1,
  // mu=|MUCON|/gamma*(1-Mcrit)*ln(2), with no copied stencil formula.
  const mucon = 1.4, gamma = 1.4, mcrit = .8, coefficient = .2 * Math.log(2);
  for (const [h0, h1] of [[.3, .9], [.8, .2], [.4, .4]]) for (const c of [-.07, .12]) {
    const s = [.2, .2 + h0, .2 + h0 + h1], speeds = s.map(x => 1.7 + .2 * x + c * x * x);
    const value = evaluate({ speeds, machSquared: [.7, 1.3], spacing: [h0, h1], mucon, gamma, mcrit });
    const correction = -coefficient * c * h1 * (h0 + h1);
    close(value.coefficient, coefficient); close(value.correction, correction);
    close(value.speed, speeds[2] + correction);
  }
});

test('a prescribed normal compression receives the first-order upstream speed bias', () => {
  // Independent perfect-gas normal-jump identities at gamma=7/5, M1=2:
  // rho2/rho1=8/3, p2/p1=9/2, M2²=1/3. The mean M² is 13/6,
  // so the Mcrit=1, MUCON=-1 limit gives mu=5/13 and qtilde2=8*q1/13.
  // This checks orientation of the local stencil; no shock is solved here.
  const gamma = 7 / 5, upstream = 2 * Math.sqrt(gamma), downstream = 3 * upstream / 8;
  assert.ok(Math.log(9 / 2) - gamma * Math.log(8 / 3) > 0);
  const p = { speeds: [upstream, upstream, downstream], machSquared: [4, 1 / 3],
    spacing: [.2, .9], mucon: -1, mcrit: 1, gamma };
  const value = evaluate(p);
  assert.equal(value.secondOrder, false); close(value.coefficient, 5 / 13);
  close(value.speed, 8 * upstream / 13);
  assert.ok(value.speed > downstream && value.speed < upstream);
  const d = linearize(p);
  close(d.apply({ speeds: [1, 0, 0] }).speed, 0);
  close(d.apply({ speeds: [0, 1, 0] }).speed, 5 / 13);
  close(d.apply({ speeds: [0, 0, 1] }).speed, 8 / 13);
  close(d.apply({ spacing: [.7, -.3] }).speed, 0);
});

test('activation uses mean Mach squared, 1-Mcrit, and the documented sonic limit', () => {
  const p = { ...base, mucon: -1.4, gamma: 1.4, mcrit: .8, machSquared: [.25, 1.75] };
  close(evaluate(p).coefficient, .2 * Math.log(2));
  // Smooth activation is already positive at the nominal critical Mach;
  // it is not a hard switch at Mcrit and does not average Mach numbers.
  assert.ok(evaluate({ ...p, machSquared: [.64, .64] }).coefficient > 0);
  const below = { ...p, mcrit: 1, machSquared: [.7, .9] };
  assert.equal(evaluate(below).coefficient, 0);
  close(linearize(below).apply({ machSquared: [.2, .3] }).coefficient, 0);
  const above = { ...p, mcrit: 1, machSquared: [1.5, 2.5] };
  close(evaluate(above).coefficient, .5);
  close(linearize(above).apply({ machSquared: [1, 1] }).coefficient, .25);
  for (const state of [below, above]) assert.throws(() => linearize(state).apply({ mcrit: 1 }));
  const sonic = { ...p, mcrit: 1, machSquared: [1, 1] };
  assert.equal(evaluate(sonic).coefficient, 0);
  assert.throws(() => linearize(sonic));
  for (const h of [1e-3, 1e-5]) {
    const right = evaluate({ ...sonic, machSquared: [1 + h, 1 + h] }).coefficient;
    close(right / h, 1 / (1 + h), 2e-11);
    assert.equal(evaluate({ ...sonic, machSquared: [1 - h, 1 - h] }).coefficient, 0);
  }
  for (const x of [0, Number.MIN_VALUE, 1e-200, 1e-12, 1e200]) {
    const tiny = { ...p, machSquared: [x, x] }, value = evaluate(tiny);
    for (const key of ['speed', 'coefficient', 'correction']) assert.ok(Number.isFinite(value[key]));
    for (const v of Object.values(linearize(tiny).apply({ machSquared: [1, 1] }))) assert.ok(Number.isFinite(v));
    if (x <= 1e-12) assert.equal(value.coefficient, 0);
  }
});

test('analytic partials match fourth-order independent differences for states, spacing and controls', t => {
  let comparisons = 0, maximumAbsoluteError = 0, maximumNormalizedError = 0;
  const tangents = [];
  for (const [key, n] of Object.entries(arrays)) for (let i = 0; i < n; i++)
    tangents.push({ [key]: Array.from({ length: n }, (_, j) => i === j ? 1 : 0) });
  for (const key of scalars) tangents.push({ [key]: 1 });
  tangents.push({ speeds: [.2, -.3, .4], machSquared: [.13, -.07], spacing: [.02, -.05],
    gamma: .11, mcrit: -.04, mucon: .12 });
  for (const mucon of [-.7, 1.25]) for (const machSquared of [[.35, .65], [.88, 1.17], [1.5, 2.1]]) {
    const p = { ...base, mucon, machSquared }, analytic = linearize(p);
    assert.deepEqual(analytic.value, evaluate(p));
    for (const tangent of tangents) {
      const expected = fourthDifference(p, tangent, 2e-5), actual = analytic.apply(tangent);
      for (const key of ['speed', 'coefficient', 'correction']) {
        const error = Math.abs(actual[key] - expected[key]);
        maximumAbsoluteError = Math.max(maximumAbsoluteError, error);
        maximumNormalizedError = Math.max(maximumNormalizedError,
          error / Math.max(1, Math.abs(actual[key]), Math.abs(expected[key])));
        comparisons++; close(actual[key], expected[key], 2e-8);
      }
    }
  }
  t.diagnostic(JSON.stringify({ comparisons, maximumAbsoluteError, maximumNormalizedError }));
});

test('geometric scale cancels from unequal spacing and its directional derivative', () => {
  for (const mucon of [-1, 1]) {
    const p = { ...base, mucon }, value = evaluate(p);
    for (const scale of [1e-8, .2, 7, 1e8]) {
      const mapped = evaluate({ ...p, spacing: p.spacing.map(h => scale * h) });
      for (const key of ['speed', 'coefficient', 'correction']) close(mapped[key], value[key]);
    }
    const d = linearize(p).apply({ spacing: p.spacing });
    for (const key of ['speed', 'coefficient', 'correction']) close(d[key], 0);
  }
});

test('zero upwinding preserves physical speed and malformed stencils are rejected without mutation', () => {
  const p = structuredClone(base), before = structuredClone(p);
  evaluate(p); linearize(p).apply({ speeds: [.1, .2, .3], spacing: [.2, -.1] });
  assert.deepEqual(p, before);
  const zero = linearize({ ...p, mucon: 0 });
  assert.equal(zero.value.speed, p.speeds[2]); assert.equal(zero.value.coefficient, 0);
  const d = zero.apply({ speeds: [.1, .2, .3], spacing: [.2, -.1] });
  close(d.speed, .3); close(d.coefficient, 0); close(d.correction, 0);
  assert.throws(() => zero.apply({ mucon: .1 }));
  for (const malformed of [{ speeds: [1, 2] }, { speeds: [1, -1, 2] }, { machSquared: [-1, 1] },
    { spacing: [0, 1] }, { spacing: [1, NaN] }, { gamma: 1 }, { mcrit: -.01 }, { mcrit: 1.01 }, { mucon: Infinity }])
    assert.throws(() => evaluate({ ...p, ...malformed }));
  for (const tangent of [{ speeds: [1, 2] }, { spacing: [0, NaN] }, { gamma: Infinity }])
    assert.throws(() => linearize(p).apply(tangent));
});
