// SPDX-License-Identifier: GPL-2.0-or-later
// Shared physical-section stencil tests. Explicit boundary choices below
// are not evidence of the inlet implementation used by MSES.
import test from 'node:test';
import assert from 'node:assert/strict';
import { evaluateStreamtubeSpeedStencil as evaluate,
  linearizeStreamtubeSpeedStencil as linearize } from '../src/euler/tests/streamtube-speed-stencil.js';

const close = (a, b, tolerance = 3e-13) => assert.ok(
  Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const keys = ['speeds', 'machSquared', 'sectionArc'];
const outputKeys = ['speeds', 'coefficients', 'corrections'];
const controls = ['mucon', 'mcrit', 'gamma'];
const base = { speeds: [1.2, 1.05, 1.31, 1.11, 1.24], machSquared: [.92, 1.05, 1.25, 1.13, .95],
  sectionArc: [0, .3, .8, 1.7, 2.1], mucon: 1.1, mcrit: .86, gamma: 1.4,
  boundary: { kind: 'upstream-history', speeds: [1.04, 1.14], machSquared: [.6, .84], sectionArc: [-.9, -.2] } };
const polynomial = (arc, c = 0) => arc.map(s => 2 + .2 * s + c * s * s);
const shift = (p, d, h) => {
  const next = structuredClone(p);
  for (const key of keys) {
    if (d[key]) next[key] = next[key].map((v, i) => v + h * d[key][i]);
    if (d.boundary?.[key]) next.boundary[key] = next.boundary[key].map((v, i) => v + h * d.boundary[key][i]);
  }
  for (const key of controls) if (d[key]) next[key] += h * d[key];
  return next;
};

test('the caller must choose upstream history or an explicitly unfiltered leading pair', () => {
  const { boundary, ...withoutBoundary } = base;
  assert.throws(() => evaluate(withoutBoundary), /explicit/);
  const value = evaluate({ ...base, boundary: { kind: 'unfiltered-first-two' } });
  assert.deepEqual(value.filtered, [false, false, true, true, true]);
  assert.deepEqual(value.speeds.slice(0, 2), base.speeds.slice(0, 2));
  assert.deepEqual(value.coefficients.slice(0, 2), [0, 0]);
  assert.deepEqual(value.corrections.slice(0, 2), [0, 0]);
  const d = linearize({ ...base, boundary: { kind: 'unfiltered-first-two' } }).apply({
    speeds: [.3, -.2, 0, 0, 0], machSquared: [1, 2, 0, 0, 0], sectionArc: [.2, -.1, 0, 0, 0] });
  assert.deepEqual(d.speeds.slice(0, 2), [.3, -.2]);
  assert.deepEqual(d.coefficients.slice(0, 2), [0, 0]);
  assert.deepEqual(d.corrections.slice(0, 2), [0, 0]);
  assert.deepEqual(evaluate(base).filtered, Array(5).fill(true));
});

test('nonuniform affine and quadratic profiles have independent polynomial solutions through the inlet', () => {
  // At mean M²=1, MUCON=gamma, Mcrit=.8: mu=.2*ln(2).
  const coefficient = .2 * Math.log(2), arc = base.boundary.sectionArc.concat(base.sectionArc);
  for (const c of [0, -.07, .12]) {
    const speeds = polynomial(arc, c);
    for (const kind of ['upstream-history', 'unfiltered-first-two']) {
      const boundary = kind === 'upstream-history' ? { kind, speeds: speeds.slice(0, 2),
        machSquared: [1, 1], sectionArc: arc.slice(0, 2) } : { kind };
      const value = evaluate({ speeds: speeds.slice(2), machSquared: Array(5).fill(1),
        sectionArc: arc.slice(2), boundary, mucon: 1.4, gamma: 1.4, mcrit: .8 });
      for (let i = 0; i < 5; i++) {
        const j = i + 2, h = arc[j] - arc[j - 1], previous = arc[j - 1] - arc[j - 2];
        const expected = kind === 'unfiltered-first-two' && i < 2 ? 0 : -coefficient * c * h * (h + previous);
        close(value.corrections[i], expected); close(value.speeds[i], speeds[j] + expected);
      }
    }
  }
});

test('a physical speed impulse has only the documented upstream stencil footprint, with no recursive tail', () => {
  const amplitude = .4, background = 2, coefficient = .2 * Math.log(2), n = 8;
  for (const mucon of [-1.4, 1.4]) {
    const speeds = Array(n).fill(background); speeds[2] += amplitude;
    const p = { speeds, machSquared: Array(n).fill(1), sectionArc: Array.from({ length: n }, (_, i) => i),
      boundary: { kind: 'upstream-history', speeds: [background, background], machSquared: [1, 1], sectionArc: [-2, -1] },
      mucon, gamma: 1.4, mcrit: .8 };
    const value = evaluate(p), second = mucon > 0;
    const weights = [1 - coefficient, second ? 2 * coefficient : coefficient, second ? -coefficient : 0];
    value.speeds.forEach((q, i) => close(q, background + amplitude * (weights[i - 2] ?? 0)));
    const delta = Array(n).fill(0); delta[2] = 1;
    linearize(p).apply({ speeds: delta }).speeds.forEach((v, i) => close(v, weights[i - 2] ?? 0));
  }
});

test('splitting the streamtube preserves every shared filtered section when physical history is transferred', () => {
  for (const mucon of [-1.1, 1.1]) {
    const p = { ...base, mucon }, complete = evaluate(p), split = 3;
    const a = { ...p }, b = { ...p, boundary: { kind: 'upstream-history' } };
    for (const key of keys) {
      a[key] = p[key].slice(0, split); b[key] = p[key].slice(split);
      b.boundary[key] = p[key].slice(split - 2, split);
    }
    const first = evaluate(a), second = evaluate(b);
    for (const key of outputKeys) assert.deepEqual(first[key].concat(second[key]), complete[key]);
  }
});

test('supplied upstream states have the exact finite directional footprint', () => {
  const p = { ...base, machSquared: Array(5).fill(1), mucon: 1.4, gamma: 1.4, mcrit: .8,
    boundary: { ...base.boundary, machSquared: [1, 1] } };
  const coefficient = .2 * Math.log(2), { apply } = linearize(p);
  const old = apply({ boundary: { speeds: [1, 0] } });
  old.speeds.forEach((v, i) => close(v, i === 0 ? -coefficient * 2 / 7 : 0));
  const previous = apply({ boundary: { speeds: [0, 1] } });
  previous.speeds.forEach((v, i) => close(v, i === 0 ? coefficient * 9 / 7 : i === 1 ? -1.5 * coefficient : 0));
  const mach = apply({ boundary: { machSquared: [0, 1] } });
  mach.coefficients.forEach((v, i) => close(v, i === 0 ? .25 : 0));
  const unusedMach = apply({ boundary: { machSquared: [1, 0] } });
  for (const key of outputKeys) unusedMach[key].forEach(v => close(v, 0));
});

test('full section, upstream-history and control derivatives match fourth-order differences', t => {
  let comparisons = 0, maximumNormalizedError = 0;
  for (const mucon of [-.8, 1.1]) for (const kind of ['upstream-history', 'unfiltered-first-two']) {
    const p = { ...base, mucon, boundary: kind === 'upstream-history' ? base.boundary : { kind } };
    const { value, apply } = linearize(p); assert.deepEqual(value, evaluate(p));
    const directions = [];
    for (const key of keys) for (let i = 0; i < p.speeds.length; i++)
      directions.push({ [key]: p.speeds.map((_, j) => i === j ? 1 : 0) });
    if (kind === 'upstream-history') for (const key of keys) for (let i = 0; i < 2; i++)
      directions.push({ boundary: { [key]: [Number(i === 0), Number(i === 1)] } });
    for (const key of controls) directions.push({ [key]: 1 });
    directions.push({ speeds: [.1, -.2, .3, -.4, .2], machSquared: [.05, -.1, .15, -.1, .2],
      sectionArc: [.1, -.04, .03, -.05, .2], gamma: .11, mcrit: -.03, mucon: .08,
      ...(kind === 'upstream-history' ? { boundary: { speeds: [.2, -.1], machSquared: [.03, -.04], sectionArc: [.03, -.1] } } : {}) });
    for (const direction of directions) {
      const h = 2e-5, values = [-2, -1, 1, 2].map(k => evaluate(shift(p, direction, k * h))), actual = apply(direction);
      for (const key of outputKeys) for (let i = 0; i < p.speeds.length; i++) {
        const expected = (values[0][key][i] - 8 * values[1][key][i] + 8 * values[2][key][i] - values[3][key][i]) / (12 * h);
        maximumNormalizedError = Math.max(maximumNormalizedError,
          Math.abs(actual[key][i] - expected) / Math.max(1, Math.abs(actual[key][i]), Math.abs(expected)));
        comparisons++; close(actual[key][i], expected, 3e-8);
      }
    }
  }
  t.diagnostic(JSON.stringify({ comparisons, maximumNormalizedError }));
});

test('arc translation and geometric scaling preserve the value and have zero directional change', () => {
  const value = evaluate(base), { apply } = linearize(base);
  for (const scale of [1e-8, .2, 7, 1e8]) {
    const p = { ...base, sectionArc: base.sectionArc.map(s => scale * (s + 3)),
      boundary: { ...base.boundary, sectionArc: base.boundary.sectionArc.map(s => scale * (s + 3)) } };
    const next = evaluate(p);
    for (const key of outputKeys) next[key].forEach((v, i) => close(v, value[key][i]));
  }
  for (const tangent of [{ sectionArc: Array(5).fill(1), boundary: { sectionArc: [1, 1] } },
    { sectionArc: base.sectionArc, boundary: { sectionArc: base.boundary.sectionArc } }])
    for (const key of outputKeys) apply(tangent)[key].forEach(v => close(v, 0));
});

test('input and history arrays are copied, including typed solver arrays', () => {
  const p = structuredClone(base), before = structuredClone(p);
  const linearization = linearize(p), reference = linearization.apply({ speeds: [.2, 0, 0, 0, 0] });
  assert.deepEqual(p, before);
  p.speeds[0] = 100; p.boundary.speeds[1] = 100; p.sectionArc[0] = -100;
  assert.deepEqual(linearization.value, evaluate(before));
  assert.deepEqual(linearization.apply({ speeds: [.2, 0, 0, 0, 0] }), reference);
  const typed = { ...before, boundary: { ...before.boundary } };
  for (const key of keys) {
    typed[key] = Float64Array.from(before[key]); typed.boundary[key] = Float64Array.from(before.boundary[key]);
  }
  assert.deepEqual(evaluate(typed), evaluate(before));
});

test('invalid arcs/history and malformed tangents fail explicitly, without hidden extrapolation', () => {
  for (const change of [{ speeds: [] }, { machSquared: [1] }, { sectionArc: [0, .3, .3, 1.7, 2.1] },
    { sectionArc: [0, .3, NaN, 1.7, 2.1] }, { gamma: 1 },
    { boundary: { kind: 'upstream-history' } },
    { boundary: { ...base.boundary, sectionArc: [-.9, 0] } },
    { boundary: { ...base.boundary, speeds: [1, -1] } },
    { boundary: { kind: 'unfiltered-first-two', speeds: [1, 1] } }])
    assert.throws(() => evaluate({ ...base, ...change }));
  const { apply } = linearize(base);
  for (const tangent of [{ speeds: [1] }, { boundary: { sectionArc: [0, NaN] } }, { mucon: Infinity }])
    assert.throws(() => apply(tangent));
  assert.throws(() => linearize({ ...base, boundary: { kind: 'unfiltered-first-two' } }).apply({ boundary: { speeds: [1, 0] } }));
  // If no filtered row exists, the explicit closure is simply qtilde=q;
  // no otherwise-unused sonic activation is evaluated or invented.
  const one = { speeds: [.5], machSquared: [1], sectionArc: [0], mcrit: 1,
    boundary: { kind: 'unfiltered-first-two' } };
  assert.deepEqual(linearize(one).apply({ speeds: [.2] }), { speeds: [.2], coefficients: [0], corrections: [0] });
  assert.throws(() => evaluate({ ...one, gamma: 1 }));
});
