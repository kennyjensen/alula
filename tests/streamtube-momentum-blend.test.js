// SPDX-License-Identifier: GPL-2.0-or-later
// Local algebra for our explicitly chosen ISMOM4-style switch. These tests
// do not compute shocks or establish conservative hybrid/airfoil solutions.
import test from 'node:test';
import assert from 'node:assert/strict';
import { evaluateStreamtubeMomentumBlend as evaluate,
  linearizeStreamtubeMomentumBlend as linearize } from '../src/euler/streamtube-momentum-blend.js';

const close = (a, b, tolerance = 3e-13) => assert.ok(
  Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const base = { states: [{ rho: 1, q: 2, p: 1 }, { rho: 1, q: 1, p: 1 }],
  transportSpeeds: [2, 1], streamwiseResidual: 10, isentropicResidual: -6, epsilonP: 1 / 16 };
const numbers = value => [value.fraction, value.residual, ...value.lossIndicators, value.compression];
const scalarKeys = ['streamwiseResidual', 'isentropicResidual', 'epsilonP'];
const stateKeys = ['rho', 'q', 'p'];
const shifted = (p, d, h) => {
  const r = structuredClone(p);
  for (let i = 0; i < 2; i++) {
    for (const key of stateKeys) r.states[i][key] += h * (d.states?.[i]?.[key] ?? 0);
    r.transportSpeeds[i] += h * (d.transportSpeeds?.[i] ?? 0);
  }
  for (const key of scalarKeys) r[key] += h * (d[key] ?? 0);
  return r;
};
const fd4 = (p, tangent, h) => {
  const a = numbers(evaluate(shifted(p, tangent, -2 * h)));
  const b = numbers(evaluate(shifted(p, tangent, -h)));
  const c = numbers(evaluate(shifted(p, tangent, h)));
  const d = numbers(evaluate(shifted(p, tangent, 2 * h)));
  return a.map((x, i) => (x - 8 * b[i] + 8 * c[i] - d[i]) / (12 * h));
};
// Analytically prescribe phi_a=2, phi_b=1, with normalized loss and
// compression coordinates independent of each other. Positive gas inputs
// exercise partials; they are not a solved perfect-gas shock state.
const controlled = (a, b, c) => {
  const p = structuredClone(base), q = Math.exp(c * Math.sqrt(p.epsilonP));
  p.states[0] = { rho: 1, q, p: q / 2 };
  p.transportSpeeds = [q + a * p.epsilonP, 1 + b * p.epsilonP];
  return p;
};

test('zero speed bias selects the exact entropy residual and Jacobian', () => {
  for (const residual of [-0, Number.MIN_VALUE, -Math.PI]) {
    const p = { ...base, streamwiseResidual: Number.MAX_VALUE, isentropicResidual: residual };
    const l = linearize(p);
    assert.equal(l.value.fraction, 0); assert.equal(l.value.residual, residual);
    // Even independent perturbations that introduce bias have zero blend
    // derivative at its compact support boundary, so no R1 contamination.
    const d = l.apply({ states: [{ rho: .2, q: -.1 }, { p: .3 }], transportSpeeds: [.4, -.2],
      streamwiseResidual: Number.MAX_VALUE, isentropicResidual: -Math.E, epsilonP: .02 });
    assert.equal(d.fraction, 0); assert.equal(d.residual, -Math.E);
  }
});

test('full compression and endpoint loss select exact momentum values and derivatives', () => {
  for (const loss of [[1, 0], [0, 1], [1.5, 2]]) {
    const p = controlled(...loss, 2);
    p.streamwiseResidual = -0; p.isentropicResidual = Number.MAX_VALUE;
    const l = linearize(p);
    assert.equal(l.value.fraction, 1); assert.equal(l.value.residual, -0);
    const d = l.apply({ states: [{ q: .2 }, { rho: -.3, p: .1 }], transportSpeeds: [-.2, .3],
      streamwiseResidual: Math.PI, isentropicResidual: Number.MAX_VALUE, epsilonP: -.01 });
    assert.equal(d.fraction, 0); assert.equal(d.residual, Math.PI);
  }
});

test('acceleration, zero compression and zero phi difference suppress the momentum blend', () => {
  for (const c of [-2, -.1, 0]) {
    const l = linearize(controlled(2, 2, c));
    assert.equal(l.value.fraction, 0); assert.equal(l.value.residual, base.isentropicResidual);
    assert.equal(l.apply({ states: [{ q: 1 }, {}] }).fraction, 0);
  }
  const p = structuredClone(base);
  p.states[1].p = .5; p.transportSpeeds = [2.2, 1.3];
  const value = evaluate(p);
  assert.equal(value.fraction, 0);
  assert.ok(value.lossIndicators.every(x => x === 0));
});

test('known midpoint weights and endpoint OR agree with exact polynomial values', () => {
  // B(1/2)=1/2 by symmetry; OR of two halves is 3/4. This distinguishes
  // the endpoint combination from averaging or multiplying its losses.
  for (const [a, b, c, fraction] of [[.5, 0, 2, .5], [.5, .5, 2, .75], [1, 1, .5, .5], [.5, .5, .5, .375]]) {
    const value = evaluate(controlled(a, b, c));
    close(value.fraction, fraction); close(value.residual, 16 * fraction - 6);
    close(value.lossIndicators[0], a / 16); close(value.lossIndicators[1], b / 16);
    close(value.compression, c / 4);
  }
  // B(1/4)=53/512 and B(3/4)=459/512; independent rational values
  // distinguish our quintic from a merely C1 cubic smoothstep.
  close(evaluate(controlled(.25, 0, 2)).fraction, 53 / 512);
  close(evaluate(controlled(.75, 0, 2)).fraction, 459 / 512);
});

test('opposite endpoint biases cannot cancel a significant positive loss', () => {
  const p = { ...base, transportSpeeds: [2 + 1 / 16, 1 - 1 / 16] };
  const value = evaluate(p);
  assert.equal((p.transportSpeeds[0] - 2) + (p.transportSpeeds[1] - 1), 0);
  assert.deepEqual(value.lossIndicators, [1 / 16, -1 / 16]);
  assert.equal(value.fraction, 1); assert.equal(value.residual, base.streamwiseResidual);
  // The bias sign still matters; negative loss alone cannot activate R1.
  assert.equal(evaluate({ ...p, transportSpeeds: [2 - 1 / 16, 1 - 1 / 16] }).fraction, 0);
});

test('all state, speed, residual and tolerance partials match fourth-order differences', t => {
  const directions = [];
  for (let i = 0; i < 2; i++) for (const key of stateKeys)
    directions.push({ states: Array.from({ length: 2 }, (_, j) => i === j ? { [key]: 1 } : {}) });
  for (let i = 0; i < 2; i++) directions.push({ transportSpeeds: [Number(i === 0), Number(i === 1)] });
  for (const key of scalarKeys) directions.push({ [key]: 1 });
  directions.push({ states: [{ rho: .12, q: -.07, p: .05 }, { rho: -.03, q: .08, p: -.13 }],
    transportSpeeds: [.04, -.09], streamwiseResidual: .11, isentropicResidual: -.17, epsilonP: .012 });
  let comparisons = 0, maximumAbsoluteError = 0, maximumNormalizedError = 0;
  for (const controls of [[.37, .69, .41], [-.2, .65, 2], [1.4, .6, .55], [-.3, -.2, .5]]) {
    const p = controlled(...controls), l = linearize(p);
    assert.deepEqual(l.value, evaluate(p));
    for (const d of directions) {
      const actual = numbers(l.apply(d)), expected = fd4(p, d, 2e-6);
      actual.forEach((v, i) => {
        const error = Math.abs(v - expected[i]);
        maximumAbsoluteError = Math.max(maximumAbsoluteError, error);
        maximumNormalizedError = Math.max(maximumNormalizedError, error / Math.max(1, Math.abs(v), Math.abs(expected[i])));
        comparisons++; close(v, expected[i], 3e-8);
      });
    }
  }
  t.diagnostic(JSON.stringify({ comparisons, maximumAbsoluteError, maximumNormalizedError }));
});

test('switch endpoints have zero first and continuous second derivatives', () => {
  // Probe across both endpoints, including the constant branches. For a
  // C2 compact gate, its second symmetric difference approaches zero there.
  for (const kind of ['loss', 'compression']) for (const knot of [0, 1]) {
    const p = kind === 'loss' ? controlled(knot, 0, 2) : controlled(2, 2, knot);
    const tangent = kind === 'loss' ? { transportSpeeds: [p.epsilonP, 0] }
      : { states: [{ q: p.states[0].q * Math.sqrt(p.epsilonP) }, {}] };
    const l = linearize(p);
    close(l.apply(tangent).fraction, 0, 1e-24);
    close(fd4(p, tangent, 2e-6)[0], 0, 2e-9);
    const seconds = [2e-4, 1e-4].map(h => (
      evaluate(shifted(p, tangent, h)).fraction - 2 * l.value.fraction
      + evaluate(shifted(p, tangent, -h)).fraction) / (h * h));
    assert.ok(Math.abs(seconds[1]) < .6 * Math.abs(seconds[0]) + 1e-7);
    assert.ok(Math.abs(seconds[1]) < .003);
  }
});

test('consistent velocity and density unit changes preserve sensors and scale pressure residuals', () => {
  const p = controlled(.31, .62, .57), original = evaluate(p);
  for (const velocity of [1e-4, 3, 1e4]) for (const density of [.2, 7]) {
    const pressure = density * velocity ** 2;
    const mapped = evaluate({ ...p, states: p.states.map(s => ({ rho: density * s.rho,
      q: velocity * s.q, p: pressure * s.p })), transportSpeeds: p.transportSpeeds.map(q => velocity * q),
      streamwiseResidual: pressure * p.streamwiseResidual, isentropicResidual: pressure * p.isentropicResidual });
    close(mapped.fraction, original.fraction, 3e-12); close(mapped.compression, original.compression, 3e-12);
    mapped.lossIndicators.forEach((v, i) => close(v, original.lossIndicators[i], 3e-12));
    close(mapped.residual / pressure, original.residual, 3e-12);
  }
  for (const [velocity, density] of [[1, 0], [0, 1]]) {
    const pressure = density + 2 * velocity;
    const d = linearize(p).apply({ states: p.states.map(s => ({ rho: density * s.rho,
      q: velocity * s.q, p: pressure * s.p })), transportSpeeds: p.transportSpeeds.map(q => velocity * q),
      streamwiseResidual: pressure * p.streamwiseResidual, isentropicResidual: pressure * p.isentropicResidual });
    close(d.fraction, 0); close(d.compression, 0); d.lossIndicators.forEach(v => close(v, 0));
    close(d.residual, pressure * original.residual);
  }
});

test('mandatory tolerance, physical domains and tangents are validated without input mutation', () => {
  const p = controlled(.3, .7, .4), original = structuredClone(p);
  const l = linearize(p), expected = l.apply({ transportSpeeds: [.1, -.2], epsilonP: .01 });
  evaluate(p); assert.deepEqual(p, original);
  p.states[0].q = 99; p.transportSpeeds[1] = 12; p.epsilonP = 10;
  l.value.fraction = 0; l.value.lossIndicators[0] = 100; l.value.compression = 30;
  assert.deepEqual(l.apply({ transportSpeeds: [.1, -.2], epsilonP: .01 }), expected);
  assert.ok(numbers(linearize(original).apply()).every(v => v === 0));
  for (const bad of [{ epsilonP: undefined }, { epsilonP: 0 }, { epsilonP: -.1 }, { epsilonP: NaN },
    { epsilonP: Infinity }, { states: [original.states[0]] }, { states: [{ rho: 1, q: 0, p: 1 }, original.states[1]] },
    { states: [{ rho: -1, q: 1, p: 1 }, original.states[1]] }, { states: [{ rho: 1, q: 1, p: NaN }, original.states[1]] },
    { transportSpeeds: [-1, 1] }, { transportSpeeds: [1] }, { streamwiseResidual: Infinity }, { isentropicResidual: NaN }]) {
    assert.throws(() => evaluate({ ...original, ...bad }), /momentum-blend/);
    assert.throws(() => linearize({ ...original, ...bad }), /momentum-blend/);
  }
  for (const d of [{ states: [{}] }, { states: [{ q: NaN }, {}] }, { states: [null, {}] },
    { transportSpeeds: [0, Infinity] }, { epsilonP: NaN }, { isentropicResidual: Infinity }])
    assert.throws(() => linearize(original).apply(d), /momentum-blend tangent/);
});
