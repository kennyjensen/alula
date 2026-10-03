import test from 'node:test';
import assert from 'node:assert/strict';
import { traceArcCurve } from '../src/inviscid/arc-curve.js';

const close = (a, b, tolerance = 1e-10) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);
const last = r => { assert.equal(r.converged, true, r.reason); return r.points.at(-1); };

test('arc tracing measures inclined uniform-flow distance and potential in both directions', () => {
  const seed = { x: .2, y: .3 }, velocity = () => ({ u: 1.3, v: -.4 }), speed = Math.sqrt(1.85);
  for (const endArc of [-2, 3]) {
    const r = traceArcCurve({ seed, initialArc: .5, endArc, initialPotential: 7, velocity, tolerance: 1e-12 });
    last(r);
    for (const p of r.points) {
      close(p.x, seed.x + 1.3 / speed * (p.arc - .5), 2e-13);
      close(p.y, seed.y - .4 / speed * (p.arc - .5), 2e-13);
      close(p.potential, 7 + speed * (p.arc - .5), 2e-13);
    }
  }
  const p = last(traceArcCurve({ seed, endArc: -10, velocity, plane: { normal: { x: 2, y: -1 }, offset: -3 }, tolerance: 1e-11 }));
  close(2 * p.x - p.y, -3, 3e-11);
  close(p.potential, 1.3 * (p.x - seed.x) - .4 * (p.y - seed.y), 2e-11);
});

test('arc tracing refines on an exact vortex and locates a curved plane event without chord projection', () => {
  const radius = .4, strength = .7, endArc = 1.8;
  const velocity = p => ({ u: -strength * p.y / (p.x * p.x + p.y * p.y), v: strength * p.x / (p.x * p.x + p.y * p.y) });
  const errors = [];
  for (const tolerance of [1e-5, 1e-7, 1e-9]) {
    const r = traceArcCurve({ seed: { x: radius, y: 0 }, endArc, velocity, tolerance, maxStep: .4 }), p = last(r);
    errors.push(Math.hypot(p.x - radius * Math.cos(endArc / radius), p.y - radius * Math.sin(endArc / radius)));
    close(p.potential, strength * endArc / radius, 20 * tolerance);
    for (const point of r.points) close(Math.hypot(point.x, point.y), radius, 10 * tolerance);
  }
  assert.ok(errors[1] < .15 * errors[0] && errors[2] < .15 * errors[1], String(errors));
  assert.ok(errors.at(-1) < 2e-8, String(errors));
  const p = last(traceArcCurve({ seed: { x: radius, y: 0 }, endArc, velocity, tolerance: 1e-11,
    plane: { normal: { x: 1, y: 0 }, offset: .1 } }));
  const angle = Math.acos(.1 / radius);
  close(p.arc, radius * angle, 2e-10); close(p.y, radius * Math.sin(angle), 2e-10); close(p.potential, strength * angle, 4e-10);
});

test('arc sampling removes the exact stagnation potential-transfer ratio while retaining harmonic identities', () => {
  const inverse = (phi, psi) => { const r = Math.hypot(phi, psi); return { x: Math.sqrt((r + phi) / 2), y: Math.sqrt((r - phi) / 2) }; };
  const velocity = p => ({ u: 2 * p.x, v: -2 * p.y }), distance = (p, q) => Math.hypot(p.x - q.x, p.y - q.y);
  const errors = [];
  for (const h of [.2, .1, .05, .025]) {
    const phi = -4 * h * h, seed = inverse(phi, .1);
    const end = last(traceArcCurve({ seed, initialPotential: phi, endArc: 2, velocity, tolerance: 1e-12,
      plane: { normal: { x: 1, y: -1 }, offset: 0 } }));
    const middle = last(traceArcCurve({ seed, initialPotential: phi, endArc: end.arc / 2, velocity, tolerance: 1e-12 }));
    for (const p of [middle, end]) { close(2 * p.x * p.y, .1, 2e-10); close(p.x * p.x - p.y * p.y, p.potential, 2e-10); }
    errors.push(Math.abs(distance(seed, middle) / distance(middle, end) - 1));
  }
  assert.ok(errors[0] < .03 && errors.at(-1) < 1e-6, String(errors));
  assert.ok(errors[2] < errors[1] && errors[3] < errors[2], String(errors));
});

test('arc tracing respects velocity units and reports stagnation, geometry and budget failures', () => {
  const base = { seed: { x: .4, y: .1 }, endArc: .7, tolerance: 1e-11 }, velocity = p => ({ u: p.x, v: -p.y });
  const p = last(traceArcCurve({ ...base, velocity }));
  const q = last(traceArcCurve({ ...base, velocity: p => ({ u: 13 * p.x, v: -13 * p.y }), potentialTolerance: 13e-11 }));
  close(p.x, q.x, 2e-11); close(p.y, q.y, 2e-11); close(13 * p.potential, q.potential, 2e-10);
  for (const [controls, reason] of [
    [{ seed: { x: 0, y: 0 } }, /stagnation/],
    [{ admissible: p => p.x < .5 }, /admissible/],
    [{ admissibleSegment: (a, b) => b.x < .5 }, /admissible/],
    [{ maxSteps: 1 }, /step limit/],
    [{ endArc: .01, plane: { normal: { x: 1, y: 0 }, offset: 2 } }, /arc limit/],
  ]) { const r = traceArcCurve({ ...base, velocity, ...controls }); assert.equal(r.converged, false); assert.match(r.reason, reason); }
  assert.throws(() => traceArcCurve({ ...base, velocity, potentialTolerance: 0 }), /controls/);
});
