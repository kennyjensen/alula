import test from 'node:test';
import assert from 'node:assert/strict';
import { tracePotentialCurve } from '../src/inviscid/potential-curve.js';

const close = (a, b, tolerance = 2e-8) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);

test('potential-coordinate tracing recovers inclined uniform flow and locates an actual plane event', () => {
  const velocity = () => ({ u: 1.3, v: -.4 }), seed = { x: .2, y: .3 }, speedSquared = 1.85;
  for (const endPotential of [-2, 3]) {
    const r = tracePotentialCurve({ seed, initialPotential: .5, endPotential, velocity, tolerance: 1e-11 });
    assert.equal(r.converged, true, r.reason);
    for (const p of r.points) { close(p.x, seed.x + 1.3 * (p.potential - .5) / speedSquared, 1e-13); close(p.y, seed.y - .4 * (p.potential - .5) / speedSquared, 1e-13); }
  }
  const r = tracePotentialCurve({ seed, endPotential: -10, velocity, plane: { normal: { x: 2, y: -1 }, offset: -3 }, tolerance: 1e-11 });
  assert.equal(r.converged, true); assert.equal(r.reason, 'plane');
  const p = r.points.at(-1); close(2 * p.x - p.y, -3, 3e-11);
  close(p.potential, 1.3 * (p.x - seed.x) - .4 * (p.y - seed.y), 2e-11);
});

test('curved potential traces preserve the exact irrotational vortex and converge under integration refinement', () => {
  const circulationScale = .7, radius = .4, endPotential = 1.8;
  const velocity = p => ({ u: -circulationScale * p.y / (p.x * p.x + p.y * p.y), v: circulationScale * p.x / (p.x * p.x + p.y * p.y) });
  const errors = [];
  for (const tolerance of [1e-5, 1e-7, 1e-9]) {
    const r = tracePotentialCurve({ seed: { x: radius, y: 0 }, endPotential, velocity, tolerance, maxStep: .4, maxSpatialStep: .4 });
    assert.equal(r.converged, true, r.reason);
    const p = r.points.at(-1), angle = endPotential / circulationScale;
    errors.push(Math.hypot(p.x - radius * Math.cos(angle), p.y - radius * Math.sin(angle)));
    for (const p of r.points) close(Math.hypot(p.x, p.y), radius, 10 * tolerance);
  }
  assert.ok(errors[1] < .15 * errors[0], String(errors)); assert.ok(errors[2] < .15 * errors[1], String(errors));
  assert.ok(errors.at(-1) < 2e-8, String(errors));
});

test('potential tracing respects a stagnation singularity, geometry rejection, and explicit incomplete traces', () => {
  const velocity = p => ({ u: -p.x, v: p.y });
  const r = tracePotentialCurve({ seed: { x: -.01, y: 0 }, endPotential: -.3, velocity, plane: { normal: { x: 1, y: 0 }, offset: -.5 }, tolerance: 1e-10 });
  assert.equal(r.converged, true, r.reason);
  close(r.points.at(-1).x, -.5, 2e-10); close(r.points.at(-1).potential, -.5 * (.25 - .0001), 5e-8);
  const stopped = tracePotentialCurve({ seed: { x: 0, y: 0 }, endPotential: 1, velocity });
  assert.equal(stopped.converged, false); assert.match(stopped.reason, /singular/);
  const limited = tracePotentialCurve({ seed: { x: -.1, y: 0 }, endPotential: -1, velocity, maxSteps: 1 });
  assert.equal(limited.converged, false); assert.match(limited.reason, /step limit/);
  const outside = tracePotentialCurve({ seed: { x: -.1, y: 0 }, endPotential: -1, velocity, admissible: p => p.x > -.2 });
  assert.equal(outside.converged, false); assert.match(outside.reason, /admissible/);
  const missed = tracePotentialCurve({ seed: { x: -.1, y: 0 }, endPotential: -.01, velocity, plane: { normal: { x: 1, y: 0 }, offset: -2 } });
  assert.equal(missed.converged, false); assert.match(missed.reason, /potential limit/);
});
