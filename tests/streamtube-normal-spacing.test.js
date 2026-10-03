import test from 'node:test';
import assert from 'node:assert/strict';
import { stagnationStreamtubeMass, createNormalMassWeights, resolveNormalMassWeights } from '../src/geometry/streamtube-normal-spacing.js';
import { tracePotentialCurve } from '../src/inviscid/potential-curve.js';
import { fitNormalStreamtubeFraction } from '../src/geometry/streamtube-spacing-fit.js';
const near = (a, b, tolerance = 2e-12) => assert.ok(Math.abs(a - b) < tolerance, `${a} differs from ${b}`);

test('normal refinement resolves a frozen 45-fold mass jump without clipping either wall target or changing total mass', () => {
  const first = .0025590586249133263, last = .01773570456912759, input = { count: 5, first, last };
  const original = createNormalMassWeights(input);
  assert.ok(original[1] / original[0] > 44);
  const r = resolveNormalMassWeights({ ...input, maximumAdjacentRatio: 3 });
  assert.equal(r.requestedTubes, 5); assert.equal(r.weights.length, 14);
  assert.equal(r.weights[0], first); assert.equal(r.weights.at(-1), last);
  near(r.weights.reduce((s, w) => s + w, 0), 1);
  for (let i = 1; i < r.weights.length; i++) assert.ok(Math.abs(Math.log(r.weights[i]) - Math.log(r.weights[i - 1])) <= Math.log(3));
  resolveNormalMassWeights({ count: 5, first: last, last: first, maximumAdjacentRatio: 3 }).weights.forEach((w, i) => near(w, r.weights.at(-1 - i), 2e-14));
  assert.deepEqual(resolveNormalMassWeights(input).weights, original);
  assert.throws(() => resolveNormalMassWeights({ ...input, maximumAdjacentRatio: 3, maximumTubes: 5 }), /budget/);
});

test('normal spacing recovers exact linear-stagnation streamlines and oblique potential cross-lines', () => {
  for (const rate of [.3, 3, 100]) for (const slope of [-3, 0, 2]) for (const sign of [-1, 1]) {
    const distance = .04, density = 1.3;
    const mass = stagnationStreamtubeMass({ strainRate: rate, distance, potentialSlope: slope, density });
    // Independent exact-flow oracle: integrate (u,v)=(K*x,-K*y) from
    // x=+/-1, y>0 on its hyperbolic streamline, to the requested potential.
    const psi = sign * mass / density, seed = { x: sign, y: mass / density / rate };
    const phi = p => .5 * rate * (p.x * p.x - p.y * p.y);
    const r = tracePotentialCurve({ seed, initialPotential: phi(seed), endPotential: slope * psi,
      velocity: p => ({ u: rate * p.x, v: -rate * p.y }), tolerance: 1e-11,
      maxStep: rate * .02, maxSpatialStep: .05, minStep: 1e-16 });
    assert.equal(r.converged, true, r.reason);
    const p = r.points.at(-1);
    near(Math.hypot(p.x, p.y), distance, 2e-8); near(rate * p.x * p.y, psi, 1e-9);
    near(phi(p), slope * psi, 1e-9);
  }
});

test('measured normal-spacing fit recovers exact curved-cylinder potential flow beyond the linear strain approximation', () => {
  // Uniform flow past the unit circle: F(z)=z+1/z. On phi=-2,
  // solve z^2-F*z+1=0 on the exterior branch, at streamfunction psi>0.
  const point = psi => {
    const modulus = Math.hypot(psi * psi, 4 * psi);
    const real = Math.sqrt(.5 * (modulus - psi * psi)), imag = Math.sqrt(.5 * (modulus + psi * psi));
    return { x: -1 - .5 * real, y: .5 * (psi + imag) };
  };
  for (const targetDistance of [.01, .15, .6, 1.1]) {
    const initialFraction = stagnationStreamtubeMass({ strainRate: 2, distance: targetDistance });
    const r = fitNormalStreamtubeFraction({ initialFraction, targetDistance, relativeTolerance: 1e-7,
      distanceAtFraction: f => { const p = point(f); return Math.hypot(p.x + 1, p.y); } });
    const p = point(r.fraction), r2 = p.x * p.x + p.y * p.y;
    near(Math.hypot(p.x + 1, p.y), targetDistance, 1e-7 * targetDistance);
    near(p.x * (1 + 1 / r2), -2); near(p.y * (1 - 1 / r2), r.fraction);
    assert.ok(r.history.length > 1); assert.ok(r.fraction < initialFraction);
    if (initialFraction < 1) assert.equal(r.history[0].fraction, initialFraction, 'Retain a valid asymptotic initial estimate.');
    else {
      near(initialFraction, 1.21);
      assert.ok(r.history[0].fraction > 0 && r.history[0].fraction < 1, 'Measure a feasible interior trial before fitting.');
      assert.ok(Math.abs(r.fraction - r.history[0].fraction) > .1, 'The interior trial must not become a clipped accepted mass.');
    }
  }
  const trials = [];
  assert.throws(() => fitNormalStreamtubeFraction({ initialFraction: 4, targetDistance: 2,
    distanceAtFraction: f => { trials.push(f); const p = point(f); return Math.hypot(p.x + 1, p.y); } }), /did not converge/);
  assert.ok(trials.length > 1 && trials.every(f => f > 0 && f < 1), 'An unattainable target remains a physical fit failure.');
  assert.ok(Math.hypot(point(1).x + 1, point(1).y) < 2, 'Even the full available streamfunction span cannot meet the requested distance.');
  assert.throws(() => fitNormalStreamtubeFraction({ initialFraction: .1, targetDistance: .1, distanceAtFraction: () => NaN }), /Invalid measured/);
  assert.throws(() => fitNormalStreamtubeFraction({ initialFraction: .1, targetDistance: 2, distanceAtFraction: () => 1 }), /not monotone/);
  assert.throws(() => fitNormalStreamtubeFraction({ initialFraction: .1, targetDistance: .1, distanceAtFraction: () => { throw new Error('body collision'); } }), /body collision/);
});

test('stagnation mass scales with strain, density and squared cell distance without a mass floor', () => {
  const controls = { strainRate: 7, distance: .01, potentialSlope: .8, density: 2 };
  const m = stagnationStreamtubeMass(controls);
  near(stagnationStreamtubeMass({ ...controls, distance: .005 }), m / 4);
  near(stagnationStreamtubeMass({ ...controls, density: 6 }), 3 * m);
  // Double geometric units at unchanged velocity: K halves, mass doubles.
  near(stagnationStreamtubeMass({ ...controls, distance: .02, strainRate: 3.5 }), 2 * m);
  assert.throws(() => stagnationStreamtubeMass({ ...controls, strainRate: 0 }), /Invalid/);
  assert.throws(() => stagnationStreamtubeMass({ ...controls, distance: 1e-300 }), /Unresolved/);
});

test('normal mass distributions preserve asymmetric endpoint constraints, positive total and reflection', () => {
  for (const count of [3, 7, 19]) for (const [first, last] of [[.0002, .001], [.05, .05], [.0002, undefined], [undefined, .03]]) {
    const w = createNormalMassWeights({ count, first, last }), reversed = createNormalMassWeights({ count, first: last, last: first });
    if (first !== undefined) assert.equal(w[0], first); if (last !== undefined) assert.equal(w.at(-1), last);
    near(w.reduce((a, b) => a + b, 0), 1); assert.ok(w.every(v => v > 0));
    w.forEach((v, i) => near(v, reversed[count - 1 - i]));
    const log = w.map(Math.log);
    if (first === undefined || last === undefined) for (let i = 1; i < count - 1; i++) near(log[i + 1] - 2 * log[i] + log[i - 1], 0);
    else for (let i = 1; i < count - 2; i++) near(log[i + 2] - 3 * log[i + 1] + 3 * log[i] - log[i - 1], 0);
  }
  for (const controls of [{ count: 2, first: .1 }, { count: 7, first: .7, last: .4 }, { count: 7, last: 0 }])
    assert.throws(() => createNormalMassWeights(controls), /Incompatible/);
});
