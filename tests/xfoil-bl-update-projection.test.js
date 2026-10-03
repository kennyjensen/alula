import test from 'node:test';
import assert from 'node:assert/strict';
import { projectXfoilDisplacement, proposeCoupledDensityNewtonWithXfoilDstrProjection } from '../scripts/validation/xfoil-bl-update-projection.js';
import { proposeCoupledDensityNewton } from '../src/euler/streamtube-density-newton.js';
import { dslim } from '../src/viscous/xfoil/xbl.js';
import { hkin } from '../src/viscous/xfoil/xblsys.js';

test('DSLIM independently reproduces the native scalar on surface, wake and finite gap states', t => {
  let count = 0, maximumRelativeError = 0;
  for (const wake of [false, true]) for (const mach of [0, .185, .6]) for (const targetHk of [.99, 1.00001, 1.3]) {
    const theta = .0013, ue = 1.3, gamma = 1.4, wakeGap = wake ? .0024 : 0;
    const factor = mach ** 2 / (1 + .5 * (gamma - 1) * mach ** 2);
    const m2 = factor * ue ** 2 / (1 - .5 * (gamma - 1) * factor * ue ** 2);
    const deltaStar = wakeGap + theta * ((1 + .113 * m2) * targetHk + .29 * m2);
    const input = { theta, deltaStar, ue, wakeGap, wake, mach, gamma }, result = projectXfoilDisplacement(input);
    const native = dslim({ hkin }, deltaStar - wakeGap, theta, ue, m2, wake ? 1.00005 : 1.02) + wakeGap;
    const error = Math.abs(native - result.deltaStar) / native; maximumRelativeError = Math.max(maximumRelativeError, error);
    assert.ok(error < 2e-15); assert.ok(result.projectedRawHk >= result.minimumHk - 4e-15);
    if (targetHk === 1.3) { assert.equal(result.active, false); assert.equal(result.deltaStar, deltaStar); }
    assert.equal(result.massDefectChange, ue * (result.deltaStar - deltaStar));
    const scaled = projectXfoilDisplacement({ ...input, theta: theta * 1000, deltaStar: deltaStar * 1000, wakeGap: wakeGap * 1000 });
    assert.ok(Math.abs(scaled.deltaStar / 1000 - result.deltaStar) < 4e-18);
    const again = projectXfoilDisplacement({ ...input, deltaStar: result.deltaStar });
    assert.ok(Math.abs(again.deltaStar - result.deltaStar) < 4e-18); count++;
  }
  t.diagnostic(JSON.stringify({ cases: count, maximumRelativeError }));
});

function stub() {
  const system = { ne: 1, n: 9, bl: { scale: .01, hasFiniteBase: true, stations: [{ id: 0, kind: 'surface' }, { id: 1, kind: 'wake' }],
    kernel: { parameters: { mach: .185, gamma: 1.4 } }, thicknesses: x => Array.from(x),
    geometry: () => ({ coordinates: [{}, { wakeGap: .002 }] }) }, euler: {
    layout: { n: 1, densityCount: 1, independentWakeBanks: true, globals: { stagnation: [] }, positions: [] }, curves: [],
    decode: () => ({ nodes: [] }), setDisplacement: x => { system.lastDisplacement = x; } } };
  return system;
}

test('research wrapper projects only displacement after the unchanged density Newton proposal', () => {
  const s = stub(), state = Float64Array.from([0, 2, 1, 1.08, 1, .04, 1, 1.5, 1]);
  const direction = Float64Array.from([.1, .3, 0, -.07, .01, .005, 0, -.295, 0]);
  const before = Array.from(state), d = Array.from(direction), ordinary = proposeCoupledDensityNewton(s, state, direction);
  const p = proposeCoupledDensityNewtonWithXfoilDstrProjection(s, state, direction);
  assert.equal(p.step, ordinary.step); assert.deepEqual(p.limiter, ordinary.limiter);
  assert.equal(p.projection.changes.length, 2); assert.equal(p.meritComparable, false);
  for (let k = 0; k < s.n; k++) if (![3, 7].includes(k)) assert.equal(p.x[k], ordinary.x[k]);
  assert.deepEqual(Array.from(state), before); assert.deepEqual(Array.from(direction), d);
  assert.deepEqual(s.lastDisplacement, before.slice(1));
  assert.equal(p.projection.changes[1].wakeGap, .2);
});

test('inactive states are identical and invalid thermal or unsupported moving-gap cases reject', () => {
  const s = stub(), x = Float64Array.from([0, 2, 1, 2, 1, .04, 1, 2, 1]), direction = new Float64Array(9);
  const p = proposeCoupledDensityNewtonWithXfoilDstrProjection(s, x, direction);
  assert.deepEqual(p.x, x); assert.equal(p.projection.active, false); assert.equal(p.meritComparable, undefined);
  s.euler.layout.independentWakeBanks = false;
  assert.throws(() => proposeCoupledDensityNewtonWithXfoilDstrProjection(s, x, direction), /independently solved/);
  assert.throws(() => projectXfoilDisplacement({ theta: .001, deltaStar: .002, ue: 100, mach: .2 }), /thermal/);
});
