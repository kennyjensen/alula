import test from 'node:test';
import assert from 'node:assert/strict';
import { projectHalfspaces } from '../src/numerics/halfspace-projection.js';
import { correctSecondOrderStep } from '../src/numerics/second-order-step.js';
import { takeDoglegStep } from '../src/numerics/dogleg.js';

const close = (a, b, tolerance = 1e-12) => assert.ok(Math.abs(a - b) <= tolerance, `${a} != ${b}`);
const curved = (h, units = [1, 1], radius = 2 * h) => {
  const initial = Float64Array.of(1e-8 * units[0], 0), direction = Float64Array.of(-.9e-8 * units[0], h * units[1]);
  const scales = Float64Array.from(units, v => 1 / v);
  const constraints = [{ value: 1e-8, lower: -.9e-8, gradient: new Map([[0, 1 / units[0]]]) }];
  const values = x => Float64Array.of(x[0] / units[0] - (x[1] / units[1]) ** 2);
  return { initial, direction, scales, constraints, values, radius };
};

test('bounded minimum-norm shifted halfspaces match independent analytic solutions', () => {
  const constraints = [{ gradient: new Map([[0, 1], [1, 1]]), lower: 2 },
    { gradient: new Map([[0, 1]]), lower: -.1 }];
  const result = projectHalfspaces(Float64Array.of(0, 0), constraints, { maximumNorm: 2 });
  assert.equal(result.converged, true); result.point.forEach(v => close(v, 1));
  assert.equal(projectHalfspaces([0, 0], constraints, { maximumNorm: 1 }).converged, false);
  assert.throws(() => projectHalfspaces([0], [{ gradient: new Map(), lower: 1 }], { maximumNorm: 1 }), /Inconsistent/);
});

test('curved-boundary correction is quadratic in step size and invariant to variable units', () => {
  for (const units of [[1, 1], [1e3, .02]]) for (const h of [.1, .05, .025]) {
    const input = curved(h, units), initial = input.initial.slice(), direction = input.direction.slice();
    const r = correctSecondOrderStep(input);
    assert.equal(r.corrected, true, r.reason); assert.equal(r.history.length, 1);
    close(r.correctionNorm, h * h); close(r.direction[0] / units[0], -.9e-8 + h * h);
    close(r.direction[1] / units[1], h);
    const candidate = initial.map((v, i) => v + r.direction[i]);
    close(input.values(candidate)[0], 1e-9); assert.ok(input.values(candidate)[0] > 0);
    assert.deepEqual(input.initial, initial); assert.deepEqual(input.direction, direction);
  }
});

test('corrected proposals retain the trust-radius bound and recheck nonlinear margins after scaling', () => {
  const input = curved(.1, [1, 1], Math.hypot(.9e-8, .1));
  const r = correctSecondOrderStep(input);
  assert.equal(r.corrected, true, r.reason); assert.ok(r.history.some(h => h.radiusFraction < 1));
  assert.ok(r.scaledStepNorm <= input.radius * (1 + 1e-14));
  assert.ok(input.values(input.initial.map((v, i) => v + r.direction[i]))[0] > 0);
  assert.ok(r.correctionNorm <= .25 * r.originalNorm);
});

test('oversized, unsupported and failed corrections cannot certify a step', () => {
  const oversized = correctSecondOrderStep(curved(1));
  assert.equal(oversized.corrected, false); assert.match(oversized.reason, /projection failed/);
  const good = curved(.1);
  const unsupported = correctSecondOrderStep({ ...good, values: () => [1] });
  assert.equal(unsupported.corrected, false); assert.match(unsupported.reason, /no supported/);
  for (const values of [() => { throw new Error('outside coordinate chart'); }, () => [NaN], () => []]) {
    const r = correctSecondOrderStep({ ...good, values }); assert.equal(r.corrected, false);
  }
  assert.throws(() => correctSecondOrderStep({ ...good, direction: Float64Array.of(-.01, .1) }), /linearly feasible/);
});

test('second-order trust acceptance rescores the complete corrected step and rolls back rejected active preparations', () => {
  for (const outcome of ['accept', 'residual failure', 'active event']) {
    const initial = Float64Array.of(1e-8, 0), before = initial.slice(), residual = x => x.map(v => 1 + v);
    const currentResidual = residual(initial); let phase = 'base', preparations = 0, corrections = 0;
    const r = takeDoglegStep({ initial, currentResidual, matrix: Float64Array.of(1, 0, 0, 1),
      newtonDirection: currentResidual.map(v => -v), radius: .1, maxTrials: 5,
      residual: x => { if (outcome === 'residual failure') throw new Error('failed after correction'); return residual(x); },
      admissible: x => x[0] > x[1] ** 2,
      linearizedConstraints: () => [{ value: 1e-8, lower: -.9e-8, gradient: new Map([[0, 1]]) }],
      constraintValues: x => { assert.equal(phase, 'base'); corrections++; return [x[0] - x[1] ** 2]; },
      activeSet: { snapshot: () => phase, restore: s => { phase = s; },
        prepare: () => { assert.equal(phase, 'base'); phase = 'candidate'; preparations++;
          return { changed: outcome === 'active event', changes: ['test event'] }; } } });
    assert.equal(r.accepted, outcome !== 'residual failure'); assert.equal(preparations, 5);
    assert.ok(corrections >= 2); assert.deepEqual(initial, before);
    if (outcome === 'residual failure') { assert.equal(phase, 'base'); assert.match(r.trials.at(-1).reason, /failed after correction/); }
    else {
      assert.equal(phase, 'candidate'); assert.equal(r.kind, 'projected-gradient-soc');
      assert.ok(r.x[0] > r.x[1] ** 2); assert.ok(r.scaledStepNorm <= .1 * (1 + 1e-14));
      if (outcome === 'active event') { assert.equal(r.actualReduction, null); assert.equal(r.ratio, null); }
      else {
        close(r.ratio, 1); assert.ok(r.actualReduction > 0);
        const expected = .5 * currentResidual.reduce((s, v, i) => s + v ** 2 - (1 + r.x[i]) ** 2, 0);
        close(r.predictedReduction, expected);
      }
    }
  }
});

test('a curvature correction that destroys descent is rejected and a smaller trust region is tried', () => {
  // At x > y² the projected tangent descends, but the O(y²) normal
  // correction opposes the large x gradient. The large step loses descent;
  // the smaller one remains feasible and reduces this independent quadratic.
  const initial = Float64Array.of(1e-8, 0), residual = x => Float64Array.of(100 + x[0], 1 + x[1]);
  const currentResidual = residual(initial);
  const result = takeDoglegStep({ initial, currentResidual, matrix: Float64Array.of(1, 0, 0, 1),
    newtonDirection: currentResidual.map(v => -v), radius: 2, maxTrials: 12, residual,
    admissible: x => x[0] > x[1] ** 2,
    linearizedConstraints: () => [{ value: initial[0], lower: -.9 * initial[0], gradient: new Map([[0, 1]]) }],
    constraintValues: x => [x[0] - x[1] ** 2] });
  assert.ok(result.trials.some(t => t.kind === 'projected-gradient-soc' && /No resolved dogleg model decrease/.test(t.reason)));
  assert.equal(result.accepted, true, result.reason); assert.equal(result.kind, 'projected-gradient-soc');
  assert.ok(result.trialRadius < 2); assert.ok(result.x[0] > result.x[1] ** 2);
  assert.ok(result.actualReduction > 0); close(result.ratio, 1, 1e-8);
});
