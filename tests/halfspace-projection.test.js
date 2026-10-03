import test from 'node:test';
import assert from 'node:assert/strict';
import { projectHalfspaces } from '../src/numerics/halfspace-projection.js';
import { createDoglegModel, takeDoglegStep } from '../src/numerics/dogleg.js';

const close = (a, b, tolerance = 2e-10) => assert.ok(Math.abs(a - b) <= tolerance, `${a} != ${b}`);
const row = (a, lower = 0) => ({ gradient: new Map(a.map((v, j) => [j, v])), lower });

test('sparse halfspace projection matches independently solved faces and intersecting faces', () => {
  // Projection of (-2,-1) onto x>=0 is (0,-1), onto x>=0,y>=x is (0,0).
  // Projection onto x+y>=-1 is (-1,0). Redundant and distant faces do not alter it.
  for (const [constraints, expected] of [
    [[row([1, 0])], [0, -1]],
    [[row([1, 0]), row([-1, 1])], [0, 0]],
    [[row([1, 1], -1), row([2, 2], -2), row([1, 0], -100)], [-1, 0]],
  ]) {
    const point = Float64Array.of(-2, -1), before = point.slice();
    const r = projectHalfspaces(point, constraints);
    assert.equal(r.converged, true); assert.deepEqual(point, before);
    r.point.forEach((v, i) => close(v, expected[i]));
    assert.ok(r.primalViolation <= 3e-11); assert.ok(r.stationarityNorm <= 3e-11);
    assert.ok(r.projectedNorm <= Math.hypot(...point) + 1e-10);
  }
});

test('projection verifies KKT convergence and rejects invalid inputs without accepting an unfinished iterate', () => {
  const r = projectHalfspaces(Float64Array.of(-2, -1), [row([1, 0]), row([-1, 1])], { maxSweeps: 1 });
  assert.equal(r.converged, false); assert.ok(r.kktResidual > 0);
  for (const constraints of [[row([1, 0], 1)], [row([NaN, 1])], [{ gradient: new Map([[2, 1]]), lower: 0 }]])
    assert.throws(() => projectHalfspaces([1, 1], constraints));
  const empty = projectHalfspaces([.1, -.2], []);
  assert.equal(empty.converged, true); assert.deepEqual(Array.from(empty.point), [.1, -.2]);
});

test('small intersecting positive margins are resolved even below the global projection tolerance', () => {
  for (const units of [1e-6, 1, 1e6]) {
    const bound = -1e-12 * units;
    const r = projectHalfspaces([-2 * units, -units], [row([1, 0], bound), row([-1, 1], bound)]);
    assert.equal(r.converged, true);
    assert.ok(r.point[0] >= bound - 64 * Number.EPSILON * Math.abs(bound));
    assert.ok(r.point[1] - r.point[0] >= bound - 64 * Number.EPSILON * Math.abs(bound));
    assert.ok(r.maximumPrimalToleranceRatio <= 1);
    close(r.point[0] / units, -1e-12, 1e-25);
    close(r.point[1] / units, -2e-12, 1e-25);
  }
});

test('projected Cauchy model removes an outward component and predicts the full residual decrease in any units', () => {
  // Identity residual has root (-1,-1), outside x>=0. At x~0 the
  // unconstrained Newton and Cauchy directions both point out. The feasible
  // model improvement is y motion; projection cannot certify a root here.
  for (const units of [[1, 1], [1e3, .02]]) {
    const matrix = Float64Array.of(1 / units[0], 0, 0, 1 / units[1]);
    const residual = Float64Array.of(1, 1), newton = Float64Array.of(-units[0], -units[1]);
    const model = createDoglegModel(matrix, residual, newton);
    const constraints = [{ gradient: new Map([[0, 1 / units[0]]]), lower: -1e-12 }];
    const p = model.proposeProjectedGradient(1, constraints);
    assert.equal(p.kind, 'projected-gradient'); assert.equal(p.projection.converged, true);
    close(p.direction[0] / units[0], -1e-12);
    close(p.direction[1] / units[1], -1 / Math.sqrt(2));
    const next = p.direction.map((v, i) => 1 + v / units[i]);
    close(p.predictedReduction, 1 - .5 * (next[0] ** 2 + next[1] ** 2));
    assert.ok(p.predictedReduction > 0); assert.ok(p.scaledStepNorm <= 1);
  }
});

test('the trust controller scores projected steps against actual equations and restores rejected active metadata', () => {
  for (const curved of [false, true]) {
    const initial = Float64Array.of(1e-8, 0), before = initial.slice(); let phase = 'base', assemblies = 0;
    const residual = x => Float64Array.from(x, v => 1 + v), currentResidual = residual(initial);
    const r = takeDoglegStep({ initial, currentResidual, matrix: Float64Array.of(1, 0, 0, 1),
      newtonDirection: currentResidual.map(v => -v), radius: 1, maxTrials: 3, residual,
      admissible: x => x[0] > (curved ? x[1] ** 2 : 0),
      activeSet: { snapshot: () => phase, restore: saved => { phase = saved; }, prepare: () => { phase = 'trial'; return { changed: false }; } },
      linearizedConstraints: x => {
        assert.equal(phase, 'base'); assert.deepEqual(x, initial); assemblies++;
        return [{ gradient: new Map([[0, 1]]), lower: -.9 * x[0] }];
      } });
    assert.equal(assemblies, 1); assert.deepEqual(initial, before); assert.equal(r.trials.length, 3);
    assert.equal(r.accepted, !curved);
    if (curved) { assert.equal(phase, 'base'); assert.match(r.trials[2].reason, /Inadmissible/); }
    else {
      assert.equal(phase, 'trial'); assert.equal(r.kind, 'projected-newton');
      close(r.ratio, 1); assert.ok(r.actualReduction > 0); assert.ok(r.x[0] > 0);
    }
  }
});
