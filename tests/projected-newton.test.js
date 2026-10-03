import test from 'node:test';
import assert from 'node:assert/strict';
import { createDoglegModel, takeDoglegStep } from '../src/numerics/dogleg.js';

const close = (a, b, tolerance = 1e-12) => assert.ok(Math.abs(a - b) <= tolerance, `${a} != ${b}`);

test('constrained dogleg extends Cauchy to the first linear face with independently solved segment geometry', () => {
  // J=[[1,1],[0,1]], F=(1,2), Newton=(1,-2), y>=-3/2.
  // The scaled Cauchy point is (-11/17,-33/34). Its segment to Newton
  // meets the face at fraction 18/35, giving (1/5,-3/2), Fnew=(-3/10,1/2).
  // This predicts more reduction than either projected Cauchy or Newton.
  for (const units of [[1, 1], [1e3, .02]]) {
    const matrix = Float64Array.of(1 / units[0], 1 / units[1], 0, 1 / units[1]);
    const r = Float64Array.of(1, 2), n = Float64Array.of(units[0], -2 * units[1]);
    const constraints = [{ gradient: new Map([[1, 1 / units[1]]]), lower: -1.5 }];
    const model = createDoglegModel(matrix, r, n);
    const cauchy = model.proposeProjectedGradient(model.cauchyNorm, constraints);
    close(cauchy.direction[0] / units[0], -11 / 17); close(cauchy.direction[1] / units[1], -33 / 34);
    const p = model.proposeConstrainedDogleg(3, constraints, cauchy);
    close(p.doglegSegment.maximumFraction, 18 / 35); close(p.doglegSegment.fraction, 18 / 35);
    close(p.direction[0] / units[0], .2); close(p.direction[1] / units[1], -1.5);
    close(p.predictedReduction, 2.33); assert.ok(p.scaledStepNorm <= 3);
    assert.ok(p.predictedReduction > cauchy.predictedReduction);
    assert.ok(p.predictedReduction > model.proposeProjectedNewton(3, constraints).predictedReduction);
    const step = takeDoglegStep({ initial: Float64Array.of(0, 0), currentResidual: r, matrix, newtonDirection: n,
      residual: ([x, y]) => [1 + x / units[0] + y / units[1], 2 + y / units[1]], radius: 3,
      admissible: ([, y]) => y / units[1] >= -1.5 - 1e-14, linearizedConstraints: () => constraints });
    assert.equal(step.accepted, true); assert.equal(step.kind, 'constrained-dogleg'); close(step.ratio, 1);
    close(step.actualReduction, 2.33);
  }
});

test('a weak projected Newton step cannot preempt a much stronger feasible Cauchy model decrease', () => {
  // J=[[1,1],[0,1]], root=(-1,-epsilon) lies outside x>=0. Projecting
  // Newton preserves only its tiny y motion, despite an exactly solvable
  // restricted quadratic: x=0, y=-(1+2*epsilon)/2. This tests merit progress,
  // not convergence of an infeasible root. Changing units preserves it.
  const epsilon = 1e-5;
  for (const units of [[1, 1], [1e3, .02]]) {
    const matrix = Float64Array.of(1 / units[0], 1 / units[1], 0, 1 / units[1]);
    const initial = Float64Array.of(0, 0), residual = ([x, y]) =>
      Float64Array.of(1 + epsilon + x / units[0] + y / units[1], epsilon + y / units[1]);
    const r = residual(initial), direction = Float64Array.of(-units[0], -epsilon * units[1]);
    const constraints = [{ gradient: new Map([[0, 1 / units[0]]]), lower: 0 }];
    const model = createDoglegModel(matrix, r, direction), newton = model.proposeProjectedNewton(1, constraints);
    assert.ok(newton.predictedReduction > 0 && newton.predictedReduction < 2e-5);
    let evaluations = 0;
    const step = takeDoglegStep({ initial, matrix, currentResidual: r, newtonDirection: direction, radius: 1,
      residual: x => { evaluations++; return residual(x); }, admissible: ([x]) => x >= 0,
      linearizedConstraints: () => constraints });
    assert.equal(step.accepted, true); assert.equal(step.kind, 'projected-gradient');
    assert.equal(evaluations, 1, 'model ranking must happen before an expensive residual evaluation');
    assert.ok(step.actualReduction > 1e4 * newton.predictedReduction);
    assert.ok(step.scaledStepNorm <= 1); close(step.ratio, 1);
    // The projected Cauchy segment need not reach the restricted quadratic
    // minimum; independently evaluate its full residual and require descent.
    const next = residual(step.x);
    close(step.actualReduction, .5 * (r[0] ** 2 + r[1] ** 2 - next[0] ** 2 - next[1] ** 2));
    assert.deepEqual(initial, Float64Array.of(0, 0));
  }
});

test('projected Newton preserves a useful coupled direction and scores its complete model in any variable units', () => {
  // J=[[1,3],[0,1]], F=(4,1), exact Newton=(-1,-1), D=(1,sqrt(10)).
  // At x>=0, projecting the full Newton step gives (0,-1), Fnew=(1,0).
  // The projected gradient at the same radius makes less model progress.
  for (const units of [[1, 1], [1e3, .02]]) {
    const matrix = Float64Array.of(1 / units[0], 3 / units[1], 0, 1 / units[1]);
    const model = createDoglegModel(matrix, Float64Array.of(4, 1), Float64Array.of(-units[0], -units[1]));
    const radius = Math.sqrt(11), constraints = [{ gradient: new Map([[0, 1 / units[0]]]), lower: 0 }];
    const p = model.proposeProjectedNewton(radius, constraints), g = model.proposeProjectedGradient(radius, constraints);
    assert.equal(p.kind, 'projected-newton'); assert.equal(p.projection.converged, true);
    close(p.direction[0], 0); close(p.direction[1] / units[1], -1);
    close(p.predictedReduction, 8); close(p.linearizedResidualNorm, 1);
    assert.ok(p.predictedReduction > g.predictedReduction); assert.ok(p.scaledStepNorm <= radius);
    const small = model.proposeProjectedNewton(.1, constraints);
    close(small.direction[0], 0); close(small.direction[1] / units[1], -.1 / Math.sqrt(11));
    const y = small.direction[1] / units[1], expected = .5 * (17 - (4 + 3 * y) ** 2 - (1 + y) ** 2);
    close(small.predictedReduction, expected);
  }
});

test('a projection that destroys Newton descent is rejected while a projected Cauchy direction remains available', () => {
  // Root (1,-1), feasible y>=0. Projected Newton=(1,0) increases F=(1,1).
  const model = createDoglegModel(Float64Array.of(1, 2, 0, 1), Float64Array.of(1, 1), Float64Array.of(1, -1));
  const constraints = [{ gradient: new Map([[1, 1]]), lower: 0 }];
  const p = model.proposeProjectedNewton(3, constraints);
  assert.equal(p.direction, null); assert.match(p.reason, /No projected model-descent/);
  assert.ok(model.proposeProjectedGradient(3, constraints).predictedReduction > 0);
});

test('an inactive constraint leaves the feasible full Newton model unchanged', () => {
  const model = createDoglegModel(Float64Array.of(1, 3, 0, 1), Float64Array.of(4, 1), Float64Array.of(-1, -1));
  const plain = model.propose(4), projected = model.proposeProjectedNewton(4, [{ gradient: new Map([[0, 1]]), lower: -100 }]);
  assert.deepEqual(projected.direction, plain.direction); close(projected.predictedReduction, plain.predictedReduction);
});

test('the controller accepts projected Newton only after checking the actual residual and admissibility', () => {
  const residual = x => Float64Array.of(4 + x[0] + 3 * x[1], 1 + x[1]);
  const initial = Float64Array.of(0, 0);
  const result = takeDoglegStep({ initial, currentResidual: residual(initial),
    matrix: Float64Array.of(1, 3, 0, 1), newtonDirection: Float64Array.of(-1, -1), radius: 4,
    residual, admissible: x => x[0] >= 0,
    linearizedConstraints: () => [{ gradient: new Map([[0, 1]]), lower: 0 }] });
  assert.equal(result.accepted, true, result.reason); assert.equal(result.kind, 'projected-newton');
  assert.ok(result.trials.some(t => /Inadmissible/.test(t.reason)));
  close(result.x[0], 0); close(result.x[1], -1);
  close(result.actualReduction, 8); close(result.ratio, 1);
  assert.deepEqual(initial, Float64Array.of(0, 0));
});

test('the controller falls back when Newton projection destroys descent and restores each rejected active preparation', () => {
  const residual = x => Float64Array.of(1 + x[0] + 2 * x[1], 1 + x[1]);
  const initial = Float64Array.of(0, 0); let phase = 'base';
  const result = takeDoglegStep({ initial, currentResidual: residual(initial),
    matrix: Float64Array.of(1, 2, 0, 1), newtonDirection: Float64Array.of(1, -1), radius: 3,
    residual, admissible: x => x[1] >= 0,
    linearizedConstraints: () => [{ gradient: new Map([[1, 1]]), lower: 0 }],
    activeSet: { snapshot: () => phase, restore: s => { phase = s; },
      prepare: () => { assert.equal(phase, 'base'); phase = 'candidate'; return { changed: false }; } } });
  assert.equal(result.accepted, true, result.reason); assert.equal(result.kind, 'projected-gradient');
  assert.ok(result.trials.some(t => t.kind === 'projected-newton' && /No projected model-descent/.test(t.reason)));
  // D=(1,sqrt(5)); ||D^-1 J^T F||²=14/5 and
  // ||J D^-2 J^T F||²=26/5, hence the Cauchy multiplier is 7/13.
  close(result.x[0], -7 / 13); close(result.x[1], 0);
  close(result.actualReduction, .5 * (1 - (6 / 13) ** 2));
  assert.equal(phase, 'candidate');
});

test('a feasible projected Newton step with poor nonlinear reduction is rejected before shrinking the region', () => {
  const residual = x => Float64Array.of(4 + x[0] + 3 * x[1] + 10 * x[1] ** 2, 1 + x[1]);
  const initial = Float64Array.of(0, 0);
  const result = takeDoglegStep({ initial, currentResidual: residual(initial),
    matrix: Float64Array.of(1, 3, 0, 1), newtonDirection: Float64Array.of(-1, -1), radius: 4,
    residual, admissible: x => x[0] >= 0,
    linearizedConstraints: () => [{ gradient: new Map([[0, 1]]), lower: 0 }] });
  assert.ok(result.trials.some(t => t.kind === 'projected-newton' && t.ratio < 0));
  assert.equal(result.accepted, true, result.reason); assert.ok(result.trialRadius < 4);
  assert.ok(result.actualReduction > 0); assert.ok(result.ratio > 1e-4); assert.ok(result.x[0] >= 0);
});
