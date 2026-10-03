import test from 'node:test';
import assert from 'node:assert/strict';
import { createDoglegModel, takeDoglegStep } from '../src/numerics/dogleg.js';
import { solveLinear, normInf } from '../src/numerics/linear.js';
import { sparseMatrix, sparseAdd } from '../src/numerics/sparse.js';

const close = (a, b, tolerance = 2e-12) => assert.ok(Math.abs(a - b) <= tolerance, `${a} != ${b}`);
const matrix = Float64Array.of(1, 1, 0, 2), residual = Float64Array.of(1, 2), newton = Float64Array.of(0, -1);

test('scaled dogleg covers Cauchy, segment and Newton branches with independently checked model decrease', () => {
  const model = createDoglegModel(matrix, residual, newton);
  close(model.scales[0], 1); close(model.scales[1], Math.sqrt(5));
  // A=J/D, A^T*r=(1,sqrt(5)), |g|²/|Ag|²=6/8.
  close(model.cauchyNorm, .75 * Math.sqrt(6));
  for (const [radius, kind] of [[1, 'gradient'], [2, 'dogleg'], [3, 'newton']]) {
    const p = model.propose(radius), [x, y] = p.direction;
    assert.equal(p.kind, kind);
    const expectedNorm = Math.hypot(x, Math.sqrt(5) * y);
    close(expectedNorm, Math.min(radius, Math.sqrt(5)));
    const next = [1 + x + y, 2 + 2 * y], prediction = .5 * (5 - next[0] ** 2 - next[1] ** 2);
    close(p.predictedReduction, prediction); assert.ok(prediction > 0);
    if (kind === 'gradient') { close(x, -1 / Math.sqrt(6)); close(y, -1 / Math.sqrt(6)); }
  }
});

test('dense/sparse dogleg and independent state-variable units give the same physical step', () => {
  const sparse = sparseMatrix([[0, 1], [1]]); sparseAdd(sparse, 0, 0, 1); sparseAdd(sparse, 0, 1, 1); sparseAdd(sparse, 1, 1, 2);
  const units = [1e3, .02], changed = matrix.map((v, i) => v / units[i % 2]);
  for (const radius of [.1, 2, 3]) {
    const a = createDoglegModel(matrix, residual, newton).propose(radius);
    const b = createDoglegModel(sparse, residual, newton).propose(radius);
    const c = createDoglegModel(changed, residual, newton.map((v, i) => units[i] * v)).propose(radius);
    a.direction.forEach((v, i) => { close(v, b.direction[i]); close(v, c.direction[i] / units[i]); });
    close(a.predictedReduction, c.predictedReduction);
  }
  assert.throws(() => createDoglegModel(Float64Array.of(0, 0, 0, 1), residual, newton), /unresolved/);
});

test('trust-region acceptance predicts a linear residual exactly and leaves rejected states unchanged', () => {
  const initial = Float64Array.of(0, 0), before = initial.slice();
  const r = takeDoglegStep({ initial, currentResidual: residual, matrix, newtonDirection: newton,
    residual: ([x, y]) => [1 + x + y, 2 + 2 * y], radius: 2 });
  assert.equal(r.accepted, true); close(r.ratio, 1); close(r.radius, 4); assert.deepEqual(initial, before);
  const rejected = takeDoglegStep({ initial, currentResidual: residual, matrix, newtonDirection: newton,
    residual: () => residual, radius: 1, maxTrials: 3 });
  assert.equal(rejected.accepted, false); assert.equal(rejected.trials.length, 3); assert.deepEqual(initial, before);
  const forbidden = takeDoglegStep({ initial, currentResidual: residual, matrix, newtonDirection: newton,
    residual: () => { throw new Error('Must not evaluate forbidden state.'); }, admissible: () => false, maxTrials: 2 });
  assert.equal(forbidden.accepted, false); assert.ok(forbidden.trials.every(t => /Inadmissible/.test(t.reason)));
});

test('repeated dogleg updates solve the Rosenbrock residual with an explicit positive-domain constraint', () => {
  const residual = ([x, y]) => Float64Array.of(10 * (y - x * x), 1 - x);
  let x = Float64Array.of(-1.2, 1), radius = 1, count = 0; const kinds = new Set();
  for (; count < 60 && normInf(residual(x)) > 1e-10; count++) {
    const r = residual(x), matrix = Float64Array.of(-20 * x[0], 10, -1, 0), direction = solveLinear(matrix, r.map(v => -v));
    const step = takeDoglegStep({ initial: x, currentResidual: r, matrix, newtonDirection: direction,
      residual, radius, admissible: v => v[1] > 0 });
    assert.equal(step.accepted, true, step.reason); assert.ok(step.actualReduction > 0); kinds.add(step.kind);
    x = step.x; radius = step.radius;
  }
  assert.ok(count < 60); assert.ok(normInf(residual(x)) < 1e-10); close(x[0], 1, 1e-10); close(x[1], 1, 1e-10);
  assert.ok(kinds.has('newton') && (kinds.has('gradient') || kinds.has('dogleg')));
});

test('roundoff-sized trust changes stop as failure without imposing an absolute residual floor', () => {
  const matrix = Float64Array.of(1, 0, 0, 1), initial = Float64Array.of(0, 0);
  const stopped = takeDoglegStep({ initial, matrix, currentResidual: Float64Array.of(1, 1),
    newtonDirection: Float64Array.of(-1, -1), radius: Number.EPSILON,
    residual: () => { throw new Error('A roundoff-sized step must stop before evaluating the flow.'); } });
  assert.equal(stopped.accepted, false); assert.match(stopped.reason, /roundoff/); assert.equal(stopped.trials.length, 0);
  const tiny = Float64Array.of(1e-14, 2e-14);
  const resolved = takeDoglegStep({ initial, matrix, currentResidual: tiny, newtonDirection: tiny.map(v => -v),
    residual: x => x.map((v, i) => v + tiny[i]) });
  assert.equal(resolved.accepted, true); assert.deepEqual(resolved.residual, Float64Array.of(0, 0));
});

test('a feasible Newton ray escapes an outward Cauchy step without exceeding the trust radius', () => {
  // A^T A = [[1,2],[2,5]]. At the positive-x boundary, -A^T F points
  // toward forbidden x<0, while the exact Newton root (1,-1) is feasible.
  // Rescaling the two state units must not change that physical decision.
  for (const units of [[1, 1], [1e3, .02]]) {
    const a = Float64Array.of(1 / units[0], 2 / units[1], 0, 1 / units[1]);
    const initial = Float64Array.of(1e-12 * units[0], 0), radius = .1;
    const residual = ([x, y]) => Float64Array.of(x / units[0] + 2 * y / units[1] + 1, y / units[1] + 1);
    const r = residual(initial), direction = Float64Array.of((1 - 1e-12) * units[0], -units[1]);
    const model = createDoglegModel(a, r, direction), primary = model.propose(radius);
    assert.equal(primary.kind, 'gradient');
    assert.ok(initial[0] + primary.direction[0] < 0);
    const step = takeDoglegStep({ initial, matrix: a, currentResidual: r, newtonDirection: direction,
      residual, radius, admissible: x => x[0] > 0 });
    assert.equal(step.accepted, true, step.reason); assert.equal(step.kind, 'newton-ray');
    close(step.trialRadius, radius); close(step.scaledStepNorm, radius);
    close(step.ratio, 1); assert.ok(step.x[0] > initial[0]);
    const fraction = radius / Math.sqrt((1 - 1e-12) ** 2 + 5);
    close(step.x[0] / units[0], 1e-12 + fraction * (1 - 1e-12));
    close(step.x[1] / units[1], -fraction);
    assert.equal(step.trials.length, 2);
    assert.match(step.trials[0].reason, /Inadmissible/);
    assert.ok(step.actualReduction > 0 && step.predictedReduction > 0);
  }
});

test('inactive halfspaces leave ordinary nonlinear dogleg steps and trust-radius reductions unchanged', () => {
  const residual = ([x, y]) => Float64Array.of(10 * (y - x * x), 1 - x);
  for (const halfspaces of [[], [{ gradient: new Map([[0, 1]]), lower: -1e6 }, { gradient: new Map([[1, 1]]), lower: -1e6 }]]) {
    let x = Float64Array.of(-1.2, 1), radius = 1, updates = 0, contractions = 0;
    while (normInf(residual(x)) > 1e-10 && updates < 60) {
      const r = residual(x), matrix = Float64Array.of(-20 * x[0], 10, -1, 0), newtonDirection = solveLinear(matrix, r.map(v => -v));
      const input = { initial: x, currentResidual: r, matrix, newtonDirection, residual, radius };
      const plain = takeDoglegStep(input);
      const constrained = takeDoglegStep({ ...input, linearizedConstraints: () => halfspaces });
      assert.equal(plain.accepted, true, plain.reason); assert.equal(constrained.accepted, true, constrained.reason);
      assert.deepEqual(constrained.x, plain.x); assert.equal(constrained.kind, plain.kind);
      assert.equal(constrained.radius, plain.radius); assert.equal(constrained.trialRadius, plain.trialRadius);
      assert.deepEqual(constrained.trials, plain.trials);
      if (plain.trialRadius < radius) contractions++;
      x = plain.x; radius = plain.radius; updates++;
    }
    assert.ok(updates < 60); assert.ok(contractions > 0); assert.ok(normInf(residual(x)) < 1e-10);
  }
});

test('a poor but admissible Newton prediction contracts the radius before an unrelated active gradient projection', () => {
  // Newton preserves x=0 and is inside x>-.1, but overshoots the curved
  // residual in y. A radius-3 gradient projection would hit the x bound
  // and accept only 0.138 merit reduction. Ordinary contraction to .75
  // yields a feasible Newton ray with 1.190 reduction instead.
  const residual = ([x, y]) => Float64Array.of(1 + x + y + 2.3 * y * y, 2 + 2 * y);
  const initial = Float64Array.of(0, 0);
  const input = { initial, currentResidual: residual(initial), matrix: Float64Array.of(1, 1, 0, 2),
    newtonDirection: Float64Array.of(0, -1), residual, radius: 3, admissible: ([x]) => x > -.1 };
  const plain = takeDoglegStep(input);
  const constrained = takeDoglegStep({ ...input, linearizedConstraints: () => [{ gradient: new Map([[0, 1]]), lower: -.09 }] });
  assert.equal(constrained.accepted, true); assert.equal(constrained.kind, 'newton-ray');
  close(constrained.trialRadius, .75); assert.ok(constrained.actualReduction > 1.19);
  assert.deepEqual(constrained.x, plain.x); assert.deepEqual(constrained.trials, plain.trials);
  assert.ok(constrained.trials[0].ratio < 0); assert.equal(constrained.trials[0].reason, undefined);
  assert.deepEqual(initial, Float64Array.of(0, 0));
});
