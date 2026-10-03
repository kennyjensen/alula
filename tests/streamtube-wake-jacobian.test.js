import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { createDoglegModel } from '../src/numerics/dogleg.js';
import { solveSparseDirect } from '../src/numerics/klu.js';
import { sparseProduct } from '../src/numerics/sparse.js';

test('the actual default wake-branch Jacobian predicts whole-system directional derivatives and quadratic Newton defects', t => {
  const seed = JSON.parse(fs.readFileSync(new URL('fixtures/default-coupled-displaced-restart.json', import.meta.url)));
  const system = createCoupledStreamtubeBody(seed.input, { ...seed.options, initialEuler: seed.initialEuler, initialBL: Float64Array.from(seed.initialBL) });
  const x = system.initial, residual = system.residual(x), active = system.bl.snapshotActive(), matrix = system.jacobian(x);
  const linear = solveSparseDirect(matrix, residual.map(v => -v), { tolerance: 1e-10 });
  const model = createDoglegModel(matrix, residual, linear.x), newton = model.proposeNewtonRay(2);
  const directions = [newton, model.proposeProjectedGradient(2, system.stepConstraints(x))], diagnostics = [];
  for (const p of directions) {
    assert.ok(p.direction); const jd = sparseProduct(matrix, p.direction), errors = [];
    for (const h of [1e-3, 3e-4, 1e-4]) {
      const samples = [-2, -1, 1, 2].map(m => system.residual(x.map((v, i) => v + m * h * p.direction[i])));
      let maximum = 0;
      for (let i = 0; i < system.n; i++) {
        const fd = (samples[0][i] - 8 * samples[1][i] + 8 * samples[2][i] - samples[3][i]) / (12 * h);
        maximum = Math.max(maximum, Math.abs(fd - jd[i]) / Math.max(1, Math.abs(fd), Math.abs(jd[i])));
      }
      assert.ok(maximum < 1e-7, `${p.kind}, h=${h}: ${maximum}`); errors.push(maximum);
    }
    diagnostics.push({ kind: p.kind, errors });
  }
  const jd = sparseProduct(matrix, newton.direction), defects = [];
  for (const fraction of [.01, .001]) {
    const f = system.residual(x.map((v, i) => v + fraction * newton.direction[i]));
    const predicted = -fraction * residual.reduce((sum, v, i) => sum + v * jd[i], 0)
      - .5 * fraction ** 2 * jd.reduce((sum, v) => sum + v * v, 0);
    const actual = .5 * residual.reduce((sum, v, i) => sum + (v - f[i]) * (v + f[i]), 0);
    assert.ok(Math.abs(actual / predicted - 1) < 1e-4);
    defects.push(Math.max(...f.map((v, i) => Math.abs(v - residual[i] - fraction * jd[i]))));
  }
  assert.ok(Math.abs(defects[0] / defects[1] - 100) < 1);
  assert.deepEqual(system.bl.snapshotActive(), active); assert.deepEqual(system.residual(x), residual);
  t.diagnostic(JSON.stringify({ unknowns: system.n, linearRelativeResidual: linear.relativeResidual, directions: diagnostics, defects }));
});

test('the retained active-Hk state assembles its full Jacobian without invalid numerical probes', () => {
  const seed = JSON.parse(fs.readFileSync(new URL('fixtures/default-coupled-active-limit.json', import.meta.url)));
  const system = createCoupledStreamtubeBody(seed.input, { ...seed.options, initialEuler: seed.initialEuler, initialBL: Float64Array.from(seed.initialBL) });
  const x = system.initial, baseline = system.residual(x), matrix = system.jacobian(x), active = system.bl.snapshotActive();
  assert.equal(system.n, 7079); assert.ok(matrix.values.every(Number.isFinite)); assert.ok(system.admissible(x));
  const station = 51, k = system.ne + 4 * station;
  for (const [column, sign] of [[k + 1, -1], [k + 2, 1], [k + 3, -1]]) {
    const h = 1e-4 * x[column], direction = new Float64Array(system.n); direction[column] = 1;
    const exact = sparseProduct(matrix, direction);
    const samples = [0, 1, 2, 3, 4].map(m => { const trial = x.slice(); trial[column] += m * sign * h; return system.residual(trial); });
    for (let row = 0; row < system.n; row++) {
      const fd = (-25 * samples[0][row] + 48 * samples[1][row] - 36 * samples[2][row] + 16 * samples[3][row] - 3 * samples[4][row]) / (12 * sign * h);
      assert.ok(Math.abs(fd - exact[row]) < 2e-6 * Math.max(1, Math.abs(fd), Math.abs(exact[row])), `${row}/${column}: ${fd} != ${exact[row]}`);
    }
  }
  assert.deepEqual(system.bl.snapshotActive(), active); assert.deepEqual(system.residual(x), baseline);
});
