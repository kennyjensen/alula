import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { createDoglegModel } from '../src/numerics/dogleg.js';
import { correctSecondOrderStep } from '../src/numerics/second-order-step.js';
import { solveSparseDirect } from '../src/numerics/klu.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';

test('coupled constraints include physical BL shape and enthalpy at every surface and wake station', () => {
  const system = createCoupledStreamtubeBody(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }));
  const x = system.initial, before = x.slice(), baseline = system.residual(x), constraints = system.stepConstraints(x);
  assert.deepEqual(system.constraintValues(x), constraints.map(c => c.value));
  assert.equal(constraints.filter(c => c.kind === 'kinematic-shape').length, system.bl.stations.length);
  assert.equal(constraints.filter(c => c.kind === 'static-enthalpy').length, system.bl.stations.length);
  assert.ok(constraints.every(c => c.value > 0));
  const station = system.bl.wakes[0].ids.at(-1), k = system.ne + 4 * station;
  for (const kind of ['kinematic-shape', 'static-enthalpy']) {
    const row = constraints.findIndex(c => c.kind === kind && c.station === station);
    for (const col of [k + 1, k + 2, k + 3]) {
      const h = 1e-4 * Math.abs(x[col]);
      const values = [-2, -1, 1, 2].map(m => { const trial = x.slice(); trial[col] += m * h; return system.constraintValues(trial)[row]; });
      const difference = (values[0] - 8 * values[1] + 8 * values[2] - values[3]) / (12 * h);
      assert.ok(Math.abs(difference - (constraints[row].gradient.get(col) ?? 0)) < 1e-8);
    }
    const bad = x.slice();
    if (kind === 'kinematic-shape') bad[k + 2] = 1.0001 * bad[k + 1];
    else bad[k + 3] = 20;
    assert.ok(system.constraintValues(bad)[row] < 0, `${kind} must reject the invalid state`);
  }
  assert.deepEqual(x, before); assert.deepEqual(system.residual(x), baseline);
});

test('a supplied coupled restart restores the valid displaced grid even when its auxiliary undisplaced grid is folded', () => {
  const seed = JSON.parse(fs.readFileSync(new URL('fixtures/default-coupled-displaced-restart.json', import.meta.url)));
  const system = createCoupledStreamtubeBody(seed.input, { ...seed.options,
    initialEuler: seed.initialEuler, initialBL: Float64Array.from(seed.initialBL) });
  const value = system.evaluate(system.initial), mesh = streamtubeMeshSnapshot({ system: system.euler, nodes: value.outer.nodes });
  assert.equal(mesh.quality.valid, true); assert.ok(system.admissible(system.initial));
  assert.deepEqual(Array.from(system.initial.subarray(system.ne)), seed.initialBL);
  for (const [key, expected] of Object.entries(seed.expectedFamilies)) assert.ok(Math.abs(value.families[key] - expected) < 1e-8);
  value.outer.nodes.forEach((g, k) => g.forEach((row, i) => row.forEach((p, j) => {
    const q = seed.initialEuler.nodes[k][i][j];
    assert.ok(Math.hypot(p.x - q.x, p.y - q.y) < 2e-13);
  })));
  system.euler.setDisplacement(system.bl.thicknesses(new Float64Array(seed.initialBL.length)));
  assert.throws(() => system.euler.evaluate(system.initial.subarray(0, system.ne)), /Folded or degenerate/);
  assert.deepEqual(system.residual(system.initial), value.residual);
});

test('retained default compressible bound admits useful Newton proposals and resolves the corner correction input', () => {
  const seed = JSON.parse(fs.readFileSync(new URL('fixtures/default-coupled-compressible-limit.json', import.meta.url)));
  const system = createCoupledStreamtubeBody(seed.input, { ...seed.options,
    initialEuler: seed.initialEuler, initialBL: Float64Array.from(seed.initialBL) });
  const x = system.initial, base = system.evaluate(x), active = system.bl.snapshotActive(), matrix = system.jacobian(x);
  const linear = solveSparseDirect(matrix, base.residual.map(v => -v), { tolerance: 1e-10 });
  const model = createDoglegModel(matrix, base.residual, linear.x), constraints = system.stepConstraints(x);
  for (const radius of [2, 8]) {
    const p = model.proposeProjectedNewton(radius, constraints);
    assert.ok(p.direction); assert.ok(p.predictedReduction > 0);
    assert.ok(p.scaledStepNorm <= radius * (1 + 64 * Number.EPSILON), `scaled norm ${p.scaledStepNorm} exceeds radius ${radius}`);
    const trial = x.map((v, i) => v + p.direction[i]);
    assert.ok(system.constraintValues(trial).every(v => v > 0)); assert.ok(system.admissible(trial));
    const next = system.residual(trial);
    const decrease = .5 * base.residual.reduce((sum, v, i) => sum + (v - next[i]) * (v + next[i]), 0);
    assert.ok(decrease / p.predictedReduction > .8);
  }
  const radius = 8, p = model.proposeProjectedGradient(radius, constraints);
  const corrected = correctSecondOrderStep({ initial: x, direction: p.direction, scales: model.scales, radius,
    constraints, values: system.constraintValues });
  assert.equal(corrected.corrected, true); assert.ok(corrected.correctionNorm < 1e-3);
  assert.ok(system.constraintValues(x.map((v, i) => v + corrected.direction[i])).every(v => v > 0));
  // This larger gradient step also crosses a material-trip interval. It is
  // a constraint check only; the event controller must still bound/transfer it.
  assert.deepEqual(system.bl.snapshotActive(), active); assert.deepEqual(system.residual(x), base.residual);
});
