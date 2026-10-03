import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';

test('projected Newton escapes the actual default BL shape bound with full residual decrease and an improved grid corner', t => {
  const seed = JSON.parse(fs.readFileSync(new URL('fixtures/default-coupled-shape-bound.json', import.meta.url)));
  const system = createCoupledStreamtubeBody(seed.input, { ...seed.options,
    initialEuler: seed.initialEuler, initialBL: Float64Array.from(seed.initialBL) });
  const initial = system.initial.slice(), active = system.bl.snapshotActive(), baseline = system.residual(initial), meshes = [];
  const result = solveCoupledStreamtubeBody(system, { stepMethod: 'dogleg', initialTrustRadius: 16, maxIterations: 1,
    onMesh: mesh => meshes.push(mesh) });
  assert.equal(result.x.length, 7079); assert.equal(result.history.length, 2); assert.equal(meshes.length, 1);
  assert.equal(result.boundaryLayer.surfaces.length, 4); assert.equal(result.boundaryLayer.wakes.length, 2);
  const h = result.history[1];
  assert.equal(h.stepKind, 'projected-newton'); assert.equal(h.trialRadius, 16);
  assert.equal(h.projection.converged, true); assert.equal(h.activeChange, undefined);
  assert.ok(h.actualReduction > 8.8); assert.ok(Math.abs(h.reductionRatio - 1) < .02);
  assert.ok(result.linearDiagnostics.maxRelativeResidual < 1e-10);
  const final = system.residual(result.x);
  const decrease = .5 * baseline.reduce((sum, v, i) => sum + (v - final[i]) * (v + final[i]), 0);
  assert.ok(Math.abs(decrease - h.actualReduction) < 1e-8);
  assert.ok(system.constraintValues(result.x).every(v => v > 0));
  assert.ok(system.admissible(result.x)); assert.equal(result.mesh.quality.valid, true);
  assert.ok(result.mesh.quality.minCornerSine > .009);
  assert.deepEqual(meshes[0].vertices, result.mesh.vertices);
  assert.equal(meshes[0].initialization.flowSolved, false); assert.equal(result.converged, false);
  assert.equal(result.reason, 'iteration limit');
  assert.deepEqual(system.bl.snapshotActive(), active); assert.deepEqual(system.initial, initial);
  t.diagnostic(JSON.stringify({ step: h, quality: result.mesh.quality, families: result.families }));
});

test('native BL roundoff cannot admit a state outside the polynomial domain used by the next coupled step', () => {
  const seed = JSON.parse(fs.readFileSync(new URL('fixtures/default-coupled-roundoff-shape.json', import.meta.url)));
  const system = createCoupledStreamtubeBody(seed.input, { ...seed.options,
    initialEuler: seed.initialEuler, initialBL: Float64Array.from(seed.initialBL) });
  const x = system.initial;
  // The original native expression rounds to raw Hk > 1, while the
  // equivalent polynomial is -5.10e-18 in scaled thickness units.
  assert.ok(system.residual(x).every(Number.isFinite));
  const failed = system.stepConstraints(x).filter(c => c.value <= 0);
  assert.equal(failed.length, 1); assert.equal(failed[0].kind, 'kinematic-shape'); assert.equal(failed[0].station, 51);
  assert.equal(system.admissible(x), false);
  assert.throws(() => solveCoupledStreamtubeBody(system, { stepMethod: 'dogleg', maxIterations: 1 }), /Invalid coupled streamtube solve controls or initial state/);
  // A represented positive margin remains allowed; this check introduces
  // no empirical Hk floor, tolerance relaxation or closure modification.
  const inside = x.slice(); inside[system.ne + 4 * 51 + 2] += 16 * Number.EPSILON * inside[system.ne + 4 * 51 + 2];
  assert.ok(system.constraintValues(inside).every(v => v > 0)); assert.ok(system.admissible(inside));
});
