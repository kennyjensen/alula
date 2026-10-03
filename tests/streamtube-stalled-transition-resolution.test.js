import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { stalledBoundaryLayerResolutionPlan, recoverCoupledTransition } from '../src/euler/streamtube-transition-recovery.js';

test('a dominant near-transition stall earns one bounded refinement with the original acceptance policy', () => {
  const fixture = JSON.parse(gunzipSync(fs.readFileSync(new URL('./fixtures/rae16x9-transition-stall.json.gz', import.meta.url))));
  const cp = fixture.checkpoint, original = structuredClone(cp);
  const source = solveCoupledStreamtubeIses(undefined, { ...cp.continuation, resume: cp, maxIterations: 0, tolerance: 1e-10 });
  source.history = fixture.history;
  const plan = stalledBoundaryLayerResolutionPlan(source);
  assert.ok(plan);
  assert.equal(plan.normalFactor, 1);
  assert.ok(plan.streamwiseSubdivisions.every(n => n === 1 || n === 2));
  assert.ok(plan.refinedNx <= plan.parentNx + 16);
  assert.equal(stalledBoundaryLayerResolutionPlan({ ...source, converged: true }), null);
  assert.equal(stalledBoundaryLayerResolutionPlan({ ...source, automaticRefinement: {} }), null);
  assert.equal(stalledBoundaryLayerResolutionPlan({ ...source, families: { ...source.families, euler: 1 } }), null);
  assert.equal(stalledBoundaryLayerResolutionPlan({ ...source, history: source.history.map(h => ({ ...h, step: 1 })) }), null);
  assert.equal(stalledBoundaryLayerResolutionPlan(source, { maxNodes: 10 }), null);
  const result = recoverCoupledTransition(source, { plan, maxIterations: 40, tolerance: 1e-10 });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.stepAcceptance, cp.continuation.stepAcceptance);
  assert.equal(result.mesh.quality.valid, true);
  assert.equal(result.mesh.cells.length, 1944);
  assert.equal(result.solverInput.mach, cp.restart.input.mach);
  assert.deepEqual(result.solverInput.upwind, cp.restart.input.upwind);
  for (const key of ['reynolds', 'ncrit', 'transitionMode'])
    assert.equal(result.conditions[key], source.conditions[key]);
  for (const value of Object.values(result.families)) assert.ok(value <= 1e-10);
  assert.deepEqual(cp, original);
});

test('a laminar surface stall away from transition converges after bounded local refinement', () => {
  const fixture = JSON.parse(gunzipSync(fs.readFileSync(new URL('./fixtures/nlr-laminar-resolution-stall.json.gz', import.meta.url))));
  const cp = fixture.checkpoint, original = structuredClone(cp);
  const source = solveCoupledStreamtubeIses(undefined, { ...cp.continuation, resume: cp, maxIterations: 0, tolerance: 1e-10 });
  source.history = fixture.history;
  const plan = stalledBoundaryLayerResolutionPlan(source);
  assert.equal(plan.reason, 'stalled-boundary-layer-resolution');
  const surface = source.boundaryLayer.surfaces.find(s => s.ids.includes(plan.dominantStation));
  assert.ok(Math.abs(surface.ids.indexOf(plan.dominantStation) - surface.transition) > 2);
  assert.equal(source.boundaryLayer.stations[plan.dominantStation].regime, 'laminar');
  assert.equal(plan.refinedNx - plan.parentNx, 9);
  assert.equal(plan.normalFactor, 1);
  assert.equal(stalledBoundaryLayerResolutionPlan(source, { maxNodes: 10 }), null);
  const result = recoverCoupledTransition(source, { plan, maxIterations: 40, tolerance: 1e-10 });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.mesh.quality.valid, true);
  assert.equal(result.stepAcceptance, cp.continuation.stepAcceptance);
  for (const key of ['mach', 'alpha']) assert.equal(result.solverInput[key], source.solverInput[key]);
  for (const key of ['reynolds', 'ncrit', 'transitionMode']) assert.equal(result.conditions[key], source.conditions[key]);
  assert.deepEqual(result.solverInput.upwind, source.solverInput.upwind);
  for (const value of Object.values(result.families)) assert.ok(value <= 1e-10);
  assert.equal(stalledBoundaryLayerResolutionPlan(result), null);
  assert.deepEqual(cp, original);
});
