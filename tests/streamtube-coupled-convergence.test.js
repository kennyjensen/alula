// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { coupledChangeDecision, coupledConvergenceSatisfied, smallCoupledChanges,
  acceptedCoupledChanges, coupledChangeSnapshot, newtonCoupledChanges, coupledResidualChanges } from '../src/euler/streamtube-coupled-convergence.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const small = () => ({ density: { count: 100, rms: 1e-5, maximum: 5e-5 } });
const families = { euler: 1e-6, boundaryLayer: 2e-5, edgeMatching: 1e-6 };
const evidence = (overrides = {}) => coupledChangeDecision({ accepted: small(), newton: small(), residual: small(),
  families, stable: true, ...overrides });

test('MSES-style stopping requires RMS and maximum limits in every family', () => {
  const decision = evidence(); assert.equal(decision.converged, true);
  assert.equal(coupledConvergenceSatisfied({ families, convergence: decision }, 1e-10), true);
  // One outlier cannot disappear into a large grid's RMS.
  assert.equal(evidence({ accepted: { grid: { count: 1e6, rms: 1e-5, maximum: .002 } } }).converged, false);
  assert.equal(evidence({ accepted: { grid: { count: 100, rms: 2e-4, maximum: 3e-4 } } }).converged, false);
});

test('damped stagnation, unresolved equations, events, and missing/nonfinite evidence fail shut', () => {
  assert.equal(evidence({ newton: { shear: { count: 10, rms: .01, maximum: .07 } } }).converged, false);
  assert.equal(evidence({ residual: { boundaryLayer: { count: 1000, rms: .01, maximum: .8 } } }).converged, false);
  assert.equal(evidence({ stable: false }).converged, false);
  for (const newton of [undefined, {}, { density: { count: 1, rms: NaN, maximum: Infinity } }])
    assert.equal(evidence({ newton }).converged, false);
  assert.equal(coupledConvergenceSatisfied({ families, convergence: evidence({ families: { ...families, euler: 0 } }) }, 1e-10), false);
  assert.equal(coupledConvergenceSatisfied({ families, convergence: { method: 'mses-changes', converged: true } }, 1e-10), false);
  assert.equal(coupledConvergenceSatisfied({ families: { euler: 1e-12, boundaryLayer: 2e-12, edgeMatching: 0 } }, 1e-10), true);
});

function fixture(length = 1, scale = .001) {
  const state = Float64Array.of(0, 123, .02, 2, 2.2, 1);
  const layers = [{ aux: .02, theta: 2 * scale, deltaStar: 2.2 * scale, wakeGap: .1 * scale, ue: 1 }];
  const system = { ne: 2, conditions: { edgeMatching: 'pressure' }, euler: {
    layout: { densityCount: 1, globalOffset: 2, positions: [{ column: 1 }], rows: [{ kind: 'streamwise' }, { kind: 'wakeGap' }] },
    conditions: { lengthScale: length, pressureScale: 50 },
  }, bl: { scale, stations: [{ id: 0, regime: 'wake' }], snapshotActive: () => [1],
    geometry: () => ({ coordinates: [{ wakeGapDerivatives: new Map([[1, scale]]) }] }) } };
  return { system, state, value: { outer: { nodes: [[[{ x: length, y: 2 * length }]]] }, layers: { states: layers } } };
}

test('maintained-grid changes are invariant to chart rebasing and coordinate units', () => {
  const run = length => {
    const { system, state, value } = fixture(length), before = coupledChangeSnapshot(system, state, value);
    const next = structuredClone(value); next.outer.nodes[0][0][0].y += length * 1e-5;
    const rebased = state.slice(); rebased[1] = -800; // Chart origin changed, physical grid barely moved.
    return acceptedCoupledChanges(system, before, rebased, next);
  };
  assert.ok(Math.abs(run(1).grid.maximum - run(100).grid.maximum) < 1e-15);
  assert.equal(smallCoupledChanges(run(100)), true);
});

test('Newton checks use physical BL scales, total wake-gap derivatives, and low-Mach pressure units', () => {
  for (const scale of [.001, .01]) {
    const { system, state, value } = fixture(1, scale);
    // Move numerical delta* and geometric base gap by equal amounts: fluid thickness stays fixed.
    const d = Float64Array.of(0, 1e-5, 0, 0, 1e-5, 0);
    const changes = newtonCoupledChanges(system, state, d, value);
    assert.equal(changes.displacement.maximum, 0); assert.equal(changes.shape.maximum, 0);
    const residual = coupledResidualChanges(system, Float64Array.of(3e-5, 0, 0, 0, 0, 0));
    assert.equal(residual.streamwise.maximum, .0015); assert.equal(smallCoupledChanges(residual), false);
  }
});

test('strict mode remains reproducible; resumed change checks require fresh evidence', () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  const controls = { edgeMatching: 'section-velocity', maxIterations: 12, tolerance: 1e-10, stepAcceptance: 'admissible' };
  const strict = solveCoupledStreamtubeIses(input, controls);
  const explicit = solveCoupledStreamtubeIses(input, { ...controls, convergence: 'residual' });
  assert.equal(strict.converged, true); assert.deepEqual(explicit.x, strict.x);
  const practical = solveCoupledStreamtubeIses(input, { ...controls, convergence: 'mses' });
  assert.equal(practical.converged, true, practical.reason);
  assert.equal(practical.reason, 'solution changes');
  assert.equal(practical.residualConverged, false);
  assert.ok(practical.history.length < strict.history.length);
  assert.ok(coupledConvergenceSatisfied(practical, 1e-10));
  assert.deepEqual(practical.checkpoint.convergence, practical.convergence);
  assert.equal(practical.checkpoint.continuation.convergence, 'mses');
  const replay = solveCoupledStreamtubeIses(undefined, { ...practical.checkpoint.continuation,
    resume: practical.checkpoint, maxIterations: 0, tolerance: 1e-10 });
  assert.equal(replay.converged, practical.residualConverged); // No stored change flag is trusted on resume.
  assert.equal(replay.convergence, undefined);
  assert.equal(practical.mesh.quality.valid, true);
  const polished = solveCoupledStreamtubeIses(undefined, { ...practical.checkpoint.continuation,
    resume: practical.checkpoint, convergence: 'residual', maxIterations: 4, tolerance: 1e-10 });
  assert.equal(polished.converged, true);
  for (let b = 0; b < practical.flow.diagnosticForces.length; b++)
    assert.ok(Math.abs(polished.flow.diagnosticForces[b].cl - practical.flow.diagnosticForces[b].cl) < 1e-6);
  assert.throws(() => solveCoupledStreamtubeIses(input, { convergence: 'bogus' }), /convergence policy/);
});
