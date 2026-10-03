import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directChannelConservation } from './oracles/streamtube.js';

test('coupled trust region closes all two-element surface/wake equations with independent conservation', t => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  const baseline = solveCoupledStreamtubeBody(createCoupledStreamtubeBody(input), { tolerance: 1e-10 });
  assert.equal(baseline.converged, true, baseline.reason);
  const system = createCoupledStreamtubeBody(input), meshes = [];
  const r = solveCoupledStreamtubeBody(system, { stepMethod: 'dogleg', initialTrustRadius: .01,
    maxIterations: 30, tolerance: 1e-10, onMesh: mesh => meshes.push(mesh) });
  assert.equal(r.converged, true, r.reason); assert.equal(r.stepMethod, 'dogleg');
  assert.equal(r.boundaryLayer.surfaces.length, 4); assert.equal(r.boundaryLayer.wakes.length, 2);
  assert.equal(r.linearDiagnostics.solves, r.history.length - 1);
  assert.ok(r.linearDiagnostics.maxRelativeResidual <= 1e-10);
  assert.ok(r.history.some(h => ['gradient', 'dogleg', 'newton-ray'].includes(h.stepKind)));
  assert.equal(r.history.at(-1).stepKind, 'newton');
  assert.equal(meshes.length, r.history.length - 1);
  for (let i = 0; i < meshes.length; i++) {
    assert.equal(meshes[i].iteration.iteration, i + 1); assert.equal(meshes[i].initialization.flowSolved, false);
    assert.equal(meshes[i].quality.valid, true);
  }
  assert.deepEqual(meshes.at(-1).vertices, r.mesh.vertices);
  // Different partial-step paths rebase the normal chart along different
  // curves. They need not end at bitwise-identical tangential grid spacing.
  // When full Newton steps fit, both controllers must agree tightly.
  const full = solveCoupledStreamtubeBody(createCoupledStreamtubeBody(input), { stepMethod: 'dogleg', tolerance: 1e-10 });
  assert.equal(full.converged, true, full.reason);
  assert.ok(full.history.slice(1).every(h => h.stepKind === 'newton'));
  full.boundaryLayer.stations.forEach((station, id) => {
    for (const key of ['ue', 'theta', 'deltaStar', 'aux']) assert.ok(Math.abs(station[key] - baseline.boundaryLayer.stations[id][key]) < 2e-9, `${id} ${key}`);
  });
  for (let g = 0; g < r.flow.nodes.length; g++) {
    const c = directChannelConservation({ nodes: r.flow.nodes[g], sections: r.flow.sections.map(row => row[g]),
      cells: r.flow.cells.map(row => row[g]) }, system.euler.conditions.gamma);
    for (const key of ['maxLocal', 'total', 'internalCancellation']) assert.ok(c[key].every(v => Math.abs(v) < 2e-9));
  }
  t.diagnostic(JSON.stringify({ unknowns: r.x.length, iterations: r.history.length - 1, families: r.families,
    differentPartialPathGridDifference: Math.max(...r.mesh.vertices.map((p, i) => Math.hypot(p.x - baseline.mesh.vertices[i].x, p.y - baseline.mesh.vertices[i].y))),
    differentPartialPathEdgeDifference: Math.max(...r.boundaryLayer.stations.map((s, i) => Math.abs(s.ue - baseline.boundaryLayer.stations[i].ue))) }));
});

test('coupled trust-region controls and an iteration limit cannot certify unfinished flow', () => {
  const system = createCoupledStreamtubeBody(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }));
  for (const controls of [{ stepMethod: 'unknown' }, { initialTrustRadius: 0 }, { initialTrustRadius: Infinity }])
    assert.throws(() => solveCoupledStreamtubeBody(system, controls), /step method or trust radius/);
  const before = system.initial.slice();
  const r = solveCoupledStreamtubeBody(system, { stepMethod: 'dogleg', maxIterations: 1, initialTrustRadius: .001, tolerance: 1e-10 });
  assert.equal(r.converged, false); assert.equal(r.status, 'unconverged');
  assert.equal(r.mesh.initialization.flowSolved, false); assert.deepEqual(system.initial, before);
  assert.equal(r.history.length, 2); assert.ok(r.history.at(-1).actualReduction > 0);
});

test('a coupled Jacobian failure retains the last valid iterate and its explicit nonconvergence', () => {
  for (const stepMethod of ['newton', 'dogleg']) {
    const system = createCoupledStreamtubeBody(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }));
    const expected = system.initial.slice();
    system.jacobian = () => { throw new Error('local derivative left its valid domain'); };
    const r = solveCoupledStreamtubeBody(system, { stepMethod });
    assert.equal(r.converged, false); assert.match(r.reason, /local derivative left its valid domain/);
    assert.deepEqual(r.x, expected); assert.equal(r.history.length, 1);
    assert.equal(r.mesh.quality.valid, true); assert.equal(r.mesh.initialization.flowSolved, false);
  }
});

test('projected coupled step escapes the retained default corner limit with full nonlinear reduction and positive geometry', t => {
  const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/default-coupled-geometry-limit.json', import.meta.url)));
  const system = createCoupledStreamtubeBody(fixture.input, { ...fixture.options,
    initialEuler: fixture.initialEuler, initialBL: Float64Array.from(fixture.initialBL) });
  const before = system.initial.slice(), active = system.bl.snapshotActive(), meshes = [];
  const r = solveCoupledStreamtubeBody(system, { stepMethod: 'dogleg', initialTrustRadius: 1, maxIterations: 1,
    onMesh: mesh => meshes.push(mesh) });
  assert.equal(r.x.length, 7079); assert.equal(r.history.length, 2);
  assert.equal(r.projectedSteps, true); assert.equal(r.history[1].stepKind, 'projected-gradient');
  assert.equal(r.history[1].trialRadius, 1); assert.equal(r.history[1].projection.converged, true);
  // The projected Cauchy proposal predicts more decrease at this state and
  // is now tried first. Still require actual reduction of the full system,
  // positive geometry, unchanged active phases and every surface/wake.
  assert.ok(r.history[1].actualReduction > 80); assert.ok(Math.abs(r.history[1].reductionRatio - 1) < .01);
  assert.ok(r.history[1].residual < r.history[0].residual);
  assert.equal(r.converged, false); assert.equal(r.mesh.quality.valid, true);
  assert.ok(r.mesh.quality.minCornerSine > 0); assert.equal(meshes.length, 1);
  assert.equal(r.boundaryLayer.surfaces.length, 4); assert.equal(r.boundaryLayer.wakes.length, 2);
  assert.deepEqual(system.bl.snapshotActive(), active); assert.deepEqual(system.initial, before);
  assert.deepEqual(meshes[0].vertices, r.mesh.vertices);
  t.diagnostic(JSON.stringify({ unknowns: r.x.length, step: r.history[1], quality: r.mesh.quality }));
});

test('a positive but unresolved retained corner is rejected before another coupled solve', () => {
  const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/default-coupled-unresolved-corner.json', import.meta.url)));
  const system = createCoupledStreamtubeBody(fixture.input, { ...fixture.options,
    initialEuler: fixture.initialEuler, initialBL: Float64Array.from(fixture.initialBL) });
  assert.ok(system.residual(system.initial).every(Number.isFinite), 'finite equation evaluation alone cannot certify usable geometry');
  assert.equal(system.admissible(system.initial), false);
  assert.throws(() => solveCoupledStreamtubeBody(system, { stepMethod: 'dogleg', maxIterations: 1 }), /Invalid coupled streamtube solve controls or initial state/);
});

test('second-order coupled corrections cross the retained curvature limit with and without a material-trip event', t => {
  const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/default-coupled-curvature-limit.json', import.meta.url)));
  for (const radius of [.25, 1]) {
    // Velocity units retain the intended trip crossing at radius 1. Direct
    // pressure residuals change the trust-region metric and event radius.
    const system = createCoupledStreamtubeBody(fixture.input, { ...fixture.options, edgeMatching: 'section-velocity',
      initialEuler: fixture.initialEuler, initialBL: Float64Array.from(fixture.initialBL) });
    const before = system.initial.slice(), meshes = [];
    const r = solveCoupledStreamtubeBody(system, { stepMethod: 'dogleg', initialTrustRadius: radius,
      maxIterations: 1, onMesh: mesh => meshes.push(mesh) });
    assert.equal(r.x.length, 7079); assert.equal(r.history.length, 2); assert.equal(meshes.length, 1);
    const h = r.history[1]; assert.equal(h.stepKind, 'projected-gradient-soc');
    assert.equal(h.correction.corrected, true); assert.equal(h.trialRadius, radius);
    assert.ok(h.correction.correctionNorm < 1e-8 * radius);
    assert.equal(r.secondOrderSteps, true); assert.equal(r.converged, false); assert.equal(r.mesh.quality.valid, true);
    assert.ok(r.mesh.quality.minCornerSine > 1e-12); assert.deepEqual(meshes[0].vertices, r.mesh.vertices);
    assert.equal(r.boundaryLayer.surfaces.length, 4); assert.equal(r.boundaryLayer.wakes.length, 2);
    assert.deepEqual(system.initial, before);
    if (radius === .25) {
      assert.ok(h.actualReduction > 80); assert.ok(h.reductionRatio > .99); assert.equal(h.activeChange, undefined);
    } else {
      assert.equal(h.activeChange, true); assert.equal(h.actualReduction, null); assert.equal(h.reductionRatio, null);
      assert.ok(h.transitionChanges.some(c => c.body === 0 && c.side === 'lower' && c.from === 6 && c.to === 5));
    }
    t.diagnostic(JSON.stringify({ radius, step: h, quality: r.mesh.quality }));
  }
});
