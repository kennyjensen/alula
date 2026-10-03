import test from 'node:test';
import assert from 'node:assert/strict';
import { createPanelStreamtubeGrid } from '../src/euler/streamtube-body-initializer.js';
import { solveInviscid } from '../src/inviscid/linear-vortex.js';
import { streamfunctionAt } from '../src/inviscid/streamfunction.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

for (const boundaryCondition of ['normal-velocity', 'streamfunction']) test(`${boundaryCondition} initialization carries the appropriate physical body level`, () => {
  const input = intrinsicBodyFixture({ elements: 1, contourPanels: 80, bodySegments: 8, tubes: 3, alpha: 2 });
  const before = structuredClone(input);
  const panel = solveInviscid({ elements: input.bodies.map(b => ({ points: b.points })), alpha: input.alpha, boundaryCondition });
  const prepared = createPanelStreamtubeGrid(input, { panelBoundaryCondition: boundaryCondition, recordPotentialCoordinates: true });
  const profile = prepared.diagnostics.profiles[0], curve = createContourCurve(input.bodies[0].points);
  const sampled = streamfunctionAt(curve.evaluate(profile.stagnationParameter).point, panel.field);
  assert.equal(profile.sampledStagnationStreamfunction, sampled);
  const expected = boundaryCondition === 'streamfunction' ? panel.diagnostics.surfaceStreamfunctions[0] : sampled;
  assert.equal(profile.streamfunction, expected);
  assert.equal(profile.stagnationStreamfunctionDefect, sampled - expected);
  assert.equal(profile.streamfunctionSource, boundaryCondition === 'streamfunction' ? 'solved body Dirichlet constant' : 'sampled stagnation point');
  assert.equal(prepared.input.captureLevels[1], expected);
  const body = prepared.input.bodies[0];
  for (let i = body.leadingIndex; i <= body.trailingIndex; i++) {
    const offset = i - body.leadingIndex;
    for (const [side, g, j] of [['lower', 0, 3], ['upper', 1, 0]]) {
      assert.equal(prepared.potentialCoordinates[g][i][j].y, expected);
      const wall = curve.branch(side, prepared.system.fractions[0][side][offset], profile.stagnationParameter).point;
      assert.deepEqual(prepared.nodes[g][i][j], wall, 'The level correction must not move the prescribed C2 wall.');
    }
  }
  if (boundaryCondition === 'streamfunction') {
    assert.equal(profile.solvedSurfaceStreamfunction, expected);
    assert.ok(Math.abs(sampled - expected) > 1e-10, 'The case must distinguish sampled and solved levels.');
  } else assert.equal(profile.solvedSurfaceStreamfunction, null);
  assert.ok(prepared.diagnostics.maxStreamfunctionDrift < 3e-8);
  assert.deepEqual(input, before);
});
