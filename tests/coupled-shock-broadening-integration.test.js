// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { solveCoupledWithShockBroadening } from '../scripts/validation/coupled-shock-broadening-driver.js';

const copy = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const assertNumericallyEqual = (actual, expected, label = 'value') => {
  if (typeof actual === 'number' && typeof expected === 'number') {
    if (Number.isInteger(actual) && Number.isInteger(expected)) {
      assert.equal(actual, expected, label);
      return;
    }
    assert.ok(Math.abs(actual - expected) <= 1e-12 * Math.max(1, Math.abs(actual), Math.abs(expected)),
      `${label}: ${actual} differs from ${expected}`);
  } else if (Array.isArray(actual) && Array.isArray(expected)) {
    assert.equal(actual.length, expected.length, `${label} length`);
    actual.forEach((value, i) => assertNumericallyEqual(value, expected[i], `${label}[${i}]`));
  } else if (actual && expected && typeof actual === 'object' && typeof expected === 'object') {
    assert.deepEqual(Object.keys(actual).sort(), Object.keys(expected).sort(), `${label} keys`);
    for (const key of Object.keys(actual)) assertNumericallyEqual(actual[key], expected[key], `${label}.${key}`);
  } else assert.deepEqual(actual, expected, label);
};

test('constant-threshold multielement composition and serialized restart match uninterrupted ordinary Newton', t => {
  const cp = copy(JSON.parse(fs.readFileSync('docs/coupled-current-profile-preparation/two-element-six-update/initial.json')).checkpoint);
  const f = cp.restart;
  // A fixed .75 target makes the documented schedule exactly constant. This
  // isolates one-step composition/restart from any changed equation policy.
  f.input = { ...f.input, streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 },
    upwind: { mucon: 1, mcrit: .75, boundary: { kind: 'unfiltered-first-two' } } };
  f.options.blThermodynamics = 'historical-common-isentrope';
  const initial = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  const value = initial.evaluate(initial.initial);
  f.initialEuler = { x: copy(initial.initial.slice(0, initial.ne)), nodes: copy(value.outer.nodes),
    undisplacedNodes: copy(value.outer.undisplacedNodes) };
  cp.families = value.families;
  const before = copy(cp), c = cp.continuation, controls = { tolerance: 1e-10, maxIterations: 2 };
  const ordinary = solveCoupledStreamtubeIses(undefined, { resume: cp, ...controls,
    iterationGeometry: c.iterationGeometry, stepAcceptance: c.stepAcceptance,
    stagnationLimiter: c.stagnationLimiter, blUpdate: c.blUpdate });
  assert.equal(ordinary.history.length - 1, 2, ordinary.reason);
  const combined = solveCoupledWithShockBroadening(cp, controls);
  assert.equal(combined.operations.acceptedUpdates, 2, combined.reason);
  assert.equal(combined.linearDiagnostics.solves, ordinary.linearDiagnostics.solves);
  assertNumericallyEqual(combined.checkpoint, copy(ordinary.checkpoint));
  assertNumericallyEqual(copy(combined.result.residual), copy(ordinary.residual));
  assert.equal(combined.converged, ordinary.converged);
  assert.equal(combined.targetRestored, true);

  const first = solveCoupledWithShockBroadening(cp, { ...controls, maxIterations: 1 });
  const wrapper = copy(first.controller), preserved = copy(wrapper);
  const second = solveCoupledWithShockBroadening(wrapper, { ...controls, maxIterations: 1 });
  assertNumericallyEqual(second.checkpoint, combined.checkpoint);
  assert.equal(second.controller.acceptedUpdates, combined.controller.acceptedUpdates);
  assert.equal(second.controller.previousDensityChange, combined.controller.previousDensityChange);
  assert.deepEqual(wrapper, preserved); assert.deepEqual(cp, before);
  assert.equal(ordinary.boundaryLayer.surfaces.length, 4); assert.equal(ordinary.boundaryLayer.wakes.length, 2);
  t.diagnostic(JSON.stringify({ unknowns: initial.n, sourceConstructors: 1, sourceEvaluations: 1,
    nativeProfileMarches: 0, flowStartupSolves: 0, totalOrdinaryUpdates: 6,
    families: ordinary.families, packedAndMaintenanceReplayEquivalent: true,
    scope: 'Composition/restart at a constant threshold within floating-point roundoff; no adaptive convergence or physical validation claim.' }));
});
