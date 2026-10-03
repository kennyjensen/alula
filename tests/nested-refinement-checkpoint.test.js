// SPDX-License-Identifier: GPL-2.0-or-later
// A continuation-transfer regression, not a converged-refinement claim.
import test from 'node:test';
import assert from 'node:assert/strict';
import { nestedRefinementCheckpoint } from '../scripts/validation/nested-refinement-checkpoint.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { refineCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-refinement.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';

const serialize = x => JSON.parse(JSON.stringify(x, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const create = f => createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
const cumulative = counts => counts.reduce((r, n) => [...r, r.at(-1) + n], [0]);
const near = (a, b, tolerance = 2e-12) => assert.ok(Math.abs(a - b) <= tolerance, `${a} != ${b}`);

test('coordinated nested refinement maps stored inlet history and replays all surfaces and wakes without initial SMOVE', t => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    streamwiseMode: 'isentropic', wakeGeometry: 'independent-banks', wakeOutlet: 'banks' };
  const controls = { reynolds: 1e6, edgeMatching: 'section-velocity', transitionMode: 'automatic',
    maxIterations: 12, tolerance: 1e-10, stepAcceptance: 'admissible' };
  const root = solveCoupledStreamtubeIses(input, controls);
  assert.equal(root.converged, true, root.reason);
  const checkpoint = serialize(root.checkpoint), originalRoot = serialize(root.checkpoint);
  // A stored target need not equal distances recaptured from the solved grid.
  // Exercise that distinction while keeping endpoints and monotonicity.
  checkpoint.continuation.fractions = checkpoint.continuation.fractions.map(row => row.map(f => f + .01 * f * (1 - f)));
  checkpoint.continuation.preferredOrdering = 'colamd'; checkpoint.continuation.pivotTolerance = 1;
  const frozen = structuredClone(checkpoint), parent = create(checkpoint.restart), before = parent.evaluate(parent.initial);
  const old = parent.euler.layout, streamwiseSubdivisions = Array(old.nx).fill(1);
  const normalSubdivisions = old.tubes.map(n => Array(n).fill(1));
  old.bodies.forEach((body, b) => {
    for (let i = Math.max(0, body.leadingIndex - 2); i < Math.min(body.trailingIndex, body.leadingIndex + 4); i++) streamwiseSubdivisions[i] = 2;
    normalSubdivisions[b][old.tubes[b] - 1] = 4; normalSubdivisions[b + 1][0] = 4;
  });
  const config = { streamwiseSubdivisions, normalSubdivisions }, mapped = refineCoupledStreamtubeBody(checkpoint.restart.input, parent, config);
  const restart = serialize({ input: mapped.input, options: mapped.options, initialEuler: mapped.initialEuler, initialBL: mapped.initialBL });
  const frozenRestart = structuredClone(restart), frozenConfig = structuredClone(config);
  const result = nestedRefinementCheckpoint(checkpoint, restart, config), child = create(result.checkpoint.restart), after = child.evaluate(child.initial);
  const retained = cumulative(streamwiseSubdivisions), normal = normalSubdivisions.map(cumulative);
  assert.deepEqual(result.diagnostics.retainedStreamwiseStations, retained);
  assert.deepEqual(result.diagnostics.retainedNormalStations, normal);
  assert.deepEqual(child.conditions, parent.conditions);
  assert.deepEqual(child.euler.conditions, parent.euler.conditions);
  assert.deepEqual(after.outer.captured, before.outer.captured);
  const capture = old.globals.capture.find(col => col !== null);
  assert.ok(Math.abs(parent.initial[capture]) > 1e-6, 'Exercise physical masses at nonzero capture.');
  for (let b = 0; b < old.elements; b++) {
    assert.equal(child.euler.layout.bodies[b].leadingIndex, retained[old.bodies[b].leadingIndex]);
    assert.equal(child.euler.layout.bodies[b].trailingIndex, retained[old.bodies[b].trailingIndex]);
    const previous = checkpoint.continuation.fractions[b], next = result.checkpoint.continuation.fractions[b];
    for (let i = 0; i < previous.length; i++) assert.equal(next[retained[i]], previous[i]);
    for (let i = 0; i < previous.length - 1; i++) for (let k = 1; k < streamwiseSubdivisions[i]; k++)
      assert.equal(next[retained[i] + k], previous[i] + (previous[i + 1] - previous[i]) * k / streamwiseSubdivisions[i]);
    for (let i = 0; i < previous.length - 1; i++) for (let k = 1; k < streamwiseSubdivisions[i]; k++) {
      const p = before.outer.nodes[b][i].at(-1), q = before.outer.nodes[b][i + 1].at(-1);
      const point = after.outer.nodes[b][retained[i] + k].at(-1), t = k / streamwiseSubdivisions[i];
      near(point.x, (1 - t) * p.x + t * q.x); near(point.y, (1 - t) * p.y + t * q.y);
    }
  }
  const { fractions: _old, ...history } = checkpoint.continuation;
  const { fractions: _new, ...newHistory } = result.checkpoint.continuation;
  assert.deepEqual(newHistory, history);
  before.outer.nodes.forEach((grid, g) => grid.forEach((row, i) => row.forEach((p, j) => {
    const q = after.outer.nodes[g][retained[i]][normal[g][j]]; near(q.x, p.x); near(q.y, p.y);
  })));
  before.outer.allocation.groups.forEach((group, g) => group.forEach((tube, j) => {
    const mass = after.outer.allocation.groups[g].slice(normal[g][j], normal[g][j + 1]).reduce((sum, q) => sum + q.massFlow, 0);
    near(mass / tube.massFlow, 1, 1e-13);
  }));
  let auxChanges = 0, phaseChanges = 0;
  for (const p of parent.bl.stations) {
    const q = child.bl.stations.find(s => s.kind === p.kind && s.body === p.body && s.side === p.side && s.i === retained[p.i]);
    assert.ok(q);
    for (const key of ['theta', 'deltaStar', 'ue']) near(after.layers.states[q.id][key], before.layers.states[p.id][key], 1e-14);
    if (after.layers.states[q.id].aux !== before.layers.states[p.id].aux) auxChanges++;
    if (q.regime !== p.regime) phaseChanges++;
  }
  assert.ok(auxChanges > 0, 'Exercise the refiner amplification reinitialization.');
  assert.ok(phaseChanges >= 4, 'Old first stations become ordinary downstream stations on all four surfaces.');
  assert.equal(result.diagnostics.auxiliaryChanges.length, auxChanges);
  assert.equal(result.diagnostics.phaseChanges.length, phaseChanges);
  const geometry = directStreamtubeVolumeGeometry(after.outer.nodes);
  assert.equal(geometry.valid, true); assert.deepEqual(geometry.concavePrimal, []);
  const saved = serialize(result.checkpoint), savedBefore = structuredClone(saved), replay = create(saved.restart);
  assert.deepEqual(replay.evaluate(replay.initial).residual, after.residual);
  const zero = solveCoupledStreamtubeIses(undefined, { ...controls, resume: saved, maxIterations: 0 });
  assert.equal(zero.initialRedistribution.resumed, true); assert.deepEqual(zero.initialRedistribution.passages, []);
  assert.equal(zero.history.length, 1); assert.equal(zero.linearDiagnostics.solves, 0);
  assert.deepEqual(zero.families, after.families); assert.deepEqual(zero.flow.nodes, after.outer.nodes);
  assert.equal(zero.boundaryLayer.surfaces.length, 4); assert.equal(zero.boundaryLayer.wakes.length, 2);
  assert.deepEqual(serialize(zero.checkpoint.continuation), saved.continuation);
  assert.equal(result.diagnostics.initialSMOVERepeated, false); assert.equal(result.diagnostics.exactSerializedReplay, true);

  const corruptPhysics = structuredClone(restart); corruptPhysics.options.reynolds *= 2;
  assert.throws(() => nestedRefinementCheckpoint(checkpoint, corruptPhysics, config), /physics|inadmissible|replay|thickness/i);
  const corruptCounts = structuredClone(config); corruptCounts.streamwiseSubdivisions[0] = 2;
  assert.throws(() => nestedRefinementCheckpoint(checkpoint, restart, corruptCounts), /subdivisions/);
  const corruptHistory = structuredClone(checkpoint); corruptHistory.continuation.fractions[0][1] = 0;
  assert.throws(() => nestedRefinementCheckpoint(corruptHistory, restart, config), /inlet fractions/);
  const missingState = structuredClone(restart); delete missingState.initialEuler;
  assert.throws(() => nestedRefinementCheckpoint(checkpoint, missingState, config), /cold initialization/);
  const corruptMass = structuredClone(restart); corruptMass.input.weights[0][0] *= 1.001;
  assert.throws(() => nestedRefinementCheckpoint(checkpoint, corruptMass, config), /weights|masses|inadmissible|geometry/i);
  assert.deepEqual(checkpoint, frozen); assert.deepEqual(restart, frozenRestart); assert.deepEqual(config, frozenConfig);
  assert.deepEqual(saved, savedBefore); assert.deepEqual(serialize(root.checkpoint), originalRoot);
  assert.deepEqual(parent.evaluate(parent.initial).residual, before.residual);
  t.diagnostic(JSON.stringify({ parentUnknowns: parent.n, unknowns: child.n, parentIterations: root.history.length - 1,
    streamwiseSubdivisions, normalSubdivisions, families: after.families, nodeError: result.diagnostics.nodeError,
    physicalBLError: result.diagnostics.physicalBLError, inletInterpolationError: result.diagnostics.inletInterpolationError,
    wakeDistanceError: result.diagnostics.wakeDistanceError, auxiliaryError: result.diagnostics.auxiliaryError, auxChanges, phaseChanges }));
});
