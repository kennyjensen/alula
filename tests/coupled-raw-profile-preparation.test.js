// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { prepareRawCoupledProposal } from '../scripts/validation/coupled-raw-profile-preparation.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { transferStreamtubeGeometry } from '../src/euler/streamtube-geometry.js';
import { relaxXfoilProfiles } from '../src/viscous/xfoil-profile-relaxation.js';

const sourcePath = 'docs/coupled-current-profile-preparation/two-element-six-update/initial.json';
const retained = JSON.parse(fs.readFileSync(sourcePath)).checkpoint;
const copy = v => JSON.parse(JSON.stringify(v, (_, x) => ArrayBuffer.isView(x) ? Array.from(x) : x));
const rawState = cp => Float64Array.from([...cp.restart.initialEuler.x, ...cp.restart.initialBL]);
// Geometry-only conversion of the retained tiny centerline fixture. All
// physical nodes/density/mass/BL fields stay fixed; the additional bank
// coordinates introduce no initialization march, Jacobian or Newton solve.
const f = retained.restart;
const old = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
const input = { ...f.input, wakeGeometry: 'independent-banks' };
const target = createStreamtubeBodySystem({ ...input, displacement: old.bl.thicknesses(old.initial.subarray(old.ne)) });
const converted = transferStreamtubeGeometry(old.euler, old.initial.subarray(0, old.ne), target);
const nodes = target.decode(converted);
const initialEuler = { x: converted, nodes: nodes.nodes, undisplacedNodes: nodes.undisplacedNodes };
const prepared = createCoupledStreamtubeBody(input, { ...f.options, initialEuler, initialBL: f.initialBL });
const value = prepared.evaluate(prepared.initial);
const source = copy({ ...retained, families: value.families, restart: { ...f, input,
  initialEuler: { x: prepared.initial.slice(0, prepared.ne), nodes: value.outer.nodes, undisplacedNodes: value.outer.undisplacedNodes } } });

test('the tiny independent-bank fixture conversion preserves physical geometry, density, mass and BL', t => {
  const before = old.euler.decode(old.initial.subarray(0, old.ne));
  const after = prepared.euler.decode(prepared.initial.subarray(0, prepared.ne));
  const difference = Math.max(...before.nodes.flatMap((g, gi) => g.flatMap((r, i) => r.map((p, j) =>
    Math.hypot(p.x - after.nodes[gi][i][j].x, p.y - after.nodes[gi][i][j].y)))));
  assert.ok(difference <= 1e-14);
  assert.deepEqual(Array.from(prepared.initial.slice(0, target.layout.densityCount)), Array.from(old.initial.slice(0, old.euler.layout.densityCount)));
  assert.deepEqual(after.captured, before.captured); assert.deepEqual(after.stagnation, before.stagnation);
  assert.deepEqual(after.allocation.groups.map(g => g.map(t => t.massFlow)), before.allocation.groups.map(g => g.map(t => t.massFlow)));
  assert.deepEqual(source.restart.initialBL, retained.restart.initialBL);
  t.diagnostic(JSON.stringify({ oldUnknowns: old.n, independentUnknowns: prepared.n, maximumNodeDifference: difference,
    families: value.families, conversionOperations: { coupledConstructors: 2, eulerConstructors: 1,
      fullResidualEvaluations: 1, nativeMarches: 0, jacobians: 0, linearSolves: 0, newtonUpdates: 0 } }));
});

const identity = input => ({ surfaces: input.surfaces.map((states, k) => ({ states, transition: input.phases[k] })),
  wake: input.wake, localConvergenceWarnings: [], messages: [] });

test('raw handoff preserves an identity preparation, all two-body profiles, source and history', () => {
  const cp = copy(source), raw = rawState(cp), before = copy(cp), rawBefore = copy(raw);
  let calls = 0;
  const result = prepareRawCoupledProposal(cp, raw, { relaxProfiles: input => { calls++; return identity(input); } });
  assert.equal(calls, 2);
  assert.equal(result.report.canonicalReplayExact, true);
  assert.equal(result.report.geometryExtension.maximumInteriorChange, 0);
  assert.equal(result.report.maximumPhysicalDisplacementChange, 0);
  assert.deepEqual(result.checkpoint.restart.initialBL, cp.restart.initialBL);
  assert.deepEqual(result.checkpoint.restart.options, cp.restart.options);
  assert.deepEqual(result.checkpoint.families, cp.families);
  assert.deepEqual(result.checkpoint.continuation, cp.continuation);
  assert.deepEqual(cp, before); assert.deepEqual(copy(raw), rawBefore);
  assert.equal(result.report.operations.genericPhaseConversions, 0);
  assert.equal(result.report.operations.jacobians, 0);
  assert.equal(result.report.operations.linearSolves, 0);
  assert.equal(result.report.operations.newtonUpdates, 0);
  assert.equal(result.records.raw.states.length, cp.restart.initialBL.length / 4);
});

test('native input uses raw physical profiles and old above-critical N before any phase conversion', () => {
  const cp = copy(source), raw = rawState(cp), ne = cp.restart.initialEuler.x.length;
  // Body 1 has old phase 2: station 1 of each surface is still laminar.
  // Locate it through the already saved BL station metadata, without a solve.
  const saved = JSON.parse(fs.readFileSync(sourcePath));
  const surface = saved.boundaryLayer.surfaces.find(s => s.body === 1 && s.side === 'upper');
  const id = surface.ids[1];
  raw[ne + 4 * id] = cp.restart.options.ncrit + .125;
  raw[0] += 1e-7; // a genuine raw Euler density increment must survive adoption
  const original = copy(raw); let sawAbove = false;
  const result = prepareRawCoupledProposal(cp, raw, { relaxProfiles: (input, parameters) => {
    if (input.surfaces.some((p, k) => p.slice(0, input.phases[k]).some(s => s.aux >= parameters.ncrit))) sawAbove = true;
    return relaxXfoilProfiles(input, parameters);
  } });
  assert.ok(sawAbove);
  assert.equal(result.records.raw.states[id].aux, cp.restart.options.ncrit + .125);
  assert.equal(result.checkpoint.restart.initialEuler.x[0], raw[0]);
  assert.deepEqual(result.records.raw.phase, cp.restart.options.transitionState);
  assert.ok(result.records['native-phase-targets'].every(t => t.from === t.to && !t.amplificationReconciliation));
  assert.equal(result.report.canonicalReplayExact, true);
  assert.deepEqual(copy(raw), original);
});

test('native phase disagreement is retained and rejected without generic conversion or any full evaluation', () => {
  const cp = copy(source), raw = rawState(cp); let failure;
  try {
    prepareRawCoupledProposal(cp, raw, { relaxProfiles: input => {
      const r = identity(input); if (input.phases[0] === 1) r.surfaces[0].transition = 2; return r;
    } });
  } catch (error) { failure = error; }
  assert.ok(failure, 'A differing native phase must not be silently replaced.');
  assert.match(failure.message, /Native phase cannot hand off unchanged/);
  assert.equal(failure.preparation.report.stage, 'native phase eligibility');
  assert.equal(failure.preparation.report.operations.explicitEvaluations, 0);
  assert.equal(failure.preparation.report.operations.admissibilityCalls, 0);
  assert.equal(failure.preparation.report.operations.genericPhaseConversions, 0);
  assert.ok(failure.preparation.records['native-phase-targets'].some(t => t.from !== t.to));
  assert.ok(failure.preparation.sourceUnchanged && failure.preparation.rawUnchanged);
});

test('native thermal failure retains raw inputs and cannot publish a candidate checkpoint', () => {
  const cp = copy(source), raw = rawState(cp), ne = cp.restart.initialEuler.x.length;
  raw[ne + 3] = 20;
  let failure; try { prepareRawCoupledProposal(cp, raw); } catch (error) { failure = error; }
  assert.ok(failure); assert.match(failure.message, /thermal|enthalpy|BL edge/i);
  assert.equal(failure.preparation.report.prepared, false);
  assert.equal(failure.preparation.report.operations.explicitEvaluations, 0);
  assert.equal(failure.preparation.records.raw.states[0].ue, 20);
  assert.equal(failure.preparation.records.prepared, undefined);
  assert.ok(failure.preparation.sourceUnchanged && failure.preparation.rawUnchanged);
});

test('injected adapters must also obey returned N/phase and shear domains', () => {
  const cp = copy(source), raw = rawState(cp); let failure;
  try { prepareRawCoupledProposal(cp, raw, { relaxProfiles: (input, parameters) => {
    const r = identity(input); r.surfaces[0].states[0].aux = parameters.ncrit; return r;
  } }); } catch (error) { failure = error; }
  assert.ok(failure); assert.match(failure.message, /output auxiliary state/);
  assert.equal(failure.preparation.report.operations.explicitEvaluations, 0);
  assert.equal(failure.preparation.report.prepared, false);
  assert.ok(failure.preparation.sourceUnchanged && failure.preparation.rawUnchanged);
});
