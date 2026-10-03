// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { initializeCoupledStreamtubeFromFlow } from '../src/euler/tests/streamtube-coupled-flow-restart.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { refineCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-refinement.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const tolerance = 1e-10;
const serial = v => JSON.parse(JSON.stringify(v, (_, value) => ArrayBuffer.isView(value) ? Array.from(value) : value));
let cached, startupError;
function accepted() {
  if (startupError) throw startupError;
  if (!cached) try {
    const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2, mach: .2 }),
      streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 },
      upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
    const result = solveCoupledStreamtubeIses(input, { edgeMatching: 'section-velocity',
      blThermodynamics: 'historical-common-isentrope', transitionMode: 'automatic',
      maxIterations: 12, tolerance, stepAcceptance: 'admissible' });
    assert.equal(result.converged, true, result.reason);
    assert.equal(result.mesh.quality.valid, true);
    assert.ok(result.x.length < 500, 'Keep the one source solve on the existing tiny fixture.');
    assert.equal(result.boundaryLayer.surfaces.length, 4);
    assert.equal(result.boundaryLayer.wakes.length, 2);
    cached = { checkpoint: serial(result.checkpoint), result };
  } catch (error) { startupError = error; throw error; }
  return { checkpoint: structuredClone(cached.checkpoint), result: cached.result };
}
function resume(checkpoint) {
  const { iterationGeometry, stepAcceptance, stagnationLimiter } = checkpoint.continuation;
  return solveCoupledStreamtubeIses(undefined, { resume: checkpoint, iterationGeometry, stepAcceptance,
    stagnationLimiter, maxIterations: 0, tolerance });
}
function physicalIdentity(a, b) {
  assert.deepEqual(a.outer.nodes, b.outer.nodes);
  assert.deepEqual(a.outer.undisplacedNodes, b.outer.undisplacedNodes);
  for (const key of ['captured', 'stagnation', 'strengths']) assert.deepEqual(a.outer[key], b.outer[key]);
  for (let i = 0; i < a.outer.sections.length; i++) for (let g = 0; g < a.outer.sections[i].length; g++)
    for (let j = 0; j < a.outer.sections[i][g].length; j++) {
      assert.equal(a.outer.sections[i][g][j].rho, b.outer.sections[i][g][j].rho);
      assert.equal(a.outer.sections[i][g][j].q, b.outer.sections[i][g][j].q);
    }
  for (let g = 0; g < a.outer.allocation.groups.length; g++) for (let j = 0; j < a.outer.allocation.groups[g].length; j++)
    assert.equal(a.outer.allocation.groups[g][j].massFlow, b.outer.allocation.groups[g][j].massFlow);
  for (let i = 0; i < a.layers.states.length; i++) for (const key of ['theta', 'deltaStar', 'ue', 'aux', 's'])
    assert.equal(a.layers.states[i][key], b.layers.states[i][key]);
}

test('explicit native Hk-floor directions preserve every physical row while omitted/exact retain the legacy schema', () => {
  const { checkpoint } = accepted(), r = checkpoint.restart;
  const make = policy => createCoupledStreamtubeBody(r.input, { ...r.options,
    initialEuler: r.initialEuler, initialBL: Float64Array.from(r.initialBL),
    ...(policy === undefined ? {} : { hkFloorLinearization: policy }) });
  const legacy = make(), exact = make('exact'), native = make('native');
  assert.equal(Object.hasOwn(legacy.conditions, 'hkFloorLinearization'), false);
  assert.equal(Object.hasOwn(legacy.bl.kernel.parameters, 'hkFloorLinearization'), false);
  assert.deepEqual(exact.conditions, legacy.conditions);
  assert.deepEqual(exact.bl.kernel.parameters, legacy.bl.kernel.parameters);
  assert.equal(native.conditions.hkFloorLinearization, 'native');
  assert.equal(native.bl.kernel.parameters.hkFloorLinearization, 'native');
  for (const system of [exact, native]) {
    assert.deepEqual(system.initial, legacy.initial);
    assert.deepEqual(system.evaluate(system.initial).residual, legacy.evaluate(legacy.initial).residual);
    assert.deepEqual(system.bl.snapshotActive(), legacy.bl.snapshotActive());
    assert.equal(system.admissible(system.initial), true);
  }
  for (const policy of [null, true, 'fortran', ''])
    assert.throws(() => make(policy), /Hk-floor linearization policy/);
});

test('ISES checkpoints inherit the explicit Hk-floor direction policy and reject silent overrides', () => {
  const { checkpoint } = accepted(), native = structuredClone(checkpoint);
  native.restart.options.hkFloorLinearization = 'native';
  const zero = resume(native);
  assert.equal(zero.conditions.hkFloorLinearization, 'native');
  assert.equal(zero.coupledOptions.hkFloorLinearization, 'native');
  assert.deepEqual(serial(zero.checkpoint), native);
  assert.deepEqual(zero.families, checkpoint.families);
  assert.equal(zero.linearDiagnostics.solves, 0);
  const { iterationGeometry, stepAcceptance, stagnationLimiter } = checkpoint.continuation;
  const controls = { iterationGeometry, stepAcceptance, stagnationLimiter, maxIterations: 0, tolerance };
  const exact = solveCoupledStreamtubeIses(undefined, { ...controls, resume: checkpoint, hkFloorLinearization: 'exact' });
  assert.deepEqual(serial(exact.checkpoint), checkpoint);
  assert.equal(Object.hasOwn(exact.conditions, 'hkFloorLinearization'), false);
  for (const [saved, policy] of [[native, 'exact'], [checkpoint, 'native'], [native, 'unknown']])
    assert.throws(() => solveCoupledStreamtubeIses(undefined, { ...controls, resume: saved,
      hkFloorLinearization: policy }), /Hk-floor linearization controls do not match/);
  const malformed = structuredClone(native); malformed.restart.options.hkFloorLinearization = 'unknown';
  assert.throws(() => resume(malformed), /Hk-floor linearization policy/);
});

test('warm entry and nested complete-state refinement retain explicit native direction provenance', () => {
  const { checkpoint } = accepted(), native = structuredClone(checkpoint);
  native.restart.options.hkFloorLinearization = 'native';
  const a = initializeCoupledStreamtubeFromFlow(.2, checkpoint);
  const b = initializeCoupledStreamtubeFromFlow(.2, native);
  assert.deepEqual(serial(b.checkpoint), native);
  assert.equal(b.system.conditions.hkFloorLinearization, 'native');
  assert.deepEqual(b.value.residual, a.value.residual); physicalIdentity(a.value, b.value);
  const options = { streamwiseFactor: 2, normalFactor: 1 };
  const refinedExact = refineCoupledStreamtubeBody(checkpoint.restart.input, a.system, options);
  const refinedNative = refineCoupledStreamtubeBody(native.restart.input, b.system, options);
  assert.ok(refinedNative.system.n < 1000, 'This transfer performs no global Newton or linear solve.');
  assert.equal(refinedNative.options.hkFloorLinearization, 'native');
  assert.equal(refinedNative.system.conditions.hkFloorLinearization, 'native');
  assert.equal(refinedNative.system.bl.kernel.parameters.hkFloorLinearization, 'native');
  assert.deepEqual(refinedNative.input, refinedExact.input);
  assert.deepEqual(refinedNative.system.initial, refinedExact.system.initial);
  assert.deepEqual(refinedNative.system.evaluate(refinedNative.system.initial).residual,
    refinedExact.system.evaluate(refinedExact.system.initial).residual);
});

test('same-Mach complete coupled transfer exactly replays all BLs, wakes, phases and ISES maintenance history', t => {
  const { checkpoint, result } = accepted(), before = structuredClone(checkpoint);
  const prepared = initializeCoupledStreamtubeFromFlow(.2, checkpoint, { tolerance });
  assert.deepEqual(serial(prepared.checkpoint), checkpoint);
  assert.deepEqual(serial(prepared.initial.subarray(prepared.system.ne)), checkpoint.restart.initialBL);
  assert.deepEqual(prepared.system.bl.snapshotActive(), checkpoint.restart.options.transitionState);
  assert.deepEqual(prepared.value.families, checkpoint.families);
  assert.deepEqual(prepared.value.residual, result.residual);
  assert.equal(prepared.diagnostics.sourceConverged, true);
  assert.equal(prepared.diagnostics.targetConverged, true);
  assert.equal(prepared.system.initialization.suppliedBL, true);
  assert.equal(prepared.system.bl.surfaces.length, 4);
  assert.equal(prepared.system.bl.wakes.length, 2);
  assert.ok(Object.values(prepared.diagnostics.operations).every(v => v === 0));
  const zero = resume(prepared.checkpoint);
  assert.equal(zero.converged, true, zero.reason);
  assert.equal(zero.initialRedistribution.resumed, true);
  assert.deepEqual(zero.initialRedistribution.passages, []);
  assert.equal(zero.linearDiagnostics.solves, 0);
  assert.deepEqual(zero.families, checkpoint.families);
  assert.deepEqual(zero.residual, prepared.value.residual);
  assert.deepEqual(zero.boundaryLayer.transitionState, checkpoint.restart.options.transitionState);
  assert.deepEqual(serial(zero.flow.nodes), checkpoint.restart.initialEuler.nodes);
  assert.deepEqual(zero.checkpoint.continuation, checkpoint.continuation);
  assert.deepEqual(checkpoint, before);
  t.diagnostic(JSON.stringify({ unknowns: prepared.system.n, sourceUpdates: result.history.length - 1,
    sourceFamilies: result.families, undisplacedRoundoff: prepared.diagnostics.sourceUndisplacedReconstructionDeparture }));
});

test('new Mach retains complete physical Euler/BL state but recomputes gas, residual and restart families', () => {
  const { checkpoint } = accepted(), before = structuredClone(checkpoint);
  const source = initializeCoupledStreamtubeFromFlow(.2, checkpoint);
  const prepared = initializeCoupledStreamtubeFromFlow(.21, checkpoint);
  assert.deepEqual(prepared.initial, source.initial);
  physicalIdentity(source.value, prepared.value);
  assert.deepEqual(prepared.system.bl.snapshotActive(), source.system.bl.snapshotActive());
  assert.deepEqual(prepared.system.bl.thicknesses(prepared.initial.subarray(prepared.system.ne)),
    source.system.bl.thicknesses(source.initial.subarray(source.system.ne)));
  assert.equal(prepared.system.conditions.reynolds, source.system.conditions.reynolds);
  assert.equal(prepared.system.conditions.ncrit, source.system.conditions.ncrit);
  assert.equal(prepared.system.conditions.mach, .21);
  assert.notEqual(prepared.system.euler.conditions.h0, source.system.euler.conditions.h0);
  assert.notEqual(prepared.system.euler.conditions.rhoTotal, source.system.euler.conditions.rhoTotal);
  assert.ok(prepared.diagnostics.targetResidual > tolerance);
  assert.equal(prepared.diagnostics.targetConverged, false);
  assert.ok(prepared.value.residual.every(Number.isFinite));
  assert.notDeepEqual(prepared.value.families, source.value.families);
  const changed = structuredClone(prepared.checkpoint);
  changed.restart.input.mach = checkpoint.restart.input.mach; changed.families = checkpoint.families;
  assert.deepEqual(changed, checkpoint, 'Only Mach and the actually recomputed residual families may change.');
  const zero = resume(prepared.checkpoint);
  assert.equal(zero.converged, false);
  assert.equal(zero.conditions.mach, .21);
  assert.equal(zero.initialRedistribution.resumed, true);
  assert.equal(zero.linearDiagnostics.solves, 0);
  assert.deepEqual(zero.families, prepared.value.families);
  assert.deepEqual(zero.residual, prepared.value.residual);
  assert.deepEqual(checkpoint, before);
});

test('returned checkpoint and caller data are detached from the prepared numerical system', () => {
  const { checkpoint } = accepted();
  const p = initializeCoupledStreamtubeFromFlow(.21, checkpoint), residual = p.value.residual.slice();
  checkpoint.restart.input.hybrid.epsilonP = NaN;
  checkpoint.restart.initialBL.fill(NaN);
  checkpoint.continuation.fractions[0].fill(NaN);
  p.checkpoint.restart.input.bodies[0].points[0].x = NaN;
  p.checkpoint.restart.initialEuler.nodes[0][0][0].x = NaN;
  p.checkpoint.restart.initialEuler.undisplacedNodes[0][0][0].x = NaN;
  p.checkpoint.restart.initialEuler.x.fill(NaN);
  p.checkpoint.restart.initialBL.fill(NaN);
  p.checkpoint.continuation.fractions[0].fill(NaN);
  assert.deepEqual(p.system.evaluate(p.initial).residual, residual);
  assert.ok(p.initial.every(Number.isFinite));
});

test('coupled warm entry rejects incomplete checkpoints, changed model assumptions and invalid targets', () => {
  const { checkpoint } = accepted(), before = structuredClone(checkpoint);
  for (const mach of [0, -1, 1, 1.2, NaN, Infinity])
    assert.throws(() => initializeCoupledStreamtubeFromFlow(mach, checkpoint), /target or tolerance/);
  for (const tol of [0, -1, NaN])
    assert.throws(() => initializeCoupledStreamtubeFromFlow(.21, checkpoint, { tolerance: tol }), /target or tolerance/);
  for (const invalid of [undefined, {}, { ...checkpoint, version: 2 }])
    assert.throws(() => initializeCoupledStreamtubeFromFlow(.21, invalid), /version-1 ISES checkpoint/);
  for (const key of ['initialBL', 'initialEuler']) {
    const bad = structuredClone(checkpoint); delete bad.restart[key];
    assert.throws(() => initializeCoupledStreamtubeFromFlow(.21, bad), /complete finite Euler\/BL/);
  }
  for (const key of ['x', 'nodes', 'undisplacedNodes']) {
    const bad = structuredClone(checkpoint); delete bad.restart.initialEuler[key];
    assert.throws(() => initializeCoupledStreamtubeFromFlow(.21, bad), /complete finite Euler\/BL/);
  }
  for (const change of [c => { c.restart.input.streamwiseMode = 'momentum'; },
    c => { delete c.restart.options.blThermodynamics; },
    c => { c.restart.options.edgeMatching = 'pressure'; },
    c => { c.restart.options.edgeMatching = 'section-velocity-distance'; },
    c => { delete c.restart.input.upwind; }, c => { c.restart.input.flowModel = 'incompressible'; }]) {
    const bad = structuredClone(checkpoint); change(bad);
    assert.throws(() => initializeCoupledStreamtubeFromFlow(.21, bad), /requires explicit compressible hybrid/);
  }
  const missingPhase = structuredClone(checkpoint); delete missingPhase.restart.options.transitionState;
  assert.throws(() => initializeCoupledStreamtubeFromFlow(.21, missingPhase), /transition-interval map/);
  assert.deepEqual(checkpoint, before);
});

test('source convergence, actual geometry and maintenance history are mandatory, even when the saved families claim success', () => {
  const { checkpoint } = accepted();
  const stopped = structuredClone(checkpoint); stopped.families.euler = 1e-3;
  assert.throws(() => initializeCoupledStreamtubeFromFlow(.21, stopped), /already satisfy the residual tolerance/);
  const stale = structuredClone(checkpoint); stale.restart.initialEuler.x[0] += .001;
  const beforeStale = structuredClone(stale);
  assert.throws(() => initializeCoupledStreamtubeFromFlow(.21, stale), /does not replay exactly/);
  assert.deepEqual(stale, beforeStale);
  const folded = structuredClone(checkpoint);
  folded.restart.initialEuler.nodes[0][2][1].x = folded.restart.initialEuler.nodes[0][4][1].x;
  const beforeFolded = structuredClone(folded);
  assert.throws(() => initializeCoupledStreamtubeFromFlow(.21, folded), /valid convex physical grid/);
  assert.deepEqual(folded, beforeFolded);
  for (const change of [c => { c.continuation.preferredOrdering = 'unknown'; },
    c => { c.continuation.lastRedistributedStagnation[0] = NaN; },
    c => { c.continuation.fractions[0][1] = -1; },
    c => { c.restart.input.geometryDomain = 'positive-simple'; }]) {
    const bad = structuredClone(checkpoint); change(bad);
    assert.throws(() => initializeCoupledStreamtubeFromFlow(.21, bad), /checkpoint (maintenance|inlet)|normalized ISES geometry/);
  }
  // Reconstructing a target directly is a frozen control, not a substitute
  // for the source-convergence and maintenance-history guards above.
  const p = initializeCoupledStreamtubeFromFlow(.21, checkpoint);
  const r = p.checkpoint.restart;
  const independent = createCoupledStreamtubeBody(r.input, { ...r.options,
    initialEuler: r.initialEuler, initialBL: Float64Array.from(r.initialBL) });
  assert.deepEqual(independent.evaluate(independent.initial).families, p.checkpoint.families);
});
