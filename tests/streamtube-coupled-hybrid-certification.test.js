// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { certifyCoupledStreamtubeHybrid } from '../src/euler/tests/streamtube-coupled-hybrid-certification.js';
import { initializeCoupledStreamtubeFromFlow } from '../src/euler/tests/streamtube-coupled-flow-restart.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const serial = v => JSON.parse(JSON.stringify(v, (_, x) => ArrayBuffer.isView(x) ? Array.from(x) : x));
const research = { epsilonP: 1e-5, upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
let source, startupError;
function accepted() {
  if (startupError) throw startupError;
  if (!source) try {
    const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2, mach: .2 }), streamwiseMode: 'isentropic' };
    source = solveCoupledStreamtubeIses(input, { edgeMatching: 'section-velocity', transitionMode: 'automatic',
      maxIterations: 12, tolerance: 1e-10, stepAcceptance: 'admissible' });
    assert.equal(source.converged, true, source.reason);
    assert.equal(source.mesh.quality.valid, true); assert.ok(source.x.length < 500);
    assert.equal(source.boundaryLayer.surfaces.length, 4); assert.equal(source.boundaryLayer.wakes.length, 2);
  } catch (error) { startupError = error; throw error; }
  return serial(source.checkpoint);
}

test('hybrid certification carries explicit native Hk-floor direction metadata without changing any physical row', () => {
  const legacy = accepted(), native = structuredClone(legacy);
  native.restart.options.hkFloorLinearization = 'native';
  const exact = certifyCoupledStreamtubeHybrid(legacy, research);
  const revised = certifyCoupledStreamtubeHybrid(native, research);
  assert.equal(revised.checkpoint.restart.options.hkFloorLinearization, 'native');
  assert.equal(revised.diagnostics.residualRowsExactlyEqual, source.x.length);
  const a = initializeCoupledStreamtubeFromFlow(.2, exact.checkpoint);
  const b = initializeCoupledStreamtubeFromFlow(.2, revised.checkpoint);
  assert.equal(b.system.bl.kernel.parameters.hkFloorLinearization, 'native');
  assert.deepEqual(b.initial, a.initial);
  assert.deepEqual(b.value.residual, a.value.residual);
  const normalized = structuredClone(revised.checkpoint); delete normalized.restart.options.hkFloorLinearization;
  assert.deepEqual(normalized, exact.checkpoint);
});

test('explicit same-Mach certification preserves every row, packed BL, geometry, phase and maintenance field', t => {
  const checkpoint = accepted(), before = structuredClone(checkpoint);
  assert.throws(() => initializeCoupledStreamtubeFromFlow(.2, checkpoint), /requires explicit compressible hybrid/);
  const certified = certifyCoupledStreamtubeHybrid(checkpoint, research);
  assert.deepEqual(Object.keys(certified), ['checkpoint', 'diagnostics']);
  assert.equal(certified.diagnostics.certified, true);
  assert.equal(certified.diagnostics.maximumSpeedBias, 0);
  assert.equal(certified.diagnostics.residualRowsExactlyEqual, source.x.length);
  assert.equal(certified.diagnostics.boundaryLayers, 4); assert.equal(certified.diagnostics.wakes, 2);
  assert.equal(certified.diagnostics.strictCompleteFlowReplay, true);
  for (const key of ['newtonUpdates', 'linearSolves', 'redistributions', 'coldInitializations', 'boundaryLayerInitializations'])
    assert.equal(certified.diagnostics.operations[key], 0);
  const changed = structuredClone(certified.checkpoint);
  delete changed.restart.input.hybrid; delete changed.restart.input.upwind;
  changed.restart.input.streamwiseMode = 'isentropic'; delete changed.restart.options.blThermodynamics;
  assert.deepEqual(changed, before, 'Only the three explicit model-control fields and BL model label may change.');
  const prepared = initializeCoupledStreamtubeFromFlow(.2, certified.checkpoint);
  assert.deepEqual(prepared.value.residual, source.residual);
  assert.deepEqual(serial(prepared.initial.subarray(prepared.system.ne)), before.restart.initialBL);
  assert.deepEqual(prepared.system.bl.snapshotActive(), before.restart.options.transitionState);
  assert.deepEqual(serial(prepared.value.outer.nodes), before.restart.initialEuler.nodes);
  assert.deepEqual(certified.checkpoint.restart.initialEuler.undisplacedNodes, before.restart.initialEuler.undisplacedNodes);
  assert.deepEqual(certified.checkpoint.continuation, before.continuation);
  assert.deepEqual(checkpoint, before);
  t.diagnostic(JSON.stringify(certified.diagnostics));
});

test('any nonzero transport-speed bias rejects certification even with the same physical seed', () => {
  const checkpoint = accepted(), before = structuredClone(checkpoint);
  assert.throws(() => certifyCoupledStreamtubeHybrid(checkpoint, {
    epsilonP: 1e-5, upwind: { ...research.upwind, mcrit: 0 } }), /rejects nonzero speed bias/);
  assert.deepEqual(checkpoint, before);
});

test('research controls and baseline source equations must be explicit and current', () => {
  const checkpoint = accepted();
  for (const options of [{}, { ...research, epsilonP: 0 }, { ...research, epsilonP: NaN }])
    assert.throws(() => certifyCoupledStreamtubeHybrid(checkpoint, options), /explicit positive epsilonP/);
  for (const upwind of [undefined, {}, { ...research.upwind, mucon: NaN }, { ...research.upwind, mcrit: 2 },
    { ...research.upwind, unused: 1 }, { ...research.upwind, boundary: { kind: 'implicit' } }])
    assert.throws(() => certifyCoupledStreamtubeHybrid(checkpoint, { ...research, upwind }), /explicit normalized/);
  for (const mutate of [c => { c.restart.input.streamwiseMode = 'momentum'; },
    c => { c.restart.input.upwind = research.upwind; }, c => { c.restart.input.hybrid = { epsilonP: 1e-5 }; },
    c => { c.restart.options.edgeMatching = 'pressure'; },
    c => { c.restart.options.blThermodynamics = 'historical-common-isentrope'; }]) {
    const changed = structuredClone(checkpoint); mutate(changed);
    assert.throws(() => certifyCoupledStreamtubeHybrid(changed, research), /requires baseline compressible isentropic/);
  }
  const stale = structuredClone(checkpoint); stale.restart.options.reynolds *= 1.01;
  assert.throws(() => certifyCoupledStreamtubeHybrid(stale, research), /residual does not replay exactly/);
});

test('incomplete, stale, unconverged and invalid-maintenance source checkpoints are rejected', () => {
  const checkpoint = accepted();
  for (const c of [undefined, {}, { ...checkpoint, version: 2 }])
    assert.throws(() => certifyCoupledStreamtubeHybrid(c, research), /complete version-1/);
  for (const key of ['nodes', 'undisplacedNodes', 'x']) {
    const c = structuredClone(checkpoint); delete c.restart.initialEuler[key];
    assert.throws(() => certifyCoupledStreamtubeHybrid(c, research), /complete finite supplied/);
  }
  const unconverged = structuredClone(checkpoint); unconverged.families.euler = .01;
  assert.throws(() => certifyCoupledStreamtubeHybrid(unconverged, research), /already satisfy/);
  const stale = structuredClone(checkpoint); stale.restart.initialEuler.x[0] += .001;
  assert.throws(() => certifyCoupledStreamtubeHybrid(stale, research), /residual does not replay exactly/);
  const map = structuredClone(checkpoint); delete map.restart.options.transitionState;
  assert.throws(() => certifyCoupledStreamtubeHybrid(map, research), /transition-interval map/);
  const history = structuredClone(checkpoint); history.continuation.fractions[0][1] = 0;
  assert.throws(() => certifyCoupledStreamtubeHybrid(history, research), /inlet fractions/);
  const normalized = structuredClone(checkpoint); normalized.restart.input.geometryDomain = 'positive-simple';
  assert.throws(() => certifyCoupledStreamtubeHybrid(normalized, research), /normalized ISES/);
});

test('returned checkpoint, diagnostics and caller inputs share no mutable state', () => {
  const checkpoint = accepted(), input = structuredClone(checkpoint), controls = structuredClone(research);
  const result = certifyCoupledStreamtubeHybrid(checkpoint, controls);
  const certified = structuredClone(result.checkpoint);
  checkpoint.restart.initialBL.fill(NaN); controls.upwind.mucon = 99;
  assert.deepEqual(result.checkpoint, certified);
  result.checkpoint.restart.initialEuler.x.fill(NaN);
  result.checkpoint.restart.initialEuler.nodes[0][0][0].x = NaN;
  result.checkpoint.restart.input.upwind.mucon = 99;
  result.checkpoint.continuation.fractions[0][0] = 99;
  assert.equal(result.diagnostics.upwind.mucon, 1);
  assert.deepEqual(accepted(), input);
});
