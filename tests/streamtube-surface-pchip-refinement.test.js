// SPDX-License-Identifier: GPL-2.0-or-later
// Small coupled transfers only; no panel solve, flow Newton, Jacobian or LU.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { twoActiveFiniteBaseWakes } from './fixtures/two-active-finite-base-wakes.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { refineCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-refinement.js';
import { nestedRefinementCheckpoint } from '../src/euler/streamtube-nested-checkpoint.js';
import { nestedRefinementCheckpoint as nodeCertificate } from '../scripts/validation/nested-refinement-checkpoint.js';
import { captureStreamtubeInletFractions } from '../src/geometry/streamtube-grid-maintenance.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';

const plain = v => JSON.parse(JSON.stringify(v, (_, x) => ArrayBuffer.isView(x) ? Array.from(x) : x));
const create = f => createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
const oldText = fs.readFileSync(new URL('../docs/surface-pchip-refinement/before/src/euler/streamtube-coupled-refinement.js', import.meta.url), 'utf8');
const oldModule = await import(`data:text/javascript;base64,${Buffer.from(oldText.replace(/from '([^']+)'/g,
  (_, p) => `from '${new URL(p, new URL('../src/euler/streamtube-coupled-refinement.js', import.meta.url))}'`)).toString('base64')}`);

function checkpoint(input, options, initialEuler, initialBL) {
  input = { ...input, stagnationMotion: 'walls-only', normalStencil: 'body-stations', geometryDomain: 'convex' };
  const p = createCoupledStreamtubeBody(input, { ...options, initialEuler, initialBL });
  const v = p.evaluate(p.initial); assert(p.admissible(p.initial));
  return plain({ version: 1, families: v.families, restart: { input,
    options: { ...options, ...(p.bl.transitionMode === 'automatic' ? { transitionState: p.bl.snapshotActive() } : {}) },
    initialEuler: { x: p.initial.slice(0, p.ne), nodes: v.outer.nodes, undisplacedNodes: v.outer.undisplacedNodes },
    initialBL: p.initial.slice(p.ne) }, continuation: {
    fractions: captureStreamtubeInletFractions(v.outer.nodes, p.euler.layout.bodies),
    lastRedistributedStagnation: v.outer.stagnation, preferredOrdering: 'amd', pivotTolerance: .001,
    iterationGeometry: 'convex', stepAcceptance: 'admissible', stagnationLimiter: 'listing' } });
}
function sharp() {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2, mach: .15, contourPanels: 40 });
  for (const [b, body] of input.bodies.entries()) {
    const old = body.surfaceFractions;
    body.surfaceFractions = { upper: old.map(f => f + (.035 + .01 * b) * f * (1 - f)),
      lower: old.map(f => f - (.025 + .01 * b) * f * (1 - f)) };
  }
  Object.assign(input, { streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5, ismom: 3 },
    upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } });
  return checkpoint(input, { reynolds: 1e6, ncrit: 9, edgeMatching: 'section-velocity', transitionMode: 'automatic',
    tripFractions: [[1, 1], [1, 1]], blThermodynamics: 'historical-common-isentrope' });
}
function transfer(cp, mode) {
  const parent = create(cp.restart), counts = parent.euler.layout.tubes.map(row => Array(row).fill(2));
  const controls = { streamwiseFactor: 2, normalSubdivisions: counts, normalInterpolation: 'streamfunction-quadratic',
    ...(mode ? { streamwiseInterpolation: mode } : {}) };
  const r = refineCoupledStreamtubeBody(cp.restart.input, parent, controls), v = r.system.evaluate(r.system.initial);
  const restart = plain({ input: r.input, options: { ...r.options,
    ...(r.system.bl.transitionMode === 'automatic' ? { transitionState: r.system.bl.snapshotActive() } : {}) },
    initialEuler: { x: r.system.initial.slice(0, r.system.ne), nodes: v.outer.nodes, undisplacedNodes: v.outer.undisplacedNodes },
    initialBL: r.system.initial.slice(r.system.ne) });
  const config = { streamwiseSubdivisions: Array(parent.euler.layout.nx).fill(2), normalSubdivisions: counts,
    ...(mode ? { streamwiseInterpolation: mode } : {}) };
  return { r, v, restart, config, parent, controls };
}

test('omitted interpolation preserves the archived coupled transfer arithmetic and schema', () => {
  const cp = sharp(), before = plain(cp), { r, v, parent, controls } = transfer(cp);
  const old = oldModule.refineCoupledStreamtubeBody(cp.restart.input, parent, controls);
  for (const key of ['input', 'options', 'initialEuler', 'initialBL', 'diagnostics']) assert.deepEqual(plain(r[key]), plain(old[key]), key);
  assert.deepEqual(r.system.initial, old.system.initial);
  assert.deepEqual(v.residual, old.system.evaluate(old.system.initial).residual);
  assert.equal(Object.hasOwn(r.input.gridSpacing, 'streamwiseInterpolation'), false);
  assert.deepEqual(cp, before);
});

test('asymmetric multielement PCHIP transfer certifies all surface BLs, physical density and exact JSON restart', t => {
  const cp = sharp(), before = plain(cp), { r, v, restart, config } = transfer(cp, 'surface-pchip');
  assert.equal(r.system.bl.surfaces.length, 4); assert.equal(r.system.bl.wakes.length, 2);
  assert.equal(r.diagnostics.streamwiseInterpolation, 'surface-pchip');
  assert(r.diagnostics.physicalDensityTransfer && r.system.admissible(r.system.initial));
  const geometry = directStreamtubeVolumeGeometry(v.outer.nodes);
  assert(geometry.valid && geometry.concavePrimal.length === 0);
  const cert = nestedRefinementCheckpoint(cp, restart, config);
  assert.deepEqual(cert, nodeCertificate(cp, restart, config));
  assert(cert.diagnostics.insertedPhysicalBLError < 1e-14);
  assert(cert.diagnostics.nodeError < 2e-12 && cert.diagnostics.massRelativeError < 1e-13);
  const replay = solveCoupledStreamtubeIses(undefined, { resume: cert.checkpoint, maxIterations: 0,
    stepAcceptance: 'admissible', tolerance: 1e-10 });
  assert.equal(replay.linearDiagnostics.solves, 0);
  assert.deepEqual(plain(replay.checkpoint), cert.checkpoint);
  assert.deepEqual(plain(replay.residual), plain(v.residual));
  for (const certify of [nestedRefinementCheckpoint, nodeCertificate]) {
    const wrong = { ...config }; delete wrong.streamwiseInterpolation;
    assert.throws(() => certify(cp, restart, wrong), /interpolation does not match/);
    const altered = structuredClone(restart);
    const station = r.system.bl.stations.find(s => s.kind === 'surface' && s.i % 2 && s.i > r.system.euler.layout.bodies[s.body].leadingIndex + 2);
    altered.initialBL[4 * station.id + 3] *= 1.0001;
    assert.throws(() => certify(cp, altered, config), /surface BL values do not follow/);
  }
  assert.deepEqual(cp, before);
  t.diagnostic(JSON.stringify({ unknowns: r.system.n, families: v.families, quality: r.diagnostics.quality,
    insertedPhysicalBLError: cert.diagnostics.insertedPhysicalBLError, massError: cert.diagnostics.massRelativeError }));
});

test('two active finite-base wakes keep fluid displacement and exact history through surface PCHIP transfer', () => {
  const f = twoActiveFiniteBaseWakes(), v = f.system.evaluate(f.x);
  const cp = checkpoint(f.input, { reynolds: 1e6, ncrit: 9, edgeMatching: 'section-velocity' },
    { x: f.x.slice(0, f.system.ne), nodes: v.outer.nodes, undisplacedNodes: v.outer.undisplacedNodes }, f.x.slice(f.system.ne));
  const { r, restart, config } = transfer(cp, 'surface-pchip');
  const cert = nestedRefinementCheckpoint(cp, restart, config);
  assert.deepEqual(cert, nodeCertificate(cp, restart, config));
  assert.equal(r.system.bl.wakes.length, 2);
  assert(cert.diagnostics.finiteWakeDisplacement.length > 0);
  assert(cert.diagnostics.finiteWakeDisplacement.every(w => w.fluidError < 1e-14));
  assert.equal(cert.diagnostics.exactSerializedReplay, true);
  assert.equal(cert.diagnostics.initialSMOVERepeated, false);
});
