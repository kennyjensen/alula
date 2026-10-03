// SPDX-License-Identifier: GPL-2.0-or-later
// Supplied tiny checkpoint only: no startup, Jacobian, factorization or Newton.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';

const plain = x => JSON.parse(JSON.stringify(x, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const sourceText = fs.readFileSync(new URL('../docs/coupled-current-profile-preparation/two-element-six-update/initial.json', import.meta.url), 'utf8');
const original = JSON.parse(sourceText).checkpoint;
let cached;
function fixture() {
  if (!cached) {
    const cp = plain(original), f = cp.restart;
    const topology = { nx: f.input.outerLower.length - 1, tubes: f.input.weights.map(row => row.length),
      bodies: f.input.bodies.map(({ leadingIndex, trailingIndex }) => ({ leadingIndex, trailingIndex })) };
    let cut = 0;
    f.input = { ...f.input, streamwiseMode: 'hybrid', upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } },
      hybrid: { epsilonP: 1e-5, ismom: 3, entropyRegions: {
        version: 1, parentCheckpointSha256: createHash('sha256').update(JSON.stringify(original)).digest('hex'), topology,
        regions: topology.bodies.map((body, b) => {
          cut += topology.tubes[b];
          return { body: b, throughRow: body.leadingIndex + 1, lowerTube: cut - 1, upperTube: cut + 1 };
        }),
      } } };
    f.options.blThermodynamics = 'historical-common-isentrope';
    // This is an explicitly re-evaluated equation-selection fixture, not a
    // claim that the original accepted equations or residuals were retained.
    const system = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
    const value = system.evaluate(system.initial);
    assert.equal(system.n, 349);
    assert.deepEqual(Array.from(system.initial), [...original.restart.initialEuler.x, ...original.restart.initialBL]);
    assert.deepEqual(system.bl.snapshotActive(), original.restart.options.transitionState);
    assert.deepEqual(plain(value.outer.nodes), original.restart.initialEuler.nodes);
    assert.deepEqual(plain(value.outer.undisplacedNodes), original.restart.initialEuler.undisplacedNodes);
    cp.families = plain(value.families);
    cached = { cp, residual: Array.from(value.residual) };
  }
  return plain(cached);
}

function replay(cp, onCheckpoint) {
  const { iterationGeometry, stepAcceptance, stagnationLimiter } = cp.continuation;
  return solveCoupledStreamtubeIses(undefined, { resume: cp, iterationGeometry, stepAcceptance, stagnationLimiter,
    maxIterations: 0, tolerance: 1e-10, onCheckpoint });
}

test('ISMOM3 bound regions survive zero-update checkpoint reporting and a JSON restart exactly', () => {
  const { cp, residual } = fixture(), before = plain(cp), published = [];
  const first = replay(cp, checkpoint => published.push(plain(checkpoint)));
  assert.equal(first.linearDiagnostics.solves, 0);
  assert.equal(first.history.length, 1);
  assert.equal(published.length, 1);
  assert.equal(first.mesh.quality.valid, true);
  assert.deepEqual(plain(first.residual), residual);
  assert.deepEqual(plain(first.checkpoint), cp);
  assert.deepEqual(published[0], cp);
  assert.deepEqual(first.solverInput.hybrid, cp.restart.input.hybrid);
  const second = replay(plain(first.checkpoint));
  assert.equal(second.linearDiagnostics.solves, 0);
  assert.deepEqual(plain(second.checkpoint), cp);
  assert.deepEqual(plain(second.residual), residual);
  assert.deepEqual(cp, before, 'A supplied checkpoint must remain unchanged.');
});

test('restart rejects stale streamwise, mass-tube and LE/TE topology before publishing any checkpoint', () => {
  const changes = [
    t => { t.nx++; },
    t => { t.tubes[2]++; },
    t => { t.bodies[0].leadingIndex--; },
    t => { t.bodies[0].trailingIndex++; },
  ];
  for (const change of changes) {
    const { cp } = fixture(); change(cp.restart.input.hybrid.entropyRegions.topology);
    const before = plain(cp); let published = 0;
    assert.throws(() => replay(cp, () => { published++; }), /topology does not match this grid/);
    assert.equal(published, 0);
    assert.deepEqual(cp, before);
  }
});

test('family maxima alone cannot certify the provenance of a changed valid equation mask', () => {
  const { cp, residual } = fixture();
  cp.restart.input.hybrid.entropyRegions.regions[0].throughRow++;
  const before = plain(cp), replayed = replay(cp);
  assert.equal(replayed.linearDiagnostics.solves, 0);
  assert.deepEqual(replayed.families, cp.families);
  const changedRows = Array.from(replayed.residual).flatMap((v, i) => v === residual[i] ? [] : [i]);
  assert.ok(changedRows.length > 0, 'Different equations can share the same maximum residual.');
  assert.deepEqual(plain(replayed.checkpoint.restart), cp.restart);
  // Region preparation must independently bind the parent hash and preserve
  // before/after full residuals; the existing family replay is not that proof.
  assert.deepEqual(cp, before);
});
