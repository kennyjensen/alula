import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';

const seed = JSON.parse(fs.readFileSync(new URL('fixtures/default-coupled-trip-event.json', import.meta.url)));
const native = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/streamtube-trip-transfer.json', import.meta.url)));
const create = () => createCoupledStreamtubeBody(seed.input, { ...seed.options,
  initialEuler: seed.initialEuler, initialBL: Float64Array.from(seed.initialBL) });
const close = (a, b) => assert.ok(Math.abs(a - b) <= 2e-11 * Math.max(Math.abs(b), 1e-20), `${a} != ${b}`);

test('moving a node across a material trip transfers only its auxiliary state, matching original Fortran in both directions', () => {
  assert.equal(seed.origin.sha256, native.provenance.sourceHash);
  const system = create(), { bl, ne } = system, x = system.initial.slice();
  const saved = bl.snapshotActive(), regimes = bl.stations.map(s => s.regime);
  const trips = bl.surfaces.map(s => [s.tripParameter, s.tripArc]);
  const baseline = system.residual(x), col = system.euler.layout.globals.stagnation[0];
  const id = bl.surfaces[0].ids[saved[0]];
  for (const [index, reference] of native.cases.entries()) {
    x[col] = system.initial[col] + reference.shift;
    const before = x.slice(), geometry = bl.geometry(x.subarray(0, ne));
    const up = bl.surfaces[0].ids[saved[0] - 1];
    for (const [station, expected] of [[up, reference.upstream], [id, reference.downstream]]) {
      close(geometry.coordinates[station].s, expected.s);
      close(bl.scale * x[ne + 4 * station + 1], expected.theta);
      close(bl.scale * x[ne + 4 * station + 2], expected.deltaStar);
      close(x[ne + 4 * station + 3], expected.ue);
    }
    close(geometry.surfaceData[0].tripS, reference.tripS);
    const event = bl.updateActive(x.subarray(ne), x.subarray(0, ne));
    assert.equal(event.changed, true); assert.equal(event.changes.length, 1);
    const change = event.changes[0];
    assert.deepEqual([change.body, change.side, change.from, change.to], index ? [0, 'upper', 5, 4] : [0, 'upper', 4, 5]);
    assert.equal(change.converted.length, 1); assert.equal(change.converted[0].id, id);
    close(x[ne + 4 * id], index ? reference.transitionShear : reference.amplification.value);
    for (let k = 0; k < x.length; k++) if (k !== ne + 4 * id) assert.equal(x[k], before[k], `unexpected change to unknown ${k}`);
    assert.deepEqual(bl.surfaces.map(s => [s.tripParameter, s.tripArc]), trips);
    assert.deepEqual(bl.snapshotActive().slice(1), saved.slice(1));
    assert.equal(system.admissible(x), true);
    const noChange = x.slice(); assert.equal(bl.updateActive(x.subarray(ne), x.subarray(0, ne)).changed, false);
    assert.deepEqual(x, noChange, 'unchanged intervals must not project existing amplification or shear');
  }
  // Roll back both state and active metadata; the same full coupled residual
  // must be recovered, including displacement and wake equations.
  bl.restoreActive(saved); x.set(system.initial);
  assert.deepEqual(bl.stations.map(s => s.regime), regimes);
  assert.deepEqual(system.residual(x), baseline);
});

test('failed material-trip conversion is atomic even after an earlier station has been prepared', () => {
  const system = create(), { bl, ne } = system, x = system.initial.slice();
  x[system.euler.layout.globals.stagnation[0]] += .004;
  assert.ok(bl.activeTargets(x.subarray(0, ne))[0].to - bl.surfaces[0].transition > 1);
  const before = x.slice(), active = bl.snapshotActive(), regimes = bl.stations.map(s => s.regime);
  let calls = 0;
  bl.kernel.transitionCheck = () => ++calls === 1 ? { transition: false, amplification: .5 }
    : { transition: true, amplification: 9 };
  assert.throws(() => bl.updateActive(x.subarray(ne), x.subarray(0, ne)), /Natural transition precedes/);
  assert.equal(calls, 2); assert.deepEqual(x, before);
  assert.deepEqual(bl.snapshotActive(), active); assert.deepEqual(bl.stations.map(s => s.regime), regimes);
  assert.throws(() => bl.restoreActive([5, 7, 22, Infinity]), /Invalid streamtube transition active set/);
  assert.deepEqual(bl.snapshotActive(), active);
});
