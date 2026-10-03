import test from 'node:test';
import assert from 'node:assert/strict';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const serialize = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const near = (a, b, tolerance = 1e-12) => {
  if (typeof a === 'number' && typeof b === 'number') { assert.ok(Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`); return; }
  if (a && typeof a === 'object') { assert.deepEqual(Object.keys(a), Object.keys(b)); for (const k of Object.keys(a)) near(a[k], b[k], tolerance); }
  else assert.deepEqual(a, b);
};
// Legacy rebasing subtracts then restores BL displacement. Test its
// coordinate roundoff relatively, without a chord-sized absolute allowance.
const sameCoordinates = (a, b) => {
  if (typeof a === 'number' && typeof b === 'number') {
    assert.ok(Math.abs(a - b) <= 4 * Number.EPSILON * Math.max(Math.abs(a), Math.abs(b)), `${a} != ${b}`);
  } else if (a && typeof a === 'object') {
    assert.deepEqual(Object.keys(a), Object.keys(b));
    for (const k of Object.keys(a)) sameCoordinates(a[k], b[k]);
  } else assert.deepEqual(a, b);
};

for (const linearOrdering of ['auto', 'station-auto', 'aligned-auto'])
test(`serialized accepted checkpoints resume the complete two-element solve with ${linearOrdering} ordering`, () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  // The comparison concerns Newton replay. The optional progress monitor
  // starts a new history window on resume and is covered separately.
  const controls = { edgeMatching: 'section-velocity', maxIterations: 12, tolerance: 1e-10,
    stepAcceptance: 'admissible', linearOrdering, iterationRecovery: false };
  const full = solveCoupledStreamtubeIses(input, controls); assert.equal(full.converged, true, full.reason);
  let saved, prefix;
  assert.throws(() => solveCoupledStreamtubeIses(input, { ...controls, onCheckpoint: (c, details) => {
    saved = serialize(c); prefix = serialize(details.history);
    if (details.history.at(-1).iteration === 3) throw new Error('simulated interruption after atomic save');
  } }), /simulated interruption/);
  assert.equal(prefix.at(-1).iteration, 3);
  const unchanged = serialize(saved), zero = solveCoupledStreamtubeIses(undefined, { ...controls, resume: saved, maxIterations: 0 });
  // Reconstructing the legacy displacement chart may change physical or
  // undisplaced coordinates by roundoff. Packed fields, residual families
  // and continuation data must still replay exactly. The explicit raw-chart
  // policy is tested bit-for-bit in streamtube-coupled-geometry-replay.test.js.
  const actual = serialize(zero.checkpoint), expected = structuredClone(saved);
  for (const key of ['nodes', 'undisplacedNodes']) {
    sameCoordinates(actual.restart.initialEuler[key], expected.restart.initialEuler[key]);
    delete actual.restart.initialEuler[key]; delete expected.restart.initialEuler[key];
  }
  assert.deepEqual(actual, expected); assert.deepEqual(zero.families, saved.families);
  assert.equal(zero.initialRedistribution.resumed, true); assert.deepEqual(zero.initialRedistribution.passages, []);
  const resumed = solveCoupledStreamtubeIses(undefined, { ...controls, resume: saved, maxIterations: 9 });
  assert.equal(resumed.converged, true, resumed.reason); assert.deepEqual(serialize(saved), unchanged);
  near(resumed.x, full.x); near(resumed.flow.nodes, full.flow.nodes); near(resumed.families, full.families);
  assert.equal(resumed.boundaryLayer.surfaces.length, 4); assert.equal(resumed.boundaryLayer.wakes.length, 2);
  near(resumed.history.slice(1).map(({ iteration, ...rest }) => rest), full.history.slice(4).map(({ iteration, ...rest }) => rest));
  assert.throws(() => solveCoupledStreamtubeIses(input, { ...controls, resume: saved }), /no separate initial state/);
  assert.throws(() => solveCoupledStreamtubeIses(undefined, { ...controls, resume: saved, iterationGeometry: 'ises-sampled' }), /controls do not match/);
  const broken = structuredClone(saved); broken.families.euler += 1;
  assert.throws(() => solveCoupledStreamtubeIses(undefined, { ...controls, resume: broken }), /does not replay exactly/);
});
