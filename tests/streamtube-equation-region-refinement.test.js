// SPDX-License-Identifier: GPL-2.0-or-later
// Pure index/mass-region controls only; no Euler/BL system or numerical solve.
import test from 'node:test';
import assert from 'node:assert/strict';
import { refineStreamtubeEquationRegions } from '../src/euler/tests/streamtube-equation-region-refinement.js';

const hash = 'a'.repeat(64), nextHash = 'b'.repeat(64);
const cumulative = counts => counts.reduce((map, n) => [...map, map.at(-1) + n], [0]);
const mappedTopology = (parent, streamwise, normal) => {
  const retained = cumulative(streamwise);
  return { nx: retained.at(-1), tubes: normal.map(row => row.reduce((a, b) => a + b, 0)),
    bodies: parent.bodies.map(b => ({ leadingIndex: retained[b.leadingIndex], trailingIndex: retained[b.trailingIndex] })) };
};
const input = (parent, streamwiseSubdivisions, normalSubdivisions, rest = {}) => ({ parent,
  child: mappedTopology(parent, streamwiseSubdivisions, normalSubdivisions), streamwiseSubdivisions, normalSubdivisions,
  parentCheckpointSha256: hash, ...rest });
const selected = (regions, i, globalTube) => regions.some(r => i <= r.throughRow && globalTube >= r.lowerTube && globalTube < r.upperTube);

test('retained RAE index counts map the inclusive nose endpoint and original wall-mass bands', () => {
  const parent = { nx: 191, tubes: [14, 14], bodies: [{ leadingIndex: 64, trailingIndex: 127 }] };
  const request = input(parent, Array(191).fill(2), [Array(13).fill(1).concat(4), [4, ...Array(13).fill(1)]]);
  const before = structuredClone(request), result = refineStreamtubeEquationRegions(request);
  assert.deepEqual(result, { version: 1, parentCheckpointSha256: hash, topology: request.child,
    regions: [{ body: 0, throughRow: 148, lowerTube: 10, upperTube: 24 }] });
  assert.deepEqual(request, before);
  assert.ok(Object.isFrozen(result)); assert.ok(Object.isFrozen(result.regions));
  assert.ok(Object.isFrozen(result.regions[0])); assert.ok(Object.isFrozen(result.topology));
  assert.ok(Object.isFrozen(result.topology.tubes)); assert.ok(Object.isFrozen(result.topology.bodies[0]));
  let inherited = 0, ordinary = 0, changes = 0;
  for (let i = 1; i < request.child.nx; i++) for (let tube = 0; tube < 34; tube++) {
    const a = selected(result.regions, i, tube), b = i <= 138 && tube >= 13 && tube < 21;
    inherited += a; ordinary += b; changes += a !== b;
  }
  assert.deepEqual({ inherited, ordinary, changes }, { inherited: 2072, ordinary: 1104, changes: 968 });
  request.child.tubes[0] = 99; assert.equal(result.topology.tubes[0], 17);
});

test('nonuniform maps preserve staggered multielement regions across a narrow shared passage', () => {
  const parent = { nx: 32, tubes: [6, 3, 7], bodies: [
    { leadingIndex: 4, trailingIndex: 16 }, { leadingIndex: 12, trailingIndex: 24 }] };
  const streamwise = Array(32).fill(1); streamwise[1] = 2; streamwise[7] = 3; streamwise[17] = 4;
  const normal = [[1, 2, 1, 3, 1, 2], [4, 1, 2], [2, 1, 3, 1, 1, 2, 1]];
  const request = input(parent, streamwise, normal), result = refineStreamtubeEquationRegions(request);
  assert.deepEqual(result.regions, [
    { body: 0, throughRow: 17, lowerTube: 3, upperTube: 19 },
    { body: 1, throughRow: 28, lowerTube: 8, upperTube: 24 },
  ]);
  // Check the union on retained cross-lines and every child tube descendant,
  // independently using parent tube membership rather than mapped bounds.
  const stations = cumulative(streamwise), descendants = [];
  normal.flat().forEach((count, parentTube) => { for (let k = 0; k < count; k++) descendants.push(parentTube); });
  for (let i = 1; i < parent.nx; i++) for (let tube = 0; tube < descendants.length; tube++) {
    const old = descendants[tube], expected = (i <= 14 && old >= 2 && old < 10) || (i <= 22 && old >= 5 && old < 13);
    assert.equal(selected(result.regions, stations[i], tube), expected, `station ${i}, tube ${tube}`);
  }
  assert.equal(selected(result.regions, 18, 3), false, 'First body region ends earlier.');
  assert.equal(selected(result.regions, 18, 8), true, 'Second body region remains active in the overlap.');
});

test('default regions clamp to real inlet/outlet and transverse extents before mapping', () => {
  const parent = { nx: 12, tubes: [2, 1, 2], bodies: [
    { leadingIndex: 1, trailingIndex: 3 }, { leadingIndex: 8, trailingIndex: 10 }] };
  const request = input(parent, Array(12).fill(2), [[2, 1], [3], [1, 2]]);
  const result = refineStreamtubeEquationRegions(request);
  assert.deepEqual(result.regions, [{ body: 0, throughRow: 22, lowerTube: 0, upperTube: 9 },
    { body: 1, throughRow: 22, lowerTube: 0, upperTube: 9 }]);
  assert.equal(selected(result.regions, 23, 4), false, 'Do not extend the inherited row to a new outlet-adjacent row.');
});

test('repeated refinement maps explicit inherited regions and replaces only parent provenance', () => {
  const parent = { nx: 24, tubes: [5, 5], bodies: [{ leadingIndex: 3, trailingIndex: 16 }] };
  const firstInput = input(parent, Array(24).fill(2), [Array(5).fill(2), Array(5).fill(2)]);
  const first = refineStreamtubeEquationRegions(firstInput), original = JSON.stringify(first);
  const secondInput = input(first.topology, Array(48).fill(2), [Array(10).fill(2), Array(10).fill(2)],
    { entropyRegions: first, parentCheckpointSha256: nextHash });
  const second = refineStreamtubeEquationRegions(secondInput);
  const direct = refineStreamtubeEquationRegions(input(parent, Array(24).fill(4), [Array(5).fill(4), Array(5).fill(4)],
    { parentCheckpointSha256: nextHash }));
  assert.deepEqual(second, direct); assert.equal(JSON.stringify(first), original);
  assert.equal(second.parentCheckpointSha256, nextHash);
  assert.deepEqual(second.regions, [{ body: 0, throughRow: 52, lowerTube: 4, upperTube: 36 }]);
  assert.notEqual(second.regions[0].throughRow, second.topology.bodies[0].leadingIndex + 10);
});

test('reject malformed subdivisions, nonnested child topology, no refinement and missing source identity', () => {
  const parent = { nx: 20, tubes: [5, 5], bodies: [{ leadingIndex: 3, trailingIndex: 12 }] };
  const valid = input(parent, Array(20).fill(2), [Array(5).fill(1), Array(5).fill(1)]);
  for (const mutate of [
    q => { delete q.parentCheckpointSha256; }, q => { q.parentCheckpointSha256 = 'A'.repeat(64); },
    q => { q.parentCheckpointSha256 = 'a'.repeat(63); }, q => { delete q.streamwiseSubdivisions; },
    q => { q.streamwiseSubdivisions[0] = 0; }, q => { q.streamwiseSubdivisions[0] = 5; },
    q => { q.streamwiseSubdivisions[0] = 1.5; }, q => { q.streamwiseSubdivisions.pop(); },
    q => { q.normalSubdivisions[0][0] = -1; }, q => { q.normalSubdivisions.pop(); },
    q => { q.normalSubdivisions[0] = undefined; q.child.tubes[0] = 10; },
    q => { delete q.normalSubdivisions[0]; q.child.tubes[0] = 10; },
    q => { q.normalSubdivisions[0].pop(); }, q => { q.child.nx++; },
    q => { q.child.tubes[0]++; }, q => { q.child.bodies[0].leadingIndex++; },
    q => { q.child.bodies[0].trailingIndex++; }, q => { q.parent.bodies[0].leadingIndex = 0; },
    q => { delete q.parent.bodies[0]; }, q => { delete q.parent.tubes[0]; },
  ]) {
    const request = structuredClone(valid); mutate(request); const before = structuredClone(request);
    assert.throws(() => refineStreamtubeEquationRegions(request)); assert.deepEqual(request, before);
  }
  assert.throws(() => refineStreamtubeEquationRegions(input(parent, Array(20).fill(1), [Array(5).fill(1), Array(5).fill(1)])), /subdivided/);
});

test('reject stale inherited topology or a region detached from its associated body cut', () => {
  const parent = { nx: 20, tubes: [5, 5], bodies: [{ leadingIndex: 3, trailingIndex: 12 }] };
  const first = refineStreamtubeEquationRegions(input(parent, Array(20).fill(2), [Array(5).fill(1), Array(5).fill(1)]));
  const valid = input(first.topology, Array(40).fill(2), [Array(5).fill(1), Array(5).fill(1)], { entropyRegions: first });
  for (const mutate of [
    r => { r.topology.nx++; }, r => { r.topology.tubes[0]++; }, r => { r.topology.bodies[0].leadingIndex++; },
    r => { r.regions[0].body = 1; }, r => { r.regions[0].throughRow = 2; },
    r => { r.regions[0].throughRow = 40; }, r => { r.regions[0].lowerTube = 5; },
    r => { r.regions[0].upperTube = 5; }, r => { r.regions[0].lowerTube = -1; },
    r => { r.parentCheckpointSha256 = 'bad'; }, r => { r.regions.push({ ...r.regions[0] }); },
  ]) {
    const request = structuredClone(valid); mutate(request.entropyRegions);
    assert.throws(() => refineStreamtubeEquationRegions(request));
  }
});
