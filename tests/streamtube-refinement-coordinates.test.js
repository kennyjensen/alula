// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { nestedIntervalMap } from '../src/geometry/nested-intervals.js';
import { createStreamtubeRefinementCoordinates } from '../src/geometry/streamtube-refinement-coordinates.js';

const close = (actual, expected, tolerance = 2e-15) => assert.ok(Math.abs(actual - expected) <= tolerance,
  `${actual} differs from ${expected}`);
const single = () => ({ stations: nestedIntervalMap(5), bodies: [{ leadingIndex: 1, trailingIndex: 3 }],
  fractions: [{ upper: [0, .1, 1], lower: [0, .9, 1] }], weights: [[1, 3], [3, 1]] });
const staggered = () => ({ stations: nestedIntervalMap(10, { subdivisions: [1, 3, 2, 4, 2, 1, 3, 2, 4, 1] }),
  bodies: [{ leadingIndex: 1, trailingIndex: 6 }, { leadingIndex: 4, trailingIndex: 9 }],
  fractions: [{ upper: [0, .04, .16, .38, .73, 1], lower: [0, .21, .49, .71, .91, 1] },
    { upper: [0, .16, .38, .62, .83, 1], lower: [0, .08, .18, .42, .72, 1] }],
  weights: [[1, 3], [2, 1, 4], [5, 2, 1, 1]], tubes: [1, 2, 3] });

test('PCHIP fractions and their source-index inverse share exact retained knots', () => {
  const input = single(), result = createStreamtubeRefinementCoordinates(input);
  assert.deepEqual(result.surfaceFractions[0].upper.filter((_, i) => i % 2 === 0), [0, .1, 1]);
  assert.deepEqual(result.surfaceCoordinates[0].upper.filter((_, i) => i % 2 === 0), [1, 2, 3]);
  // Three-knot increasing PCHIP has upper slopes [0, .18, 1.3].
  close(result.surfaceFractions[0].upper[1], .0275);
  close(result.surfaceFractions[0].upper[3], .41);
  close(result.surfaceCoordinates[0].upper[1], 1.275);
  close(result.surfaceCoordinates[0].upper[3], 2 + .31 / .9);
  assert.notEqual(result.surfaceCoordinates[0].upper[1], result.surfaceCoordinates[0].lower[1]);
  for (const side of ['upper', 'lower']) for (let i = 1; i < 4; i += 2) {
    const k = (i - 1) / 2, tau = result.surfaceCoordinates[0][side][i] - 1 - k;
    close(input.fractions[0][side][k] + tau * (input.fractions[0][side][k + 1] - input.fractions[0][side][k]),
      result.surfaceFractions[0][side][i]);
  }
});

test('staggered bodies retain every parent plane and keep all inserted rails inside their parent interval', () => {
  const input = staggered(), result = createStreamtubeRefinementCoordinates(input), { stations } = input;
  for (let parent = 0; parent < stations.retained.length; parent++) {
    for (const group of result.nodeCoordinatesByGroup)
      assert.ok(group[stations.retained[parent]].every(value => value === parent));
    if (parent === stations.counts.length) continue;
    for (let i = stations.retained[parent] + 1; i < stations.retained[parent + 1]; i++)
      for (const group of result.nodeCoordinatesByGroup)
        assert.ok(group[i].every(value => value > parent && value < parent + 1));
  }
  for (const [b, body] of input.bodies.entries()) for (const side of ['upper', 'lower']) {
    const origin = stations.retained[body.leadingIndex];
    assert.equal(result.surfaceFractions[b][side].length, stations.retained[body.trailingIndex] - origin + 1);
    for (let parent = body.leadingIndex; parent <= body.trailingIndex; parent++) {
      const child = stations.retained[parent] - origin;
      assert.equal(result.surfaceFractions[b][side][child], input.fractions[b][side][parent - body.leadingIndex]);
      assert.equal(result.surfaceCoordinates[b][side][child], parent);
    }
  }
});

test('interiors use cumulative target mass and independent wall maps; inactive cuts remain shared', () => {
  const input = staggered(), result = createStreamtubeRefinementCoordinates(input), { stations } = input;
  for (let i = 0; i < stations.coordinates.length; i++) {
    assert.equal(result.nodeCoordinatesByGroup[0][i][0], stations.coordinates[i]);
    assert.equal(result.nodeCoordinatesByGroup.at(-1)[i].at(-1), stations.coordinates[i]);
    for (const [b, body] of input.bodies.entries()) {
      const lowerSide = result.nodeCoordinatesByGroup[b][i].at(-1);
      const upperSide = result.nodeCoordinatesByGroup[b + 1][i][0];
      if (i < stations.retained[body.leadingIndex] || i > stations.retained[body.trailingIndex]) {
        assert.equal(lowerSide, stations.coordinates[i]); assert.equal(upperSide, lowerSide);
      } else {
        const local = i - stations.retained[body.leadingIndex];
        assert.equal(lowerSide, result.surfaceCoordinates[b].lower[local]);
        assert.equal(upperSide, result.surfaceCoordinates[b].upper[local]);
      }
    }
  }
  const child = stations.retained[4] + 1, shared = result.nodeCoordinatesByGroup[1][child];
  assert.notEqual(shared[0], shared.at(-1));
  close(shared[1], shared[0] + (shared.at(-1) - shared[0]) * 2 / 7);
  close(shared[2], shared[0] + (shared.at(-1) - shared[0]) * 3 / 7);
});

test('four-corner section coordinates are finite, increasing, and inside their parent interval', () => {
  const { stations, ...rest } = staggered();
  const result = createStreamtubeRefinementCoordinates({ stations, ...rest });
  for (const group of result.nodeCoordinatesByGroup) for (let j = 0; j < group[0].length - 1; j++) {
    let previous = -Infinity;
    for (let i = 0; i < group.length - 1; i++) {
      const center = .25 * (group[i][j] + group[i][j + 1] + group[i + 1][j] + group[i + 1][j + 1]);
      const parent = Math.floor(stations.coordinates[i]);
      assert.ok(Number.isFinite(center) && center > previous && center > parent && center < parent + 1);
      previous = center;
    }
  }
});

test('outputs do not mutate or alias input arrays, wall sides, groups, or station rows', () => {
  const input = single(), before = JSON.stringify(input), result = createStreamtubeRefinementCoordinates(input);
  result.surfaceFractions[0].upper[0] = 9;
  result.surfaceCoordinates[0].upper[0] = 8;
  result.nodeCoordinatesByGroup[0][0][0] = 7;
  assert.equal(JSON.stringify(input), before);
  assert.equal(result.surfaceFractions[0].lower[0], 0);
  assert.equal(result.surfaceCoordinates[0].lower[0], 1);
  assert.equal(result.nodeCoordinatesByGroup[1][0][0], 0);
  assert.equal(result.nodeCoordinatesByGroup[0][1][0], .5);
});

test('rejects malformed maps, sparse arrays, invalid fraction domains, and unresolved weights', () => {
  const mutations = [
    x => { delete x.stations.counts[1]; }, x => { delete x.stations.retained[1]; },
    x => { delete x.stations.coordinates[1]; }, x => { x.stations.coordinates[1] = .6; },
    x => { x.stations.retained[2]++; }, x => { x.bodies[0].trailingIndex = 6; },
    x => { delete x.bodies[0]; }, x => { delete x.fractions[0]; },
    x => { delete x.fractions[0].upper[1]; }, x => { x.fractions[0].upper[1] = 0; },
    x => { x.fractions[0].upper[2] = 1.1; }, x => { x.fractions[0].lower.pop(); },
    x => { delete x.weights[0]; }, x => { delete x.weights[0][1]; },
    x => { x.weights[0][0] = 0; }, x => { x.weights[0][0] = Infinity; },
    x => { x.weights[0] = [1, Number.MIN_VALUE]; }, x => { x.weights.pop(); },
    x => { x.tubes = [1, , 2]; },
  ];
  for (const mutate of mutations) { const input = single(); mutate(input); assert.throws(() => createStreamtubeRefinementCoordinates(input), /Invalid surface refinement coordinates/); }
});
