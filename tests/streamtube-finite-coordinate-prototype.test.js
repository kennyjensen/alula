// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { invertSampledStreamtubeCoordinate } from '../scripts/validation/streamtube-finite-coordinate-prototype.js';
import { streamtubeGridConvexity } from '../src/geometry/streamtube-convex-step.js';

const fixture = () => ({ nodes: [0, 1, 1.01, 3.6].map(x => [{ x, y: 0 }, { x, y: 1 }]),
  baseline: [0, 1, 2, 3], corrections: [[0, 0], [1, 0], [1, 0], [0, 0]] });

test('finite inverse preserves ordered nonuniform streamlines when the centered movement reverses them', () => {
  const { nodes, baseline, corrections } = fixture(), before = structuredClone(nodes), scale = .1;
  const listed = structuredClone(nodes);
  for (let i = 1; i < 3; i++) {
    const span = baseline[i + 1] - baseline[i - 1] + scale * (corrections[i + 1][0] - corrections[i - 1][0]);
    listed[i][0].x -= scale * corrections[i][0] * (nodes[i + 1][0].x - nodes[i - 1][0].x) / span;
  }
  assert.equal(streamtubeGridConvexity([listed]).valid, false);
  const { nodes: mapped, labelsByRow } = invertSampledStreamtubeCoordinate(nodes, baseline, corrections, scale);
  assert(labelsByRow.every(row => row.every((v, i) => !i || v > row[i - 1])));
  assert(mapped.every((row, i) => !i || row[0].x > mapped[i - 1][0].x));
  assert.equal(streamtubeGridConvexity([mapped]).valid, true);
  assert.deepEqual(mapped[0], nodes[0]); assert.deepEqual(mapped.at(-1), nodes.at(-1));
  assert.deepEqual(mapped.map(row => row[1]), nodes.map(row => row[1]));
  assert.deepEqual(nodes, before);
});

test('zero correction and zero scale retain exact coordinates and return detached state', () => {
  const { nodes, baseline, corrections } = fixture(); nodes[0][0].x = -0;
  for (const [c, scale] of [[corrections, 0], [corrections.map(row => row.map(() => 0)), 1]]) {
    const mapped = invertSampledStreamtubeCoordinate(nodes, baseline, c, scale).nodes;
    assert.deepEqual(mapped, nodes); mapped[1][0].x = 100; assert.equal(nodes[1][0].x, 1);
  }
});

test('nonmonotone corrected coordinates, changing endpoints and nonfinite inputs remain failures', () => {
  const { nodes, baseline, corrections } = fixture();
  assert.throws(() => invertSampledStreamtubeCoordinate(nodes, baseline, corrections, 1),
    e => e.stage === 'coordinate monotonicity' && e.details.i === 3);
  const shifted = structuredClone(corrections); shifted[0][0] = .01;
  assert.throws(() => invertSampledStreamtubeCoordinate(nodes, baseline, shifted, .1), /Endpoints/);
  const invalid = structuredClone(corrections); invalid[1][0] = NaN;
  assert.throws(() => invertSampledStreamtubeCoordinate(nodes, baseline, invalid, .1), /Nonfinite/);
  assert.throws(() => invertSampledStreamtubeCoordinate(nodes, [0, 1, 1, 3], corrections, .1), /baseline/);
});

test('finite inverse has the documented one-sided small-step derivative', () => {
  const { nodes, baseline, corrections } = fixture(), h = 1e-7;
  const mapped = invertSampledStreamtubeCoordinate(nodes, baseline, corrections, h).nodes;
  const inverseDerivative = -(nodes[2][0].x - nodes[1][0].x) / (baseline[2] - baseline[1]);
  const listedDerivative = -(nodes[3][0].x - nodes[1][0].x) / (baseline[3] - baseline[1]);
  const fd = (mapped[2][0].x - nodes[2][0].x) / h;
  assert(Math.abs(fd - inverseDerivative) < 2e-8);
  assert(Math.abs(fd - listedDerivative) > 1);
});
