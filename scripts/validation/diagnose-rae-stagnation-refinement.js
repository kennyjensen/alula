// SPDX-License-Identifier: GPL-2.0-or-later
// Exact local cold-grid discriminator. Prints evidence; it does not certify
// an Euler solution or change the app's interpolation policy.
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { refineStreamtubeMassCoordinates } from '../../src/geometry/streamtube-mass-interpolation.js';
import { streamtubeCellGeometry } from '../../src/euler/streamtube-cell.js';
import { streamtubeGridConvexity } from '../../src/geometry/streamtube-convex-step.js';

const source = JSON.parse(fs.readFileSync(new URL('../../docs/solver-reliability/rae-inviscid-audit/fixtures/rae128-cold-stagnation-parent.json', import.meta.url)));
const before = structuredClone(source), { grid, masses, counts, lowerStagnation } = source;
const childMasses = masses.flatMap((m, j) => Array(counts[j]).fill(m / counts[j]));
const results = [];
for (const method of ['quadratic-with-stagnation', 'quadratic-fit', 'linear']) {
  const nodes = method === 'linear' ? grid.map(row => [row[0], ...counts.flatMap((n, j) =>
    Array.from({ length: n }, (_, k) => {
      const a = row[j], b = row[j + 1], t = (k + 1) / n;
      return t === 1 ? b : { x: a.x + t * (b.x - a.x), y: a.y + t * (b.y - a.y) };
    }))]) : refineStreamtubeMassCoordinates(grid, masses, counts,
    method === 'quadratic-with-stagnation' ? { lowerStagnation } : {}).nodes;
  let maximumMassFlux = 0, bottleneck;
  for (let i = 1; i < nodes.length - 1; i++) for (let j = 0; j < childMasses.length; j++) {
    const geometry = streamtubeCellGeometry([i - 1, i, i + 1].map(k => nodes[k][j]),
      [i - 1, i, i + 1].map(k => nodes[k][j + 1]));
    geometry.normalAreas.forEach((area, k) => {
      const flux = childMasses[j] / area;
      if (flux > maximumMassFlux) { maximumMassFlux = flux; bottleneck = { i: i - 1 + k, tube: j }; }
    });
  }
  nodes.forEach((row, i) => {
    let column = 0;
    assert.deepEqual(row[column], grid[i][0]);
    counts.forEach((n, j) => { column += n; assert.deepEqual(row[column], grid[i][j + 1]); });
  });
  assert.deepEqual(source, before);
  results.push({ method, maximumMassFlux, bottleneck, quality: streamtubeGridConvexity([nodes]),
    parentNodesUnchanged: true, childMasses });
}
console.log(JSON.stringify({ source: source.source, results }, null, 2));
