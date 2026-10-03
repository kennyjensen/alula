// SPDX-License-Identifier: GPL-2.0-or-later
// Explicit inherited ISMOM3 region control. This maps logical/mass intervals,
// not Cartesian geometry, and does not change any default refinement policy.
import { nestedIntervalMap } from '../../geometry/nested-intervals.js';
import { ISMOM3_LEADING_REGION, normalizeStreamtubeEquationRegions,
  validateStreamtubeEquationRegionTopology } from '../streamtube-equation-selection.js';

function topology(value) {
  if (!value || !Number.isSafeInteger(value.nx) || value.nx < 4
    || !Array.isArray(value.bodies) || value.bodies.length === 0
    || !Array.isArray(value.tubes) || value.tubes.length !== value.bodies.length + 1
    || Array.from(value.tubes).some(n => !Number.isSafeInteger(n) || n < 1)
    || !Number.isSafeInteger(value.tubes.reduce((sum, n) => sum + n, 0))
    || Array.from(value.bodies).some(b => !b || !Number.isSafeInteger(b.leadingIndex) || !Number.isSafeInteger(b.trailingIndex)
      || b.leadingIndex < 1 || b.leadingIndex >= b.trailingIndex || b.trailingIndex >= value.nx))
    throw new Error('Invalid equation-region refinement topology.');
  return { nx: value.nx, tubes: [...value.tubes], bodies: value.bodies.map(({ leadingIndex, trailingIndex }) =>
    ({ leadingIndex, trailingIndex })) };
}

export function refineStreamtubeEquationRegions({ parent, child, entropyRegions,
  streamwiseSubdivisions, normalSubdivisions, parentCheckpointSha256 } = {}) {
  if (typeof parentCheckpointSha256 !== 'string' || !/^[0-9a-f]{64}$/.test(parentCheckpointSha256))
    throw new Error('Equation-region refinement requires the parent checkpoint SHA-256.');
  const before = topology(parent), after = topology(child);
  if (before.bodies.length !== after.bodies.length || !Array.isArray(streamwiseSubdivisions)
    || !Array.isArray(normalSubdivisions) || normalSubdivisions.length !== before.tubes.length
    || Array.from(normalSubdivisions).some(counts => !Array.isArray(counts)))
    throw new Error('Supply the actual nested equation-region subdivision maps.');
  const streamwise = nestedIntervalMap(before.nx, { subdivisions: streamwiseSubdivisions });
  const normal = normalSubdivisions.map((counts, g) => nestedIntervalMap(before.tubes[g], { subdivisions: counts }));
  if (streamwise.retained.at(-1) !== after.nx
    || normal.some((map, g) => map.retained.at(-1) !== after.tubes[g])
    || before.bodies.some((body, b) => streamwise.retained[body.leadingIndex] !== after.bodies[b].leadingIndex
      || streamwise.retained[body.trailingIndex] !== after.bodies[b].trailingIndex))
    throw new Error('Equation-region maps do not describe the child topology and retained body stations.');
  if (streamwise.counts.every(n => n === 1) && normal.every(map => map.counts.every(n => n === 1)))
    throw new Error('Equation-region refinement requires at least one subdivided interval.');

  let source;
  if (entropyRegions === undefined) {
    const total = before.tubes.reduce((sum, n) => sum + n, 0);
    let cut = 0;
    source = normalizeStreamtubeEquationRegions({ version: 1, parentCheckpointSha256, topology: before,
      regions: before.bodies.map((body, b) => {
        cut += before.tubes[b];
        return { body: b, throughRow: Math.min(before.nx - 1, body.leadingIndex + ISMOM3_LEADING_REGION.downstreamCells),
          lowerTube: Math.max(0, cut - ISMOM3_LEADING_REGION.transverseCells),
          upperTube: Math.min(total, cut + ISMOM3_LEADING_REGION.transverseCells) };
      }) });
  } else source = normalizeStreamtubeEquationRegions(entropyRegions);
  validateStreamtubeEquationRegionTopology(source, before);

  // Build the retained global tube-boundary map from EVERY passage. A body
  // region can span several narrow passages; adjacency-only mapping is wrong.
  const global = [];
  let oldOffset = 0, newOffset = 0;
  for (let g = 0; g < before.tubes.length; g++) {
    normal[g].retained.forEach((boundary, j) => { global[oldOffset + j] = newOffset + boundary; });
    oldOffset += before.tubes[g]; newOffset += after.tubes[g];
  }
  const result = normalizeStreamtubeEquationRegions({ version: 1, parentCheckpointSha256, topology: after,
    regions: source.regions.map(({ body, throughRow, lowerTube, upperTube }) => ({ body,
      // Map the INCLUSIVE final parent row to its retained child cross-line.
      // Never regenerate LE+10 or the four-slot width on the child grid.
      throughRow: streamwise.retained[throughRow], lowerTube: global[lowerTube], upperTube: global[upperTube] })) });
  validateStreamtubeEquationRegionTopology(result, after);
  return result;
}
