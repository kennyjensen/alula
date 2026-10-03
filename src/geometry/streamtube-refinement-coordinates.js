// SPDX-License-Identifier: GPL-2.0-or-later
import { createMonotoneCubicMap } from '../numerics/monotone-cubic.js';

const dense = row => Array.isArray(row)
  && Array.from({ length: row.length }, (_, i) => Object.hasOwn(row, i)).every(Boolean);
const invalid = detail => { throw new Error(`Invalid surface refinement coordinates: ${detail}.`); };

// This explicit refinement correspondence changes inserted surface spacing only.
// Surface arrays span child LE..TE inclusive: array index is
// childI - stations.retained[body.leadingIndex]. Their coordinates are GLOBAL
// parent station indices, not local surface indices. Group node coordinates are
// [group][childI][childJ], with j increasing from the lower to the upper bank.
// No geometry, gas state, or governing equation is evaluated here.
export function createStreamtubeRefinementCoordinates({ stations, bodies, fractions, weights, tubes } = {}) {
  if (!stations || !dense(stations.counts) || !stations.counts.length
    || stations.counts.some(n => !Number.isInteger(n) || n < 1 || n > 4)) invalid('station subdivisions');
  const { counts, retained, coordinates } = stations, nx = counts.length;
  if (!dense(retained) || retained.length !== nx + 1 || retained[0] !== 0
    || retained.some((v, i) => !Number.isInteger(v) || i && v !== retained[i - 1] + counts[i - 1])
    || !dense(coordinates) || coordinates.length !== retained[nx] + 1) invalid('nested station map');
  const uniform = counts.every(n => n === counts[0]) ? counts[0] : null;
  const parentAt = Array(coordinates.length), fractionAt = Array(coordinates.length);
  for (let parent = 0; parent < nx; parent++) for (let k = 0; k < counts[parent]; k++) {
    const i = retained[parent] + k;
    const expected = uniform === null ? parent + k / counts[parent] : i / uniform;
    if (coordinates[i] !== expected) invalid('station coordinates do not match subdivisions');
    parentAt[i] = parent; fractionAt[i] = k / counts[parent];
  }
  if (coordinates.at(-1) !== nx) invalid('final station coordinate');
  parentAt[parentAt.length - 1] = nx; fractionAt[fractionAt.length - 1] = 0;
  if (!dense(bodies) || bodies.some(b => !b || !Number.isInteger(b.leadingIndex)
    || !Number.isInteger(b.trailingIndex) || b.leadingIndex < 0 || b.trailingIndex > nx
    || b.leadingIndex >= b.trailingIndex)) invalid('body station bounds');
  if (!dense(fractions) || fractions.length !== bodies.length) invalid('surface fraction count');
  for (const [b, body] of bodies.entries()) for (const side of ['upper', 'lower']) {
    const row = fractions[b]?.[side];
    if (!dense(row) || row.length !== body.trailingIndex - body.leadingIndex + 1
      || row[0] !== 0 || row.at(-1) !== 1
      || row.some((v, i) => !Number.isFinite(v) || i && v <= row[i - 1])) invalid('surface fractions');
  }
  if (!dense(weights) || weights.length !== bodies.length + 1
    || weights.some(row => !dense(row) || !row.length
      || row.some(v => !Number.isFinite(v) || !(v > 0)))) invalid('positive target passage weights');
  if (tubes !== undefined && (!dense(tubes) || tubes.length !== weights.length
    || tubes.some(n => !Number.isInteger(n) || n < 1))) invalid('optional parent tube counts');
  const eta = weights.map(row => {
    const total = row.reduce((sum, value) => sum + value, 0);
    if (!Number.isFinite(total)) invalid('unresolved target weight sum');
    const result = [0]; let partial = 0;
    for (let j = 0; j < row.length - 1; j++) {
      partial += row[j]; const value = partial / total;
      if (!(value > result.at(-1) && value < 1)) invalid('unresolved target weight fraction');
      result.push(value);
    }
    result.push(1); return result;
  });
  const surfaceFractions = [], surfaceCoordinates = [];
  for (const [b, body] of bodies.entries()) {
    const leading = body.leadingIndex, trailing = body.trailingIndex;
    const first = retained[leading], last = retained[trailing];
    const nextFractions = {}, nextCoordinates = {};
    for (const side of ['upper', 'lower']) {
      const old = fractions[b][side];
      const map = createMonotoneCubicMap(old.map((_, k) => k), old, { derivatives: 'pchip' });
      const values = [], sources = [];
      for (let i = first; i <= last; i++) {
        const parent = parentAt[i], k = parent - leading;
        if (fractionAt[i] === 0) {
          // Bypass all interpolation at a retained station, including LE/TE.
          values.push(old[k]); sources.push(parent); continue;
        }
        const local = uniform === null ? k + fractionAt[i] : (i - first) / uniform;
        const value = map.value(local);
        const source = parent + (value - old[k]) / (old[k + 1] - old[k]);
        if (!(value > old[k] && value < old[k + 1] && source > parent && source < parent + 1)
          || !Number.isFinite(source) || !(value > values.at(-1) && source > sources.at(-1)))
          invalid('unresolved inserted surface coordinate');
        values.push(value); sources.push(source);
      }
      nextFractions[side] = values; nextCoordinates[side] = sources;
    }
    surfaceFractions.push(nextFractions); surfaceCoordinates.push(nextCoordinates);
  }
  const bankCoordinate = (b, side, i) => {
    const body = bodies[b];
    if (!body || i < retained[body.leadingIndex] || i > retained[body.trailingIndex]) return coordinates[i];
    return surfaceCoordinates[b][side][i - retained[body.leadingIndex]];
  };
  const nodeCoordinatesByGroup = eta.map((row, g) => coordinates.map((ordinary, i) => {
    if (fractionAt[i] === 0) return row.map(() => ordinary);
    const lower = bankCoordinate(g - 1, 'upper', i), upper = bankCoordinate(g, 'lower', i);
    return row.map((value, j) => {
      const source = j === 0 ? lower : j === row.length - 1 ? upper
        : lower === upper ? lower : lower + (upper - lower) * value;
      if (!(source > parentAt[i] && source < parentAt[i] + 1)) invalid('unresolved blended node coordinate');
      return source;
    });
  }));
  for (const group of nodeCoordinatesByGroup) for (let i = 1; i < group.length; i++)
    if (group[i].some((value, j) => !(value > group[i - 1][j]))) invalid('non-increasing node rail');
  return { surfaceFractions, surfaceCoordinates, nodeCoordinatesByGroup };
}
