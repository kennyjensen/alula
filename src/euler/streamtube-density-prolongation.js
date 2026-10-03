// SPDX-License-Identifier: GPL-2.0-or-later
import { streamtubeCellGeometry as cellGeometry } from './streamtube-cell.js';
// Physical state interpolation. Physical log-density interpolation is an initial guess, not an
// entropy equation solve. Mass and h0 are exact; entropy/branch errors reported.
export function prolongPhysicalLogDensity({ source, target, h0, gamma = 1.4 }) {
  const { logDensity, masses, normalAreas } = source ?? {}, { nodeCoordinates, nodeCoordinatesByGroup,
    subdivisions, masses: childMasses, normalAreas: childAreas } = target ?? {};
  const nx = logDensity?.length, ng = masses?.length;
  if (!(Number.isFinite(h0) && h0 > 0 && Number.isFinite(gamma) && gamma > 1)
    || !Number.isInteger(nx) || nx < 1 || !Number.isInteger(ng) || ng < 1
    || masses.some(row => !Array.isArray(row) || !row.length || row.some(m => !Number.isFinite(m) || m <= 0)))
    throw new Error('Invalid physical source density conditions.');
  const shape = (rows, counts) => Array.isArray(rows) && rows.every(row => Array.isArray(row) && row.length === ng
    && row.every((group, g) => Array.isArray(group) && group.length === counts[g] && group.every(Number.isFinite)));
  const counts = masses.map(row => row.length);
  if (!shape(logDensity, counts) || !shape(normalAreas, counts) || normalAreas.length !== nx)
    throw new Error('Invalid source section arrays.');
  if (!Array.isArray(nodeCoordinates) || nodeCoordinates.length < 2 || nodeCoordinates[0] !== 0 || nodeCoordinates.at(-1) !== nx
    || nodeCoordinates.some((x, i) => !Number.isFinite(x) || i > 0 && x <= nodeCoordinates[i - 1])
    || !Array.isArray(subdivisions) || subdivisions.length !== ng || subdivisions.some((row, g) => !Array.isArray(row)
      || row.length !== counts[g] || row.some(n => !Number.isInteger(n) || n < 1)))
    throw new Error('Invalid density section refinement map.');
  const childCounts = subdivisions.map(row => row.reduce((a, b) => a + b, 0));
  if (nodeCoordinatesByGroup !== undefined) {
    if (!Array.isArray(nodeCoordinatesByGroup) || nodeCoordinatesByGroup.length !== ng)
      throw new Error('Invalid per-node density source-coordinate field.');
    for (let g = 0; g < ng; g++) {
      const grid = nodeCoordinatesByGroup[g];
      if (!Array.isArray(grid) || grid.length !== nodeCoordinates.length)
        throw new Error('Invalid per-node density source-coordinate field.');
      for (let i = 0; i < grid.length; i++) {
        const row = grid[i];
        if (!Array.isArray(row) || row.length !== childCounts[g] + 1)
          throw new Error('Invalid per-node density source-coordinate field.');
        for (let j = 0; j < row.length; j++) {
          const s = row[j];
          if (!Number.isFinite(s) || s < 0 || s > nx || i === 0 && s !== 0
            || i === grid.length - 1 && s !== nx || i > 0 && !(s > grid[i - 1][j]))
            throw new Error('Invalid per-node density source-coordinate field.');
        }
      }
    }
  }
  if (!Array.isArray(childMasses) || childMasses.length !== ng || childMasses.some((row, g) => !Array.isArray(row)
    || row.length !== childCounts[g] || row.some(m => !Number.isFinite(m) || m <= 0))
    || !shape(childAreas, childCounts) || childAreas.length !== nodeCoordinates.length - 1)
    throw new Error('Invalid target density section arrays.');
  let maximumPartitionError = 0;
  const centers = rows => rows.map(row => { let left = 0; return row.map(m => { const right = left + m, center = (left + right) / 2;
    if (!(right > left)) throw new Error('Unresolved mass coordinate.'); left = right; return center; }); });
  for (let g = 0; g < ng; g++) {
    let k = 0;
    for (let j = 0; j < masses[g].length; j++) {
      const sum = childMasses[g].slice(k, k + subdivisions[g][j]).reduce((a, b) => a + b, 0);
      const error = Math.abs(sum - masses[g][j]); maximumPartitionError = Math.max(maximumPartitionError, error);
      if (error > 64 * Number.EPSILON * Math.max(sum, masses[g][j])) throw new Error('Refinement changed a parent tube mass.');
      k += subdivisions[g][j];
    }
  }
  const parentMassCenters = centers(masses), childMassCenters = centers(childMasses);
  const parentStreamwise = Array.from({ length: nx }, (_, i) => i + .5);
  const childStreamwise = nodeCoordinatesByGroup === undefined
    ? nodeCoordinates.slice(1).map((x, i) => (x + nodeCoordinates[i]) / 2)
    : nodeCoordinates.slice(1).map((_, i) => nodeCoordinatesByGroup.map((grid, g) =>
      Array.from({ length: childCounts[g] }, (_, j) => {
        const a = grid[i][j], b = grid[i][j + 1], c = grid[i + 1][j], d = grid[i + 1][j + 1];
        // Equal coordinates on each station plane retain the original
        // (next + previous)/2 arithmetic, including nonbinary values.
        return a === b && c === d ? (c + a) / 2 : (a + b + c + d) / 4;
      })));
  const bracket = (positions, x) => {
    if (x <= positions[0]) return { a: 0, b: 0, fraction: 0 };
    if (x >= positions.at(-1)) return { a: positions.length - 1, b: positions.length - 1, fraction: 0 };
    let a = 0; while (positions[a + 1] < x) a++;
    if (positions[a] === x) return { a, b: a, fraction: 0 };
    if (positions[a + 1] === x) return { a: a + 1, b: a + 1, fraction: 0 };
    return { a, b: a + 1, fraction: (x - positions[a]) / (positions[a + 1] - positions[a]) };
  };
  const mix = (a, b, t) => t === 0 || a === b ? a : t === 1 ? b : a + t * (b - a);
  const gas = (logRho, mass, area) => {
    const rho = Math.exp(logRho), q = mass / (rho * area), h = h0 - q * q / 2, p = (gamma - 1) / gamma * rho * h;
    if (![rho, q, h, p, area].every(v => Number.isFinite(v) && v > 0)) throw new Error('Transferred physical state is thermally or geometrically inadmissible.');
    return { rho, q, h, p, machSquared: q * q / ((gamma - 1) * h), entropyCoordinate: Math.log(h) / (gamma - 1) - logRho };
  };
  const sourceGas = logDensity.map((row, i) => row.map((group, g) => group.map((v, j) => gas(v, masses[g][j], normalAreas[i][g][j]))));
  const interpolate = (field, sb, mb, g) => mix(mix(field[sb.a][g][mb.a], field[sb.a][g][mb.b], mb.fraction),
    mix(field[sb.b][g][mb.a], field[sb.b][g][mb.b], mb.fraction), sb.fraction);
  const sourceEntropy = sourceGas.map(row => row.map(group => group.map(s => s.entropyCoordinate)));
  const sourceMach = sourceGas.map(row => row.map(group => group.map(s => s.machSquared)));
  const diagnostics = { method: 'physical-log-density-on-section-and-mass-centers', initialGuessOnly: true,
    commonIsentropeInversions: 0, maximumPartitionError, maximumEntropyInterpolationDeparture: 0,
    minimumStaticEnthalpy: Infinity, minimumPressure: Infinity, minimumDensity: Infinity, maximumMach: 0,
    interpolatedMachBranchChanges: 0, collocatedValues: 0, collocatedValuesExact: true,
    endpointExtension: 'constant outside the first/last parent section or mass center',
    branchMeaning: 'Compare with bounded interpolated parent Mach squared; an initial-guess diagnostic, not a physical phase constraint.' };
  const result = childStreamwise.map((coordinate, i) => masses.map((_, g) => childMassCenters[g].map((m, j) => {
    const s = nodeCoordinatesByGroup === undefined ? coordinate : coordinate[g][j];
    const sb = bracket(parentStreamwise, s), mb = bracket(parentMassCenters[g], m), value = interpolate(logDensity, sb, mb, g);
    const state = gas(value, childMasses[g][j], childAreas[i][g][j]);
    const entropyReference = interpolate(sourceEntropy, sb, mb, g), machReference = interpolate(sourceMach, sb, mb, g);
    diagnostics.maximumEntropyInterpolationDeparture = Math.max(diagnostics.maximumEntropyInterpolationDeparture, Math.abs(state.entropyCoordinate - entropyReference));
    diagnostics.minimumStaticEnthalpy = Math.min(diagnostics.minimumStaticEnthalpy, state.h);
    diagnostics.minimumPressure = Math.min(diagnostics.minimumPressure, state.p);
    diagnostics.minimumDensity = Math.min(diagnostics.minimumDensity, state.rho);
    diagnostics.maximumMach = Math.max(diagnostics.maximumMach, Math.sqrt(state.machSquared));
    diagnostics.interpolatedMachBranchChanges += (state.machSquared > 1) !== (machReference > 1) ? 1 : 0;
    if (sb.a === sb.b && mb.a === mb.b && s === parentStreamwise[sb.a] && m === parentMassCenters[g][mb.a]) {
      diagnostics.collocatedValues++; diagnostics.collocatedValuesExact &&= value === logDensity[sb.a][g][mb.a];
    }
    return value;
  })));
  return { logDensity: result, diagnostics, sourceGas, nodeCoordinates: nodeCoordinates.slice(),
    ...(nodeCoordinatesByGroup === undefined ? {} : { nodeCoordinatesByGroup: nodeCoordinatesByGroup.map(grid => grid.map(row => row.slice())) }),
    sourceSectionCoordinates: parentStreamwise, targetSectionCoordinates: childStreamwise,
    sourceMassCenters: parentMassCenters, targetMassCenters: childMassCenters };
}

// Layout adapter: no source or target residual evaluation. The source flow has
// already been evaluated by the refiner; target normal areas are geometry only.
export function prolongStreamtubeDensities(source, sourceState, sourceFlow, target, targetState, { nodeCoordinates, nodeCoordinatesByGroup, subdivisions }) {
  for (const key of ['mach','gamma','h0','lengthScale','massScale'])
    if (source.conditions[key] !== target.conditions[key]) throw new Error(`Density refinement changed ${key}.`);
  const physical = target.decode(targetState), { nx, tubes } = target.layout;
  const areas = Array.from({ length: nx }, () => tubes.map(n => Array(n)));
  for (let g = 0; g < tubes.length; g++) for (let i = 1; i < nx; i++) for (let j = 0; j < tubes[g]; j++) {
    const geometry = cellGeometry([physical.nodes[g][i-1][j],physical.nodes[g][i][j],physical.nodes[g][i+1][j]],
      [physical.nodes[g][i-1][j+1],physical.nodes[g][i][j+1],physical.nodes[g][i+1][j+1]]);
    if (i === 1) areas[0][g][j] = geometry.normalAreas[0]; areas[i][g][j] = geometry.normalAreas[1];
  }
  const old = source.layout;
  const result = prolongPhysicalLogDensity({ h0: target.conditions.h0, gamma: target.conditions.gamma,
    source: {
      logDensity: Array.from({ length: old.nx }, (_, i) => old.tubes.map((n,g) => Array.from({length:n},(_,j) => sourceState[old.densityIndex(i,g,j)]))),
      masses: sourceFlow.allocation.groups.map(row => row.map(t => t.massFlow)),
      normalAreas: Array.from({ length: old.nx }, (_, i) => old.tubes.map((n,g) => Array.from({length:n},(_,j) => sourceFlow.cells[Math.max(0,i-1)][g][j].geometry.normalAreas[i === 0 ? 0 : 1]))),
    }, target: { nodeCoordinates, ...(nodeCoordinatesByGroup === undefined ? {} : { nodeCoordinatesByGroup }),
      subdivisions, masses: physical.allocation.groups.map(row => row.map(t => t.massFlow)), normalAreas: areas } });
  const state = targetState.slice();
  result.logDensity.forEach((row,i) => row.forEach((group,g) => group.forEach((value,j) => { state[target.layout.densityIndex(i,g,j)] = value; })));
  return { state, diagnostics: result.diagnostics };
}
