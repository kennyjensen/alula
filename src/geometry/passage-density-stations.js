// SPDX-License-Identifier: GPL-2.0-or-later
// Common normalized physical-arc progression in each connected fluid block,
// with simultaneous endpoint matching along body outlines. Original wall
// requests supply density; no-wall blocks retain a reference cut progression.
import { distributeFacingDensity } from './facing-density-stations.js';
import { distributeSurfaceByDensity } from './surface-density-stations.js';
import { fitLogRatios } from '../numerics/log-ratio-fit.js';
import { fitConstrainedStations } from '../numerics/constrained-stations.js';
import { fitConstrainedStationPositions } from '../numerics/constrained-station-positions.js';
import { createSpacingEnvelope } from './spacing-envelope.js';

export function passageStationTargets({ components, retainEndpointResolution = false }) {
  if (typeof retainEndpointResolution !== 'boolean') throw new Error('Endpoint resolution retention must be boolean.');
  if (!Array.isArray(components) || !components.length || components.some(c => !c || !(c.start < c.end)
    || !Number.isInteger(c.intervals) || c.intervals < 3 || !Array.isArray(c.boundaries)
    || c.boundaries.some(b => !['outer', 'wall', 'cut'].includes(b.kind) || !Number.isFinite(b.length) || !(b.length > 0))))
    throw new Error('Passage station fitting requires ordered blocks, positive lengths and at least three intervals.');
  const natural = components.map(c => {
    const walls = c.boundaries.filter(b => b.kind === 'wall');
    if (walls.length) return distributeFacingDensity({ profiles: walls.map(b => b.requested), ranges: walls.map(b => b.range), intervals: c.intervals });
    const cuts = c.boundaries.filter(b => b.kind === 'cut');
    if (!cuts.length || cuts.some(b => !Array.isArray(b.progress) || b.progress.length !== c.intervals + 1
      || b.progress[0] !== 0 || b.progress.at(-1) !== 1 || b.progress.some((u, i) => !Number.isFinite(u) || i && !(u > b.progress[i - 1]))))
      throw new Error('No-wall blocks require ordered reference progress on their shared cuts.');
    return { progress: cuts[0].progress.map((_, i) => cuts.reduce((sum, b) => sum + b.progress[i], 0) / cuts.length),
      intervals: c.intervals, referenceOnly: true, countPolicy: 'fixed inlet/wake count; arithmetic mean of reference cut progress' };
  });
  const bodyIds = [...new Set(components.flatMap(c => c.boundaries.filter(b => b.kind !== 'outer').map(b => b.body)))];
  if (bodyIds.some(b => !Number.isInteger(b) || b < 0)) throw new Error('Invalid passage boundary body index.');
  const connections = [], seen = new Set();
  for (const body of bodyIds) for (const side of ['upper', 'lower']) {
    const outline = components.flatMap((c, component) => c.boundaries
      .filter(b => b.body === body && (b.kind === 'cut' || b.kind === 'wall' && b.side === side))
      .map(boundary => ({ component, start: c.start, end: c.end, boundary }))).sort((a, b) => a.start - b.start);
    for (let k = 1; k < outline.length; k++) {
      const a = outline[k - 1], b = outline[k];
      if (a.end !== b.start) throw new Error('Disconnected outline in passage station fitting.');
      const key = `${body}:${a.component}:${b.component}:${a.boundary.kind}:${b.boundary.kind}`;
      if (seen.has(key)) continue; seen.add(key);
      connections.push({ left: 2 * a.component + 1, right: 2 * b.component,
        ratio: a.boundary.length / b.boundary.length, body,
        side: a.boundary.kind === 'cut' && b.boundary.kind === 'cut' ? 'shared cut' : side, rank: a.end });
    }
  }
  const maximumValues = retainEndpointResolution ? components.flatMap(c => {
    const walls = c.boundaries.filter(b => b.kind === 'wall');
    if (!walls.length) return [Infinity, Infinity];
    const envelope = createSpacingEnvelope({ profiles: walls.map(b => b.requested.map(s => (s - b.range[0]) / b.length)), start: 0, end: 1 });
    return [envelope.inverse(Math.min(1, envelope.total)), 1 - envelope.inverse(Math.max(0, envelope.total - 1))];
  }) : undefined;
  const endpointFit = fitLogRatios({ natural: natural.flatMap(c => [c.progress[1], 1 - c.progress.at(-2)]), connections, maximumValues });
  return { natural, connections, endpointFit };
}

export function reconcilePassageStations({ components, endpointMethod = 'vinokur', maximumGrowth = 1.5, retainEndpointResolution = false, projectionMethod = 'dykstra', projectionObjective = 'relative-interval' }) {
  if (!['vinokur', 'growth-constrained'].includes(endpointMethod)) throw new Error('Unknown passage endpoint placement.');
  if (!['relative-interval', 'station-position'].includes(projectionObjective)) throw new Error('Unknown passage station objective.');
  const { natural, connections, endpointFit } = passageStationTargets({ components, retainEndpointResolution });
  const fitted = components.map((c, k) => {
    const endpointSpacings = { firstSpacing: endpointFit.values[2 * k], lastSpacing: endpointFit.values[2 * k + 1] };
    const walls = c.boundaries.filter(b => b.kind === 'wall');
    try {
      if (endpointMethod === 'growth-constrained') {
        const project = projectionObjective === 'station-position' ? fitConstrainedStationPositions : fitConstrainedStations;
        const projection = project({ positions: natural[k].progress, ...endpointSpacings, maximumGrowth, method: projectionMethod });
        return { progress: projection.positions, intervals: c.intervals, projection,
          originalDensityReevaluationRequired: true, endpointSpacings, referenceOnly: !walls.length };
      }
      if (walls.length) return distributeFacingDensity({ profiles: walls.map(b => b.requested), ranges: walls.map(b => b.range), intervals: c.intervals, endpointSpacings });
      const fit = distributeSurfaceByDensity({ requested: natural[k].progress,
        anchors: [{ index: 0, position: 0 }, { index: c.intervals, position: 1 }], endpointSpacings: [endpointSpacings] });
      return { progress: fit.positions, intervals: c.intervals, referenceOnly: true, endpointSpacings };
    } catch (error) { throw new Error(`Passage ${c.start} to ${c.end}: ${error.message}`, { cause: error }); }
  });
  return { natural, fitted, connections, endpointFit, endpointMethod,
    ...(projectionObjective === 'station-position' ? { projectionObjective } : {}),
    ...(retainEndpointResolution ? { retainEndpointResolution } : {}), exactModernMsetLaw: false,
    sharedProgress: true, metric: 'wall contour arc; supplied cut arc; outer boundaries excluded from matching' };
}
