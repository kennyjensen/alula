// SPDX-License-Identifier: GPL-2.0-or-later
// Apply the passage-local spacing reconstruction to exact wall curves and
// shared panel-flow cuts. Farfield stations and actual block hits stay fixed.
import { createPassageSpacingComponents } from '../geometry/passage-spacing-components.js';
import { reconcilePassageStations } from '../geometry/passage-density-stations.js';
import { createSpacingEnvelope } from '../geometry/spacing-envelope.js';
import { planPassageCounts } from '../geometry/passage-count-planner.js';
import { streamtubeCutReferencePoints } from './streamtube-cut-reference.js';

export function matchPassageDensityGuides({ input, profiles, guides, outer, blocks, ranks, surfaceMaps, resolvedDemands, sampleCutX, countPlanning }) {
  const referenceGuides = structuredClone(guides);
  const components = createPassageSpacingComponents(blocks).map(c => {
    const firstIndex = ranks.indexOf(c.start), lastIndex = ranks.indexOf(c.end);
    if (firstIndex < 0 || lastIndex <= firstIndex) throw new Error('Missing passage density block endpoint.');
    return { ...c, firstIndex, lastIndex, intervals: lastIndex - firstIndex, boundaries: c.boundaries.map(b => {
      if (b.kind === 'wall') {
        const descriptor = surfaceMaps[b.body][b.side], a = descriptor.rows.find(r => r.rank === c.start), z = descriptor.rows.find(r => r.rank === c.end);
        if (!a || !z) throw new Error('A wall spacing block must end at actual wall hits.');
        return { ...b, length: z.value - a.value, range: [a.value, z.value], requested: resolvedDemands[b.body][b.side].map(descriptor.at) };
      }
      const points = b.kind === 'cut'
        ? streamtubeCutReferencePoints(input.bodies[b.body], referenceGuides[b.body], b.end, firstIndex, lastIndex)
        : outer[b.side === 'lower' ? 0 : 1].slice(firstIndex, lastIndex + 1);
      const positions = [0];
      for (let i = 1; i < points.length; i++) positions.push(positions.at(-1) + Math.hypot(points[i].x - points[i - 1].x, points[i].y - points[i - 1].y));
      const length = positions.at(-1);
      return { ...b, length, points, positions, progress: positions.map(s => s / length) };
    }) };
  });
  const planned = countPlanning ? planPassageCounts({ components, ...countPlanning }) : null;
  const countPlan = planned ? { counts: planned.counts, indices: planned.indices, total: planned.total, history: planned.history,
    allOriginalResolutionSatisfied: planned.allOriginalResolutionSatisfied, globalMinimumProven: planned.globalMinimumProven,
    minimumTotalForBounds: planned.minimumTotalForBounds, tieBreak: planned.tieBreak, scope: planned.scope } : null;
  if (planned && planned.counts.some((n, i) => n !== countPlanning.minimumCounts[i])) return { rebuildRequired: true, countPlan };
  const fitted = planned?.fitted ?? reconcilePassageStations({ components, endpointMethod: 'growth-constrained' });
  const achievedDensity = components.map((c, k) => {
    const walls = c.boundaries.filter(b => b.kind === 'wall');
    if (!walls.length) return { referenceOnly: true };
    const envelope = createSpacingEnvelope({ profiles: walls.map(b => b.requested.map(s => (s - b.range[0]) / b.length)), start: 0, end: 1 });
    const q = fitted.fitted[k].progress.map(envelope.metric);
    const maximumDensityInterval = Math.max(...q.slice(1).map((v, i) => v - q[i]));
    return { maximumDensityInterval, resolutionSatisfied: maximumDensityInterval <= 1 + 1e-10 };
  });
  for (let k = 0; k < components.length; k++) {
    const c = components[k], progress = fitted.fitted[k].progress;
    for (const b of c.boundaries) {
      if (b.kind === 'outer') continue;
      const profile = profiles[b.body], body = input.bodies[b.body];
      if (b.kind === 'wall') {
        const descriptor = surfaceMaps[b.body][b.side];
        progress.forEach((u, j) => {
          const i = c.firstIndex + j;
          // Exact physical hit samples from the base placement take precedence.
          if (!j || j === progress.length - 1) return;
          const fraction = descriptor.fractionAtPosition(b.range[0] + b.length * u);
          body.surfaceFractions[b.side][i - body.leadingIndex] = fraction;
          const v = profile.curve.branch(b.side, fraction, profile.stag);
          guides[b.body][b.side][i] = { ...v.point, potential: profile.phase(v.parameter) };
        });
      } else {
        const jump = profile.phase(0) - profile.phase(profile.curve.length), wake = b.end === 'wake';
        progress.forEach((u, j) => {
          if (!j || j === progress.length - 1) return;
          const s = b.length * u, m = b.positions.findIndex(p => p > s) - 1;
          if (m < 0 || m >= b.points.length - 1) throw new Error('Cut arc station is outside its reference segment.');
          const t = (s - b.positions[m]) / (b.positions[m + 1] - b.positions[m]);
          const xx = b.points[m].x + t * (b.points[m + 1].x - b.points[m].x);
          let point;
          try { point = sampleCutX(profile[b.end], xx); }
          catch (error) { error.guideRequest = { body: b.body, end: b.end, x: xx, component: k, station: c.firstIndex + j,
            progress: u, referencePoints: b.points, positions: b.positions }; throw error; }
          const i = c.firstIndex + j;
          guides[b.body].lower[i] = { ...point, potential: point.potential + (wake ? profile.phase(profile.curve.length) : 0) };
          guides[b.body].upper[i] = { ...guides[b.body].lower[i], potential: guides[b.body].lower[i].potential + (wake ? jump : 0) };
        });
      }
    }
  }
  return { ...fitted, ...(countPlan ? { countPlan } : {}), achievedDensity, components: components.map(c => ({ ...c, boundaries: c.boundaries.map(({ points, positions, ...b }) => b) })),
    cutMetric: 'reference polygonal arc mapped to x; positions resampled on the original panel trace',
    outerStationsChanged: false, physicalHitsChanged: false, physicalAcceptance: false };
}
