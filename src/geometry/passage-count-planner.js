// SPDX-License-Identifier: GPL-2.0-or-later
// Bounded boundary-only count planning. Counts increase from the supplied
// minimum schedule. Success requires measured density and growth constraints;
// no flow solve, hidden geometry repair or global nonlinear optimality claim.
import { allocateOrderedIntervalCounts } from '../numerics/ordered-interval-counts.js';
import { stationGrowthFeasibility } from '../numerics/constrained-stations.js';
import { passageStationTargets, reconcilePassageStations } from './passage-density-stations.js';
import { createSpacingEnvelope } from './spacing-envelope.js';

export function planPassageCounts(input) {
  try { return planPassageCountsAttempt(input); }
  catch (error) {
    // Retry only a measured density failure, from the original minimum counts.
    // Invalid topology, budgets, feasibility, and numerical failures retain
    // their original rejection. All successful original plans are unchanged.
    if (error.code !== 'passage-density-unresolved') throw error;
    const originalAttempt = { code: error.code, reason: error.message, history: error.history };
    try {
      const result = planPassageCountsAttempt(input, 'station-position');
      return { ...result, placementRecovery: { attempted: true, originalAttempt,
        objective: 'station-position', initialCounts: [...input.minimumCounts],
        physicalConstraintsChanged: false, stationAccuracyIndependentlyCertified: true } };
    } catch (recoveryError) {
      recoveryError.placementRecovery = { attempted: true, originalAttempt, objective: 'station-position' };
      throw recoveryError;
    }
  }
}

function planPassageCountsAttempt({ components, anchors, minimumCounts, fixedCounts = [], maximumTotal = 2000, maximumGrowth = 1.5, maximumPasses = 32 }, projectionObjective = 'relative-interval') {
  if (!Array.isArray(components) || !components.length || !Array.isArray(anchors) || anchors.length !== minimumCounts?.length + 1
    || anchors.some((v, i) => !Number.isFinite(v) || i && !(v > anchors[i - 1]))
    || !Number.isInteger(maximumPasses) || maximumPasses < 1 || !Number.isFinite(maximumGrowth) || maximumGrowth <= 1)
    throw new Error('Invalid passage count planning inputs.');
  const spans = components.map(c => ({ from: anchors.indexOf(c.start), to: anchors.indexOf(c.end) }));
  if (spans.some(s => s.from < 0 || s.to <= s.from)) throw new Error('Passage count bounds must match global anchors.');
  const envelopes = components.map(c => {
    const walls = c.boundaries.filter(b => b.kind === 'wall');
    return walls.length ? createSpacingEnvelope({ profiles: walls.map(b => b.requested.map(s => (s - b.range[0]) / b.length)), start: 0, end: 1 }) : null;
  });
  const requirements = spans.map((s, k) => ({ ...s, minimum: Math.max(3, Math.ceil((envelopes[k]?.total ?? 0) - 1e-12)) }));
  const history = [];
  for (let pass = 0; pass < maximumPasses; pass++) {
    const allocation = allocateOrderedIntervalCounts({ minimumCounts, requirements, fixedCounts, maximumTotal });
    const current = components.map((c, k) => {
      const intervals = allocation.indices[spans[k].to] - allocation.indices[spans[k].from];
      return { ...c, intervals, boundaries: c.boundaries.map(b => {
        if (b.kind !== 'cut' || !b.progress || b.progress.length === intervals + 1) return b;
        const old = b.progress, n = old.length - 1;
        return { ...b, progress: Array.from({ length: intervals + 1 }, (_, i) => {
          if (!i) return 0; if (i === intervals) return 1;
          const s = n * i / intervals, a = Math.min(n - 1, Math.floor(s));
          return old[a] + (s - a) * (old[a + 1] - old[a]);
        }) };
      }) };
    });
    const targets = passageStationTargets({ components: current, retainEndpointResolution: true });
    const deficits = [];
    for (let k = 0; k < current.length; k++) {
      const intervals = current[k].intervals, firstSpacing = targets.endpointFit.values[2 * k], lastSpacing = targets.endpointFit.values[2 * k + 1];
      const bounds = stationGrowthFeasibility({ intervals, firstSpacing, lastSpacing, maximumGrowth });
      if (bounds.feasible) continue;
      let required = intervals + 1;
      while (required <= maximumTotal && !stationGrowthFeasibility({ intervals: required, firstSpacing, lastSpacing, maximumGrowth }).feasible) required++;
      // End targets are fitted again after allocation. If no larger count
      // can fit these frozen targets, take one increment and re-evaluate;
      // this is not an infeasibility proof for the simultaneous problem.
      if (required > maximumTotal) required = intervals + 1;
      deficits.push({ component: k, reason: 'endpoint growth feasibility', required, ...bounds });
    }
    let fitted, density;
    if (!deficits.length) {
      fitted = reconcilePassageStations({ components: current, endpointMethod: 'growth-constrained', maximumGrowth, retainEndpointResolution: true, projectionMethod: 'primal-dual', projectionObjective });
      density = current.map((c, k) => {
        const envelope = envelopes[k]; if (!envelope) return { referenceOnly: true, resolutionSatisfied: true };
        const q = fitted.fitted[k].progress.map(envelope.metric), increments = q.slice(1).map((v, i) => v - q[i]);
        const maximumDensityInterval = Math.max(...increments), resolutionSatisfied = maximumDensityInterval <= 1 + 1e-10;
        if (!resolutionSatisfied) {
          // Splitting the current metric intervals this many times would
          // resolve them. A subsequent endpoint refit is measured afresh;
          // this proposal is not a lower bound on the joint optimum.
          const subdivided = increments.reduce((sum, v) => sum + Math.max(1, Math.ceil(v - 1e-10)), 0);
          deficits.push({ component: k, reason: 'original density unresolved', required: Math.max(c.intervals + 1, subdivided), maximumDensityInterval });
        }
        return { maximumDensityInterval, resolutionSatisfied };
      });
    }
    history.push({ pass, counts: allocation.counts, total: allocation.total, deficits, density });
    if (!deficits.length) return { ...allocation, components: current, fitted, density, history, maximumGrowth,
      retainEndpointResolution: true, allOriginalResolutionSatisfied: true, globalMinimumProven: false, exactModernMsetLaw: false,
      scope: 'Automatic counts for saved boundary geometry and density/growth constraints; not complete mesh or physical validation.' };
    for (const d of deficits) requirements[d.component].minimum = Math.max(requirements[d.component].minimum, d.required);
  }
  const error = new Error(`Passage count planning did not satisfy its measured constraints in ${maximumPasses} passes.`);
  if (history.at(-1)?.deficits.length && history.at(-1).deficits.every(d => d.reason === 'original density unresolved'))
    error.code = 'passage-density-unresolved';
  error.history = history; throw error;
}
