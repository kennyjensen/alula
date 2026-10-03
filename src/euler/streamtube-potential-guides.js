// SPDX-License-Identifier: GPL-2.0-or-later
// Panel LE/TE blocks first, then surface spacing. A logical block rank is
// not a velocity potential; each fluid passage retains its own branch.
import { createPotentialCrosslines } from '../geometry/streamtube-crosslines.js';
import { createPanelPotentialBlocks } from '../geometry/panel-potential-blocks.js';
import { resolveSurfaceTurns } from '../geometry/surface-turn-spacing.js';
import { createContourArc } from '../geometry/contour-arc.js';
import { createMonotoneCubicMap } from '../numerics/monotone-cubic.js';
import { createOuterStreamtubeSpacing } from '../geometry/outer-streamtube-spacing.js';
import { reconcileFacingArcSlopes } from '../geometry/facing-arc-spacing.js';
import { reconcileFacingOutlineSlopes } from '../geometry/facing-outline-spacing.js';
import { distributeSurfaceByDensity } from '../geometry/surface-density-stations.js';
import { reconcileBodyEndpointStations } from '../geometry/body-endpoint-stations.js';
import { matchPassageDensityGuides } from './passage-density-guides.js';
import { streamtubeCutEdge, streamtubeCutReferencePoints } from './streamtube-cut-reference.js';

const interpolate = (rows, options) => {
  const map = createMonotoneCubicMap(rows.map(r => r.rank), rows.map(r => r.value), options);
  const inverse = value => {
    if (!Number.isFinite(value) || value < rows[0].value || value > rows.at(-1).value)
      throw new Error('Requested surface value is outside its panel blocks.');
    const anchor = rows.find(r => r.value === value);
    if (anchor) return anchor.rank;
    let lo = rows[0].rank, hi = rows.at(-1).rank;
    for (let k = 0; k < 64; k++) {
      const mid = lo + .5 * (hi - lo);
      if (mid === lo || mid === hi) break;
      if (map.value(mid) < value) lo = mid; else hi = mid;
    }
    return lo + .5 * (hi - lo);
  };
  return { ...map, inverse };
};

export function matchPanelPotentialGuides({ input, profiles, surfaceFractions, outerPaths, sample, sampleCut, sampleCutX, growth, maxIntervals, outerSpread,
  resolveTurns = false, surfaceMetric = 'potential', cutStationSpacing = 'potential', surfaceStationPlacement = 'common-rank', passageIntervalCounts,
  passageCountPlanning = false, passageCountHistory = [] }) {
  if (typeof passageCountPlanning !== 'boolean' || passageCountPlanning && surfaceStationPlacement !== 'passage-density')
    throw new Error('Automatic passage counts require passage-density placement.');
  const localDensity = ['local-density', 'local-density-joined', 'passage-density'].includes(surfaceStationPlacement);
  if (!['common-rank', 'local-density', 'local-density-joined', 'passage-density'].includes(surfaceStationPlacement)
    || localDensity && surfaceMetric !== 'arc'
    || ['local-density-joined', 'passage-density'].includes(surfaceStationPlacement) && cutStationSpacing !== 'physical-x')
    throw new Error('Local surface density placement requires contour-arc maps.');
  if (!['potential', 'physical-x'].includes(cutStationSpacing)
    || cutStationSpacing === 'physical-x' && typeof sampleCutX !== 'function')
    throw new Error('Physical cut spacing requires an x-plane streamline sampler.');
  if (cutStationSpacing === 'physical-x' && surfaceMetric !== 'arc')
    throw new Error('Physical cut spacing requires contour-arc surface maps.');
  const trailing = profiles.map(p => ({ upper: p.phase(0), lower: p.phase(p.curve.length) }));
  const blocks = createPanelPotentialBlocks({
    bodies: profiles.map((p, b) => ({ leading: p.phiStag, trailing: trailing[b],
      inlet: p.upstream[0].potential, outletIncrement: p.wake.at(-1).potential })),
    outer: outerPaths.map(path => ({ inlet: path[0].potential, outlet: path.at(-1).potential })),
  });
  const inverseSurface = (p, side, phi) => {
    const end = trailing[profiles.indexOf(p)][side];
    if (phi < p.phiStag || phi > end) throw new Error('Surface intersection lies outside its monotone potential branch.');
    if (phi === p.phiStag) return 0;
    if (phi === end) return 1;
    let lo = 0, hi = 1;
    for (let k = 0; k < 56; k++) {
      const mid = .5 * (lo + hi);
      if (p.phase(p.curve.branch(side, mid, p.stag).parameter) < phi) lo = mid; else hi = mid;
    }
    return .5 * (lo + hi);
  };
  const known = (b, side) => blocks.events.flatMap(event => {
    const potential = event.targets[side][b];
    return potential === null ? [] : [{ rank: blocks.rank[event.id], potential }];
  }).sort((a, b) => a.rank - b.rank);
  const surfaceMaps = profiles.map((p, b) => {
    const arc = createContourArc(p.curve), origin = arc.at(p.stag);
    const leadingRank = blocks.rank[`${b}:LE`], trailingRank = blocks.rank[`${b}:TE`];
    return Object.fromEntries(['upper', 'lower'].map(side => {
      const at = fraction => surfaceMetric === 'arc'
        ? (side === 'upper' ? -1 : 1) * (arc.at(p.curve.branch(side, fraction, p.stag).parameter) - origin)
        : p.phase(p.curve.branch(side, fraction, p.stag).parameter);
      // A solid body shields its opposite surface from an incoming block.
      // Interpolate only through actual hits, not invented global potentials.
      const rows = known(b, side).filter(r => r.rank >= leadingRank && r.rank <= trailingRank).map(r => ({
        rank: r.rank, value: at(inverseSurface(p, side, r.potential)), potential: r.potential,
      }));
      const descriptor = { at, map: interpolate(rows), rows };
      descriptor.fractionAtPosition = target => {
        let lo = 0, hi = 1;
        for (let k = 0; k < 56; k++) { const mid = .5 * (lo + hi); if (at(mid) < target) lo = mid; else hi = mid; }
        return .5 * (lo + hi);
      };
      descriptor.fraction = rank => {
        const anchor = rows.find(r => r.rank === rank);
        if (anchor) return inverseSurface(p, side, anchor.potential);
        const target = descriptor.map.value(rank);
        if (surfaceMetric !== 'arc') return inverseSurface(p, side, target);
        return descriptor.fractionAtPosition(target);
      };
      return [side, descriptor];
    }));
  });
  // A common rank monitor alone cannot change the correspondence between
  // fixed arc maps. Reconcile the actual facing maps first; opposite sides
  // of one solid body are not a facing pair. Every physical hit stays fixed.
  const facingSpacing = [];
  let physicalCutAnchors, facingOutlineSpacing;
  if (cutStationSpacing === 'physical-x') {
    // Build the complete physical outline before reconciling its spacing.
    // Cut x and wall arc are both increasing lengths, joined at LE/TE. The
    // cut metric remains x, not exact curved-streamline arclength.
    physicalCutAnchors = profiles.map((p, b) => {
      const upperKnown = known(b, 'upper'), le = blocks.rank[`${b}:LE`], te = blocks.rank[`${b}:TE`];
      return Object.fromEntries(['upstream', 'wake'].map(end => {
        const incoming = end === 'upstream', edgeRank = incoming ? le : te;
        const edge = streamtubeCutEdge(p, input.bodies[b], end), path = p[end];
        const targets = incoming
          ? [{ rank: blocks.inletRank, potential: path[0].potential }, ...upperKnown.filter(r => r.rank <= le)]
          : [...upperKnown.filter(r => r.rank >= te).map(r => ({ rank: r.rank, potential: r.potential - trailing[b].upper })),
            { rank: blocks.outletRank, potential: path.at(-1).potential }];
        // Pin body edges analytically. Do not integrate a trace through the
        // stagnation point; its final stored point is only a finite seed.
        const anchors = targets.map(r => ({ ...r, point: r.rank === edgeRank
          ? { ...edge, potential: r.potential } : sampleCut(path, r.potential) }));
        return [end, { anchors, edge }];
      }));
    });
    const descriptors = surfaceMaps.flatMap(sides => [sides.upper, sides.lower]);
    const outlines = descriptors.map((descriptor, m) => {
      const cuts = physicalCutAnchors[Math.floor(m / 2)], wallLength = descriptor.at(1);
      return interpolate([
        ...cuts.upstream.anchors.slice(0, -1).map(a => ({ rank: a.rank, value: a.point.x - cuts.upstream.edge.x })),
        ...descriptor.rows,
        ...cuts.wake.anchors.slice(1).map(a => ({ rank: a.rank, value: wallLength + (a.point.x - cuts.wake.edge.x) })),
      ]);
    });
    const shared = profiles.flatMap((_, b) => ['upstream', 'wake'].flatMap(end =>
      physicalCutAnchors[b][end].anchors.map(a => ({ rank: a.rank, maps: [2 * b, 2 * b + 1] }))));
    facingOutlineSpacing = reconcileFacingOutlineSlopes({ maps: outlines,
      pairs: profiles.slice(1).map((_, g) => [2 * g, 2 * (g + 1) + 1]), shared });
    facingOutlineSpacing.metric = 'cut physical x joined to wall contour arc';
    facingOutlineSpacing.outlines = outlines.map(({ knots, values }) => ({ knots, values }));
    descriptors.forEach((descriptor, m) => {
      const outline = outlines[m];
      descriptor.map = interpolate(outline.knots.map((rank, i) => ({ rank, value: outline.values[i] })),
        { derivatives: 'prescribed', slopes: facingOutlineSpacing.slopes[m] });
    });
  } else if (surfaceMetric === 'arc') for (let g = 1; g < profiles.length; g++) {
    const first = surfaceMaps[g - 1].upper, second = surfaceMaps[g].lower;
    const fitted = reconcileFacingArcSlopes(first.map, second.map);
    if (fitted.knots.length) {
      first.map = interpolate(first.rows, { derivatives: 'prescribed', slopes: fitted.first });
      second.map = interpolate(second.rows, { derivatives: 'prescribed', slopes: fitted.second });
    }
    facingSpacing.push({ region: g, ...fitted });
  }
  const turnResolution = [], resolvedDemands = profiles.map(() => ({}));
  const demands = profiles.flatMap((p, b) => ['upper', 'lower'].map(side => {
    const refined = resolveTurns ? resolveSurfaceTurns({ curve: p.curve, side, stagnation: p.stag, fractions: surfaceFractions[b][side] }) : null;
    if (refined) turnResolution.push({ body: b, side, addedPoints: refined.addedPoints,
      maximumTangentCone: refined.maximumTangentCone, maxTurn: refined.maxTurn });
    const { at, map } = surfaceMaps[b][side];
    resolvedDemands[b][side] = refined?.fractions ?? surfaceFractions[b][side];
    return resolvedDemands[b][side].map((f, i, row) => !i ? blocks.rank[`${b}:LE`]
      : i === row.length - 1 ? blocks.rank[`${b}:TE`] : map.inverse(at(f)));
  }));
  // The common envelope combines resolution demands on the reconciled
  // maps. It does not itself alter their physical correspondence.
  const crosslines = createPotentialCrosslines({ profiles: demands, start: blocks.inletRank, end: blocks.outletRank, growth, maxIntervals,
    ...(passageIntervalCounts ? { blockIntervals: passageIntervalCounts } : {}),
    ...(resolveTurns && input.gridSpacing ? { inletIntervals: input.gridSpacing.inlet.intervals, outletIntervals: input.gridSpacing.outlet.intervals } : {}) });
  const { x } = crosslines, cutSpacing = [], surfaceDensity = [];
  const guides = profiles.map((p, b) => {
    const body = input.bodies[b], le = blocks.rank[`${b}:LE`], te = blocks.rank[`${b}:TE`];
    body.leadingIndex = x.indexOf(le); body.trailingIndex = x.indexOf(te);
    if (body.leadingIndex < 1 || body.trailingIndex <= body.leadingIndex) throw new Error('Missing ordered panel block anchor.');
    body.surfaceFractions = Object.fromEntries(['upper', 'lower'].map(side => {
      const descriptor = surfaceMaps[b][side];
      if (localDensity) {
        const distribution = distributeSurfaceByDensity({ requested: resolvedDemands[b][side].map(descriptor.at),
          anchors: descriptor.rows.map(r => ({ index: x.indexOf(r.rank), position: r.value })) });
        surfaceDensity.push({ body: b, side, ...distribution });
        return [side, distribution.positions.map((s, i, row) => {
          if (!i) return 0; if (i === row.length - 1) return 1;
          const hit = descriptor.rows.find(r => x.indexOf(r.rank) === i + body.leadingIndex);
          return hit ? descriptor.fraction(hit.rank) : descriptor.fractionAtPosition(s);
        })];
      }
      return [side, x.slice(body.leadingIndex, body.trailingIndex + 1).map((rank, i, row) =>
        !i ? 0 : i === row.length - 1 ? 1 : descriptor.fraction(rank))];
    }));
    // Sample upstream/wake cuts once. The jump changes only the potential
    // branch; it must never open a gap between their two geometric copies.
    const upperKnown = known(b, 'upper');
    const upstreamMap = interpolate([{ rank: blocks.inletRank, value: p.upstream[0].potential },
      ...upperKnown.filter(r => r.rank <= le).map(r => ({ rank: r.rank, value: r.potential }))]);
    const wakeMap = interpolate([...upperKnown.filter(r => r.rank >= te).map(r => ({ rank: r.rank, value: r.potential - trailing[b].upper })),
      { rank: blocks.outletRank, value: p.wake.at(-1).potential }]);
    const physicalCuts = {};
    if (cutStationSpacing === 'physical-x') {
      const arc = createContourArc(p.curve);
      for (const end of ['upstream', 'wake']) {
        const incoming = end === 'upstream';
        const { edge, anchors: cutAnchors } = physicalCutAnchors[b][end];
        const path = p[end], seed = incoming ? path.at(-1) : path[0];
        const direction = incoming ? p.geometricStagnationConnector?.anchorDirection : null;
        const edgeDirectionX = direction ? Math.abs(direction.x) / Math.hypot(direction.x, direction.y)
          : Math.abs(seed.x - edge.x) / Math.hypot(seed.x - edge.x, seed.y - edge.y);
        const surfaceArcSpacing = Math.min(...['upper', 'lower'].map(side => {
          const f = body.surfaceFractions[side], k = incoming ? 0 : f.length - 2;
          return Math.abs(arc.at(p.curve.branch(side, f[k + 1], p.stag).parameter)
            - arc.at(p.curve.branch(side, f[k], p.stag).parameter));
        }));
        const anchors = cutAnchors.map(a => ({ ...a, index: x.indexOf(a.rank) }));
        const firstIndex = anchors[0].index, lastIndex = anchors.at(-1).index;
        const descriptor = surfaceMaps[b].upper, wallOffset = incoming ? 0 : descriptor.at(1);
        const stations = x.slice(firstIndex, lastIndex + 1).map(rank =>
          cutAnchors.find(a => a.rank === rank)?.point.x ?? edge.x + (descriptor.map.value(rank) - wallOffset));
        if (stations.some((xx, k) => !Number.isFinite(xx) || k && xx <= stations[k - 1]))
          throw new Error('Reconciled physical cut stations are not resolved and increasing.');
        const points = stations.map((xx, k) => anchors.find(a => a.index === k + firstIndex)?.point ?? sampleCutX(path, xx));
        physicalCuts[end] = { firstIndex, points };
        const neighbor = incoming ? points.at(-2) : points[1];
        cutSpacing.push({ body: b, end, x: stations, firstIndex, metric: 'physical x',
          distribution: 'reconciled full-outline cubic sampled at common rank stations',
          blocks: anchors.slice(1).map((a, k) => {
            const previous = anchors[k], begin = previous.index - firstIndex, finish = a.index - firstIndex;
            const intervals = stations.slice(begin + 1, finish + 1).map((xx, j) => xx - stations[begin + j]);
            return { fromIndex: previous.index, toIndex: a.index, intervals: finish - begin,
              length: a.point.x - previous.point.x, firstSpacing: intervals[0], lastSpacing: intervals.at(-1),
              maximumAdjacentRatio: Math.max(1, ...intervals.slice(1).map((h, j) => Math.max(h / intervals[j], intervals[j] / h))) };
          }), surfaceArcSpacing, edgeDirectionX,
          edgeDirectionSource: direction ? 'geometric incoming connector tangent at the exact material anchor'
            : 'wall-to-trace-seed secant; x projection, not streamline arclength',
          referenceEdgeXSpacing: surfaceArcSpacing * edgeDirectionX,
          achievedEdgeChord: Math.hypot(neighbor.x - edge.x, neighbor.y - edge.y),
          exactDiscreteJoinSpacing: false,
          anchorIndices: anchors.map(a => a.index) });
      }
    }
    const rows = { upper: [], lower: [] };
    x.forEach((rank, i) => {
      let cut, increment;
      if (i < body.leadingIndex) cut = physicalCuts.upstream?.points[i - physicalCuts.upstream.firstIndex]
        ?? sampleCut(p.upstream, upstreamMap.value(rank));
      if (i > body.trailingIndex) {
        if (physicalCuts.wake) { cut = physicalCuts.wake.points[i - physicalCuts.wake.firstIndex]; increment = cut.potential; }
        else { increment = wakeMap.value(rank); cut = sampleCut(p.wake, increment); }
      }
      for (const side of ['upper', 'lower']) {
        if (cut) rows[side].push({ ...cut, potential: increment === undefined ? cut.potential : increment + trailing[b][side] });
        else {
          const v = p.curve.branch(side, body.surfaceFractions[side][i - body.leadingIndex], p.stag);
          rows[side].push({ ...v.point, potential: p.phase(v.parameter) });
        }
      }
    });
    return rows;
  });
  const endpointSpacing = [];
  if (surfaceStationPlacement === 'local-density-joined') for (const [b, p] of profiles.entries()) {
    const body = input.bodies[b], components = {}, cutReferences = {};
    for (const side of ['upper', 'lower']) {
      const d = surfaceDensity.find(d => d.body === b && d.side === side);
      components[side] = { positions: d.positions, requested: resolvedDemands[b][side].map(surfaceMaps[b][side].at),
        anchors: d.anchors.map(({ index, position }) => ({ index, position })) };
    }
    for (const end of ['upstream', 'wake']) {
      const cut = cutSpacing.find(c => c.body === b && c.end === end), firstIndex = cut.anchorIndices[0];
      const points = streamtubeCutReferencePoints(body, guides[b], end, firstIndex, cut.anchorIndices.at(-1)), positions = [0];
      for (let i = 1; i < points.length; i++) positions.push(positions.at(-1) + Math.hypot(points[i].x - points[i - 1].x, points[i].y - points[i - 1].y));
      components[end] = { positions, anchors: cut.anchorIndices.map(index => ({ index, position: positions[index - firstIndex] })) };
      cutReferences[end] = { firstIndex, points, positions, cut };
    }
    const fitted = reconcileBodyEndpointStations(components);
    for (const side of ['upper', 'lower']) {
      const descriptor = surfaceMaps[b][side], distribution = fitted[side];
      body.surfaceFractions[side] = distribution.positions.map((s, k, row) => {
        if (!k) return 0; if (k === row.length - 1) return 1;
        const hit = descriptor.rows.find(r => x.indexOf(r.rank) === k + body.leadingIndex);
        return hit ? descriptor.fraction(hit.rank) : descriptor.fractionAtPosition(s);
      });
      for (let k = 0; k < distribution.positions.length; k++) {
        const v = p.curve.branch(side, body.surfaceFractions[side][k], p.stag);
        guides[b][side][body.leadingIndex + k] = { ...v.point, potential: p.phase(v.parameter) };
      }
    }
    for (const end of ['upstream', 'wake']) {
      const { points, positions, firstIndex, cut } = cutReferences[end], distribution = fitted[end];
      const sampleArc = s => {
        const exact = positions.indexOf(s); if (exact >= 0) return { ...points[exact] };
        const k = positions.findIndex(v => v > s) - 1, t = (s - positions[k]) / (positions[k + 1] - positions[k]);
        const xx = points[k].x + t * (points[k + 1].x - points[k].x);
        const sampled = sampleCutX(p[end], xx);
        return { ...sampled, potential: sampled.potential + (end === 'wake' ? trailing[b].lower : 0) };
      };
      const sampled = distribution.positions.map(sampleArc);
      for (let k = 0; k < sampled.length; k++) {
        const i = firstIndex + k;
        // LE/TE wall points above are exact; cuts contain both endpoints.
        if (i === body.leadingIndex || i === body.trailingIndex) continue;
        guides[b].lower[i] = sampled[k];
        guides[b].upper[i] = { ...sampled[k], potential: sampled[k].potential + (end === 'wake' ? trailing[b].upper - trailing[b].lower : 0) };
      }
      cut.preEndpointFit = { x: cut.x, blocks: cut.blocks, achievedEdgeChord: cut.achievedEdgeChord,
        surfaceArcSpacing: cut.surfaceArcSpacing, referenceEdgeXSpacing: cut.referenceEdgeXSpacing };
      cut.x = sampled.map(point => point.x);
      cut.distribution = 'two-ended Vinokur in local cut-density coordinate; saved cut chord arc mapped to traced physical x';
      cut.blocks = distribution.blocks;
      cut.metric = 'cut polyline arc for endpoint fit; exported stations are physical x';
      cut.surfaceArcSpacing = fitted.joins.find(j => j.edge === (end === 'upstream' ? 'LE' : 'TE')).target;
      cut.referenceEdgeXSpacing = cut.surfaceArcSpacing * cut.edgeDirectionX;
      const edge = end === 'upstream' ? sampled.at(-1) : sampled[0], neighbor = end === 'upstream' ? sampled.at(-2) : sampled[1];
      cut.achievedEdgeChord = Math.hypot(edge.x - neighbor.x, edge.y - neighbor.y);
      cut.exactDiscreteJoinSpacing = false; // Curved resampling is checked separately.
    }
    endpointSpacing.push({ body: b, ...fitted });
  }
  const outerSpacing = [];
  const outer = outerPaths.map((path, side) => {
    const rows = [{ rank: blocks.inletRank, value: path[0].potential },
      ...blocks.events.flatMap(event => event.targets.outer[side] === null ? [] : [{ rank: blocks.rank[event.id], value: event.targets.outer[side] }]),
      { rank: blocks.outletRank, value: path.at(-1).potential }].sort((a, b) => a.rank - b.rank);
    // Shielded logical stations have interpolated targets, not additional
    // isopotentials. Keep those targets independent of the requested count.
    // Published one-ended stretching matches the body-region join slope
    // while grading the outer inlet/outlet blocks toward zero end curvature.
    const spacing = createOuterStreamtubeSpacing({ ranks: x, anchors: crosslines.anchors,
      targets: rows.map(r => ({ rank: r.rank, potential: r.value })), spread: outerSpread, endSpacing: 'vinokur' });
    outerSpacing.push(spacing.diagnostics);
    return spacing.potentials.map(potential => sample(path, potential));
  });
  const passageDensity = surfaceStationPlacement === 'passage-density'
    ? matchPassageDensityGuides({ input, profiles, guides, outer, blocks, ranks: x, surfaceMaps, resolvedDemands, sampleCutX,
      ...(passageCountPlanning ? { countPlanning: { anchors: crosslines.anchors, minimumCounts: crosslines.blocks.map(b => b.intervals),
        // Automatic-mode counts are minima. Fine endpoint resolution can
        // require more inlet/wake cells to obey the same growth bound.
        // The generic planner still supports explicitly fixed counts.
        maximumTotal: maxIntervals } } : {}) }) : null;
  if (passageDensity?.rebuildRequired) {
    if (passageCountHistory.length >= 8) throw new Error('Automatic passage counts did not stabilize on resampled boundary geometry.');
    // Reuse the same panel paths, original wall requests and physical hits.
    // Only guide station counts/placement are rebuilt; no new flow tracing.
    return matchPanelPotentialGuides({ input, profiles, surfaceFractions, outerPaths, sample, sampleCut, sampleCutX, growth, maxIntervals, outerSpread,
      resolveTurns, surfaceMetric, cutStationSpacing, surfaceStationPlacement, passageIntervalCounts: passageDensity.countPlan.counts,
      passageCountPlanning, passageCountHistory: [...passageCountHistory, passageDensity.countPlan] });
  }
  if (passageDensity?.countPlan) passageDensity.countPlan.geometryReplans = passageCountHistory;
  if (passageDensity) for (const cut of cutSpacing) {
    cut.prePassageFit = { x: cut.x, blocks: cut.blocks, achievedEdgeChord: cut.achievedEdgeChord,
      surfaceArcSpacing: cut.surfaceArcSpacing, referenceEdgeXSpacing: cut.referenceEdgeXSpacing };
    const lastIndex = cut.anchorIndices.at(-1), points = guides[cut.body].lower.slice(cut.firstIndex, lastIndex + 1);
    cut.x = points.map(p => p.x);
    const a = cut.end === 'upstream' ? points.at(-1) : points[0], b = cut.end === 'upstream' ? points.at(-2) : points[1];
    cut.achievedEdgeChord = Math.hypot(a.x - b.x, a.y - b.y);
    cut.surfaceArcSpacing = Math.min(...['upper', 'lower'].map(side => {
      const f = input.bodies[cut.body].surfaceFractions[side], k = cut.end === 'upstream' ? 0 : f.length - 2;
      return surfaceMaps[cut.body][side].at(f[k + 1]) - surfaceMaps[cut.body][side].at(f[k]);
    }));
    cut.referenceEdgeXSpacing = cut.surfaceArcSpacing * cut.edgeDirectionX;
    cut.blocks = []; // The prior full-outline cubic is no longer the placement map.
    cut.distribution = 'passage-local density with convex endpoint/growth projection';
    cut.metric = passageDensity.cutMetric;
  }
  const point = p => ({ x: p.x, y: p.y });
  input.outerLower = outer[0].map(point); input.outerUpper = outer[1].map(point);
  input.cutPaths = guides.map(rows => rows.lower.map((p, i) => ({ x: .5 * (p.x + rows.upper[i].x), y: .5 * (p.y + rows.upper[i].y) })));
  return { guides, outer, diagnostics: { coordinate: 'ordered panel potential blocks', surfaceMetric, cutStationSpacing, cutSpacing, outerSpread,
    ...(localDensity ? { surfaceStationPlacement, surfaceDensity,
      ...(passageDensity ? { passageDensity, surfaceDensityStage: 'before passage fit', facingOutlineSpacingStage: 'count construction before passage fit' } : {}),
      ...(endpointSpacing.length ? { endpointSpacing, surfaceDensityStage: 'before joint endpoint fit',
        facingOutlineSpacingStage: 'count construction before local density and joint endpoint placement' } : {}) } : {}),
    outerSpreading: 'C1 monotone through LE/TE anchors', outerSpacing, facingSpacing, facingOutlineSpacing, turnResolution, ...crosslines, panelBlocks: blocks,
    bodies: profiles.map((p, b) => ({ leading: p.phiStag, trailing: .5 * (trailing[b].upper + trailing[b].lower),
      offsets: { upper: .5 * (trailing[b].upper - trailing[b].lower), lower: .5 * (trailing[b].lower - trailing[b].upper) },
      leadingIndex: input.bodies[b].leadingIndex, trailingIndex: input.bodies[b].trailingIndex })) } };
}
