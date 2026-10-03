// SPDX-License-Identifier: GPL-2.0-or-later
// Compression reporting and MSES dissipation broadening. Compression
// candidates are observational; they do not choose governing equations or
// certify a shock or bifurcation.
export function auditStreamtubeShocks(flow, { gamma = 1.4, limit = 8 } = {}) {
  const candidates = [], flaggedCandidates = [], nodes = flow.nodes;
  let candidateCount = 0, weakMomentumCandidates = 0, maximumSkewDegrees = 0;
  let supersonicSections = 0, supersonicBoundarySections = 0;
  for (let i = 0; i < flow.sections.length; i++) for (let g = 0; g < flow.sections[i].length; g++)
    for (let j = 0; j < flow.sections[i][g].length; j++) {
      if (flow.sections[i][g][j].machSquared <= 1) continue;
      supersonicSections++;
      if (i === 0 || i === flow.sections.length - 1 || g === 0 && j === 0
        || g === flow.sections[i].length - 1 && j === flow.sections[i][g].length - 1) supersonicBoundarySections++;
    }
  const mid = (a, b) => ({ x: .5 * (a.x + b.x), y: .5 * (a.y + b.y) });
  const logP0 = s => Math.log(s.p) + gamma / (gamma - 1) * Math.log1p(.5 * (gamma - 1) * s.machSquared);
  for (let k = 0; k < flow.cells.length; k++) for (let g = 0; g < flow.cells[k].length; g++)
    for (let j = 0; j < flow.cells[k][g].length; j++) {
      const cell = flow.cells[k][g][j], [a, b] = cell.states;
      if (!(a.machSquared >= .9 ** 2 && b.q < a.q && b.p > a.p)) continue;
      candidateCount++;
      const i = k + 1, row = nodes[g][i], left = mid(nodes[g][i - 1][j], nodes[g][i - 1][j + 1]);
      const right = mid(nodes[g][i + 1][j], nodes[g][i + 1][j + 1]);
      const sx = right.x - left.x, sy = right.y - left.y, cx = row[j + 1].x - row[j].x, cy = row[j + 1].y - row[j].y;
      const along = Math.hypot(sx, sy), across = Math.hypot(cx, cy);
      const skewDegrees = Math.asin(Math.min(1, Math.abs(sx * cx + sy * cy) / (along * across))) * 180 / Math.PI;
      maximumSkewDegrees = Math.max(maximumSkewDegrees, skewDegrees);
      const blend = flow.hybridCells?.[k][g][j], fraction = blend?.fraction ?? null;
      const sonicCrossing = a.machSquared > 1 && b.machSquared < 1;
      const weakMomentum = sonicCrossing && fraction !== null && fraction < .5;
      if (weakMomentum) weakMomentumCandidates++;
      const pressure = .5 * (a.p + b.p), center = mid(row[j], row[j + 1]);
      candidates.push({ i, group: g, tube: j, ...center, upstreamMach: Math.sqrt(a.machSquared),
        downstreamMach: Math.sqrt(b.machSquared), sonicCrossing, weakMomentum,
        pressureRise: (b.p - a.p) / pressure, momentumFraction: fraction,
        momentumResidualOverPressure: cell.streamwiseResidual / pressure,
        entropyResidualOverPressure: cell.isentropicResidual / pressure,
        logTotalPressureChange: logP0(b) - logP0(a), lossIndicators: blend?.lossIndicators,
        skewDegrees, streamwiseSpacing: .5 * along, alongToAcross: .5 * along / across });
      if (weakMomentum && flaggedCandidates.length < limit) flaggedCandidates.push(candidates.at(-1));
      candidates.sort((a, b) => b.pressureRise - a.pressureRise);
      if (candidates.length > limit) candidates.length = limit;
    }
  return { interpretation: 'Compression candidates; not a shock certificate or bifurcation diagnosis',
    candidateCount, weakMomentumCandidates, maximumSkewDegrees, supersonicSections, supersonicBoundarySections, candidates, flaggedCandidates };
}

// MSES 3.05 manual section 1.2.7, preferred sixth-order near-root schedule.
// The previous ACCEPTED fractional density change is frozen for this Newton
// linearization. MUCON and the governing ISMOM selection are unchanged.
export function msesTemporaryMcrit(targetMcrit, previousDensityChange) {
  if (!Number.isFinite(targetMcrit) || targetMcrit < .75 || targetMcrit > 1
    || !Number.isFinite(previousDensityChange) || previousDensityChange < 0)
    throw new Error('Invalid MSES dissipation-enhancement controls.');
  const epsilon = .15, d = previousDensityChange;
  if (d === 0) return targetMcrit;
  if (d > 10) return .75; // exp(-r*r) has already underflowed at this bound.
  const r = (d / epsilon) * (d * d / (d * d + (epsilon / 4) ** 2));
  return .75 + (targetMcrit - .75) * Math.exp(-r * r);
}

// Compare scheduled and hypothetical undamped broadening without changing
// the MSES accepted-update law or applying an unbounded Newton correction.
export function auditDissipationDamping(targetMcrit, undampedDensityChange, step) {
  if (!Number.isFinite(step) || step < 0 || step > 1) throw new Error('Invalid accepted Newton fraction.');
  const acceptedDensityChange = step * undampedDensityChange;
  const acceptedMcrit = msesTemporaryMcrit(targetMcrit, acceptedDensityChange);
  const undampedMcrit = msesTemporaryMcrit(targetMcrit, undampedDensityChange);
  return { acceptedDensityChange, undampedDensityChange, step, acceptedMcrit, undampedMcrit,
    dampingSuppressesBroadening: undampedMcrit < targetMcrit - .01 && acceptedMcrit > targetMcrit - .001 };
}
