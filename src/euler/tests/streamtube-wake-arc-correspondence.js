// SPDX-License-Identifier: GPL-2.0-or-later
// Optional recovery initialization, not an equation or a Newton constraint.
// Resample the current mean-wake polyline at the immutable reference-cut arc
// fractions, retain original-index normal gaps, then pair equal bank advances.
// Mean TE/outlet centers and both solid TE points stay fixed. Individual outlet
// bank points may move. Full geometry, gas and coupled replay are caller gates.
import { streamtubeWakeGap } from '../streamtube-wake-geometry.js';

export function initializeStreamtubeWakeArcCorrespondence({ nodes, layout, massFractions, referenceCutPaths } = {}) {
  const { nx, elements, tubes, bodies, independentWakeBanks } = layout ?? {};
  if (!Number.isInteger(nx) || nx < 2 || !Number.isInteger(elements) || elements < 1
    || independentWakeBanks !== true || !Array.isArray(tubes) || tubes.length !== elements + 1
    || !tubes.every(v => Number.isInteger(v) && v > 0) || !Array.isArray(bodies) || bodies.length !== elements
    || bodies.some(b => !Number.isInteger(b?.leadingIndex) || !Number.isInteger(b?.trailingIndex)
      || b.leadingIndex < 0 || b.leadingIndex >= b.trailingIndex || b.trailingIndex >= nx))
    throw new Error('Invalid independent-bank wake arc layout.');
  const finitePoint = p => Number.isFinite(p?.x) && Number.isFinite(p?.y);
  if (!Array.isArray(nodes) || nodes.length !== tubes.length || nodes.some((group, g) =>
    !Array.isArray(group) || group.length !== nx + 1 || group.some(row => !Array.isArray(row)
      || row.length !== tubes[g] + 1 || row.some(p => !finitePoint(p)))))
    throw new Error('Invalid wake arc coordinates.');
  if (!Array.isArray(massFractions) || massFractions.length !== tubes.length || massFractions.some((row, g) =>
    !Array.isArray(row) || row.length !== tubes[g] + 1 || row[0] !== 0 || row.at(-1) !== 1
    || row.some((v, j) => !Number.isFinite(v) || j > 0 && v <= row[j - 1])))
    throw new Error('Supply strictly increasing cumulative passage mass fractions from zero to one.');
  if (!Array.isArray(referenceCutPaths) || referenceCutPaths.length !== elements || referenceCutPaths.some(row =>
    !Array.isArray(row) || row.length !== nx + 1 || row.some(p => !finitePoint(p))))
    throw new Error('Supply one finite reference cut path per body with the original station count.');
  let coordinateScale = 0;
  for (const group of nodes) for (const row of group) for (const p of row)
    coordinateScale = Math.max(coordinateScale, Math.abs(p.x), Math.abs(p.y));
  const roundoffTolerance = 128 * Number.EPSILON * coordinateScale;
  const fail = (body, station, message) => {
    const error = new Error(`Wake arc correspondence body ${body}, station ${station}: ${message}`);
    Object.assign(error, { code: 'WAKE_ARC_CORRESPONDENCE_GEOMETRY', body, station }); throw error;
  };
  const difference = (p, q) => ({ x: p.x - q.x, y: p.y - q.y });
  const dot = (p, q) => p.x * q.x + p.y * q.y;
  const distance = (p, q) => Math.hypot(p.x - q.x, p.y - q.y);
  const center = (grid, b, i) => ({ x: .5 * (grid[b][i].at(-1).x + grid[b + 1][i][0].x),
    y: .5 * (grid[b][i].at(-1).y + grid[b + 1][i][0].y) });
  const separation = (grid, b, i) => difference(grid[b + 1][i][0], grid[b][i].at(-1));
  const unit = (v, b, i) => { const d = Math.hypot(v.x, v.y);
    if (!(d > 0) || !Number.isFinite(d)) fail(b, i, 'Degenerate wake-center frame.');
    return { x: v.x / d, y: v.y / d }; };
  const arc = (points, b, te, kind) => {
    const values = [0];
    for (let k = 1; k < points.length; k++) {
      const d = distance(points[k - 1], points[k]), next = values.at(-1) + d;
      if (!(d > 0) || !Number.isFinite(next) || !(next > values.at(-1)))
        fail(b, te + k, `Degenerate or unresolved ${kind} arc interval.`);
      values.push(next);
    }
    return values;
  };
  const gap = (grid, b, i) => {
    try {
      const ids = [i - 1, i, Math.min(nx, i + 1)];
      const { apply, ...result } = streamtubeWakeGap(ids.map(k => grid[b][k].at(-1)), ids.map(k => grid[b + 1][k][0]));
      if (!Number.isFinite(result.gap) || result.gap < -roundoffTolerance) fail(b, i, 'Invalid normal wake gap.');
      return result;
    } catch (error) {
      if (error.code === 'WAKE_ARC_CORRESPONDENCE_GEOMETRY') throw error;
      fail(b, i, error.message);
    }
  };
  const moved = structuredClone(nodes).map(group => group.map(row => row.map(p => ({ ...p }))));
  const diagnostics = { initialGuessOnly: true, equationsChanged: false, geometryAcceptanceRequired: true,
    coordinateRule: 'reference-arc-fractions-and-equal-center-edge-projected-bank-advance',
    referenceScope: 'Original-index gap; current polyline; mean TE/outlet fixed, individual outlet banks free.',
    roundoffTolerance, maximumGapChange: 0, maximumCenterChange: 0, maximumCenterSamplingError: 0,
    maximumNodeDisplacement: 0, maximumTangentialMismatch: 0, minimumFrameProjection: Infinity,
    minimumBankForwardAdvance: Infinity, maximumForwardAdvanceMismatch: 0,
    minimumInputGap: Infinity, minimumOutputGap: Infinity, bodies: [] };
  const prescribed = [];
  for (let b = 0; b < elements; b++) {
    const te = bodies[b].trailingIndex;
    const current = Array.from({ length: nx - te + 1 }, (_, k) => center(nodes, b, te + k));
    const currentArc = arc(current, b, te, 'current'), referenceArc = arc(referenceCutPaths[b].slice(te), b, te, 'reference');
    const fractions = referenceArc.map(v => v / referenceArc.at(-1));
    const sampled = [], sampling = [];
    let segment = 0;
    for (let k = 0; k < fractions.length; k++) {
      const wanted = fractions[k] * currentArc.at(-1);
      while (segment + 1 < currentArc.length - 1 && currentArc[segment + 1] < wanted) segment++;
      const weight = (wanted - currentArc[segment]) / (currentArc[segment + 1] - currentArc[segment]);
      if (!(weight >= 0 && weight <= 1) || !Number.isFinite(weight)) fail(b, te + k, 'Invalid arc interpolation.');
      const a = current[segment], c = current[segment + 1];
      sampled.push(k === 0 ? { ...current[0] } : k === fractions.length - 1 ? { ...current.at(-1) }
        : { x: (1 - weight) * a.x + weight * c.x, y: (1 - weight) * a.y + weight * c.y });
      sampling.push({ station: te + k, fraction: fractions[k], sourceSegment: te + segment, weight });
    }
    const firstDirection = unit(difference(sampled[1], sampled[0]), b, te + 1);
    let previous = separation(nodes, b, te);
    const teGap = previous.y * firstDirection.x - previous.x * firstDirection.y;
    if (!Number.isFinite(teGap) || teGap < -roundoffTolerance) fail(b, te, 'Negative normal TE gap in resampled frame.');
    const body = { body: b, trailingIndex: te, stations: nx - te,
      referenceArcLength: referenceArc.at(-1), originalCenterArcLength: currentArc.at(-1),
      sampledCenterArcLength: arc(sampled, b, te, 'sampled').at(-1),
      firstCenterAdvanceBefore: currentArc[1], firstCenterAdvanceAfter: distance(sampled[0], sampled[1]),
      teNormalGap: teGap, sampling, firstBefore: gap(nodes, b, te + 1), maximumBankDisplacement: 0 };
    prescribed[b] = [];
    for (let k = 1; k < sampled.length; k++) {
      const i = te + k, e = unit(difference(sampled[k], sampled[k - 1]), b, i);
      const t = unit(difference(sampled[Math.min(k + 1, sampled.length - 1)], sampled[k - 1]), b, i);
      const n = { x: -t.y, y: t.x }, before = gap(nodes, b, i), denominator = dot(e, t);
      if (!(denominator > 64 * Number.EPSILON)) fail(b, i, 'Singular or backward resampled pairing frame.');
      const tau = (dot(e, previous) - before.gap * dot(e, n)) / denominator;
      const d = { x: tau * t.x + before.gap * n.x, y: tau * t.y + before.gap * n.y };
      if (![d.x, d.y].every(Number.isFinite)) fail(b, i, 'Nonfinite paired separation.');
      moved[b][i][tubes[b]] = { x: sampled[k].x - .5 * d.x, y: sampled[k].y - .5 * d.y };
      moved[b + 1][i][0] = { x: sampled[k].x + .5 * d.x, y: sampled[k].y + .5 * d.y };
      previous = d; prescribed[b][k] = before.gap;
      diagnostics.minimumInputGap = Math.min(diagnostics.minimumInputGap, before.gap);
      diagnostics.minimumFrameProjection = Math.min(diagnostics.minimumFrameProjection, denominator);
      for (const [g, j] of [[b, tubes[b]], [b + 1, 0]])
        body.maximumBankDisplacement = Math.max(body.maximumBankDisplacement, distance(moved[g][i][j], nodes[g][i][j]));
    }
    body.outletBankMovement = { lower: distance(moved[b][nx].at(-1), nodes[b][nx].at(-1)),
      upper: distance(moved[b + 1][nx][0], nodes[b + 1][nx][0]) };
    for (let k = 1; k < sampled.length; k++) {
      const i = te + k, c = center(moved, b, i), e = unit(difference(c, center(moved, b, i - 1)), b, i), after = gap(moved, b, i);
      const ds = distance(c, center(moved, b, i - 1));
      diagnostics.maximumGapChange = Math.max(diagnostics.maximumGapChange, Math.abs(after.gap - prescribed[b][k]));
      diagnostics.maximumCenterChange = Math.max(diagnostics.maximumCenterChange, distance(c, current[k]));
      diagnostics.maximumCenterSamplingError = Math.max(diagnostics.maximumCenterSamplingError, distance(c, sampled[k]));
      diagnostics.maximumTangentialMismatch = Math.max(diagnostics.maximumTangentialMismatch,
        Math.abs(dot(e, difference(separation(moved, b, i), separation(moved, b, i - 1)))));
      diagnostics.minimumOutputGap = Math.min(diagnostics.minimumOutputGap, after.gap);
      for (const [g, j] of [[b, tubes[b]], [b + 1, 0]]) {
        const forward = dot(e, difference(moved[g][i][j], moved[g][i - 1][j]));
        diagnostics.minimumBankForwardAdvance = Math.min(diagnostics.minimumBankForwardAdvance, forward);
        diagnostics.maximumForwardAdvanceMismatch = Math.max(diagnostics.maximumForwardAdvanceMismatch, Math.abs(forward - ds));
      }
    }
    body.firstAfter = gap(moved, b, te + 1); diagnostics.bodies.push(body);
  }
  // A shared passage may receive distinct corrections on both banks. Apply
  // their physical-mass interpolation together, exactly once per interior node.
  for (let g = 0; g < tubes.length; g++) for (let i = 0; i <= nx; i++) for (let j = 1; j < tubes[g]; j++) {
    const eta = massFractions[g][j];
    for (const key of ['x', 'y']) {
      const delta = (1 - eta) * (moved[g][i][0][key] - nodes[g][i][0][key])
        + eta * (moved[g][i].at(-1)[key] - nodes[g][i].at(-1)[key]);
      if (delta !== 0) moved[g][i][j][key] += delta;
    }
  }
  for (let g = 0; g < tubes.length; g++) for (let i = 0; i <= nx; i++) for (let j = 0; j <= tubes[g]; j++) {
    const p = moved[g][i][j];
    if (!finitePoint(p)) throw new Error('Nonfinite wake arc correspondence result.');
    diagnostics.maximumNodeDisplacement = Math.max(diagnostics.maximumNodeDisplacement, distance(p, nodes[g][i][j]));
  }
  return { nodes: moved, diagnostics };
}
