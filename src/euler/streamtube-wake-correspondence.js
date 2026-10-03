// SPDX-License-Identifier: GPL-2.0-or-later
// One-time coordinate initializer reconstructed from the documented research
// control. It is not an ISES listing or an iterated physical wake constraint.
// The caller must check the resulting grid and replay the unchanged equations.
// Keep each wake center C and its centered-normal gap g. For center-edge unit
// tangent e, impose e dot (D_i-D_(i-1)) = 0, where D is the bank separation.
// Writing D_i = tau_i*t_i + g_i*n_i gives
// tau_i = (e dot D_(i-1) - g_i*(e dot n_i))/(e dot t_i).
// Both bank intervals then have the center's forward projected advance. This
// coordinate condition does not establish convexity or physical admissibility.
import { streamtubeWakeGap } from './streamtube-wake-geometry.js';

export function initializeStreamtubeWakeCorrespondence({ nodes, layout, massFractions } = {}) {
  const { nx, elements, tubes, bodies, independentWakeBanks } = layout ?? {};
  if (!Number.isInteger(nx) || nx < 2 || !Number.isInteger(elements) || elements < 1
    || independentWakeBanks !== true || !Array.isArray(tubes) || tubes.length !== elements + 1
    || !tubes.every(n => Number.isInteger(n) && n > 0)
    || !Array.isArray(bodies) || bodies.length !== elements || bodies.some(b =>
      !Number.isInteger(b?.leadingIndex) || !Number.isInteger(b?.trailingIndex)
      || b.leadingIndex < 0 || b.leadingIndex >= b.trailingIndex || b.trailingIndex >= nx))
    throw new Error('Invalid independent-bank wake correspondence layout.');
  if (!Array.isArray(nodes) || nodes.length !== tubes.length || nodes.some((group, g) =>
    !Array.isArray(group) || group.length !== nx + 1 || group.some(row =>
      !Array.isArray(row) || row.length !== tubes[g] + 1 || row.some(p =>
        !Number.isFinite(p?.x) || !Number.isFinite(p?.y)))))
    throw new Error('Invalid wake correspondence coordinates.');
  // These are cumulative physical tube-mass fractions, not equal-index weights.
  if (!Array.isArray(massFractions) || massFractions.length !== tubes.length
    || massFractions.some((row, g) => !Array.isArray(row) || row.length !== tubes[g] + 1
      || row[0] !== 0 || row.at(-1) !== 1 || row.some((v, j) =>
        !Number.isFinite(v) || (j > 0 && v <= row[j - 1]))))
    throw new Error('Supply strictly increasing cumulative passage mass fractions from zero to one.');

  let coordinateScale = 0;
  for (const group of nodes) for (const row of group) for (const p of row)
    coordinateScale = Math.max(coordinateScale, Math.abs(p.x), Math.abs(p.y));
  // Roundoff only: retain zero/slightly negative computed gaps without clipping.
  // Grid/physical acceptance remains with the caller, including positive gaps
  // when required by the coupled case. No division by a wake gap is performed.
  const roundoffTolerance = 128 * Number.EPSILON * coordinateScale;
  const fail = (body, station, message) => {
    const error = new Error(`Wake correspondence body ${body}, station ${station}: ${message}`);
    error.code = 'WAKE_CORRESPONDENCE_GEOMETRY'; error.body = body; error.station = station;
    throw error;
  };
  const center = (grid, b, i) => ({
    x: .5 * (grid[b][i].at(-1).x + grid[b + 1][i][0].x),
    y: .5 * (grid[b][i].at(-1).y + grid[b + 1][i][0].y),
  });
  const separation = (grid, b, i) => ({
    x: grid[b + 1][i][0].x - grid[b][i].at(-1).x,
    y: grid[b + 1][i][0].y - grid[b][i].at(-1).y,
  });
  const gap = (grid, b, i) => {
    try {
      const indices = [i - 1, i, Math.min(nx, i + 1)];
      const { apply, ...value } = streamtubeWakeGap(indices.map(k => grid[b][k].at(-1)),
        indices.map(k => grid[b + 1][k][0]));
      if (![value.gap, value.tangentialOffset].every(Number.isFinite)) fail(b, i, 'Nonfinite bank separation.');
      if (value.gap < -roundoffTolerance) fail(b, i, 'Negative normal wake gap.');
      return value;
    } catch (error) {
      if (error.code === 'WAKE_CORRESPONDENCE_GEOMETRY') throw error;
      fail(b, i, error.message);
    }
  };
  // Independent banks can coincide and even share an input point at zero gap.
  // Split point aliases before adding opposite bank corrections.
  const moved = structuredClone(nodes).map(group => group.map(row => row.map(p => ({ ...p })))), diagnostics = {
    initialGuessOnly: true, equationsChanged: false, geometryAcceptanceRequired: true,
    coordinateRule: 'equal-center-edge-projected-bank-advance',
    tangentialMismatchDefinition: 'center-edge projection of successive bank-separation difference',
    roundoffTolerance, maximumGapChange: 0, maximumCenterChange: 0,
    maximumTangentialMismatch: 0, maximumNodeDisplacement: 0,
    minimumFrameProjection: Infinity, minimumBankForwardAdvance: Infinity,
    maximumForwardAdvanceMismatch: 0,
    minimumInputGap: Infinity, minimumOutputGap: Infinity, bodies: [],
  };
  for (let b = 0; b < elements; b++) {
    const te = bodies[b].trailingIndex, a = center(nodes, b, te), c = center(nodes, b, te + 1);
    const length = Math.hypot(c.x - a.x, c.y - a.y);
    if (!(length > 0) || !Number.isFinite(length)) fail(b, te, 'Degenerate first wake-center segment.');
    const tangent = { x: (c.x - a.x) / length, y: (c.y - a.y) / length };
    const lower = nodes[b][te].at(-1), upper = nodes[b + 1][te][0];
    const tau = (upper.x - lower.x) * tangent.x + (upper.y - lower.y) * tangent.y;
    const teGap = -(upper.x - lower.x) * tangent.y + (upper.y - lower.y) * tangent.x;
    if (![tau, teGap].every(Number.isFinite)) fail(b, te, 'Nonfinite TE bank separation.');
    if (teGap < -roundoffTolerance) fail(b, te, 'Negative normal TE gap.');
    const body = { body: b, trailingIndex: te, tangentialOffset: tau, teNormalGap: teGap,
      stations: nx - te, maximumBankDisplacement: 0, firstBefore: gap(nodes, b, te + 1) };
    let previous = separation(nodes, b, te);
    for (let i = te + 1; i <= nx; i++) {
      // All frames and gaps come from the original centers; corrections to the
      // other boundary of a shared passage cannot change this body's pairing.
      const before = gap(nodes, b, i), left = center(nodes, b, i - 1), right = center(nodes, b, i);
      const advance = Math.hypot(right.x - left.x, right.y - left.y);
      if (!(advance > 0) || !Number.isFinite(advance)) fail(b, i, 'Degenerate wake-center segment.');
      const ex = (right.x - left.x) / advance, ey = (right.y - left.y) / advance;
      const denominator = ex * before.tangent.x + ey * before.tangent.y;
      // Reject unresolved or backward frames, without changing the prescribed
      // gap, substituting a tangent, or manufacturing a minimum forward step.
      if (!(denominator > 64 * Number.EPSILON)) fail(b, i, 'Singular or backward wake pairing frame.');
      const pairedTau = (ex * previous.x + ey * previous.y
        - before.gap * (ex * before.normal.x + ey * before.normal.y)) / denominator;
      if (!Number.isFinite(pairedTau)) fail(b, i, 'Nonfinite paired bank separation.');
      const delta = pairedTau - before.tangentialOffset;
      diagnostics.minimumInputGap = Math.min(diagnostics.minimumInputGap, before.gap);
      diagnostics.minimumFrameProjection = Math.min(diagnostics.minimumFrameProjection, denominator);
      for (const [g, j, sign] of [[b, tubes[b], -1], [b + 1, 0, 1]]) {
        const dx = .5 * sign * delta * before.tangent.x, dy = .5 * sign * delta * before.tangent.y;
        // Avoid changing signed zeros in an exact no-op.
        if (dx !== 0) moved[g][i][j].x += dx;
        if (dy !== 0) moved[g][i][j].y += dy;
        body.maximumBankDisplacement = Math.max(body.maximumBankDisplacement, Math.hypot(dx, dy));
      }
      previous = separation(moved, b, i);
    }
    diagnostics.bodies.push(body);
  }
  for (let g = 0; g < tubes.length; g++) for (let j = 1; j < tubes[g]; j++) {
    const eta = massFractions[g][j];
    for (let i = 0; i <= nx; i++) for (const key of ['x', 'y']) {
      const delta = (1 - eta) * (moved[g][i][0][key] - nodes[g][i][0][key])
        + eta * (moved[g][i].at(-1)[key] - nodes[g][i].at(-1)[key]);
      if (delta !== 0) moved[g][i][j][key] += delta;
    }
  }
  for (let g = 0; g < tubes.length; g++) for (let i = 0; i <= nx; i++) for (let j = 0; j <= tubes[g]; j++) {
    const p = moved[g][i][j], q = nodes[g][i][j];
    if (!Number.isFinite(p.x) || !Number.isFinite(p.y)) throw new Error('Nonfinite wake correspondence result.');
    diagnostics.maximumNodeDisplacement = Math.max(diagnostics.maximumNodeDisplacement, Math.hypot(p.x - q.x, p.y - q.y));
  }
  for (let b = 0; b < elements; b++) {
    const body = diagnostics.bodies[b];
    for (let i = bodies[b].trailingIndex + 1; i <= nx; i++) {
      const before = gap(nodes, b, i), after = gap(moved, b, i), a = center(nodes, b, i), c = center(moved, b, i);
      diagnostics.maximumGapChange = Math.max(diagnostics.maximumGapChange, Math.abs(after.gap - before.gap));
      diagnostics.maximumCenterChange = Math.max(diagnostics.maximumCenterChange, Math.hypot(a.x - c.x, a.y - c.y));
      const left = center(nodes, b, i - 1), advance = Math.hypot(a.x - left.x, a.y - left.y);
      const ex = (a.x - left.x) / advance, ey = (a.y - left.y) / advance;
      const previous = separation(moved, b, i - 1), current = separation(moved, b, i);
      diagnostics.maximumTangentialMismatch = Math.max(diagnostics.maximumTangentialMismatch,
        Math.abs(ex * (current.x - previous.x) + ey * (current.y - previous.y)));
      for (const [g, j] of [[b, tubes[b]], [b + 1, 0]]) {
        const forward = ex * (moved[g][i][j].x - moved[g][i - 1][j].x)
          + ey * (moved[g][i][j].y - moved[g][i - 1][j].y);
        diagnostics.minimumBankForwardAdvance = Math.min(diagnostics.minimumBankForwardAdvance, forward);
        diagnostics.maximumForwardAdvanceMismatch = Math.max(diagnostics.maximumForwardAdvanceMismatch, Math.abs(forward - advance));
      }
      diagnostics.minimumOutputGap = Math.min(diagnostics.minimumOutputGap, after.gap);
    }
    body.firstAfter = gap(moved, b, bodies[b].trailingIndex + 1);
  }
  return { nodes: moved, diagnostics };
}
