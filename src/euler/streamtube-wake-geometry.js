// SPDX-License-Identifier: GPL-2.0-or-later

// Nonlinear normal-gap condition for independently moving wake banks.
// D = (upper-lower) dot n, with n normal to the centered mean-bank secant.
// No division by D: the value and derivative remain regular at zero gap.
export function streamtubeWakeGap(lower, upper) {
  if (![lower, upper].every(row => Array.isArray(row) && row.length === 3
    && row.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y)))) throw new Error('Invalid wake-gap stencil.');
  const sx = .5 * ((lower[2].x - lower[0].x) + (upper[2].x - upper[0].x));
  const sy = .5 * ((lower[2].y - lower[0].y) + (upper[2].y - upper[0].y));
  const length = Math.hypot(sx, sy);
  if (!(length > 0) || !Number.isFinite(length)) throw new Error('Degenerate wake centerline secant.');
  const tangent = { x: sx / length, y: sy / length }, normal = { x: -sy / length, y: sx / length };
  const dx = upper[1].x - lower[1].x, dy = upper[1].y - lower[1].y;
  return { gap: dx * normal.x + dy * normal.y, tangentialOffset: dx * tangent.x + dy * tangent.y, normal, tangent, length,
    apply: (dlower, dupper) => {
      const dsx = .5 * ((dlower[2].x - dlower[0].x) + (dupper[2].x - dupper[0].x));
      const dsy = .5 * ((dlower[2].y - dlower[0].y) + (dupper[2].y - dupper[0].y));
      const along = tangent.x * dsx + tangent.y * dsy;
      const dnx = -(dsy - tangent.y * along) / length, dny = (dsx - tangent.x * along) / length;
      return (dupper[1].x - dlower[1].x) * normal.x + (dupper[1].y - dlower[1].y) * normal.y + dx * dnx + dy * dny;
    } };
}

// Initial-guess geometry after a complete BL profile change. This is not a
// residual/Jacobian operator or an update to the independent wake chart.
// Widths are total physical lengths from bl.thicknesses(), including any
// finite-base contribution already contained in the packed wake delta-star.
// This preserves a gap residual; it does not impose a solution or certify cells.

export function incrementIndependentWakeWidths({ layout, nodes, beforeWidths, afterWidths }) {
  const require = (ok, message) => { if (!ok) throw new Error(message); };
  require(layout?.independentWakeBanks === true && Number.isInteger(layout.nx) && layout.nx >= 2
    && Array.isArray(layout.bodies) && layout.bodies.length > 0 && Array.isArray(layout.tubes)
    && layout.tubes.length === layout.bodies.length + 1,
  'Wake-width increments require an independent-bank layout.');
  const { nx, bodies, tubes } = layout;
  require(tubes.every(n => Number.isInteger(n) && n >= 1) && Array.isArray(nodes) && nodes.length === tubes.length
    && nodes.every((group, g) => Array.isArray(group) && group.length === nx + 1
      && group.every(row => Array.isArray(row) && row.length === tubes[g] + 1
        && row.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y)))),
  'Wake-width increments require matching finite physical coordinates.');
  require([beforeWidths, afterWidths].every(widths => Array.isArray(widths) && widths.length === bodies.length),
    'Supply previous and next total wake widths for every body.');
  for (const [b, body] of bodies.entries()) {
    require(Number.isInteger(body.trailingIndex) && body.trailingIndex >= 0 && body.trailingIndex < nx,
      'Invalid trailing-edge index for wake-width increments.');
    require([beforeWidths[b], afterWidths[b]].every(row => Array.isArray(row) && row.length === nx - body.trailingIndex
      && row.every(d => Number.isFinite(d) && d >= 0)),
    'Wake-width increments require nonnegative finite physical total widths.');
  }
  const diagnostics = { active: false, changedNodes: 0, maximumOffset: 0,
    initialGuessOnly: true, physicalTotalWidth: true, equationsChanged: false,
    zeroesGapResidual: false, geometryAdmissibilityChecked: false };
  if (afterWidths.every((row, b) => row.every((d, k) => d === beforeWidths[b][k]))) return { nodes, diagnostics };
  const moved = nodes.map(group => group.map(row => row.map(p => ({ ...p }))));
  for (const [body, spec] of bodies.entries()) for (let i = spec.trailingIndex + 1; i <= nx; i++) {
    const k = i - spec.trailingIndex - 1, delta = afterWidths[body][k] - beforeWidths[body][k];
    if (delta === 0) continue;
    // Match the actual nonlinear gap stencil and its arithmetic. Read every
    // normal from the unmodified grid, including at TE+1 and at the outlet.
    const indices = [i - 1, i, Math.min(i + 1, nx)];
    const lower = indices.map(j => nodes[body][j].at(-1)), upper = indices.map(j => nodes[body + 1][j][0]);
    const { normal } = streamtubeWakeGap(lower, upper);
    for (const [g, j, sign] of [[body, tubes[body], -1], [body + 1, 0, 1]]) {
      const p = nodes[g][i][j];
      moved[g][i][j] = { ...p, x: p.x + .5 * sign * delta * normal.x, y: p.y + .5 * sign * delta * normal.y };
      require(Number.isFinite(moved[g][i][j].x) && Number.isFinite(moved[g][i][j].y), 'Nonfinite wake-width increment.');
      diagnostics.changedNodes++;
    }
    diagnostics.maximumOffset = Math.max(diagnostics.maximumOffset, .5 * Math.abs(delta));
  }
  diagnostics.active = diagnostics.changedNodes > 0;
  return { nodes: moved, diagnostics };
}
