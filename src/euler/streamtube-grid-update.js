// SPDX-License-Identifier: GPL-2.0-or-later
// Geometric predictor for a common Newton step. The caller must recompute
// the complete physical state (including curved walls and BL displacement)
// at the smaller step; interpolated coordinates are not an Euler update.
import { limitStreamtubeGridStep } from '../geometry/streamtube-convex-step.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';
import { initializeStreamtubeWakeCorrespondence } from './streamtube-wake-correspondence.js';

// The displayed mesh merges shared cuts/endpoints. Check that exact
// connectivity and predicate too, including its floating-point arithmetic.
export function requireConvexPublishedGrid(system, nodes) {
  const quality = streamtubeMeshSnapshot({ system, nodes }).quality;
  if (!quality.valid) {
    const cell = quality.invalidCells[0], { nx, tubes } = system.layout;
    let local = cell, group = 0;
    while (group < tubes.length - 1 && local >= nx * tubes[group]) local -= nx * tubes[group++];
    const i = Math.floor(local / tubes[group]), tube = local % tubes[group];
    throw Object.assign(new Error(`Grid update would publish a nonconvex cell at passage ${group}, interval ${i}, tube ${tube}.`), {
      code: 'streamtube-grid-nonconvex', diagnostics: { cell: { group, i, tube }, quality },
    });
  }
  return quality;
}

export function requireConvexGridUpdate(previous, proposed) {
  const limit = limitStreamtubeGridStep(previous, proposed);
  if (limit.limited) {
    const cell = limit.limiter;
    throw Object.assign(new Error(`Grid update would lose convexity at passage ${cell.group}, interval ${cell.i}, tube ${cell.tube}, corner ${cell.corner}.`), {
      code: 'streamtube-grid-step', stepFraction: limit.step,
      diagnostics: { cell: { group: cell.group, i: cell.i, tube: cell.tube, corner: cell.corner },
        stepFraction: limit.step, boundaryStep: cell.boundaryStep },
    });
  }
  return limit;
}

// A displacement update can reverse one first-wake bank before its paired
// bank. Try a coordinate correction on this uncommitted trial, keeping its
// mean wake and normal gap. This is a grid globalization operation, not an
// extra wake equation. Preserve ordinary backtracking when it can resolve the
// contact within its search budget; changing that path can prevent convergence
// even if the attempted wake correction itself is rejected. The caller still
// owns physical and residual acceptance.
export function prepareConvexWakeGridUpdate(previous, proposed, layout, allocation,
  minimumUncorrectedFraction = 2 ** -12) {
  if (!Number.isFinite(minimumUncorrectedFraction) || minimumUncorrectedFraction <= 0 || minimumUncorrectedFraction > 1)
    throw new Error('Invalid minimum uncorrected grid-step fraction.');
  let rejected;
  try { requireConvexGridUpdate(previous, proposed); return { nodes: proposed }; }
  catch (error) { rejected = error; }
  const cell = rejected.diagnostics?.cell;
  const firstWake = layout.independentWakeBanks && rejected.code === 'streamtube-grid-step' && cell
    && layout.bodies.some((body, b) => cell.i === body.trailingIndex
      && (cell.group === b && cell.tube === layout.tubes[b] - 1 || cell.group === b + 1 && cell.tube === 0));
  if (!firstWake || rejected.stepFraction >= minimumUncorrectedFraction) throw rejected;

  const massFractions = allocation.groups.map(group => {
    const total = group.reduce((sum, tube) => sum + tube.massFlow, 0), fractions = [0];
    for (const tube of group) fractions.push(fractions.at(-1) + tube.massFlow / total);
    fractions[fractions.length - 1] = 1;
    return fractions;
  });
  let pairing, pairingFailure;
  try { pairing = initializeStreamtubeWakeCorrespondence({ nodes: proposed, layout, massFractions }); }
  catch (error) {
    if (error.code !== 'WAKE_CORRESPONDENCE_GEOMETRY') throw error;
    pairingFailure = { code: error.code, message: error.message };
  }
  if (pairing) for (const fraction of [1 / 16, 1 / 8, 1 / 4, 1 / 2, 1]) {
    const nodes = proposed.map((group, g) => group.map((row, i) => row.map((p, j) => ({
      x: p.x + fraction * (pairing.nodes[g][i][j].x - p.x),
      y: p.y + fraction * (pairing.nodes[g][i][j].y - p.y),
    }))));
    try { requireConvexGridUpdate(previous, nodes); }
    catch (error) {
      if (error.code !== 'streamtube-grid-step') throw error;
      continue;
    }
    return { nodes, correction: { method: 'paired-wake-trial', fraction, cell,
      equationsChanged: false, physicalAcceptanceRequired: true } };
  }
  // The raw-grid boundary is not a bound on the corrected-coordinate path:
  // a full correction can move the contact to the opposite wake bank. Halve
  // Newton and reconstruct/recheck the entire trial instead of jumping to
  // the raw boundary (which can be smaller by many orders of magnitude).
  throw Object.assign(new Error(rejected.message), {
    code: rejected.code, stepFraction: .5,
    diagnostics: { ...rejected.diagnostics, rawGridStepFraction: rejected.stepFraction,
      coordinateRepairBacktrack: .5, ...(pairingFailure ? { pairingFailure } : {}) },
  });
}
