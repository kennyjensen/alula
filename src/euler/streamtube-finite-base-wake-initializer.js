// SPDX-License-Identifier: GPL-2.0-or-later
// Initializer-only recovery. The solid TE vector supplies a declared
// tangential pairing; prescribed native normal width and both endpoints stay.
import { createStreamtubeBodySystem } from './streamtube-body.js';
import { createStreamtubeDisplacement } from './streamtube-displacement.js';
import { streamtubeGridConvexity } from '../geometry/streamtube-convex-step.js';

// A panel trace is not a boundary of the prescribed constant-width base
// wake. Opening its banks can overtake a closely spaced interior trace far
// downstream, including at the outlet that SLOR subsequently holds fixed.
// Reconstruct only such wake interiors as a Coons strip: retain the TE
// crossline and both longitudinal boundaries, grading its departure from
// linear normal spacing to zero at the outlet. This is initial geometry,
// not a change to the wake width or to any flow equation.
export function recoverFiniteBaseWakeInterior({ input, system, nodes }, { admissibleNode } = {}) {
  if (!system.inviscidBaseWake || input.displacement !== undefined
    || (input.wakeGeometry ?? 'centerline') !== 'centerline') return null;
  const originalQuality = streamtubeGridConvexity(nodes);
  if (originalQuality.valid) return null;
  const groups = new Map();
  for (const cell of originalQuality.invalidCells) {
    const body = input.bodies.findIndex((b, k) => b.trailingEdge?.kind === 'finite-base'
      && cell.i > b.trailingIndex && (cell.group === k || cell.group === k + 1));
    if (body < 0) return null;
    // Wake opening can overtake interior traces as well as the first bank
    // tube. The other boundary may still be a solid surface: preserve it,
    // and let the original fluid-region and convexity gates judge the strip.
    for (const g of [body, body + 1]) if (nodes[g])
      groups.set(g, Math.min(groups.get(g) ?? Infinity, input.bodies[body].trailingIndex));
  }
  if (typeof admissibleNode !== 'function') throw new Error('Finite-base wake initialization requires the original solid-region guard.');
  const next = structuredClone(nodes);
  const report = { attempted: true, accepted: false, method: 'wake-transfinite-interior', originalQuality,
    groups: [...groups].map(([group, firstCrossline]) => ({ group, firstCrossline })),
    longitudinalBoundariesChanged: false, upstreamCrosslineChanged: false,
    outletInteriorRebuilt: true, panelTracesPreserved: false };
  for (const [g, start] of groups) {
    const rows = nodes[g], n = rows[0].length - 1, weights = input.weights[g];
    const total = weights.reduce((sum, w) => sum + w, 0), eta = [0], arc = [0];
    for (const w of weights) eta.push(eta.at(-1) + w / total);
    eta[n] = 1;
    for (let i = start + 1; i < rows.length; i++) {
      const length = j => Math.hypot(rows[i][j].x - rows[i - 1][j].x, rows[i][j].y - rows[i - 1][j].y);
      arc.push(arc.at(-1) + .5 * (length(0) + length(n)));
    }
    for (let i = start + 1; i < rows.length; i++) for (let j = 1; j < n; j++) {
      const fraction = eta[j], remaining = 1 - arc[i - start] / arc.at(-1);
      const p = next[g][i][j];
      for (const key of ['x', 'y']) p[key] = (1 - fraction) * rows[i][0][key] + fraction * rows[i][n][key]
        + remaining * (rows[start][j][key] - (1 - fraction) * rows[start][0][key] - fraction * rows[start][n][key]);
      if (!admissibleNode(p)) return { nodes, report: { ...report, reason: 'Wake interior leaves the original fluid region.', node: { g, i, j } } };
    }
  }
  const quality = streamtubeGridConvexity(next);
  if (!quality.valid) return { nodes, report: { ...report, quality, reason: 'Wake interpolation did not produce a convex initial grid.' } };
  return { nodes: next, report: { ...report, accepted: true, quality, fluidRegionChecked: true } };
}

export function finiteBaseWakeRecoveryPlan(input, system, nodes) {
  if (!system.inviscidBaseWake || input.displacement !== undefined
    || (input.wakeGeometry ?? 'centerline') !== 'centerline') return null;
  const quality = streamtubeGridConvexity(nodes);
  if (quality.valid) return null;
  const bodies = new Set();
  for (const cell of quality.invalidCells) {
    const body = input.bodies.findIndex((b, k) => b.trailingEdge?.kind === 'finite-base'
      && b.wakeTangentialReference === undefined && b.trailingIndex === cell.i
      && (cell.group === k && cell.tube === system.layout.tubes[k] - 1
        || cell.group === k + 1 && cell.tube === 0));
    if (body < 0) {
      const downstream = input.bodies.some((b, k) => b.trailingEdge?.kind === 'finite-base'
        && cell.i > b.trailingIndex && (cell.group === k || cell.group === k + 1));
      if (!downstream) return null;
    } else bodies.add(body);
  }
  if (!bodies.size) return null;
  return { code: 'finite-base-first-wake-tangency', bodies: [...bodies], originalQuality: quality };
}

export function recoverFiniteBaseWakeTangency({ input, system, nodes }, { admissibleNode } = {}) {
  const plan = finiteBaseWakeRecoveryPlan(input, system, nodes);
  if (!plan) return null;
  if (typeof admissibleNode !== 'function') throw new Error('Finite-base wake initialization requires the original solid-region guard.');
  const nextInput = structuredClone(input);
  for (const b of plan.bodies) nextInput.bodies[b].wakeTangentialReference = 'material-te';
  const candidate = createStreamtubeBodySystem(nextInput);
  const make = source => createStreamtubeDisplacement({ layout: source.layout, curves: source.curves,
    fractions: source.fractions, thicknesses: source.displacement });
  const bare = make(system).restore(nodes, system.initialStagnation, system.originalNodes, 1e-10 * system.conditions.lengthScale);
  const physical = make(candidate).apply(bare, candidate.initialStagnation).nodes;
  let next = nodes.map((group, g) => {
    const weights = input.weights[g], total = weights.reduce((s, v) => s + v, 0), eta = [0];
    weights.forEach(w => eta.push(eta.at(-1) + w / total)); eta[weights.length] = 1;
    return group.map((row, i) => {
      const a = physical[g][i][0], z = physical[g][i].at(-1);
      const dl = { x: a.x - row[0].x, y: a.y - row[0].y };
      const du = { x: z.x - row.at(-1).x, y: z.y - row.at(-1).y };
      return row.map((p, j) => j === 0 ? a : j === row.length - 1 ? z
        : { x: p.x + (1 - eta[j]) * dl.x + eta[j] * du.x, y: p.y + (1 - eta[j]) * dl.y + eta[j] * du.y });
    });
  });
  const fail = (message, details = {}) => { throw Object.assign(new Error(message), {
    code: 'finite-base-wake-tangency-rejected', diagnostics: { ...plan, ...details } }); };
  for (let g = 0; g < next.length; g++) for (let i = 0; i < next[g].length; i++) for (let j = 0; j < next[g][i].length; j++) {
    const a = nodes[g][i][j], p = next[g][i][j];
    if (system.layout.nodes[g][i][j].kind === 'wall') {
      if (a.x !== p.x || a.y !== p.y) fail('Finite-base wake initialization changed a solid boundary point.', { node: { g, i, j } });
    } else if (!admissibleNode(p)) fail('Finite-base wake initialization left the original admissible fluid region.', { node: { g, i, j } });
  }
  let quality = streamtubeGridConvexity(next), interiorRecovery;
  if (!quality.valid) {
    const interior = recoverFiniteBaseWakeInterior({ input: nextInput, system: candidate, nodes: next }, { admissibleNode });
    if (interior?.report.accepted) {
      next = interior.nodes; quality = streamtubeGridConvexity(next); interiorRecovery = interior.report;
    }
  }
  if (!quality.valid) fail('Finite-base material tangency did not produce a convex initial grid.', { quality });
  const initial = candidate.adoptGeometry(candidate.initial, next), decoded = candidate.decode(initial).nodes;
  let replayError = 0;
  for (let g = 0; g < next.length; g++) for (let i = 0; i < next[g].length; i++) for (let j = 0; j < next[g][i].length; j++)
    replayError = Math.max(replayError, Math.hypot(next[g][i][j].x - decoded[g][i][j].x, next[g][i][j].y - decoded[g][i][j].y));
  if (!(replayError <= 64 * Number.EPSILON * system.conditions.lengthScale)) fail('Finite-base material wake chart did not retain its initialized nodes.', { replayError });
  return { input: nextInput, system: candidate, initial, nodes: decoded, report: { attempted: true, accepted: true,
    ...plan, quality, replayError, ...(interiorRecovery ? { interiorRecovery } : {}), chart: 'material-te', solidTeVector: 'fixed upper solid TE minus lower solid TE',
    normalGapChanged: false, surfaceNodesChanged: false, fluidRegionChecked: true,
    scope: 'Initial finite-base centerline tangential correspondence; subsequent coupled independent banks keep the existing normal-gap equation.' } };
}
