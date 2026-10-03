// SPDX-License-Identifier: GPL-2.0-or-later
// Inlet tangential adjustment (UPDATE pp.262–263) and one-pass interior
// DEKINK (p.287), reconstructed from their mathematical operations. An inlet
// tangent movement that would fold the grid may use a finite curve-following
// reparameterization when that complete alternative stays convex.
import { limitStreamtubeGridStep, interpolateStreamtubeGridNodes } from './streamtube-convex-step.js';
const distance = (a, b) => Math.hypot(b.x - a.x, b.y - a.y);
const clone = nodes => nodes.map(group => group.map(row => row.map(p => ({ ...p }))));
const validate = nodes => {
  if (!Array.isArray(nodes) || !nodes.length || !nodes.every(group => Array.isArray(group) && group.length >= 3
    && group.every(row => Array.isArray(row) && row.length === group[0].length && row.length >= 2
      && row.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y))))) throw new Error('Invalid grid maintenance coordinates.');
};
const inletArc = (nodes, body, leadingIndex) => {
  if (!Number.isInteger(leadingIndex) || leadingIndex < 1 || leadingIndex >= nodes[body].length
    || nodes[body + 1]?.length !== nodes[body].length) throw new Error('Invalid inlet cut range.');
  const arc = [0];
  for (let i = 1; i <= leadingIndex; i++) {
    const length = distance(nodes[body][i - 1].at(-1), nodes[body][i].at(-1));
    if (!(length > 0)) throw new Error('Degenerate inlet cut segment.');
    arc.push(arc.at(-1) + length);
  }
  for (let i = 0; i <= leadingIndex; i++) if (distance(nodes[body][i].at(-1), nodes[body + 1][i][0]) > 1e-12 * arc.at(-1))
    throw new Error('Disconnected inlet cut banks.');
  return arc;
};

export function captureStreamtubeInletFractions(nodes, bodies) {
  validate(nodes);
  if (!Array.isArray(bodies) || bodies.length !== nodes.length - 1) throw new Error('Invalid inlet body count.');
  return bodies.map((b, body) => { const arc = inletArc(nodes, body, b.leadingIndex); return arc.map(v => v / arc.at(-1)); });
}

function inletCoordinates(nodes, bodies, fractions) {
  validate(nodes);
  if (!Array.isArray(bodies) || bodies.length !== nodes.length - 1 || !Array.isArray(fractions) || fractions.length !== bodies.length)
    throw new Error('Invalid inlet adjustment data.');
  return bodies.map((b, body) => {
    const arc = inletArc(nodes, body, b.leadingIndex), target = fractions[body];
    if (!Array.isArray(target) || target.length !== arc.length || target[0] !== 0 || target.at(-1) !== 1
      || !target.every(Number.isFinite) || target.some((v, i) => i && v <= target[i - 1])) throw new Error('Invalid stored inlet fractions.');
    return { body, leadingIndex: b.leadingIndex, arc, target };
  });
}

function inletArcError(nodes, coordinates) {
  let maximum = 0;
  for (const { body, leadingIndex, target } of coordinates) {
    const arc = inletArc(nodes, body, leadingIndex);
    for (let i = 1; i < leadingIndex; i++) maximum = Math.max(maximum, Math.abs(arc.at(-1) * target[i] - arc[i]));
  }
  return maximum;
}

function sampleInletCoordinates(nodes, coordinates) {
  const moved = clone(nodes); let maxDisplacement = 0, maxArcErrorBefore = 0;
  for (const { body, leadingIndex, arc, target } of coordinates) {
    let segment = 0;
    for (let i = 1; i < leadingIndex; i++) {
      const desired = arc.at(-1) * target[i];
      while (segment + 1 < leadingIndex && arc[segment + 1] <= desired) segment++;
      const fraction = (desired - arc[segment]) / (arc[segment + 1] - arc[segment]);
      const a = nodes[body][segment].at(-1), c = nodes[body][segment + 1].at(-1);
      const q = { x: a.x + fraction * (c.x - a.x), y: a.y + fraction * (c.y - a.y) };
      moved[body][i][moved[body][i].length - 1] = q;
      moved[body + 1][i][0] = q;
      maxDisplacement = Math.max(maxDisplacement, distance(nodes[body][i].at(-1), q));
      maxArcErrorBefore = Math.max(maxArcErrorBefore, Math.abs(desired - arc[i]));
    }
  }
  return { nodes: moved, maxDisplacement, maxArcErrorBefore,
    maxArcErrorAfter: inletArcError(moved, coordinates), reparameterization: { method: 'polyline-arclength' } };
}

// This is an unaccepted geometric proposal, including when the source is an
// invalid raw Newton trial. Its caller must check the complete movement from
// the last accepted grid, physical admissibility and residual decrease.
export function reparameterizeStreamtubeInlets(nodes, bodies, fractions) {
  return sampleInletCoordinates(nodes, inletCoordinates(nodes, bodies, fractions));
}

export function adjustStreamtubeInlets(nodes, bodies, fractions, { preserveConvexity = false } = {}) {
  const coordinates = inletCoordinates(nodes, bodies, fractions);
  let moved = clone(nodes), maxDisplacement = 0, maxArcErrorBefore = 0;
  for (const { body, leadingIndex, arc, target } of coordinates) {
    for (let i = 1; i < leadingIndex; i++) {
      const a = nodes[body][i - 1].at(-1), p = nodes[body][i].at(-1), c = nodes[body][i + 1].at(-1);
      const length = distance(a, c), shift = arc.at(-1) * target[i] - arc[i];
      if (!(length > 0)) throw new Error('Degenerate inlet adjustment tangent.');
      const q = { x: p.x + shift * (c.x - a.x) / length, y: p.y + shift * (c.y - a.y) / length };
      moved[body][i][moved[body][i].length - 1] = q; moved[body + 1][i][0] = q;
      maxDisplacement = Math.max(maxDisplacement, Math.abs(shift)); maxArcErrorBefore = Math.max(maxArcErrorBefore, Math.abs(shift));
    }
  }
  let limit = preserveConvexity ? limitStreamtubeGridStep(nodes, moved) : null, reparameterization;
  if (limit?.limited) {
    // Keep successful UPDATE tangent movements unchanged. Near a curved cut,
    // that approximation can leave the curve and flatten an adjacent cell.
    // Test the same stored fractions on the original polyline before clipping
    // the correction toward a degenerate grid. Failed alternatives retain the
    // original bounded movement; no physical or residual check is bypassed.
    const alongCurve = sampleInletCoordinates(nodes, coordinates).nodes;
    const candidateLimit = limitStreamtubeGridStep(nodes, alongCurve);
    if (!candidateLimit.limited) {
      reparameterization = { method: 'polyline-arclength', replacedTangentScale: limit.step };
      moved = alongCurve; limit = candidateLimit;
    }
  }
  if (limit?.limited) {
    moved = interpolateStreamtubeGridNodes(nodes, moved, limit.step);
    maxDisplacement *= limit.step;
  }
  return { nodes: moved, maxDisplacement, maxArcErrorBefore, maxArcErrorAfter: inletArcError(moved, coordinates),
    ...(reparameterization ? { reparameterization } : {}),
    ...(limit?.limited ? { convexity: { scale: limit.step, limiter: limit.limiter } } : {}) };
}

export function dekinkStreamtubeInteriors(nodes, { preserveConvexity = false } = {}) {
  validate(nodes);
  const repairs = []; let moved = clone(nodes), maxDisplacement = 0;
  nodes.forEach((group, g) => {
    for (let j = 1; j < group[0].length - 1; j++) for (let i = 1; i < group.length - 1; i++) {
      const a = group[i - 1][j], b = group[i][j], c = group[i + 1][j];
      const dot = (b.x - a.x) * (c.x - b.x) + (b.y - a.y) * (c.y - b.y);
      if (dot < 0) {
        const p = { x: .5 * (a.x + c.x), y: .5 * (a.y + c.y) };
        moved[g][i][j] = p; maxDisplacement = Math.max(maxDisplacement, distance(p, b)); repairs.push({ group: g, i, j, dot });
      }
    }
  });
  const limit = preserveConvexity ? limitStreamtubeGridStep(nodes, moved) : null;
  if (limit?.limited) {
    moved = interpolateStreamtubeGridNodes(nodes, moved, limit.step);
    maxDisplacement *= limit.step;
  }
  return { nodes: moved, repairs, maxDisplacement,
    ...(limit?.limited ? { convexity: { scale: limit.step, limiter: limit.limiter } } : {}) };
}
