// SPDX-License-Identifier: GPL-2.0-or-later
import { solveSparseDirectAlignedMany } from '../numerics/klu.js';
import { projectHalfspaces } from '../numerics/halfspace-projection.js';

const cross = (a, b) => a.x * b.y - a.y * b.x;
const subtract = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });
const field = (nodes, fn) => nodes.map((group, g) => group.map((row, i) => row.map((p, j) => fn(p, g, i, j))));

function corners(nodes, visit) {
  for (let g = 0; g < nodes.length; g++) for (let i = 0; i < nodes[g].length - 1; i++)
    for (let j = 0; j < nodes[g][i].length - 1; j++) {
      const indices = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]];
      for (let corner = 0; corner < 4; corner++) {
        const picked = [0, 1, 2].map(k => indices[(corner + k) % 4]);
        const [a, b, c] = picked.map(([ii, jj]) => nodes[g][ii][jj]);
        const u = subtract(b, a), v = subtract(c, b);
        const derivative = movement => {
          const [da, db, dc] = picked.map(([ii, jj]) => movement[g][ii][jj]);
          return cross(subtract(db, da), v) + cross(u, subtract(dc, db));
        };
        visit({ g, i, j, corner, picked, value: cross(u, v), derivative });
      }
    }
}

// Add a few tangential coordinate freedoms when the normal-coordinate
// Newton direction points through a convexity boundary. For a geometric
// basis B, solve J z = R_geometry B. The combined motion B - geometry(z)
// then changes the grid while cancelling its first-order equation defect.
// A small convex projection chooses a combination of these freedoms. This
// is only a proposal: the caller checks the complete nonlinear gas state,
// every grid corner, and maintained residual decrease before accepting it.
export function tangentialNewton(system, state, value, matrix, direction, order) {
  const nodes = value.nodes, epsilon = 1e-7;
  const geometry = system.geometryDerivatives(state);
  const motion = d => field(nodes, (_, g, i, j) => {
    let x = 0, y = 0;
    for (const [column, derivative] of geometry[g][i][j]) {
      x += derivative.x * d[column]; y += derivative.y * d[column];
    }
    return { x, y };
  });
  const newton = motion(direction), candidates = [];
  corners(nodes, corner => {
    const d = corner.derivative(newton);
    if (d < 0 && corner.value / -d < .02) candidates.push({ ...corner, ratio: corner.value / -d });
  });
  candidates.sort((a, b) => a.ratio - b.ratio);
  const points = new Map();
  for (const c of candidates) for (const [i, j] of c.picked) {
    const node = system.layout.nodes[c.g][i][j];
    if (i === 0 || i === system.layout.nx || node.column === null) continue;
    if (!points.has(node.column) && points.size < 16) points.set(node.column, { group: c.g, i, j, node });
  }
  if (!points.size) return null;
  const seeds = [];
  for (const { group, i, j, node } of points.values()) {
    const delta = field(nodes, () => ({ x: 0, y: 0 }));
    const prev = nodes[group][i - 1][j], next = nodes[group][i + 1][j], p = nodes[group][i][j];
    const dx = next.x - prev.x, dy = next.y - prev.y, length = Math.hypot(dx, dy);
    const scale = .25 * Math.min(Math.hypot(p.x - prev.x, p.y - prev.y), Math.hypot(p.x - next.x, p.y - next.y));
    const shift = { x: scale * dx / length, y: scale * dy / length };
    delta[group][i][j] = shift;
    // The two banks of an inviscid cut have one coordinate unknown.
    if (node.kind === 'cut' && !node.side) {
      delta[node.body][i][system.layout.tubes[node.body]] = shift;
      delta[node.body + 1][i][0] = shift;
    }
    let perturbed;
    try {
      const x = system.adoptGeometry(state, field(nodes, (p, g, i, j) => ({
        x: p.x + epsilon * delta[g][i][j].x, y: p.y + epsilon * delta[g][i][j].y })));
      perturbed = system.evaluate(x);
    } finally { system.adoptGeometry(state, nodes); }
    seeds.push({ delta, rhs: value.residual.map((r, k) => (perturbed.residual[k] - r) / epsilon), cell: { group, i, j } });
  }
  const linear = solveSparseDirectAlignedMany(matrix, seeds.map(s => s.rhs), order);
  const bases = seeds.map(({ delta, cell }, k) => {
    const response = linear[k].x, normal = motion(response);
    return { delta, cell, response, gauge: field(delta, (p, g, i, j) => ({
      x: p.x - normal[g][i][j].x, y: p.y - normal[g][i][j].y })) };
  });
  const constraints = [];
  corners(nodes, c => constraints.push({ value: c.value, newton: c.derivative(newton),
    gradient: new Map(bases.map((b, k) => [k, c.derivative(b.gauge)])) }));
  for (const step of [.1, .03, .01, .003]) {
    // Move away from an imminent contact, rather than remaining tangent to
    // a nearly flat cell. Other corners retain 10% of their initial area at
    // this trial scale. The exact nonlinear limiter still owns acceptance.
    const halfspaces = constraints.map(c => ({ gradient: c.gradient,
      lower: Math.max(-.9 * c.value / step,
        c.newton < 0 && c.value < -.01 * c.newton ? -.25 * c.newton : -Infinity) - c.newton }));
    let projection;
    try { projection = projectHalfspaces(new Float64Array(bases.length), halfspaces,
      { maximumNorm: 1e3, maxSweeps: 1000 }); }
    catch { continue; }
    if (!projection.converged) continue;
    const weights = projection.point;
    const corrected = Float64Array.from(direction, (v, k) => v - bases.reduce((sum, b, m) => sum + weights[m] * b.response[k], 0));
    const delta = field(nodes, (_, g, i, j) => ({
      x: bases.reduce((sum, b, m) => sum + weights[m] * b.delta[g][i][j].x, 0),
      y: bases.reduce((sum, b, m) => sum + weights[m] * b.delta[g][i][j].y, 0) }));
    return { direction: corrected, linear,
      mix: (trial, a) => field(trial, (p, g, i, j) => ({ x: p.x + a * delta[g][i][j].x, y: p.y + a * delta[g][i][j].y })),
      diagnostics: { method: 'compensated-tangential-coordinates', cells: bases.map(b => b.cell), weights: [...weights], step } };
  }
  return { linear, diagnostics: { method: 'compensated-tangential-coordinates', accepted: false, reason: 'No feasible bounded coordinate projection.' } };
}
