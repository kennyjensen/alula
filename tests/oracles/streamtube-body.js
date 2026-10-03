// Direct physical-boundary flux integration. This file does not import the
// intrinsic body residual, pressure reconstruction, or force integration.
import { directChannelConservation } from './streamtube.js';

const midpoint = (a, b) => ({ x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 });
const difference = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });

export function directBodyConservation(result, bodies, { gamma, lengthScale, flowModel = 'compressible' }) {
  const enthalpyFactor = flowModel === 'incompressible' ? 1 : gamma / (gamma - 1);
  const groups = result.nodes.length, nx = result.nodes[0].length - 1;
  const outer = [0, 0, 0, 0], cutTraction = [0, 0];
  const blocks = result.nodes.map((nodes, g) => directChannelConservation({ nodes,
    sections: result.sections.map(row => row[g]), cells: result.cells.map(row => row[g]) }, gamma, { flowModel }));
  // Cross-sections at the inlet/outlet of the union of the dual volumes.
  for (let g = 0; g < groups; g++) for (let j = 0; j < result.nodes[g][0].length - 1; j++) {
    for (const [i, sign] of [[0, -1], [nx - 1, 1]]) {
      const nodes = result.nodes[g], lower = midpoint(nodes[i][j], nodes[i + 1][j]);
      const upper = midpoint(nodes[i][j + 1], nodes[i + 1][j + 1]);
      const area = difference(upper, lower), normal = { x: sign * area.y, y: -sign * area.x };
      const direction = difference(midpoint(nodes[i + 1][j], nodes[i + 1][j + 1]), midpoint(nodes[i][j], nodes[i][j + 1]));
      const state = result.sections[i][g][j], length = Math.hypot(direction.x, direction.y);
      const vx = state.q * direction.x / length, vy = state.q * direction.y / length;
      const mass = state.rho * (vx * normal.x + vy * normal.y);
      outer[0] += mass; outer[1] += mass * vx + state.p * normal.x; outer[2] += mass * vy + state.p * normal.y;
      outer[3] += mass * (enthalpyFactor * state.p / state.rho + .5 * (vx * vx + vy * vy));
    }
  }
  // Only the two actual outer streamlines; artificial body cuts are internal.
  for (const [g, sign] of [[0, 1], [groups - 1, -1]]) {
    const j = g === 0 ? 0 : result.nodes[g][0].length - 1;
    for (let i = 1; i < nx; i++) {
      const edge = difference(result.nodes[g][i + 1][j], result.nodes[g][i - 1][j]);
      const pressure = g === 0 ? result.cells[i - 1][g][0].interfacePressure.lower : result.cells[i - 1][g].at(-1).interfacePressure.upper;
      outer[1] += sign * pressure * edge.y / 2; outer[2] -= sign * pressure * edge.x / 2;
    }
  }
  // Include the cut half-segments beside LE/TE: their cancellation requires
  // the independent Kutta rows, not just pressure continuity on free nodes.
  for (let b = 0; b < bodies.length; b++) for (let k = 0; k < nx; k++) {
    if (k >= bodies[b].leadingIndex && k < bodies[b].trailingIndex) continue;
    const edge = difference(result.nodes[b][k + 1].at(-1), result.nodes[b][k].at(-1));
    for (const i of [k, k + 1]) if (i > 0 && i < nx) {
      const jump = result.cells[i - 1][b + 1][0].interfacePressure.lower - result.cells[i - 1][b].at(-1).interfacePressure.upper;
      cutTraction[0] += jump * edge.y / 2; cutTraction[1] -= jump * edge.x / 2;
    }
  }
  const bodyForce = result.diagnosticForces.reduce((force, c) => [force[0] + .5 * lengthScale * c.cx, force[1] + .5 * lengthScale * c.cy], [0, 0]);
  const balance = [outer[0], outer[1] + bodyForce[0] + cutTraction[0], outer[2] + bodyForce[1] + cutTraction[1], outer[3]];
  return { order: ['mass', 'xMomentum', 'yMomentum', 'totalEnthalpy'], blocks, outer, cutTraction, bodyForce, balance };
}
