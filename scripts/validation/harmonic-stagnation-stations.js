// SPDX-License-Identifier: GPL-2.0-or-later
// Exact W=z² flow with harmonic station coordinate X=x-y, rather than phi.
// x=(X+sqrt(X²+2psi))/2, y=(-X+sqrt(X²+2psi))/2.
// On the two wall rays X is signed arc from stagnation. For equal X
// increments, edge lengths lie between h/sqrt(2) and h. On either wall,
// tangential station shift equals the normal gap: no sqrt(psi)/psi shear.
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid } from '../../src/geometry/elliptic-streamtube-grid.js';

export const harmonicStagnationPoint = (station, psi) => {
  if (!Number.isFinite(station) || !Number.isFinite(psi) || psi < 0) throw new Error('Invalid harmonic stagnation coordinates.');
  const root = Math.hypot(station, Math.sqrt(2 * psi));
  // Rationalize the small root to avoid cancellation near either wall.
  return station >= 0 ? { x: .5 * (root + station), y: psi ? psi / (root + station) : 0 }
    : { x: psi ? psi / (root - station) : 0, y: .5 * (root - station) };
};

export function runHarmonicStagnationStations({ nt = 8, stretch = 0, offset = 0 } = {}) {
  const nx = 2 * nt, eta = Array.from({ length: nt + 1 }, (_, j) => stretch ? Math.expm1(stretch * j / nt) / Math.expm1(stretch) : j / nt);
  const nodes = Array.from({ length: nx + 1 }, (_, i) => eta.map(e => harmonicStagnationPoint(2 * i / nx - 1, offset + e)));
  const massFlows = eta.slice(1).map((e, j) => e - eta[j]);
  const system = createEllipticStreamtubeGrid({ nodes, massFlows, discretization: 'giles-1985' });
  const started = performance.now(), result = smoothEllipticStreamtubeGrid(system, { maxSweeps: 600, tolerance: 1e-10, omega: 1.3 });
  let stationError = 0, massError = 0, maximumAdjacentRatio = 1, boundaryMovement = 0, nearWallShiftOverGap = 0;
  result.nodes.forEach((row, i) => row.forEach((p, j) => {
    const station = 2 * i / nx - 1;
    stationError = Math.max(stationError, Math.abs(p.x - p.y - station));
    massError = Math.max(massError, Math.abs(2 * p.x * p.y - offset - eta[j]));
    if (!i || i === nx || !j || j === nt) boundaryMovement = Math.max(boundaryMovement, Math.hypot(p.x - nodes[i][j].x, p.y - nodes[i][j].y));
    if (!offset && j === 1 && Math.abs(station) >= .25) nearWallShiftOverGap = Math.max(nearWallShiftOverGap,
      station > 0 ? Math.abs(p.x - station) / p.y : Math.abs(p.y + station) / p.x);
    if (i && i < nx) {
      const a = result.nodes[i - 1][j], b = result.nodes[i + 1][j], first = Math.hypot(p.x - a.x, p.y - a.y), last = Math.hypot(p.x - b.x, p.y - b.y);
      maximumAdjacentRatio = Math.max(maximumAdjacentRatio, first / last, last / first);
    }
  }));
  return { nx, nt, stretch, offset, converged: result.converged, reason: result.reason, quality: result.quality,
    sweeps: result.history.length - 1, residual: result.history.at(-1).residual, boundaryMovement, stationError, massError,
    maximumAdjacentRatio, nearWallShiftOverGap, seconds: (performance.now() - started) / 1000 };
}
