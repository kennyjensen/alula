// SPDX-License-Identifier: GPL-2.0-or-later
// Exact psi=2*x*y flow in the first quadrant, with a lower cut/wall corner.
// Only the known corner station uses the mean of adjacent F or P controls.
import { harmonicStagnationPoint } from './harmonic-stagnation-stations.js';
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid } from '../../src/geometry/elliptic-streamtube-grid.js';
import { smoothOrthogonalBoundaryGrid } from '../../src/geometry/tests/orthogonal-boundary-slor.js';

export function runOrthogonalStagnationControl({ nt = 8, maxSweeps = 600, controlUpdate = 'direct', sourceForm = 'metric-stretch', streamwiseSourceDiscretization = 'centered' } = {}) {
  if (!Number.isInteger(nt) || nt < 4 || nt > 48 || nt % 2 || !Number.isInteger(maxSweeps) || maxSweeps < 0)
    throw new Error('Invalid bounded corner-control experiment.');
  if (!['direct', 'damped'].includes(controlUpdate)) throw new Error('Unknown corner-control iteration.');
  if (!['metric-stretch', 'poisson'].includes(sourceForm)) throw new Error('Unknown corner-control source form.');
  const nx = 2 * nt, xi = Array.from({ length: nx + 1 }, (_, i) => i / nx);
  const eta = Array.from({ length: nt + 1 }, (_, j) => Math.expm1(3 * j / nt) / Math.expm1(3));
  const initial = xi.map(u => eta.map(e => harmonicStagnationPoint(2 * u - 1, e)));
  const massFlows = eta.slice(1).map((e, j) => e - eta[j]);
  // Fix these rates during refinement, using the same documented 0.45
  // per first inward mass interval conversion on the eight-tube baseline.
  const decay = { lower: .45 * Math.expm1(3) / Math.expm1(3 / 8), upper: .45 / (1 - Math.expm1(21 / 8) / Math.expm1(3)) };
  const controls = { sourceForm, background: initial.map(row => row.map(() => 0)), sides: ['lower'], decay, corners: { lower: [nt] } };
  const options = { nodes: initial, massFlows, streamwiseCoordinates: xi, discretization: 'giles-1985', streamwiseSourceDiscretization, orthogonalBoundaryControl: controls };
  const started = performance.now(), result = controlUpdate === 'damped'
    ? smoothOrthogonalBoundaryGrid(options, { maxSweeps, omega: 1.3, tolerance: 1e-10 })
    : smoothEllipticStreamtubeGrid(createEllipticStreamtubeGrid(options), { maxSweeps, omega: 1.3, tolerance: 1e-10 });
  let massError = 0, maximumMovement = 0, maximumWallDerivativeShear = 0, firstRowShear = 0, maximumAdjacentRatio = 1;
  let worstMass = null;
  result.nodes.forEach((row, i) => row.forEach((p, j) => {
    const error = Math.abs(2 * p.x * p.y - eta[j]);
    if (error > massError) { massError = error; worstMass = { i, j, x: p.x, y: p.y }; }
    maximumMovement = Math.max(maximumMovement, Math.hypot(p.x - initial[i][j].x, p.y - initial[i][j].y));
    if ((!i || i === nx || !j || j === nt) && (p.x !== initial[i][j].x || p.y !== initial[i][j].y))
      throw new Error('The corner experiment moved a fixed boundary.');
    if (i && i < nx) {
      const a = result.nodes[i - 1][j], b = result.nodes[i + 1][j], left = Math.hypot(p.x - a.x, p.y - a.y), right = Math.hypot(p.x - b.x, p.y - b.y);
      maximumAdjacentRatio = Math.max(maximumAdjacentRatio, left / right, right / left);
    }
  }));
  // The boundary angle is undefined at the corner. Measure it separately
  // on the fixed regular-wall window, |signed wall arc| >= 0.25.
  for (let i = 1; i < nx; i++) {
    const station = 2 * xi[i] - 1; if (Math.abs(station) < .25) continue;
    const a = eta[1], b = eta[2], row = result.nodes[i];
    const derivative = Object.fromEntries(['x', 'y'].map(key => [key,
      b * (row[1][key] - row[0][key]) / (a * (b - a)) - a * (row[2][key] - row[0][key]) / (b * (b - a))]));
    maximumWallDerivativeShear = Math.max(maximumWallDerivativeShear,
      station > 0 ? Math.abs(derivative.x / derivative.y) : Math.abs(derivative.y / derivative.x));
    firstRowShear = Math.max(firstRowShear, station > 0 ? Math.abs((row[1].x - station) / row[1].y) : Math.abs((row[1].y + station) / row[1].x));
  }
  return { nt, nx, controlUpdate, sourceForm, streamwiseSourceDiscretization, converged: result.converged, reason: result.reason, sweeps: result.sweeps ?? result.history.length - 1,
    residual: result.history.at(-1).residual, massError, worstMass, maximumWallDerivativeShear, firstRowShear,
    maximumMovement, maximumAdjacentRatio, quality: result.quality, history: result.history, nodes: result.nodes, initial, xi, eta,
    controls: { decay, cornerStations: controls.corners, sourceExtension: 'matched exponential tails', transverseSource: 0,
      boundaryAngleWindow: 'absolute signed wall arc at least 0.25, excluding fixed inlet/outlet', omega: 1.3, maxSweeps },
    seconds: (performance.now() - started) / 1000, physicalAcceptance: false, guiDefaultChanged: false, boundariesExactlyFixed: true };
}
