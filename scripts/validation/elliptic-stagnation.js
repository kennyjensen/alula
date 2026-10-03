// SPDX-License-Identifier: GPL-2.0-or-later
// Exact incompressible stagnation flow: phi=x*x-y*y, psi=2*x*y.
// Both xi=(phi+1)/2 and eta=psi-offset are harmonic in physical space.
// The inverse z=sqrt(phi+i*psi) is singular at the boundary stagnation point.
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid } from '../../src/geometry/elliptic-streamtube-grid.js';

export function stagnationPoint(phi, psi) {
  const radius = Math.hypot(phi, psi);
  // Use the product 2*x*y=psi to avoid subtracting nearly equal numbers.
  if (phi >= 0) {
    const x = Math.sqrt((radius + phi) / 2);
    return { x, y: x ? psi / (2 * x) : 0 };
  }
  const y = Math.sqrt((radius - phi) / 2);
  return { x: psi / (2 * y), y };
}

// Analytic evaluation of the original stencil at phi=0, psi=h, with
// uniform delta-phi=r*h and delta-psi=h. This does not call grid operators.
// The factor four converts phi derivatives to xi=(phi+1)/2 derivatives.
export function stagnationCenterOracle(h, r = 1) {
  const radius = Math.hypot(1, r);
  const coefficient = .5 * ((Math.sqrt(radius + 1) - Math.SQRT2) / (r * r)
    + (1 - Math.SQRT2) / (radius + 1));
  return { alpha: 1 / (2 * h), gamma: 2 / (h * (radius + 1)), beta: 0,
    residualX: 4 * coefficient / h ** 2.5, residualY: 4 * coefficient / h ** 2.5 };
}

export function runStagnationCase({ nt, offset = 0, stretch = 0 }) {
  const started = performance.now(), nx = 2 * nt;
  const xi = Array.from({ length: nx + 1 }, (_, i) => i / nx);
  // A fixed smooth stretch under refinement; adjacent cell ratios approach
  // one. Keeping a geometric adjacent ratio fixed is a different experiment.
  const eta = Array.from({ length: nt + 1 }, (_, j) => stretch
    ? Math.expm1(stretch * j / nt) / Math.expm1(stretch) : j / nt);
  const massFlows = eta.slice(1).map((v, j) => v - eta[j]);
  const exact = xi.map(u => eta.map(v => stagnationPoint(2 * u - 1, v + offset)));
  const system = createEllipticStreamtubeGrid({ nodes: exact, massFlows,
    streamwiseCoordinates: xi, discretization: 'giles-1985' });
  // Start from the exact continuum map: movement then measures the discrete
  // equation's error, without a separate initialization error or panel field.
  const initialResidual = system.residuals(exact).residual;
  const result = smoothEllipticStreamtubeGrid(system, { maxSweeps: 600, tolerance: 1e-9, omega: 1.3 });
  let absolute = 0, localIntervals = 0, awayFromStagnation = 0, coordinate = 0, changedBoundaries = 0;
  let absolutePeak = null, localPeak = null;
  result.nodes.forEach((row, i) => row.forEach((point, j) => {
    const reference = exact[i][j], boundary = !i || i === nx || !j || j === nt;
    if (boundary) {
      if (point.x !== reference.x || point.y !== reference.y) changedBoundaries++;
      return;
    }
    const phi = 2 * xi[i] - 1, psi = eta[j] + offset;
    const error = Math.abs(2 * point.x * point.y - psi);
    const width = Math.min(massFlows[j - 1], massFlows[j]);
    const location = { i, j, phi, psi, width, error, localIntervals: error / width };
    if (error > absolute) { absolute = error; absolutePeak = location; }
    if (error / width > localIntervals) { localIntervals = error / width; localPeak = location; }
    if (Math.hypot(phi, psi) >= .25) awayFromStagnation = Math.max(awayFromStagnation, error);
    coordinate = Math.max(coordinate, Math.hypot(point.x - reference.x, point.y - reference.y));
  }));
  return { nx, nt, offset, stretch, converged: result.converged, reason: result.reason,
    sweeps: result.history.length - 1, initialResidual, residual: result.history.at(-1).residual,
    invalidCells: result.quality.invalidCells.length, changedBoundaries,
    errors: { absolute, localIntervals, awayFromStagnation, coordinate }, absolutePeak, localPeak,
    seconds: (performance.now() - started) / 1000 };
}

export function stagnationRefinement(coarse, fine) {
  return Object.fromEntries(['absolute', 'localIntervals', 'awayFromStagnation', 'coordinate']
    .map(key => [key, Math.log2(coarse.errors[key] / fine.errors[key])]));
}
