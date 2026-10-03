// SPDX-License-Identifier: GPL-2.0-or-later
// Bounded exact-flow experiment. Not a public initializer or a recovered
// MSET control law. Alternate boundary reconstruction with one unchanged
// Giles/Thomas SLOR sweep, retaining Q=0 and fixed boundary coordinates.
import { createBoundaryStretchControl } from '../../src/geometry/boundary-stretch-control.js';
import { createOrthogonalBoundaryControl } from '../../src/geometry/orthogonal-boundary-control.js';
import { createEllipticStreamtubeGrid } from '../../src/geometry/elliptic-streamtube-grid.js';

export function runOrthogonalCylinderControl({ nt = 6, maxSweeps = 600, controlRelaxation = .3, method = 'implicit', exponential = false, decayRates } = {}) {
  if (!Number.isInteger(nt) || nt < 3 || nt > 48 || !Number.isInteger(maxSweeps) || maxSweeps < 0
    || !Number.isFinite(controlRelaxation) || !(controlRelaxation > 0 && controlRelaxation <= 1)
    || !['explicit', 'implicit'].includes(method) || typeof exponential !== 'boolean' || exponential && method !== 'implicit'
    || decayRates !== undefined && !exponential) throw new Error('Invalid bounded angle-control experiment.');
  const started = performance.now(), nx = 2 * nt;
  const xi = Array.from({ length: nx + 1 }, (_, i) => i / nx);
  const eta = Array.from({ length: nt + 1 }, (_, j) => Math.expm1(3 * j / nt) / Math.expm1(3));
  const exactPoint = (u, e) => {
    const theta = Math.PI * (.75 - .5 * u), q = .64 * e / Math.sin(theta), radius = .5 * (q + Math.hypot(q, 2));
    return { x: radius * Math.cos(theta), y: radius * Math.sin(theta) };
  };
  const initial = xi.map(u => eta.map(e => exactPoint(u, e))), massFlows = eta.slice(1).map((e, j) => e - eta[j]);
  const background = createBoundaryStretchControl({ nodes: initial, xi, eta });
  // GRAPE's default 0.45 per computational interval, expressed here in
  // first inward mass intervals. This nonuniform-coordinate extension is
  // an explicit reconstruction choice; it is not an MSET default.
  const decay = exponential ? decayRates ?? { lower: .45 / eta[1], upper: .45 / (1 - eta[nt - 1]) } : undefined;
  let nodes = structuredClone(initial), boundary = [...background.lower], converged = false, reason = 'sweep limit';
  const history = [], boundaryControlHistory = [];
  const makeSystem = () => createEllipticStreamtubeGrid({ nodes, massFlows, streamwiseCoordinates: xi,
    discretization: 'giles-1985', ...(method === 'implicit' ? { orthogonalBoundaryControl: { background: background.values, sides: ['lower'], decay } }
      : { streamwiseStretch: boundary.map((f, i) => eta.map(e => (1 - e) * f + e * background.upper[i])) }) });
  for (let iteration = 0; iteration <= maxSweeps; iteration++) {
    try {
      const target = createOrthogonalBoundaryControl({ nodes, xi, eta, sides: ['lower'] }).lower;
      if (method === 'implicit') boundary = boundary.map((f, i) => target[i]?.stretch ?? f);
      let controlResidual = 0;
      for (let i = 1; i < nx; i++) controlResidual = Math.max(controlResidual, Math.abs(target[i].stretch - boundary[i]) / Math.max(1, Math.abs(target[i].stretch)));
      const system = makeSystem(), residual = system.residuals(nodes).residual, quality = system.quality(nodes);
      history.push({ iteration, residual, controlResidual, minimumCornerSine: quality.minCornerSine, invalidCells: quality.invalidCells.length });
      if (!quality.valid) { reason = 'folded grid'; break; }
      if (residual <= 1e-10 && controlResidual <= 1e-9) { converged = true; reason = 'converged'; break; }
      if (iteration === maxSweeps) break;
      // Sorenson Eq. (17) form, with the documented defaults 0.3 and 1.0.
      // Applying it to normalized metric-scaled F (rather than GRAPE's
      // P/Q amplitudes) is an explicit choice in this experiment.
      if (method === 'explicit') boundary = boundary.map((f, i) => !i || i === nx ? f : f + Math.sign(target[i].stretch - f)
        * Math.min(controlRelaxation * Math.abs(target[i].stretch - f), Math.max(Math.abs(f), 1)));
      boundaryControlHistory.push(boundary.slice());
      nodes = makeSystem().sweep(nodes, 1.3).nodes;
    } catch (error) { reason = error.message; break; }
  }
  let massError = 0, wallShiftOverGap = 0, maximumMovement = 0, maximumWallDerivativeShear = 0;
  nodes.forEach((row, i) => row.forEach((p, j) => {
    const r2 = p.x * p.x + p.y * p.y;
    massError = Math.max(massError, Math.abs(p.y * (1 - 1 / r2) - .64 * eta[j]));
    maximumMovement = Math.max(maximumMovement, Math.hypot(p.x - initial[i][j].x, p.y - initial[i][j].y));
    if (j === 1) wallShiftOverGap = Math.max(wallShiftOverGap,
      Math.abs(Math.atan2(p.y, p.x) - Math.PI * (.75 - .5 * xi[i])) / (Math.sqrt(r2) - 1));
    if ((!i || i === nx || !j || j === nt) && (p.x !== initial[i][j].x || p.y !== initial[i][j].y))
      throw new Error('The angle-control experiment moved a fixed boundary.');
  }));
  for (let i = 1; i < nx; i++) {
    // Independent one-sided quadratic derivative of the actual grid;
    // not the orthogonal derivative assumed in the control construction.
    const a = eta[1], b = eta[2], t = Math.PI * (.75 - .5 * xi[i]);
    const derivative = Object.fromEntries(['x', 'y'].map(key => [key,
      b * (nodes[i][1][key] - nodes[i][0][key]) / (a * (b - a))
      - a * (nodes[i][2][key] - nodes[i][0][key]) / (b * (b - a))]));
    maximumWallDerivativeShear = Math.max(maximumWallDerivativeShear,
      Math.abs(derivative.x * Math.sin(t) - derivative.y * Math.cos(t))
      / Math.abs(derivative.x * Math.cos(t) + derivative.y * Math.sin(t)));
  }
  return { nt, nx, method, converged, reason, sweeps: history.at(-1)?.iteration ?? 0, history, boundaryControlHistory,
    massError, wallShiftOverGap, maximumWallDerivativeShear, maximumMovement,
    quality: makeSystem().quality(nodes), initial, nodes, xi, eta, background, boundary,
    controls: { ...(method === 'explicit' ? { controlRelaxation, changeLimit: 1 } : { controlElimination: 'analytic with same-line sensitivity' }),
      omega: 1.3, maxSweeps, sourceExtension: decay ? { method: 'matched exponential tails', decay,
        normalization: decayRates ? 'prescribed rates in normalized mass' : '0.45 per first inward mass interval' } : 'linear in normalized mass',
      boundaryAngle: 'lower only, 90 degrees', upper: 'fixed background stretch', transverseSource: 0 },
    boundariesExactlyFixed: true, seconds: (performance.now() - started) / 1000, physicalAcceptance: false, exactMsetSpacingLaw: false };
}
