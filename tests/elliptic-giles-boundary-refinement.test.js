import test from 'node:test';
import assert from 'node:assert/strict';
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';

// z(w) = w + epsilon*sin(pi*(w-i)), on the closed unit square.
// |z'(w)-1| <= epsilon*pi*cosh(pi) < 1. Integrating this derivative
// bound along a segment proves injectivity and a positive Jacobian on the
// whole square. The top is horizontal with x_v=0, while the bottom curves.
const epsilon = .02;
const derivativePerturbationBound = epsilon * Math.PI * Math.cosh(Math.PI);
function exactMap(u, v) {
  const a = Math.PI * u, b = Math.PI * (v - 1);
  return { x: u + epsilon * Math.sin(a) * Math.cosh(b),
    y: v + epsilon * Math.cos(a) * Math.sinh(b) };
}

// Independent inversion of the analytic map: ordinary two-real-variable
// complex Newton, no grid metric/residual or finite-difference helper.
function inverseMap(point, initialU, initialV) {
  let u = initialU, v = initialV;
  for (let iteration = 0; iteration <= 16; iteration++) {
    const p = exactMap(u, v), rx = p.x - point.x, ry = p.y - point.y;
    const residual = Math.hypot(rx, ry);
    if (residual <= 2e-14) {
      assert.ok(u >= -1e-12 && u <= 1 + 1e-12 && v >= -1e-12 && v <= 1 + 1e-12,
        `Inverse left the prescribed conformal chart: ${u},${v}`);
      // The same global derivative bound on a slightly enlarged rectangle
      // bounds inverse-coordinate uncertainty from the checked residual.
      const lowerLipschitz = 1 - epsilon * Math.PI * Math.cosh(Math.PI * (1 + 1e-12));
      return { u, v, residual, coordinateUncertainty: residual / lowerLipschitz, iteration };
    }
    assert.ok(iteration < 16, `Analytic inverse did not resolve its residual: ${residual}`);
    const a = 1 + epsilon * Math.PI * Math.cos(Math.PI * u) * Math.cosh(Math.PI * (v - 1));
    const b = epsilon * Math.PI * Math.sin(Math.PI * u) * Math.sinh(Math.PI * (v - 1));
    const denominator = a * a + b * b;
    assert.ok(denominator > 0 && Number.isFinite(denominator));
    u -= (a * rx - b * ry) / denominator;
    v -= (b * rx + a * ry) / denominator;
  }
  throw new Error('Unreachable analytic-inverse iteration limit.');
}

test('literal Giles upper x-copy boundary converges at first order on a regular analytic curved domain', t => {
  const start = performance.now();
  assert.ok(derivativePerturbationBound < 1);
  const cases = [];
  for (const [nx, nt] of [[8, 4], [16, 8], [32, 16]]) {
    const caseStart = performance.now(), h = 1 / nt;
    const xi = Array.from({ length: nx + 1 }, (_, i) => i / nx);
    const exact = xi.map(u => Array.from({ length: nt + 1 }, (_, j) => exactMap(u, j * h)));
    const system = createEllipticStreamtubeGrid({ nodes: exact, massFlows: Array(nt).fill(h),
      streamwiseCoordinates: xi, discretization: 'giles-1985', boundaryConditions: { upper: 'giles-vertical' } });
    const initial = system.residuals(exact);
    const initialBoundaryResidual = initial.boundaryRows.find(row => row.i === nx / 2 && row.j === nt).x;
    // At u=1/2, the exact one-sided difference is known without evaluating
    // any production grid derivative. It is O(h), although x_v is zero.
    const expectedBoundaryResidual = epsilon * (1 - Math.cosh(Math.PI * h)) / (h * system.lengthScale);
    assert.ok(Math.abs(initialBoundaryResidual - expectedBoundaryResidual) < 64 * Number.EPSILON / (h * system.lengthScale));
    assert.ok(initialBoundaryResidual < 0);

    const result = smoothEllipticStreamtubeGrid(system, { maxSweeps: 1200, omega: 1.3, tolerance: 1e-10 });
    const final = result.history.at(-1);
    t.diagnostic(JSON.stringify({ nx, nt, converged: result.converged, reason: result.reason,
      sweeps: result.history.length - 1, residual: final.residual, seconds: (performance.now() - caseStart) / 1000 }));
    assert.equal(result.converged, true, result.reason);
    assert.equal(result.quality.valid, true);
    assert.ok(final.residual <= 1e-10);

    const errors = { coordinate: 0, x: 0, y: 0, xi: 0, eta: 0, upperX: 0 };
    let inverseResidual = 0, inverseCoordinateUncertainty = 0, inverseIterations = 0, minimumCornerDeterminant = Infinity;
    result.nodes.forEach((row, i) => row.forEach((point, j) => {
      const reference = exact[i][j];
      if (i === 0 || i === nx || j === 0) assert.deepEqual(point, reference);
      if (j === nt) {
        assert.equal(point.y, 1);
        if (i > 0 && i < nx) assert.equal(point.x, result.nodes[i][nt - 1].x);
      }
      const ex = Math.abs(point.x - reference.x), ey = Math.abs(point.y - reference.y);
      errors.x = Math.max(errors.x, ex); errors.y = Math.max(errors.y, ey);
      errors.coordinate = Math.max(errors.coordinate, Math.hypot(ex, ey));
      if (j === nt) errors.upperX = Math.max(errors.upperX, ex);
      const inverse = inverseMap(point, xi[i], j * h);
      errors.xi = Math.max(errors.xi, Math.abs(inverse.u - xi[i]));
      errors.eta = Math.max(errors.eta, Math.abs(inverse.v - j * h));
      inverseResidual = Math.max(inverseResidual, inverse.residual);
      inverseCoordinateUncertainty = Math.max(inverseCoordinateUncertainty, inverse.coordinateUncertainty);
      inverseIterations = Math.max(inverseIterations, inverse.iteration);
    }));
    // Independently check all four corner determinants of every Q1 cell.
    for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) {
      const p = [result.nodes[i][j], result.nodes[i + 1][j], result.nodes[i + 1][j + 1], result.nodes[i][j + 1]];
      for (let k = 0; k < 4; k++) {
        const a = p[k], b = p[(k + 1) % 4], c = p[(k + 2) % 4];
        const determinant = (b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x);
        minimumCornerDeterminant = Math.min(minimumCornerDeterminant, determinant);
        assert.ok(determinant > 0, `Nonpositive independent corner determinant at ${i},${j},${k}`);
      }
    }
    assert.ok(inverseCoordinateUncertainty < 1e-12);
    const record = { nx, nt, sweeps: result.history.length - 1,
      initialBoundaryResidual, expectedBoundaryResidual, interiorResidual: final.interiorResidual,
      boundaryResidual: final.boundaryResidual, errors, inverseResidual, inverseCoordinateUncertainty,
      inverseIterations, minimumCornerDeterminant, seconds: (performance.now() - caseStart) / 1000 };
    cases.push(record); t.diagnostic(JSON.stringify(record));
  }
  const orders = cases.slice(1).map((fine, k) => ({
    from: [cases[k].nx, cases[k].nt], to: [fine.nx, fine.nt],
    boundaryResidual: Math.log2(Math.abs(cases[k].initialBoundaryResidual / fine.initialBoundaryResidual)),
    ...Object.fromEntries(Object.keys(fine.errors).map(key => [key, Math.log2(cases[k].errors[key] / fine.errors[key])])),
  }));
  t.diagnostic(JSON.stringify({ derivativePerturbationBound, orders, totalSeconds: (performance.now() - start) / 1000,
    scope: 'Analytic mixed-boundary stencil verification; no airfoil, panel, FE or browser solve, and no airfoil acceptance.' }));
  for (const [k, order] of orders.entries()) {
    // The exact residual expansion is -epsilon*pi^2*h/(2*lengthScale)
    // plus O(h^3). The x-copy boundary drives a first-order map/xi error;
    // y/eta can have smaller leading coefficients and need only decrease.
    assert.ok(order.boundaryResidual > .9 && order.boundaryResidual < 1.2);
    for (const key of ['coordinate', 'xi', 'upperX'])
      assert.ok(order[key] > .5 && order[key] < 1.6, `${key} order outside broad first-order range: ${order[key]}`);
    for (const key of Object.keys(cases[k].errors))
      assert.ok(cases[k + 1].errors[key] < cases[k].errors[key], `${key} error did not decrease`);
  }
});
