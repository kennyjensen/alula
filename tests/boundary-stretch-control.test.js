import test from 'node:test';
import assert from 'node:assert/strict';
import { createBoundaryStretchControl } from '../src/geometry/boundary-stretch-control.js';
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { solveLinear } from '../src/numerics/linear.js';
const close = (a, b, tolerance = 1e-12) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);
const masses = eta => eta.slice(1).map((v, j) => v - eta[j]);

test('arc form preserves constant spacing through a right-angle stagnation turn on nonuniform xi', () => {
  const eta = [0, .2, 1];
  for (const h of [.2, .1, .05]) {
    const xi = [0, .5 - h, .5, .5 + 1.5 * h, 1];
    const nodes = xi.map(u => eta.map(e => ({ x: e + Math.max(0, u - .5), y: e + Math.max(0, .5 - u) })));
    const vector = createBoundaryStretchControl({ nodes, xi, eta });
    const arc = createBoundaryStretchControl({ nodes, xi, eta, metric: 'polygon-arc' });
    // Exact s=u has F=-s''/s'=0. The secant vector expression instead
    // contains the corner angle and diverges as this stencil is refined.
    close(vector.lower[2], -2 * (.5 * h) / (h * h + (1.5 * h) ** 2));
    close(arc.lower[2], 0); close(arc.upper[2], 0);
  }
  const xi = [0, .25, .5, .75, 1], scalar = u => u + .3 * (u - .5) ** 2;
  const nodes = xi.map(u => eta.map(e => ({ x: e + Math.max(0, scalar(u) - .5), y: e + Math.max(0, .5 - scalar(u)) })));
  const arc = createBoundaryStretchControl({ nodes, xi, eta, metric: 'polygon-arc' });
  close(arc.lower[2], -.6); close(arc.upper[2], -.6);
});

test('boundary stretch matches exact polynomial derivatives and interpolates in the mass coordinate', () => {
  const xi = [0, .08, .2, .37, .6, .8, 1], eta = [0, .13, .7, 1];
  const nodes = xi.map(x => eta.map(e => ({ x: x + (.3 + .4 * e) * x * x, y: 2 * e })));
  for (const discretization of ['quadratic', 'giles-1985']) {
    const control = createBoundaryStretchControl({ nodes, xi, eta, discretization });
    for (let i = 1; i < xi.length - 1; i++) {
      const exact = c => -2 * c / (1 + 2 * c * xi[i] + (discretization === 'giles-1985' ? c * (xi[i + 1] - 2 * xi[i] + xi[i - 1]) : 0));
      close(control.lower[i], exact(.3)); close(control.upper[i], exact(.7));
      eta.forEach((e, j) => close(control.values[i][j], (1 - e) * exact(.3) + e * exact(.7)));
    }
  }
  const a = createBoundaryStretchControl({ nodes, xi, eta });
  const rotate = p => ({ x: 4 + 3 * (.8 * p.x - .6 * p.y), y: -2 + 3 * (.6 * p.x + .8 * p.y) });
  const b = createBoundaryStretchControl({ nodes: nodes.map(row => row.map(rotate)), xi, eta });
  a.values.forEach((row, i) => row.forEach((v, j) => close(v, b.values[i][j], 2e-12)));
  assert.throws(() => createBoundaryStretchControl({ nodes: nodes.map((row, i) => i === 2 ? nodes[0] : row), xi, eta }), /Singular/);
});

test('metric-scaled streamwise control agrees with an independent frozen dense line and owns its input', () => {
  const xi = [0, .08, .2, .37, .6, .8, 1], eta = [0, .31, 1], n = 2 * (xi.length - 2);
  const nodes = xi.map(x => eta.map(e => ({ x: 1.4 * x + .2 * e, y: .3 * x + e })));
  const stretch = xi.map(x => eta.map(e => .4 + x + .2 * e));
  const args = { nodes, massFlows: masses(eta), streamwiseCoordinates: xi, discretization: 'giles-1985' };
  const system = createEllipticStreamtubeGrid({ ...args, streamwiseStretch: stretch });
  const coefficients = system.metrics(nodes), flat = g => system.residuals(g, coefficients).rows.flatMap(row => [row.x, row.y]);
  const base = flat(nodes), matrix = new Float64Array(n * n), h = 1e-5;
  system.residuals(nodes).rows.forEach(({ i, x, y }) => {
    // The affine metric has r_eta=(.2,1), r_xi=(1.4,.3).
    const expected = 1.04 * stretch[i][1] / ((1.04 + 2.05) * system.lengthScale);
    close(x, expected * 1.4); close(y, expected * .3);
  });
  for (let col = 0; col < n; col++) {
    const i = 1 + Math.floor(col / 2), key = col % 2 ? 'y' : 'x';
    const plus = structuredClone(nodes), minus = structuredClone(nodes); plus[i][1][key] += h; minus[i][1][key] -= h;
    const a = flat(plus), b = flat(minus);
    for (let row = 0; row < n; row++) matrix[row * n + col] = (a[row] - b[row]) / (2 * h);
  }
  const expected = solveLinear(matrix, base.map(v => -v)), actual = system.sweep(nodes, 1).nodes;
  for (let col = 0; col < n; col++) { const i = 1 + Math.floor(col / 2), key = col % 2 ? 'y' : 'x'; close(actual[i][1][key] - nodes[i][1][key], expected[col], 2e-11); }
  const before = system.residuals(nodes); stretch[2][1] = 100;
  assert.deepEqual(system.residuals(nodes), before);
  const zero = nodes.map(row => row.map(() => 0));
  const baseline = createEllipticStreamtubeGrid(args), unchanged = createEllipticStreamtubeGrid({ ...args, streamwiseStretch: zero });
  assert.deepEqual(baseline.sweep(nodes), unchanged.sweep(nodes));
  assert.throws(() => createEllipticStreamtubeGrid({ ...args, streamwiseStretch: zero, streamwiseSource: zero }), /cannot be combined/);
});

test('metric-scaled control refines to an independently specified stretched rectangle with harmonic mass', t => {
  const f = x => x + .2 * Math.sin(2 * Math.PI * x) / (2 * Math.PI);
  const control = x => 2 * Math.PI * .2 * Math.sin(2 * Math.PI * x) / (1 + .2 * Math.cos(2 * Math.PI * x));
  let previous = Infinity; const evidence = [];
  for (const nx of [8, 16, 32]) {
    const nt = 6, xi = Array.from({ length: nx + 1 }, (_, i) => i / nx), eta = Array.from({ length: nt + 1 }, (_, j) => j / nt);
    const nodes = xi.map(x => eta.map(e => ({ x: f(x) + .01 * Math.sin(Math.PI * x) * Math.sin(Math.PI * e), y: e })));
    const system = createEllipticStreamtubeGrid({ nodes, massFlows: masses(eta), streamwiseStretch: xi.map(x => eta.map(() => control(x))), discretization: 'giles-1985' });
    const result = smoothEllipticStreamtubeGrid(system, { maxSweeps: 600, tolerance: 1e-11, omega: 1.3 });
    assert.equal(result.converged, true, result.reason);
    let error = 0; result.nodes.forEach((row, i) => row.forEach((p, j) => { error = Math.max(error, Math.abs(p.x - f(xi[i]))); close(p.y, eta[j], 1e-10); }));
    assert.ok(error < previous * .3, `${previous} -> ${error}`); previous = error;
    evidence.push({ nx, error, sweeps: result.history.length - 1 });
  }
  assert.ok(previous < 2e-4); t.diagnostic(JSON.stringify(evidence));
});

test('boundary-derived streamwise control retains valid cylinder grids and refines the independent physical mass field', t => {
  const point = (xi, psi) => {
    const theta = Math.PI * (.75 - .5 * xi), q = psi / Math.sin(theta), radius = .5 * (q + Math.hypot(q, 2));
    return { x: radius * Math.cos(theta), y: radius * Math.sin(theta) };
  };
  let previous = Infinity; const evidence = [];
  for (const nt of [6, 12, 24]) {
    const nx = 2 * nt, xi = Array.from({ length: nx + 1 }, (_, i) => i / nx);
    const eta = Array.from({ length: nt + 1 }, (_, j) => Math.expm1(3 * j / nt) / Math.expm1(3));
    const nodes = xi.map(x => eta.map(e => point(x, .64 * e)));
    const control = createBoundaryStretchControl({ nodes, xi, eta });
    const system = createEllipticStreamtubeGrid({ nodes, massFlows: masses(eta), streamwiseCoordinates: xi,
      streamwiseStretch: control.values, discretization: 'giles-1985' });
    const result = smoothEllipticStreamtubeGrid(system, { maxSweeps: 600, tolerance: 1e-10, omega: 1.3 });
    assert.equal(result.converged, true, result.reason);
    let error = 0, wallShiftOverGap = 0;
    result.nodes.forEach((row, i) => row.forEach((p, j) => {
      const r2 = p.x * p.x + p.y * p.y;
      assert.ok(r2 >= 1 - 1e-12);
      error = Math.max(error, Math.abs(p.y * (1 - 1 / r2) - .64 * eta[j]));
      if (j === 1) wallShiftOverGap = Math.max(wallShiftOverGap, Math.abs(Math.atan2(p.y, p.x) - Math.PI * (.75 - .5 * xi[i])) / (Math.sqrt(r2) - 1));
    }));
    assert.ok(error < .4 * previous, `${previous} -> ${error}`); previous = error;
    assert.ok(wallShiftOverGap < 1, String(wallShiftOverGap));
    evidence.push({ nx, nt, error, wallShiftOverGap, minimumCornerSine: result.quality.minCornerSine });
  }
  t.diagnostic(JSON.stringify(evidence));
});
