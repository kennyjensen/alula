import test from 'node:test';
import assert from 'node:assert/strict';
import { harmonicMapCoefficients, createDiscreteHarmonicMapControl } from '../src/geometry/tests/harmonic-map-control.js';
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { solveLinear } from '../src/numerics/linear.js';
const uniform = n => Array.from({ length: n + 1 }, (_, i) => i / n);
const masses = eta => eta.slice(1).map((v, j) => v - eta[j]);
const close = (a, b, tol = 2e-12) => assert.ok(Math.abs(a - b) < tol, `${a} versus ${b}`);

for (const discretization of ['giles-1985', 'quadratic']) test(`a polynomial harmonic map is an exact stationary grid with ${discretization} differences`, () => {
  const xi = discretization === 'quadratic' ? [0, .08, .21, .38, .6, .84, 1] : uniform(8);
  const eta = discretization === 'quadratic' ? [0, .1, .35, .6, 1] : uniform(4);
  const s = (u, v) => u + .2 * u * u + .1 * u * v + .04 * v * v;
  const values = xi.map(u => eta.map(v => s(u, v))), nodes = values.map(row => row.map((x, j) => ({ x, y: eta[j] })));
  const coefficients = createDiscreteHarmonicMapControl({ values, xi, eta, discretization });
  for (let i = 1; i < xi.length - 1; i++) for (let j = 1; j < eta.length - 1; j++) {
    const d = 1 + .4 * xi[i] + .1 * eta[j];
    close(coefficients[i][j].xiXi, -.4 / d); close(coefficients[i][j].xiEta, -.1 / d); close(coefficients[i][j].etaEta, -.08 / d);
  }
  const input = { nodes, massFlows: masses(eta), streamwiseCoordinates: xi, discretization };
  const controlled = createEllipticStreamtubeGrid({ ...input, harmonicMapControl: coefficients }), plain = createEllipticStreamtubeGrid(input);
  assert.ok(controlled.residuals(nodes).residual < 2e-13);
  assert.ok(plain.residuals(nodes).residual > .01);
  assert.equal(controlled.coordinateEquations.eta, 'Laplace');
  assert.equal(controlled.coordinateEquations.exactMsetSpacingLaw, false);
});

for (const coefficientMode of ['analytic', 'sampled']) test(`mapped Giles SLOR with ${coefficientMode} coefficients converges to an independent conformal flow`, t => {
  const k = .7, a = .15, b = .08, c = .025;
  const S = (u, v) => u + a * Math.sin(2 * Math.PI * u) / (2 * Math.PI) + b * u * v + c * v * v;
  const Su = (u, v) => 1 + a * Math.cos(2 * Math.PI * u) + b * v;
  const point = (u, v) => ({ x: Math.exp(k * S(u, v)) * Math.cos(k * v), y: Math.exp(k * S(u, v)) * Math.sin(k * v) });
  let previous = { xi: Infinity, eta: Infinity }; const evidence = [];
  for (const n of [8, 16, 32]) {
    const xi = uniform(n), eta = uniform(n / 2), exact = xi.map(u => eta.map(v => point(u, v)));
    const nodes = exact.map((row, i) => row.map((p, j) => {
      if (!i || i === n || !j || j === n / 2) return { ...p };
      const d = .004 * Math.sin(Math.PI * xi[i]) * Math.sin(Math.PI * eta[j]);
      return { x: p.x + d, y: p.y - .7 * d };
    }));
    const harmonicMapControl = coefficientMode === 'sampled'
      ? createDiscreteHarmonicMapControl({ values: xi.map(u => eta.map(v => S(u, v))), xi, eta, discretization: 'giles-1985' })
      : xi.map(u => eta.map(v => harmonicMapCoefficients({ xi: Su(u, v),
        xiXi: -2 * Math.PI * a * Math.sin(2 * Math.PI * u), xiEta: b, etaEta: 2 * c })));
    const result = smoothEllipticStreamtubeGrid(createEllipticStreamtubeGrid({ nodes, massFlows: masses(eta),
      streamwiseCoordinates: xi, harmonicMapControl, discretization: 'giles-1985' }), { omega: 1.1, tolerance: 1e-10, maxSweeps: 600 });
    assert.equal(result.converged, true, result.reason); assert.equal(result.quality.valid, true);
    const error = { xi: 0, eta: 0 };
    result.nodes.forEach((row, i) => row.forEach((p, j) => {
      const v = Math.atan2(p.y, p.x) / k, target = Math.log(Math.hypot(p.x, p.y)) / k;
      let u = xi[i]; for (let k = 0; k < 12; k++) u -= (S(u, v) - target) / Su(u, v);
      error.xi = Math.max(error.xi, Math.abs(u - xi[i])); error.eta = Math.max(error.eta, Math.abs(v - eta[j]));
      if (!i || i === n || !j || j === n / 2) assert.deepEqual(p, exact[i][j]);
    }));
    for (const key of ['xi', 'eta']) assert.ok(error[key] < .4 * previous[key], `${key}: ${error[key]} after ${previous[key]}`);
    previous = error; evidence.push({ n, ...error, sweeps: result.history.length - 1 });
  }
  assert.ok(previous.xi < 2e-4); assert.ok(previous.eta < 2e-4); t.diagnostic(JSON.stringify(evidence));
});

test('metric-dependent control uses the same scalar Thomas sweep as an independent frozen dense system', () => {
  const xi = uniform(7), eta = [0, .35, 1], nodes = xi.map(u => eta.map(v => ({ x: u + .2 * v + .03 * Math.sin(Math.PI * u) * Math.sin(Math.PI * v), y: v + .1 * u })));
  const control = xi.map(u => eta.map(v => ({ xiXi: .3 * Math.cos(u), xiEta: -.1, etaEta: .2 * v })));
  const input = { nodes, massFlows: masses(eta), streamwiseCoordinates: xi, harmonicMapControl: control, discretization: 'giles-1985' };
  const system = createEllipticStreamtubeGrid(input), frozen = system.metrics(nodes);
  const residual = q => system.residuals(q, frozen).rows.flatMap(r => [r.x, r.y]);
  const r = residual(nodes), n = r.length, matrix = new Float64Array(n * n), h = 1e-5;
  for (let col = 0; col < n; col++) {
    const i = 1 + Math.floor(col / 2), key = col % 2 ? 'y' : 'x', p = structuredClone(nodes), m = structuredClone(nodes);
    p[i][1][key] += h; m[i][1][key] -= h;
    const a = residual(p), b = residual(m); for (let row = 0; row < n; row++) matrix[row * n + col] = (a[row] - b[row]) / (2 * h);
  }
  const expected = solveLinear(matrix, r.map(v => -v)), swept = system.sweep(nodes, 1).nodes;
  expected.forEach((v, col) => { const i = 1 + Math.floor(col / 2), key = col % 2 ? 'y' : 'x'; close(swept[i][1][key] - nodes[i][1][key], v, 1e-10); });
  // Computational coefficients are unchanged by physical scaling/rotation.
  const map = p => ({ x: 2 + 3 * (.8 * p.x - .6 * p.y), y: -1 + 3 * (.6 * p.x + .8 * p.y) });
  const transformed = createEllipticStreamtubeGrid({ ...input, nodes: nodes.map(row => row.map(map)) }).sweep(nodes.map(row => row.map(map)), 1).nodes;
  swept.forEach((row, i) => row.forEach((p, j) => { const q = map(p); close(transformed[i][j].x, q.x); close(transformed[i][j].y, q.y); }));
  control[3][1].xiXi = 1e10;
  assert.deepEqual(system.sweep(nodes, 1).nodes, swept, 'the constructor must own its fixed map coefficients');
});

test('zero map coefficients exactly recover the existing Laplace equations and invalid maps are rejected', () => {
  const xi = uniform(4), eta = uniform(2), values = xi.map(u => eta.map(v => 4 * u + v));
  const nodes = xi.map(u => eta.map(v => ({ x: u, y: v }))), massFlows = masses(eta);
  const zero = createDiscreteHarmonicMapControl({ values, xi, eta });
  const a = createEllipticStreamtubeGrid({ nodes, massFlows }), b = createEllipticStreamtubeGrid({ nodes, massFlows, harmonicMapControl: zero });
  assert.deepEqual(a.residuals(nodes), b.residuals(nodes)); assert.deepEqual(a.sweep(nodes), b.sweep(nodes));
  assert.deepEqual(a.coordinateEquations, b.coordinateEquations);
  assert.throws(() => harmonicMapCoefficients({ xi: 0, xiXi: 1, xiEta: 0, etaEta: 0 }), /positive/);
  assert.throws(() => createDiscreteHarmonicMapControl({ values: values.toReversed(), xi, eta }), /ordered/);
  assert.throws(() => createEllipticStreamtubeGrid({ nodes, massFlows, harmonicMapControl: zero, streamwiseSource: values }), /cannot be combined/);
  assert.throws(() => createEllipticStreamtubeGrid({ nodes, massFlows, harmonicMapControl: zero.slice(1) }), /must match/);
});
