import test from 'node:test';
import assert from 'node:assert/strict';
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { solveLinear } from '../src/numerics/linear.js';

const close = (a, b, tolerance) => assert.ok(Math.abs(a - b) <= tolerance, `${a} != ${b} (tolerance ${tolerance})`);
const masses = eta => eta.slice(1).map((v, j) => v - eta[j]);
const boundary = (actual, expected) => actual.forEach((row, i) => row.forEach((p, j) => {
  if (!i || i === actual.length - 1 || !j || j === row.length - 1) assert.deepEqual(p, expected[i][j]);
}));

// Independent manufactured physical coordinates. For z=exp(k*(f(xi)+i*eta)),
// conformal transformation gives Δz xi = -f''/(f'^3 |F'|²), Δz eta=0.
// P is prescribed at logical nodes; it is not recomputed at moving nodes.
const kappa = .8, stretch = .2;
const f = xi => xi + stretch * Math.sin(2 * Math.PI * xi) / (2 * Math.PI);
const fp = xi => 1 + stretch * Math.cos(2 * Math.PI * xi);
const fpp = xi => -2 * Math.PI * stretch * Math.sin(2 * Math.PI * xi);
const exactPoint = (xi, eta) => ({ x: Math.exp(kappa * f(xi)) * Math.cos(kappa * eta), y: Math.exp(kappa * f(xi)) * Math.sin(kappa * eta) });
const exactSource = xi => -fpp(xi) / (fp(xi) ** 3 * kappa ** 2 * Math.exp(2 * kappa * f(xi)));
const inverseXi = point => {
  const target = Math.log(Math.hypot(point.x, point.y)) / kappa;
  let xi = target;
  for (let i = 0; i < 12; i++) xi -= (f(xi) - target) / fp(xi);
  return xi;
};
function patch(nx, nt, perturb = .006) {
  const eta = Array.from({ length: nt + 1 }, (_, j) => Math.expm1(.7 * j / nt) / Math.expm1(.7));
  const exact = Array.from({ length: nx + 1 }, (_, i) => eta.map(e => exactPoint(i / nx, e)));
  const nodes = exact.map((row, i) => row.map((p, j) => {
    if (!i || i === nx || !j || j === nt) return { ...p };
    const d = perturb * Math.sin(Math.PI * i / nx) * Math.sin(Math.PI * eta[j]);
    return { x: p.x + d, y: p.y - .7 * d };
  }));
  return { nodes, massFlows: masses(eta), streamwiseSource: exact.map((row, i) => row.map(() => exactSource(i / nx))), exact, eta };
}

test('prescribed streamwise Poisson control refines to independent curved coordinates with harmonic mass eta', t => {
  let previousXi = Infinity, previousEta = Infinity; const evidence = [];
  for (const [nx, nt] of [[8, 4], [16, 8], [32, 16]]) {
    const input = patch(nx, nt), before = structuredClone(input);
    const system = createEllipticStreamtubeGrid(input);
    const result = smoothEllipticStreamtubeGrid(system, { omega: 1.1, tolerance: 1e-10, maxSweeps: 600 });
    assert.equal(result.converged, true, result.reason); assert.equal(result.quality.valid, true);
    assert.equal(result.coordinateEquations.xi, 'Poisson with prescribed fixed P');
    assert.equal(result.coordinateEquations.eta, 'Laplace');
    assert.equal(result.coordinateEquations.exactMsetSpacingLaw, false);
    let xiError = 0, etaError = 0;
    result.nodes.forEach((row, i) => row.forEach((p, j) => {
      xiError = Math.max(xiError, Math.abs(inverseXi(p) - i / nx));
      etaError = Math.max(etaError, Math.abs(Math.atan2(p.y, p.x) / kappa - input.eta[j]));
    }));
    assert.ok(xiError < .4 * previousXi, `xi refinement ${xiError} after ${previousXi}`);
    assert.ok(etaError < .4 * previousEta, `eta refinement ${etaError} after ${previousEta}`);
    previousXi = xiError; previousEta = etaError;
    boundary(result.nodes, input.exact); assert.deepEqual(input, before);
    evidence.push({ nx, nt, sweeps: result.history.length - 1, xiError, etaError });
  }
  assert.ok(previousXi < 5e-4); assert.ok(previousEta < 2e-4);
  t.diagnostic(JSON.stringify(evidence));
});

for (const xi of [undefined, [0, .08, .2, .36, .55, .7, .88, 1]]) test(`source residual sign and frozen dense line oracle on ${xi ? 'nonuniform' : 'uniform'} xi`, () => {
  const nx = 7, eta = [0, .28, 1];
  const nodes = Array.from({ length: nx + 1 }, (_, i) => eta.map(e => {
    const u = xi?.[i] ?? i / nx, d = .025 * Math.sin(Math.PI * u) * Math.sin(Math.PI * e);
    return { x: 1.4 * u + .3 * e + d, y: .2 * u + e - .6 * d };
  }));
  const source = nodes.map((row, i) => row.map((_, j) => 1.1 * Math.sin(2 * Math.PI * i / nx + eta[j])));
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: masses(eta), streamwiseSource: source, streamwiseCoordinates: xi });
  const baseline = createEllipticStreamtubeGrid({ nodes, massFlows: masses(eta), streamwiseCoordinates: xi });
  const firstXi = i => {
    if (!xi) return [-nx / 2, 0, nx / 2];
    const a = xi[i] - xi[i - 1], b = xi[i + 1] - xi[i];
    return [-b / (a * (a + b)), (b - a) / (a * b), a / (b * (a + b))];
  };
  const frozen = system.metrics(nodes), zero = baseline.residuals(nodes), controlled = system.residuals(nodes, frozen);
  // This term is derived from Δxi=P, separately from the production residual.
  // A wrong sign or use of J instead of J² must fail before the dense oracle.
  controlled.rows.forEach((row, n) => {
    const { i, j } = row, m = frozen[i][j];
    for (const key of ['x', 'y']) {
      const weights = firstXi(i), dxi = weights.reduce((sum, weight, k) => sum + weight * nodes[i - 1 + k][j][key], 0);
      const expected = m.jacobian ** 2 * source[i][j] * dxi / ((m.alpha + m.gamma) * system.lengthScale);
      close(row[key] - zero.rows[n][key], expected, 5e-15);
    }
  });
  const flat = grid => system.residuals(grid, frozen).rows.flatMap(row => [row.x, row.y]);
  const zeroFrozen = baseline.metrics(nodes), flatZero = grid => baseline.residuals(grid, zeroFrozen).rows.flatMap(row => [row.x, row.y]);
  const base = flat(nodes), n = base.length, matrix = new Float64Array(n * n), h = 1e-5;
  for (let col = 0; col < n; col++) {
    const i = 1 + Math.floor(col / 2), key = col % 2 ? 'y' : 'x';
    const plus = structuredClone(nodes), minus = structuredClone(nodes);
    plus[i][1][key] += h; minus[i][1][key] -= h;
    const a = flat(plus), b = flat(minus);
    for (let row = 0; row < n; row++) matrix[row * n + col] = (a[row] - b[row]) / (2 * h);
    // Nonuniform xi gives the source a central coefficient as well as
    // the two off-diagonals. Its omission must fail independently.
    const az = flatZero(plus), bz = flatZero(minus), m = frozen[i][1];
    const sourceDiagonal = m.jacobian ** 2 * source[i][1] * firstXi(i)[1] / ((m.alpha + m.gamma) * system.lengthScale);
    close(matrix[col * n + col] - (az[col] - bz[col]) / (2 * h), sourceDiagonal, 2e-9);
  }
  const expected = solveLinear(matrix, base.map(value => -value)), next = system.sweep(nodes, 1).nodes;
  for (let col = 0; col < n; col++) {
    const i = 1 + Math.floor(col / 2), key = col % 2 ? 'y' : 'x';
    close(next[i][1][key] - nodes[i][1][key], expected[col], 4e-11);
  }
  assert.ok(Math.max(...flat(next).map(Math.abs)) < 5e-13);
  boundary(next, nodes);
});

test('Poisson control scales as inverse physical length squared and does not alias caller inputs', () => {
  const input = patch(8, 5), before = structuredClone(input), system = createEllipticStreamtubeGrid(input);
  const scale = 3.7, angle = .6, c = Math.cos(angle), s = Math.sin(angle);
  const map = p => ({ x: 4 + scale * (c * p.x - s * p.y), y: -3 + scale * (s * p.x + c * p.y) });
  const mapped = createEllipticStreamtubeGrid({ nodes: input.nodes.map(row => row.map(map)),
    massFlows: input.massFlows.map(m => m * 7), streamwiseSource: input.streamwiseSource.map(row => row.map(p => p / scale ** 2)) });
  const a = system.sweep(system.initial, 1.1), b = mapped.sweep(mapped.initial, 1.1);
  a.nodes.forEach((row, i) => row.forEach((p, j) => {
    const expected = map(p); close(b.nodes[i][j].x, expected.x, 3e-13); close(b.nodes[i][j].y, expected.y, 3e-13);
  }));
  // Residual vectors rotate; their componentwise infinity norm need not.
  const originalRows = system.residuals(system.initial).rows, mappedRows = mapped.residuals(mapped.initial).rows;
  originalRows.forEach((row, i) => {
    close(mappedRows[i].x, c * row.x - s * row.y, 2e-12);
    close(mappedRows[i].y, s * row.x + c * row.y, 2e-12);
  });
  assert.deepEqual(input, before);
  const originalResidual = system.residuals(system.initial);
  // Constructor must own all state used later, including the fixed P array.
  input.nodes[2][2].x += .1; input.massFlows[1] *= 2; input.streamwiseSource[2][2] = 1e8;
  assert.deepEqual(system.residuals(system.initial), originalResidual);
  assert.deepEqual(system.sweep(system.initial, 1.1), a);
  assert.deepEqual(system.initial, before.nodes);
});

test('zero control recovers the Laplace baseline exactly and invalid source arrays are rejected', () => {
  const input = patch(5, 3), zero = input.nodes.map(row => row.map(() => 0));
  const baseline = createEllipticStreamtubeGrid({ nodes: input.nodes, massFlows: input.massFlows });
  const controlled = createEllipticStreamtubeGrid({ ...input, streamwiseSource: zero });
  assert.deepEqual(controlled.residuals(input.nodes), baseline.residuals(input.nodes));
  assert.deepEqual(controlled.sweep(input.nodes), baseline.sweep(input.nodes));
  assert.equal(controlled.coordinateEquations.xi, 'Laplace');
  for (const invalid of [null, [], zero.slice(1), zero.map(row => row.slice(1)),
    zero.map((row, i) => row.map((p, j) => i === 2 && j === 1 ? NaN : p)),
    zero.map((row, i) => row.map((p, j) => i === 2 && j === 1 ? Infinity : p))]) {
    assert.throws(() => createEllipticStreamtubeGrid({ ...input, streamwiseSource: invalid }), /Poisson source/);
  }
});
