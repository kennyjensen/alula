import test from 'node:test';
import assert from 'node:assert/strict';
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { solveLinear } from '../src/numerics/linear.js';

const close = (a, b, tol) => assert.ok(Math.abs(a - b) < tol, `${a} != ${b}`);
const masses = eta => eta.slice(1).map((v, j) => v - eta[j]);
const affine = (nx, eta, perturb = 0) => Array.from({ length: nx + 1 }, (_, i) => eta.map((v, j) => {
  const u = i / nx, d = perturb * Math.sin(Math.PI * u) * Math.sin(Math.PI * v);
  return { x: 2 * u + .4 * v + d, y: .3 * u + v - .6 * d };
}));
const assertBoundary = (actual, expected) => actual.forEach((row, i) => row.forEach((p, j) => {
  if (!i || i === actual.length - 1 || !j || j === row.length - 1) assert.deepEqual(p, expected[i][j]);
}));

test('inverse Laplace SLOR preserves and recovers a sheared affine grid on nonuniform mass coordinates', () => {
  const eta = [0, .07, .19, .43, .72, 1], exact = affine(12, eta), nodes = affine(12, eta, .035);
  // Keep the prescribed boundary bit-for-bit; sin(pi) is not exactly zero.
  nodes.forEach((row, i) => row.forEach((p, j) => {
    if (!i || i === 12 || !j || j === eta.length - 1) row[j] = { ...exact[i][j] };
  }));
  const before = structuredClone(nodes), system = createEllipticStreamtubeGrid({ nodes, massFlows: masses(eta) });
  const exactSystem = createEllipticStreamtubeGrid({ nodes: exact, massFlows: masses(eta) });
  assert.ok(exactSystem.residuals(exact).residual < 2e-13);
  const r = smoothEllipticStreamtubeGrid(system, { tolerance: 1e-11 });
  assert.equal(r.converged, true, r.reason); assert.equal(r.quality.valid, true);
  r.nodes.forEach((row, i) => row.forEach((p, j) => {
    close(p.x, exact[i][j].x, 1e-11); close(p.y, exact[i][j].y, 1e-11);
  }));
  assertBoundary(r.nodes, before); assert.deepEqual(nodes, before); assert.deepEqual(system.initial, before);
});

test('a frozen nonuniform SLOR line agrees with an independent dense finite-difference solve', () => {
  const eta = [0, .27, 1], nodes = affine(7, eta, .08);
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: masses(eta) }), frozen = system.metrics(nodes);
  const flatten = p => system.residuals(p, frozen).rows.flatMap(r => [r.x, r.y]);
  const base = flatten(nodes), n = base.length, matrix = new Float64Array(n * n), h = 1e-5;
  for (let k = 0; k < n; k++) {
    const i = 1 + Math.floor(k / 2), key = k % 2 ? 'y' : 'x';
    const plus = structuredClone(nodes), minus = structuredClone(nodes);
    plus[i][1][key] += h; minus[i][1][key] -= h;
    const a = flatten(plus), b = flatten(minus);
    for (let row = 0; row < n; row++) matrix[row * n + k] = (a[row] - b[row]) / (2 * h);
  }
  const delta = solveLinear(matrix, base.map(v => -v)), next = system.sweep(nodes, 1).nodes;
  for (let k = 0; k < n; k++) {
    const i = 1 + Math.floor(k / 2), key = k % 2 ? 'y' : 'x';
    close(next[i][1][key] - nodes[i][1][key], delta[k], 2e-11);
  }
  assert.ok(Math.max(...flatten(next).map(Math.abs)) < 1e-12);
  assertBoundary(next, nodes);
});

test('all geometric coefficients remain frozen across a complete overrelaxed line sweep', () => {
  const eta = [0, .13, .38, .69, 1], nodes = affine(6, eta, .04), omega = 1.3;
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: masses(eta) }), frozen = system.metrics(nodes);
  const locations = [];
  for (let i = 1; i < 6; i++) for (let j = 1; j < eta.length - 1; j++) for (const key of ['x', 'y']) locations.push({ i, j, key });
  const flat = p => system.residuals(p, frozen).rows.flatMap(r => [r.x, r.y]), base = flat(nodes), n = base.length;
  const matrix = new Float64Array(n * n), delta = new Float64Array(n), h = 1e-5;
  locations.forEach(({ i, j, key }, col) => {
    const a = structuredClone(nodes), b = structuredClone(nodes); a[i][j][key] += h; b[i][j][key] -= h;
    const plus = flat(a), minus = flat(b);
    for (let row = 0; row < n; row++) matrix[row * n + col] = (plus[row] - minus[row]) / (2 * h);
  });
  for (let j = 1; j < eta.length - 1; j++) {
    const indices = locations.flatMap((p, k) => p.j === j ? [k] : []);
    const a = indices.flatMap(row => indices.map(col => matrix[row * n + col]));
    const b = indices.map(row => -base[row] - delta.reduce((sum, v, col) => sum + matrix[row * n + col] * v, 0));
    const line = solveLinear(a, b); indices.forEach((k, p) => { delta[k] += omega * line[p]; });
  }
  const r = system.sweep(nodes, omega);
  locations.forEach(({ i, j, key }, k) => close(r.nodes[i][j][key] - nodes[i][j][key], delta[k], 3e-11));
});

test('mass-coordinate SLOR refines to the exact harmonic coordinates and streamlines of a circular vortex', t => {
  let previous = Infinity; const evidence = [];
  for (const [nx, nt] of [[8, 4], [16, 8], [32, 16]]) {
    const eta = Array.from({ length: nt + 1 }, (_, j) => Math.expm1(.8 * j / nt) / Math.expm1(.8));
    const exact = Array.from({ length: nx + 1 }, (_, i) => eta.map(e => {
      const theta = -.4 + .8 * i / nx, r = 2 * Math.exp(-.5 * e);
      return { x: r * Math.cos(theta), y: r * Math.sin(theta) };
    }));
    const nodes = exact.map((row, i) => row.map((p, j) => !i || i === nx || !j || j === nt ? { ...p }
      : { x: (1 - eta[j]) * row[0].x + eta[j] * row[nt].x, y: (1 - eta[j]) * row[0].y + eta[j] * row[nt].y }));
    const system = createEllipticStreamtubeGrid({ nodes, massFlows: masses(eta) });
    const result = smoothEllipticStreamtubeGrid(system, { tolerance: 1e-10, maxSweeps: 600 });
    assert.equal(result.converged, true, result.reason); assertBoundary(result.nodes, nodes);
    let error = 0;
    result.nodes.forEach((row, i) => row.forEach((p, j) => {
      // Independent physical harmonic functions: xi=(theta+.4)/.8,
      // eta=-2 log(r/2). Neither invokes the numerical grid residual.
      error = Math.max(error, Math.abs((Math.atan2(p.y, p.x) + .4) / .8 - i / nx),
        Math.abs(-2 * Math.log(Math.hypot(p.x, p.y) / 2) - eta[j]));
    }));
    assert.ok(error < .3 * previous, `harmonic-coordinate refinement: ${error} after ${previous}`);
    previous = error;
    evidence.push({ nx, nt, sweeps: result.history.length - 1, residual: result.history.at(-1).residual, error, minCornerSine: result.quality.minCornerSine });
  }
  assert.ok(previous < 1e-4); t.diagnostic(JSON.stringify(evidence));
});

test('SLOR preserves rotation, translation, length and total-mass scaling, and reports unsolved or invalid grids', () => {
  const eta = [0, .12, .35, .6, 1], nodes = affine(8, eta, .03), before = structuredClone(nodes);
  const angle = .7, scale = 2.3, c = Math.cos(angle), s = Math.sin(angle);
  const map = p => ({ x: 4 + scale * (p.x * c - p.y * s), y: -3 + scale * (p.x * s + p.y * c) });
  const original = createEllipticStreamtubeGrid({ nodes, massFlows: masses(eta) });
  const mapped = createEllipticStreamtubeGrid({ nodes: nodes.map(row => row.map(map)), massFlows: masses(eta).map(m => 7 * m) });
  const a = original.sweep(nodes).nodes, b = mapped.sweep(mapped.initial).nodes;
  a.forEach((row, i) => row.forEach((p, j) => { const q = map(p); close(q.x, b[i][j].x, 2e-13); close(q.y, b[i][j].y, 2e-13); }));
  assert.deepEqual(nodes, before);
  const stopped = smoothEllipticStreamtubeGrid(original, { maxSweeps: 0 });
  assert.equal(stopped.converged, false); assert.equal(stopped.reason, 'sweep limit'); assert.deepEqual(stopped.nodes, nodes);
  const reversed = createEllipticStreamtubeGrid({ nodes: affine(5, eta).map(row => row.map(p => ({ x: -p.x, y: p.y }))), massFlows: masses(eta) });
  const invalid = smoothEllipticStreamtubeGrid(reversed);
  assert.equal(invalid.converged, false); assert.equal(invalid.reason, 'folded converged grid');
  assert.throws(() => original.sweep(nodes, 2), /omega/);
  assert.throws(() => createEllipticStreamtubeGrid({ nodes, massFlows: [1, 0, 1, 1] }), /Invalid/);
  const movedBoundary = structuredClone(nodes); movedBoundary[0][1].x += .01;
  assert.throws(() => original.sweep(movedBoundary), /boundary/);
});
