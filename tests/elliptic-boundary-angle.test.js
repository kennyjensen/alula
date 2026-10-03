import test from 'node:test';
import assert from 'node:assert/strict';
import { createEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { createOrthogonalBoundaryFeedback } from '../src/geometry/orthogonal-boundary-control.js';
import { solveLinear } from '../src/numerics/linear.js';

test('implicit boundary-angle line sweep agrees with an independent dense finite-difference Jacobian', () => {
  const xi = [0, .09, .28, .52, .79, 1], eta = [0, .17, .52, 1], nx = xi.length - 1, nt = eta.length - 1;
  const nodes = xi.map(u => eta.map(e => ({ x: u + .035 * Math.sin(Math.PI * u) * e * (1 - e), y: e + .03 * u })));
  const background = xi.map(u => eta.map(e => .2 * u * (1 + e))), config = { background, sides: ['lower', 'upper'] };
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: eta.slice(1).map((e, j) => e - eta[j]),
    streamwiseCoordinates: xi, discretization: 'giles-1985', orthogonalBoundaryControl: config });
  const feedback = createOrthogonalBoundaryFeedback({ nodes, xi, eta, ...config });
  const expected = structuredClone(nodes), omega = .8, count = 2 * (nx - 1), h = 1e-6;
  for (let j = 1; j < nt; j++) {
    const frozen = system.metrics(expected);
    // Recompute F on each perturbed state, holding alpha,beta,gamma fixed.
    // This is a dense oracle for the Picard metrics/Newton control split.
    const residual = trial => {
      const coefficients = frozen.map((row, i) => row.map((c, k) => !c ? c : { ...c,
        streamwiseDrift: c.alpha * feedback.evaluate(trial, i, k).stretch }));
      return system.residuals(trial, coefficients).rows.filter(row => row.j === j).flatMap(row => [row.x, row.y]);
    };
    const base = residual(expected), matrix = new Float64Array(count * count);
    for (let c = 0; c < count; c++) {
      const i = 1 + Math.floor(c / 2), key = c % 2 ? 'y' : 'x';
      const plus = structuredClone(expected), minus = structuredClone(expected); plus[i][j][key] += h; minus[i][j][key] -= h;
      const a = residual(plus), b = residual(minus);
      for (let r = 0; r < count; r++) matrix[r * count + c] = (a[r] - b[r]) / (2 * h);
    }
    const correction = solveLinear(matrix, base.map(v => -v));
    correction.forEach((v, c) => { expected[1 + Math.floor(c / 2)][j][c % 2 ? 'y' : 'x'] += omega * v; });
  }
  const actual = system.sweep(nodes, omega).nodes;
  actual.forEach((row, i) => row.forEach((p, j) => {
    for (const key of ['x', 'y']) assert.ok(Math.abs(p[key] - expected[i][j][key]) < 5e-10, `${i},${j},${key}`);
    if (!i || i === nx || !j || j === nt) assert.deepEqual(p, nodes[i][j]);
  }));
  const before = system.residuals(nodes); background[2][1] = 90; config.sides.length = 0;
  assert.deepEqual(system.residuals(nodes), before);
});

test('orthogonal affine grid is an exact fixed point and incompatible controls are rejected', () => {
  const nodes = Array.from({ length: 6 }, (_, i) => Array.from({ length: 5 }, (_, j) => ({ x: i + 2, y: j - 1 })));
  const background = nodes.map(row => row.map(() => 0));
  const args = { nodes, massFlows: [1, 1, 1, 1], orthogonalBoundaryControl: { background }, discretization: 'giles-1985' };
  const system = createEllipticStreamtubeGrid(args);
  assert.ok(system.residuals(nodes).residual < 1e-12);
  system.sweep(nodes).nodes.forEach((row, i) => row.forEach((p, j) => {
    assert.ok(Math.abs(p.x - nodes[i][j].x) < 1e-12); assert.ok(Math.abs(p.y - nodes[i][j].y) < 1e-12);
  }));
  assert.throws(() => createEllipticStreamtubeGrid({ ...args, streamwiseStretch: background }), /cannot be combined/);
  assert.throws(() => createEllipticStreamtubeGrid({ ...args, boundaryConditions: { lower: 'giles-vertical' } }), /fixed boundaries/);
  assert.throws(() => createEllipticStreamtubeGrid({ ...args, discretization: 'quadratic' }), /Giles discretization/);
});
