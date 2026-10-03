import test from 'node:test';
import assert from 'node:assert/strict';
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid, gilesStreamwiseCoordinates } from '../src/geometry/elliptic-streamtube-grid.js';
import { solveLinear } from '../src/numerics/linear.js';

const close = (a, b, tolerance = 2e-11) => assert.ok(Math.abs(a - b) < tolerance, `${a} differs from ${b}`);
const masses = eta => eta.slice(1).map((e, j) => e - eta[j]);
const specimen = () => {
  const xi = [0, .08, .23, .48, .77, 1], eta = [0, .09, .28, .65, 1];
  const exact = xi.map(u => eta.map(v => ({ x: 2 * u + .4 * v, y: .3 * u + v })));
  const nodes = exact.map((row, i) => row.map((p, j) => !i || i === xi.length - 1 || !j || j === eta.length - 1 ? { ...p }
    : { x: p.x + .025 * Math.sin(Math.PI * xi[i]) * Math.sin(Math.PI * eta[j]),
      y: p.y - .012 * Math.sin(Math.PI * xi[i]) * Math.sin(Math.PI * eta[j]) }));
  return { xi, eta, nodes, exact };
};

test('original ISET station construction follows fourth-root segment spacing and geometric scaling', () => {
  const path = [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: 1, y: 16 }, { x: -80, y: 16 }];
  const actual = gilesStreamwiseCoordinates(path), expected = [0, 1 / 6, 1 / 2, 1];
  actual.forEach((v, i) => close(v, expected[i], 1e-15));
  const moved = path.map(p => ({ x: 20 - 7 * p.y, y: -31 + 7 * p.x }));
  gilesStreamwiseCoordinates(moved).forEach((v, i) => close(v, actual[i], 1e-15));
  assert.throws(() => gilesStreamwiseCoordinates([path[0], path[0], path[1]]), /nonzero/);
});

// Independent mathematical transcription of the centered-secant / edge-slope
// stencil in Giles's printed ELLIP (pp.187–188). It does not call production
// metric, derivative, residual or line-assembly helpers. Dense perturbation
// columns below test the actual tridiagonal assembly and update order.
function frozenOperator(nodes, xi, eta, i, j) {
  const hm = xi[i] - xi[i - 1], hp = xi[i + 1] - xi[i], km = eta[j] - eta[j - 1], kp = eta[j + 1] - eta[j];
  const hx = (hm + hp) / 2, hy = (km + kp) / 2;
  const dxi = key => (nodes[i + 1][j][key] - nodes[i - 1][j][key]) / (2 * hx);
  const deta = key => (nodes[i][j + 1][key] - nodes[i][j - 1][key]) / (2 * hy);
  const alpha = deta('x') ** 2 + deta('y') ** 2, beta = deta('x') * dxi('x') + deta('y') * dxi('y');
  const gamma = dxi('x') ** 2 + dxi('y') ** 2;
  return (grid, key) => {
    const center = grid[i][j][key];
    const xx = ((grid[i + 1][j][key] - center) / hp - (center - grid[i - 1][j][key]) / hm) / hx;
    const yy = ((grid[i][j + 1][key] - center) / kp - (center - grid[i][j - 1][key]) / km) / hy;
    const xy = (grid[i + 1][j + 1][key] - grid[i - 1][j + 1][key] - grid[i + 1][j - 1][key] + grid[i - 1][j - 1][key]) / (4 * hx * hy);
    return alpha * xx - 2 * beta * xy + gamma * yy;
  };
}

test('Giles-mode nonuniform line sweeps match independently assembled dense equations with metrics updated between lines', () => {
  const { nodes, xi, eta } = specimen(), expected = structuredClone(nodes), omega = 1.3, n = xi.length - 2;
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: masses(eta), streamwiseCoordinates: xi, discretization: 'giles-1985' });
  for (let j = 1; j < eta.length - 1; j++) {
    const operators = Array.from({ length: n }, (_, k) => frozenOperator(expected, xi, eta, k + 1, j));
    const matrix = new Float64Array(n * n), step = 1e-5;
    for (let col = 0; col < n; col++) {
      const plus = structuredClone(expected), minus = structuredClone(expected);
      plus[col + 1][j].x += step; minus[col + 1][j].x -= step;
      for (let row = 0; row < n; row++) matrix[n * row + col] = (operators[row](plus, 'x') - operators[row](minus, 'x')) / (2 * step);
    }
    for (const key of ['x', 'y']) {
      const correction = solveLinear(matrix, operators.map(at => -at(expected, key)));
      for (let k = 0; k < n; k++) expected[k + 1][j][key] += omega * correction[k];
    }
  }
  const actual = system.sweep(nodes, omega).nodes;
  actual.forEach((row, i) => row.forEach((p, j) => { close(p.x, expected[i][j].x); close(p.y, expected[i][j].y); }));
  assert.equal(system.coordinateEquations.metricUpdate, 'each line');
});

test('Giles-mode SLOR recovers exact sheared affine coordinates with independently stretched xi and mass labels', () => {
  const { nodes, xi, eta, exact } = specimen();
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: masses(eta), streamwiseCoordinates: xi, discretization: 'giles-1985' });
  const before = structuredClone(system.xi); xi[1] = .02;
  assert.deepEqual(system.xi, before);
  const result = smoothEllipticStreamtubeGrid(system, { tolerance: 1e-11 });
  assert.equal(result.converged, true, result.reason); assert.equal(result.quality.valid, true);
  result.nodes.forEach((row, i) => row.forEach((p, j) => {
    close(p.x, exact[i][j].x); close(p.y, exact[i][j].y);
    if (!i || i === nodes.length - 1 || !j || j === eta.length - 1)
      assert.deepEqual(p, nodes[i][j]);
  }));
  assert.deepEqual(result.nodes.at(-1), nodes.at(-1));
  assert.throws(() => createEllipticStreamtubeGrid({ nodes, massFlows: masses(eta), streamwiseCoordinates: [0, .2, .2, .5, .8, 1] }), /strictly/);
});
