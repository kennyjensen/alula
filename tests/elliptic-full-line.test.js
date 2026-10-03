import test from 'node:test';
import assert from 'node:assert/strict';
import { createEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { createOrthogonalBoundaryControl } from '../src/geometry/orthogonal-boundary-control.js';
import { harmonicStagnationPoint } from '../scripts/validation/harmonic-stagnation-stations.js';
import { solveLinear } from '../src/numerics/linear.js';

for (const sourceForm of ['poisson', 'metric-stretch']) for (const streamwiseSourceDiscretization of ['centered', 'grape-1980'])
  for (const damped of [false, true]) test(`full line Jacobian: ${sourceForm}, ${streamwiseSourceDiscretization}, ${damped ? 'damped' : 'direct'} matches nonlinear finite differences`, () => {
    const xi = [0, .17, .35, .5, .61, .83, 1], eta = [0, .025, .1, .4, 1], nx = xi.length - 1, nt = eta.length - 1;
    const nodes = xi.map(u => eta.map(e => harmonicStagnationPoint(2 * u - 1 + .15 * (2 * u - 1) ** 2, e)));
    const field = sourceForm === 'poisson' ? 'poisson' : 'stretch';
    const config = { sourceForm, sides: ['lower'], corners: { lower: [3] }, decay: { lower: 10, upper: 2 }, background: nodes.map(row => row.map(() => 0)) };
    if (damped) config.previous = { lower: createOrthogonalBoundaryControl({ nodes, xi, eta, ...config }).lower.map(q => q ? q[field] - .1 : 0) };
    const options = { nodes, massFlows: eta.slice(1).map((e, j) => e - eta[j]), streamwiseCoordinates: xi,
      discretization: 'giles-1985', streamwiseSourceDiscretization, orthogonalBoundaryControl: config, lineLinearization: 'full-metrics' };
    const system = createEllipticStreamtubeGrid(options), expected = structuredClone(nodes), n = 2 * (nx - 1), h = 3e-7, omega = .15;
    for (let j = 1; j < nt; j++) {
      const residual = trial => {
        // Fully nonlinear metrics and boundary sources, including J²*P.
        // Undo diagnostic normalization so its changing denominator does
        // not contribute to this Jacobian of the actual inverse PDE.
        const metrics = system.metrics(trial);
        return system.residuals(trial, metrics).rows.filter(row => row.j === j).flatMap(row => {
          const c = metrics[row.i][j], scale = (c.alpha + c.gamma) * system.lengthScale;
          return [row.x * scale, row.y * scale];
        });
      };
      const base = residual(expected), matrix = new Float64Array(n * n);
      for (let c = 0; c < n; c++) {
        const i = 1 + Math.floor(c / 2), key = c % 2 ? 'y' : 'x';
        const plus = structuredClone(expected), minus = structuredClone(expected); plus[i][j][key] += h; minus[i][j][key] -= h;
        const a = residual(plus), b = residual(minus);
        for (let r = 0; r < n; r++) matrix[r * n + c] = (a[r] - b[r]) / (2 * h);
      }
      const delta = solveLinear(matrix, base.map(v => -v));
      delta.forEach((v, c) => { expected[1 + Math.floor(c / 2)][j][c % 2 ? 'y' : 'x'] += omega * v; });
    }
    const actual = system.sweep(nodes, omega).nodes;
    actual.forEach((row, i) => row.forEach((p, j) => {
      for (const key of ['x', 'y']) assert.ok(Math.abs(p[key] - expected[i][j][key]) < 3e-9, `(${i},${j}) ${key}: ${p[key]} != ${expected[i][j][key]}`);
      if (!i || i === nx || !j || j === nt) assert.deepEqual(p, nodes[i][j]);
    }));
    assert.equal(system.coordinateEquations.lineLinearization, 'full-metrics');
    const frozen = createEllipticStreamtubeGrid({ ...options, lineLinearization: 'frozen-metrics' }).sweep(nodes, omega).nodes;
    const difference = Math.max(...actual.flatMap((row, i) => row.map((p, j) => Math.hypot(p.x - frozen[i][j].x, p.y - frozen[i][j].y))));
    assert.ok(difference > 1e-5, `This case must distinguish metric derivatives: ${difference}`);
  });

test('full line metrics reject unsupported discretizations and missing boundary-control equations', () => {
  const nodes = Array.from({ length: 4 }, (_, i) => Array.from({ length: 4 }, (_, j) => ({ x: i / 3, y: j / 3 })));
  const options = { nodes, massFlows: [1, 1, 1], lineLinearization: 'full-metrics' };
  assert.throws(() => createEllipticStreamtubeGrid(options), /Full line metrics/);
  assert.throws(() => createEllipticStreamtubeGrid({ ...options, discretization: 'giles-1985' }), /Full line metrics/);
});
