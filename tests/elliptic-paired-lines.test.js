import test from 'node:test';
import assert from 'node:assert/strict';
import { createEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { createOrthogonalBoundaryControl } from '../src/geometry/orthogonal-boundary-control.js';
import { harmonicStagnationPoint } from '../scripts/validation/harmonic-stagnation-stations.js';
import { solveLinear } from '../src/numerics/linear.js';

const check = ({ sourceForm, streamwiseSourceDiscretization, damped, corner = false, overlap = false, partial = false }) => {
  const xi = [0, .17, .35, .5, .61, .83, 1], eta = overlap ? [0, .05, .21, .7, 1] : [0, .04, .13, .27, .43, .75, .9, 1];
  const nx = xi.length - 1, nt = eta.length - 1;
  const nodes = xi.map(u => eta.map(e => corner ? harmonicStagnationPoint(2 * u - 1 + .15 * (2 * u - 1) ** 2, e)
    : { x: u + .1 * u ** 2 + .08 * Math.sin(Math.PI * u) * e * (1 - e), y: e + .06 * Math.sin(Math.PI * u) * Math.sin(Math.PI * e) }));
  const config = { sourceForm, sides: corner ? ['lower'] : ['lower', 'upper'], ...(corner ? { corners: { lower: [3] } } : {}),
    decay: { lower: 3, upper: 4 }, background: nodes.map(row => row.map(() => 0)) };
  if (partial) config.activeStations = { lower: xi.map((_, i) => i === 2 || i === 4), upper: xi.map((_, i) => i === 3 || i === 4) };
  if (damped) {
    const field = sourceForm === 'poisson' ? 'poisson' : 'stretch', control = createOrthogonalBoundaryControl({ nodes, xi, eta, ...config });
    config.previous = Object.fromEntries(config.sides.map(side => [side, control[side].map(q => q ? q[field] - .1 : 0)]));
  }
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: eta.slice(1).map((e, j) => e - eta[j]), streamwiseCoordinates: xi,
    discretization: 'giles-1985', streamwiseSourceDiscretization, lineLinearization: 'full-metrics', lineGrouping: 'boundary-pairs', orthogonalBoundaryControl: config });
  const groups = system.coordinateEquations.rowGroups, expectedGroups = overlap ? [[1, 2, 3]] : corner ? [[1, 2], [3], [4], [5], [6]] : [[1, 2], [3], [4], [5, 6]];
  assert.deepEqual(groups, expectedGroups);
  const expected = structuredClone(nodes), h = 3e-7, omega = .1;
  for (const rows of groups) {
    const size = 2 * rows.length, n = size * (nx - 1), analytic = system.linearizeRows(expected, rows);
    const residual = trial => {
      const metrics = system.metrics(trial), samples = system.residuals(trial, metrics).rows;
      return samples.filter(row => rows.includes(row.j)).flatMap(row => {
        const m = metrics[row.i][row.j], scale = (m.alpha + m.gamma) * system.lengthScale;
        return [row.x * scale, row.y * scale];
      });
    };
    const base = residual(expected), matrix = new Float64Array(n * n);
    for (let c = 0; c < n; c++) {
      const i = 1 + Math.floor(c / size), j = rows[Math.floor(c % size / 2)], key = c % 2 ? 'y' : 'x';
      const plus = structuredClone(expected), minus = structuredClone(expected); plus[i][j][key] += h; minus[i][j][key] -= h;
      const a = residual(plus), b = residual(minus);
      for (let r = 0; r < n; r++) {
        const numeric = (a[r] - b[r]) / (2 * h), blockRow = Math.floor(r / size), blockColumn = Math.floor(c / size), offset = blockColumn - blockRow;
        const exact = Math.abs(offset) > 1 ? 0 : [analytic.lower, analytic.diagonal, analytic.upper][offset + 1][blockRow][(r % size) * size + c % size];
        assert.ok(Math.abs(numeric - exact) < 2e-6 * Math.max(1, Math.abs(numeric)), `row ${r}, column ${c}: ${exact} != ${numeric}`);
        matrix[r * n + c] = numeric;
      }
    }
    const delta = solveLinear(matrix, base.map(v => -v));
    for (let c = 0; c < n; c++) expected[1 + Math.floor(c / size)][rows[Math.floor(c % size / 2)]][c % 2 ? 'y' : 'x'] += omega * delta[c];
  }
  system.sweep(nodes, omega).nodes.forEach((row, i) => row.forEach((p, j) => {
    for (const key of ['x', 'y']) assert.ok(Math.abs(p[key] - expected[i][j][key]) < 2e-8, `(${i},${j}) ${key}`);
    if (!i || i === nx || !j || j === nt) assert.deepEqual(p, nodes[i][j]);
  }));
};

for (const sourceForm of ['poisson', 'metric-stretch']) for (const streamwiseSourceDiscretization of ['centered', 'grape-1980'])
  test(`paired boundary rows ${sourceForm}/${streamwiseSourceDiscretization} match a full nonlinear dense Jacobian`, () => check({ sourceForm, streamwiseSourceDiscretization }));
for (const streamwiseSourceDiscretization of ['centered', 'grape-1980'])
  test(`paired damped P rows ${streamwiseSourceDiscretization} include source and metric coupling`, () => check({ sourceForm: 'poisson', streamwiseSourceDiscretization, damped: true }));
test('paired rows include neighboring-station corner-source derivatives', () => check({ sourceForm: 'poisson', streamwiseSourceDiscretization: 'grape-1980', corner: true }));
test('overlapping boundary pairs merge and retain both boundary controls', () => check({ sourceForm: 'poisson', streamwiseSourceDiscretization: 'centered', overlap: true }));
for (const sourceForm of ['poisson', 'metric-stretch'])
  test(`paired rows with partial ${sourceForm} boundary support match independent dense derivatives`, () =>
    check({ sourceForm, streamwiseSourceDiscretization: 'centered', partial: true }));
