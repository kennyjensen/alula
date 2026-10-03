import test from 'node:test';
import assert from 'node:assert/strict';
import { createOrthogonalBoundaryControl, createOrthogonalBoundaryFeedback } from '../src/geometry/orthogonal-boundary-control.js';
import { createEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { harmonicStagnationPoint } from '../scripts/validation/harmonic-stagnation-stations.js';
import { solveLinear } from '../src/numerics/linear.js';
const close = (a, b, tolerance = 2e-10) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);

test('boundary P uses the physical tangential metric and scales as inverse length squared', () => {
  const xi = [0, .1, .3, .6, .8, 1], eta = [0, .1, .4, .7, 1];
  const nodes = xi.map(u => eta.map(e => ({ x: u + .3 * u * u, y: 2 * e + .15 * u })));
  const f = createOrthogonalBoundaryControl({ nodes, xi, eta });
  const p = createOrthogonalBoundaryControl({ nodes, xi, eta, sourceForm: 'poisson' });
  for (const side of ['lower', 'upper']) for (let i = 1; i < xi.length - 1; i++) {
    const tangentX = 1 + .6 * xi[i] + .3 * (xi[i + 1] - 2 * xi[i] + xi[i - 1]);
    const gamma = tangentX ** 2 + .15 ** 2;
    close(p[side][i].poisson, f[side][i].stretch / gamma);
    for (let j = 0; j < 2; j++) for (const key of ['x', 'y'])
      close(p[side][i].interiorDerivatives[j][key], f[side][i].interiorDerivatives[j][key] / gamma);
  }
  const transform = p => ({ x: 4 + 3 * (.8 * p.x - .6 * p.y), y: -2 + 3 * (.6 * p.x + .8 * p.y) });
  const scaled = createOrthogonalBoundaryControl({ nodes: nodes.map(row => row.map(transform)), xi, eta, sourceForm: 'poisson' });
  for (const side of ['lower', 'upper']) for (let i = 1; i < xi.length - 1; i++) close(scaled[side][i].poisson, p[side][i].poisson / 9);
});

test('corner P averages physical sources from neighbors with different tangential metrics', () => {
  const xi = [0, .12, .31, .5, .66, .85, 1], eta = [0, .03, .12, .4, 1];
  const nodes = xi.map(u => eta.map(e => harmonicStagnationPoint(2 * (u - .5) + .8 * (u - .5) ** 2, e)));
  const control = createOrthogonalBoundaryControl({ nodes, xi, eta, sides: ['lower'], corners: { lower: [3] }, sourceForm: 'poisson' });
  close(control.lower[3].poisson, .5 * (control.lower[2].poisson + control.lower[4].poisson));
  assert.notEqual(control.lower[2].gamma, control.lower[4].gamma);
  assert.equal(control.lower[3].stretch, undefined);
});

test('interior inverse-Poisson drift is J squared times P on a sheared affine grid', () => {
  const xi = [0, .1, .35, .7, 1], eta = [0, .12, .43, 1];
  const nodes = xi.map(u => eta.map(e => ({ x: 1.4 * u + .2 * e, y: .3 * u + e })));
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: eta.slice(1).map((e, j) => e - eta[j]), streamwiseCoordinates: xi,
    discretization: 'giles-1985', orthogonalBoundaryControl: { sourceForm: 'poisson', background: nodes.map(row => row.map(() => 0)) } });
  // Exact metrics: alpha=1.04, beta=.58, gamma=2.05, J=1.34.
  // On either straight boundary P = +/- m*beta/J². Extending P linearly
  // gives this closed drift; using alpha*P would fail on this shear.
  const mL = 2 * (1 / eta[1] + 1 / eta[2]), mU = 2 * (1 / (1 - eta[2]) + 1 / (1 - eta[1]));
  for (const row of system.residuals(nodes).rows) {
    const drift = .58 * (-(1 - eta[row.j]) * mL + eta[row.j] * mU) / ((1.04 + 2.05) * system.lengthScale);
    close(row.x, 1.4 * drift); close(row.y, .3 * drift);
  }
  assert.equal(system.coordinateEquations.sourceUnits, 'inverse physical length squared');
});

for (const streamwiseSourceDiscretization of ['centered', 'grape-1980']) for (const damped of [false, true]) test(`${streamwiseSourceDiscretization} ${damped ? 'damped' : 'direct'} Poisson-source sweep matches an independent dense Jacobian`, () => {
  const xi = [0, .2, .4, .5, .6, .8, 1], eta = [0, .025, .1, .4, 1], nx = xi.length - 1, nt = eta.length - 1;
  const nodes = xi.map(u => eta.map(e => harmonicStagnationPoint(2 * u - 1 + (streamwiseSourceDiscretization === 'grape-1980' ? .15 * (2 * u - 1) ** 2 : 0), e)));
  const config = { sourceForm: 'poisson', sides: ['lower'], corners: { lower: [3] }, decay: { lower: 10, upper: 2 }, background: nodes.map(row => row.map(() => 0)) };
  if (damped) config.previous = { lower: createOrthogonalBoundaryControl({ nodes, xi, eta, ...config }).lower.map(q => q ? q.poisson - .1 : 0) };
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: eta.slice(1).map((e, j) => e - eta[j]), streamwiseCoordinates: xi,
    discretization: 'giles-1985', streamwiseSourceDiscretization, orthogonalBoundaryControl: config });
  const feedback = createOrthogonalBoundaryFeedback({ nodes, xi, eta, ...config });
  const expected = structuredClone(nodes), n = 2 * (nx - 1), h = 1e-6, omega = .7;
  for (let j = 1; j < nt; j++) {
    const frozen = system.metrics(expected);
    const residual = trial => {
      const coefficients = frozen.map((row, i) => row.map((c, k) => !c ? c : { ...c,
        streamwiseDrift: c.jacobian ** 2 * feedback.evaluate(trial, i, k).poisson }));
      return system.residuals(trial, coefficients).rows.filter(row => row.j === j).flatMap(row => [row.x, row.y]);
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
  system.sweep(nodes, omega).nodes.forEach((row, i) => row.forEach((p, j) => { close(p.x, expected[i][j].x, 2e-9); close(p.y, expected[i][j].y, 2e-9); }));
});
