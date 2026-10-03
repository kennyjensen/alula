import test from 'node:test';
import assert from 'node:assert/strict';
import { createOrthogonalBoundaryControl, createOrthogonalBoundaryFeedback, relaxBoundaryStretch } from '../src/geometry/orthogonal-boundary-control.js';
import { createEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { harmonicStagnationPoint } from '../scripts/validation/harmonic-stagnation-stations.js';
import { solveLinear } from '../src/numerics/linear.js';
const close = (a, b, tolerance = 2e-10) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);
const fixture = () => {
  const xi = [0, .2, .4, .5, .6, .8, 1], eta = [0, .025, .1, .4, 1];
  return { xi, eta, nodes: xi.map(u => eta.map(e => harmonicStagnationPoint(2 * u - 1, e))),
    sides: ['lower'], corners: { lower: [3] } };
};

test('a marked stagnation corner averages neighboring controls without asserting a corner normal', () => {
  const args = fixture(), result = createOrthogonalBoundaryControl(args);
  assert.equal(result.lower[3].corner, true);
  close(result.lower[3].stretch, .5 * (result.lower[2].stretch + result.lower[4].stretch));
  close(result.lower[3].stretch, 0);
  assert.equal(result.lower[3].normal, undefined);
  assert.deepEqual(result.corners, { lower: [3] });
  for (const corners of [{ lower: [1] }, { lower: [5] }, { lower: [2, 3] }, { upper: [3] }])
    assert.throws(() => createOrthogonalBoundaryControl({ ...args, corners }), /isolated stations/);
});

test('averaged corner sensitivities act on neighboring stations and agree with independent differences', () => {
  const args = fixture(), background = args.nodes.map(row => row.map(() => 0));
  const control = createOrthogonalBoundaryFeedback({ ...args, background, decay: { lower: 10, upper: 2 } });
  const h = 1e-6;
  for (const j of [1, 2]) {
    const baseline = control.evaluate(args.nodes, 3, j);
    assert.deepEqual(baseline.derivative, { x: 0, y: 0 });
    for (const [side, station] of [2, 4].entries()) for (const key of ['x', 'y']) {
      const plus = structuredClone(args.nodes), minus = structuredClone(args.nodes); plus[station][j][key] += h; minus[station][j][key] -= h;
      close(baseline.neighborDerivatives[side][key],
        (control.evaluate(plus, 3, j).stretch - control.evaluate(minus, 3, j).stretch) / (2 * h), 2e-6);
    }
  }
  const before = control.evaluate(args.nodes, 3, 1); args.corners.lower.length = 0;
  assert.deepEqual(control.evaluate(args.nodes, 3, 1), before);
});

for (const damped of [false, true]) test(`${damped ? 'Damped' : 'Direct'} corner coupling matches a dense sweep Jacobian`, () => {
  const args = fixture(), { xi, eta, nodes } = args, nx = xi.length - 1, nt = eta.length - 1;
  const config = { background: nodes.map(row => row.map(() => 0)), sides: args.sides, corners: args.corners, decay: { lower: 10, upper: 2 } };
  if (damped) config.previous = { lower: createOrthogonalBoundaryControl(args).lower.map(q => q ? q.stretch - .1 : 0) };
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: eta.slice(1).map((e, j) => e - eta[j]), streamwiseCoordinates: xi,
    discretization: 'giles-1985', orthogonalBoundaryControl: config });
  const feedback = createOrthogonalBoundaryFeedback({ ...args, ...config }), expected = structuredClone(nodes);
  const n = 2 * (nx - 1), h = 1e-6, omega = .7;
  for (let j = 1; j < nt; j++) {
    const fixed = system.metrics(expected);
    const residual = trial => {
      const coefficients = fixed.map((row, i) => row.map((c, k) => !c ? c : { ...c,
        streamwiseDrift: c.alpha * feedback.evaluate(trial, i, k).stretch }));
      return system.residuals(trial, coefficients).rows.filter(row => row.j === j).flatMap(row => [row.x, row.y]);
    };
    const initialResidual = residual(expected), matrix = new Float64Array(n * n);
    for (let c = 0; c < n; c++) {
      const i = 1 + Math.floor(c / 2), key = c % 2 ? 'y' : 'x';
      const plus = structuredClone(expected), minus = structuredClone(expected); plus[i][j][key] += h; minus[i][j][key] -= h;
      const a = residual(plus), b = residual(minus);
      for (let r = 0; r < n; r++) matrix[r * n + c] = (a[r] - b[r]) / (2 * h);
    }
    const delta = solveLinear(matrix, initialResidual.map(v => -v));
    delta.forEach((d, c) => { expected[1 + Math.floor(c / 2)][j][c % 2 ? 'y' : 'x'] += omega * d; });
  }
  const actual = system.sweep(nodes, omega).nodes;
  actual.forEach((row, i) => row.forEach((p, j) => {
    close(p.x, expected[i][j].x, 2e-9); close(p.y, expected[i][j].y, 2e-9);
    if (!i || i === nx || !j || j === nt) assert.deepEqual(p, nodes[i][j]);
  }));
});

test('source update uses the published limiter and averages corners after limiting their neighbors', () => {
  for (const [requested, previous, value, sensitivity] of [[8, 1, 2, 0], [-4, -2, -2.6, .3], [.5, 0, .15, .3]]) {
    const result = relaxBoundaryStretch({ requested, previous }); close(result.value, value); close(result.sensitivity, sensitivity);
  }
  assert.throws(() => relaxBoundaryStretch({ requested: 1, previous: 0, relaxation: 0 }), /Invalid/);
  const args = fixture(), raw = createOrthogonalBoundaryControl(args), previous = { lower: args.xi.map(() => 0) };
  previous.lower[2] = 100;
  const changed = createOrthogonalBoundaryControl({ ...args, previous });
  const expected = [2, 4].map(i => relaxBoundaryStretch({ requested: raw.lower[i].stretch, previous: previous.lower[i] }).value);
  close(changed.lower[3].stretch, .5 * (expected[0] + expected[1]));
  assert.notEqual(changed.lower[3].stretch, relaxBoundaryStretch({ requested: raw.lower[3].stretch, previous: 0 }).value);
  const limited = changed.lower[4]; assert.equal(limited.limited, true);
  assert.ok(limited.interiorDerivatives.every(d => d.x === 0 && d.y === 0));
  assert.throws(() => createOrthogonalBoundaryControl({ ...args, previous: {} }), /match every selected/);
});
