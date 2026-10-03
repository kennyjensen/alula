import test from 'node:test';
import assert from 'node:assert/strict';
import { reconstructOrthogonalBoundary, createOrthogonalBoundaryControl, createOrthogonalBoundaryFeedback } from '../src/geometry/orthogonal-boundary-control.js';
const close = (a, b, tolerance = 2e-11) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);

test('feedback cache follows in-place boundary and interior edits, corners, and rejected trials', () => {
  const xi = [0, .15, .3, .5, .7, .85, 1], eta = [0, .1, .25, .55, .8, .9, 1];
  const nodes = xi.map(x => eta.map(y => ({ x, y })));
  for (const sourceForm of ['poisson', 'metric-stretch']) for (const corners of [{}, { lower: [3], upper: [3] }]) {
    const config = { nodes, xi, eta, sourceForm, corners, background: nodes.map(row => row.map(() => 0)) };
    const cached = createOrthogonalBoundaryFeedback(config);
    const check = grid => {
      const fresh = createOrthogonalBoundaryFeedback({ ...config, nodes: grid });
      for (let j = 1; j < eta.length - 1; j++) for (let k = 1; k < eta.length - 1; k++)
        assert.deepEqual(cached.evaluate(grid, 3, j, k), fresh.evaluate(grid, 3, j, k));
    };
    const grid = structuredClone(nodes); check(grid);
    for (const i of [1, 2, 3, 4, 5]) for (const j of [0, 1, 2, 4, 5, 6]) for (const key of ['x', 'y']) {
      const before = grid[i][j][key]; grid[i][j][key] += 1e-5; check(grid);
      grid[i][j][key] = before; check(grid);
    }
    check(structuredClone(grid));
  }
});

test('nonuniform Hermite reconstruction recovers a cubic normal curve and the harmonic-mass boundary PDE', () => {
  for (const [a, b] of [[.1, .2], [.04, .13], [.12, .17]]) for (const curvature of [0, -1.4, .6]) {
    const lambda = 1.3, speed = 2, tangentialSecond = .7, normalSecond = -curvature * lambda ** 2;
    const point = { x: .4, y: -.3 }, tangent = { x: speed, y: 0 }, secondTangent = { x: .8, y: curvature * speed ** 2 };
    const sample = r => ({ x: point.x + .5 * tangentialSecond * r ** 2 + .43 * r ** 3,
      y: point.y + lambda * r + .5 * normalSecond * r ** 2 - .21 * r ** 3 });
    const result = reconstructOrthogonalBoundary({ point, tangent, secondTangent, firstInterior: sample(a), secondInterior: sample(b), firstDistance: a, secondDistance: b });
    close(result.normalSpeed, lambda); close(result.normalSecond.x, tangentialSecond); close(result.normalSecond.y, normalSecond);
    close(result.stretch, -.8 / speed - speed * tangentialSecond / lambda ** 2);
    for (const key of ['x', 'y']) close(lambda ** 2 * secondTangent[key] + speed ** 2 * result.normalSecond[key]
      + lambda ** 2 * result.stretch * tangent[key], 0, 2e-10);
    // Rotation, translation and physical scaling must not change F.
    const vector = p => ({ x: 3 * (.8 * p.x - .6 * p.y), y: 3 * (.6 * p.x + .8 * p.y) });
    const transform = p => { const q = vector(p); return { x: q.x + 4, y: q.y - 2 }; };
    const changed = reconstructOrthogonalBoundary({ point: transform(point), tangent: vector(tangent), secondTangent: vector(secondTangent),
      firstInterior: transform(sample(a)), secondInterior: transform(sample(b)), firstDistance: a, secondDistance: b });
    close(changed.stretch, result.stretch); close(changed.normalSpeed, 3 * lambda);
  }
});

test('positive branch is continuous at zero curvature and rejects incompatible or singular boundary data', () => {
  const args = { point: { x: 0, y: 0 }, tangent: { x: 1, y: 0 }, secondTangent: { x: 0, y: 0 },
    firstInterior: { x: 0, y: .1 }, secondInterior: { x: 0, y: .2 }, firstDistance: .1, secondDistance: .2 };
  for (const k of [-1e-12, 0, 1e-12]) close(reconstructOrthogonalBoundary({ ...args, secondTangent: { x: 0, y: k } }).normalSpeed, 1);
  assert.throws(() => reconstructOrthogonalBoundary({ ...args, tangent: { x: 0, y: 0 } }), /Singular/);
  assert.throws(() => reconstructOrthogonalBoundary({ ...args, firstDistance: .2 }), /ordered/);
  assert.throws(() => reconstructOrthogonalBoundary({ ...args, normalSign: -1 }), /positive/);
  // m=30, Kn=30: curvature=7.5 is the double root, 8 has no real root.
  for (const k of [7.5, 8]) assert.throws(() => reconstructOrthogonalBoundary({ ...args, secondTangent: { x: 0, y: k } }), /resolved positive/);
});

test('exact cylinder wall normal speed converges under independent physical mass refinement', t => {
  const evidence = [];
  for (const theta of [Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4]) {
    const normal = { x: Math.cos(theta), y: Math.sin(theta) }, tangent = { x: Math.sin(theta), y: -Math.cos(theta) };
    const point = psi => { const q = psi / Math.sin(theta), r = .5 * (q + Math.hypot(q, 2)); return { x: r * normal.x, y: r * normal.y }; };
    const exact = 1 / (2 * Math.sin(theta)); let previous = Infinity;
    for (const h of [.1, .05, .025]) {
      const control = reconstructOrthogonalBoundary({ point: normal, tangent, secondTangent: { x: -normal.x, y: -normal.y },
        firstInterior: point(h), secondInterior: point(2.4 * h), firstDistance: h, secondDistance: 2.4 * h });
      const error = Math.abs(control.normalSpeed - exact);
      assert.ok(error < previous * .2, `${previous} -> ${error}`); previous = error;
      close(control.stretch, 0); evidence.push({ theta, h, error });
    }
  }
  t.diagnostic(JSON.stringify(evidence));
});

test('grid wrapper uses inward normals on both boundaries and reports the responsible station on failure', () => {
  const xi = [0, .1, .4, .65, 1], eta = [0, .1, .43, .8, 1];
  const nodes = xi.map(u => eta.map(e => ({ x: u + .3 * u * u, y: 2 * e })));
  const result = createOrthogonalBoundaryControl({ nodes, xi, eta, discretization: 'quadratic' });
  for (const side of ['lower', 'upper']) {
    assert.equal(result[side][0], null); assert.equal(result[side].at(-1), null);
    for (let i = 1; i < xi.length - 1; i++) {
      close(result[side][i].normalSpeed, 2); close(result[side][i].stretch, -.6 / (1 + .6 * xi[i]));
      close(result[side][i].normal.y, side === 'lower' ? 1 : -1);
    }
  }
  const invalid = structuredClone(nodes); invalid[2][1].y = -1;
  assert.throws(() => createOrthogonalBoundaryControl({ nodes: invalid, xi, eta }), /lower boundary station 2/);
  assert.throws(() => createOrthogonalBoundaryControl({ nodes, xi, eta: [0, .1, 1] }), /two interior rows/);
  assert.equal(createOrthogonalBoundaryControl({ nodes, xi, eta, sides: ['lower'] }).upper, undefined);
});

test('implicit control sensitivities include the selected normal-speed root and match centered differences', () => {
  const args = { point: { x: .1, y: .2 }, tangent: { x: 1.4, y: .3 }, secondTangent: { x: .6, y: -.5 },
    firstInterior: { x: .095, y: .26 }, secondInterior: { x: .083, y: .35 }, firstDistance: .08, secondDistance: .2 };
  const result = reconstructOrthogonalBoundary(args), h = 1e-6;
  ['firstInterior', 'secondInterior'].forEach((name, row) => {
    for (const key of ['x', 'y']) {
      const plus = structuredClone(args), minus = structuredClone(args); plus[name][key] += h; minus[name][key] -= h;
      const numeric = (reconstructOrthogonalBoundary(plus).stretch - reconstructOrthogonalBoundary(minus).stretch) / (2 * h);
      close(result.interiorDerivatives[row][key], numeric, 2e-6);
    }
  });
  const xi = [0, .15, .4, .7, 1], eta = [0, .12, .55, 1];
  const nodes = xi.map(u => eta.map(e => ({ x: u + .02 * Math.sin(Math.PI * u) * e * (1 - e), y: e })));
  const background = xi.map(u => eta.map(e => .2 * u * (1 + e)));
  const feedback = createOrthogonalBoundaryFeedback({ nodes, xi, eta, background });
  const before = feedback.evaluate(nodes, 2, 1); background[2][1] = 80;
  assert.deepEqual(feedback.evaluate(nodes, 2, 1), before);
  for (let j = 1; j < eta.length - 1; j++) for (const key of ['x', 'y']) {
    const plus = structuredClone(nodes), minus = structuredClone(nodes); plus[2][j][key] += h; minus[2][j][key] -= h;
    close(feedback.evaluate(nodes, 2, j).derivative[key],
      (feedback.evaluate(plus, 2, j).stretch - feedback.evaluate(minus, 2, j).stretch) / (2 * h), 2e-6);
  }
});

test('exponential boundary tails match both prescribed values and recover the linear extension limit', () => {
  const xi = [0, .2, .5, .8, 1], eta = [0, .1, .4, .75, 1];
  const nodes = xi.map(u => eta.map(e => ({ x: u + .03 * Math.sin(Math.PI * u) * e * (1 - e), y: e })));
  const background = xi.map(u => eta.map(e => .2 * u * (1 + e))), i = 2;
  const target = createOrthogonalBoundaryControl({ nodes, xi, eta });
  const decay = { lower: 3, upper: 1.4 }, args = { nodes, xi, eta, background };
  const feedback = createOrthogonalBoundaryFeedback({ ...args, decay });
  const low = target.lower[i].stretch - background[i][0], high = target.upper[i].stretch - background[i].at(-1);
  // Independently match amplitudes of exp(-3*eta), exp(-1.4*(1-eta)).
  const p = (low - Math.exp(-1.4) * high) / (1 - Math.exp(-4.4)), r = high - Math.exp(-3) * p;
  eta.forEach((e, j) => close(feedback.evaluate(nodes, i, j).stretch,
    background[i][j] + p * Math.exp(-3 * e) + r * Math.exp(-1.4 * (1 - e))));
  const before = feedback.evaluate(nodes, i, 1); decay.lower = 90;
  assert.deepEqual(feedback.evaluate(nodes, i, 1), before);
  const linear = createOrthogonalBoundaryFeedback(args), limit = createOrthogonalBoundaryFeedback({ ...args, decay: { lower: 1e-12, upper: 2e-12 } });
  eta.forEach((e, j) => close(limit.evaluate(nodes, i, j).stretch, linear.evaluate(nodes, i, j).stretch, 1e-10));
  assert.throws(() => createOrthogonalBoundaryFeedback({ ...args, decay: { lower: 0, upper: 1 } }), /positive/);
});
