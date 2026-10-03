import test from 'node:test';
import assert from 'node:assert/strict';
import { createCoonsGridGeometry } from '../src/geometry/tests/coons-grid-geometry.js';
import { solveCurvedHarmonicReference } from '../src/geometry/tests/curved-harmonic-reference.js';
import { conformalPolynomialFixture, cylinderHarmonicFixture } from './fixtures/curved-harmonic-fields.js';

const closePoint = (a, b, tolerance = 2e-12) => {
  assert.ok(Math.hypot(a.x - b.x, a.y - b.y) < tolerance, `${JSON.stringify(a)} != ${JSON.stringify(b)}`);
};

test('a single curved Coons cell exactly reproduces a conformal quadratic map and its derivatives', () => {
  const data = conformalPolynomialFixture({ nx: 1, nt: 1 });
  const geometry = createCoonsGridGeometry(data);
  for (const s of [0, .17, .5, .83, 1]) for (const t of [0, .19, .5, .81, 1]) {
    const expected = data.geometry.at(0, 0, s, t), actual = geometry.at(0, 0, s, t);
    closePoint(actual.point, expected.point); closePoint(actual.ds, expected.ds); closePoint(actual.dt, expected.dt);
  }
});

test('curved cell derivatives agree with finite differences and neighboring cells share complete edges', () => {
  const data = conformalPolynomialFixture(), geometry = createCoonsGridGeometry(data), h = 1e-6;
  for (let i = 0; i < data.nx; i++) for (let j = 0; j < data.nt; j++) {
    const s = .27, t = .43, actual = geometry.at(i, j, s, t);
    const a = geometry.at(i, j, s + h, t).point, b = geometry.at(i, j, s - h, t).point;
    const c = geometry.at(i, j, s, t + h).point, d = geometry.at(i, j, s, t - h).point;
    closePoint(actual.ds, { x: (a.x - b.x) / (2 * h), y: (a.y - b.y) / (2 * h) }, 3e-10);
    closePoint(actual.dt, { x: (c.x - d.x) / (2 * h), y: (c.y - d.y) / (2 * h) }, 3e-10);
    for (const q of [0, .31, .77, 1]) {
      if (i + 1 < data.nx) closePoint(geometry.at(i, j, 1, q).point, geometry.at(i + 1, j, 0, q).point);
      if (j + 1 < data.nt) closePoint(geometry.at(i, j, q, 1).point, geometry.at(i, j + 1, q, 0).point);
    }
  }
});

test('the curved reference reproduces exact harmonic eta on a fully analytic conformal geometry', () => {
  const data = conformalPolynomialFixture({ nx: 3, nt: 3 });
  for (const refinement of [1, 2]) {
    // The analytic map is supplied on every cell here; the production Coons
    // geometry curves only boundary cells and has a separate recovery test.
    const result = solveCurvedHarmonicReference(data, { refinement });
    assert.ok(result.maximumTubeIntervals < 1e-11, `Conformal eta discrepancy ${result.maximumTubeIntervals}`);
    result.values.forEach((row, i) => row.forEach((eta, j) => {
      assert.ok(Math.abs(eta - data.psiAt(data.nodes[i][j])) < 1e-12);
    }));
  }
});

test('cylinder boundary callbacks retain exact samples and differentiate arcs and potential-plane edges', () => {
  const data = cylinderHarmonicFixture(), { boundary, nodes, nx, nt } = data;
  for (let i = 0; i <= nx; i++) {
    assert.deepEqual(boundary.bottom(i).point, nodes[i][0]); assert.deepEqual(boundary.top(i).point, nodes[i][nt]);
  }
  for (let j = 0; j <= nt; j++) {
    assert.deepEqual(boundary.left(j).point, nodes[0][j]); assert.deepEqual(boundary.right(j).point, nodes[nx][j]);
  }
  const h = 1e-6;
  for (const [edge, k] of [['bottom', .37], ['bottom', 3.37], ['bottom', 8.37], ['bottom', 10.37],
    ['top', 3.37], ['left', 2.37], ['right', 4.37]]) {
    const actual = boundary[edge](k), a = boundary[edge](k + h).point, b = boundary[edge](k - h).point;
    closePoint(actual.derivative, { x: (a.x - b.x) / (2 * h), y: (a.y - b.y) / (2 * h) }, 3e-9);
    if (edge === 'bottom') assert.ok(Math.abs(data.psiAt(actual.point)) < 1e-14);
  }
  // Stagnation joins are vertices. Angle-based arc derivatives stay finite
  // there without evaluating the singular inverse-potential derivative.
  for (const k of [nx / 6, 5 * nx / 6]) {
    const sample = boundary.bottom(k); assert.ok(Number.isFinite(sample.derivative.x) && Number.isFinite(sample.derivative.y));
  }
  const geometry = createCoonsGridGeometry(data);
  for (const k of [nx / 6, 5 * nx / 6]) {
    closePoint(geometry.at(k - 1, 0, 1, 0).ds, boundary.bottom(k, { interval: k - 1 }).derivative);
    closePoint(geometry.at(k, 0, 0, 0).ds, boundary.bottom(k, { interval: k }).derivative);
    assert.ok(Math.hypot(geometry.at(k - 1, 0, 1, 0).ds.x - geometry.at(k, 0, 0, 0).ds.x,
      geometry.at(k - 1, 0, 1, 0).ds.y - geometry.at(k, 0, 0, 0).ds.y) > .1,
    'The cut/arc corner must retain its distinct one-sided tangents.');
  }
});

test('every curved-boundary callback receives its one-sided parent interval, including constructor knots', () => {
  const data = cylinderHarmonicFixture();
  const boundary = Object.fromEntries(Object.entries(data.boundary).map(([side, callback]) => [side, (k, options) => {
    assert.ok(Number.isInteger(options.interval), `${side} interval missing at ${k}`);
    assert.ok(k >= options.interval && k <= options.interval + 1);
    return callback(k, options);
  }]));
  const geometry = createCoonsGridGeometry({ ...data, boundary });
  for (let i = 0; i < data.nx; i++) for (const s of [0, .4, 1]) geometry.at(i, 0, s, .3);
  assert.throws(() => createCoonsGridGeometry({ ...data, boundary: { ...boundary,
    bottom: (k, options) => {
      const q = boundary.bottom(k, options);
      if (k === 2 && options.interval === 1) q.point.y += .01;
      return q;
    },
  } }), /one-sided curve does not match/);
});

test('curved-boundary cylinder reference refines toward the exact physical streamfunction', t => {
  const data = cylinderHarmonicFixture(), geometry = createCoonsGridGeometry(data), evidence = [];
  let previous = Infinity;
  for (const refinement of [1, 2, 4]) {
    const result = solveCurvedHarmonicReference({ ...data, geometry }, { refinement });
    let error = 0;
    result.values.forEach((row, i) => row.forEach((eta, j) => {
      const exactEta = data.psiAt(data.nodes[i][j]) / data.height;
      error = Math.max(error, Math.abs(eta - exactEta));
    }));
    assert.ok(error < .5 * previous, `Cylinder reference error ${error} after ${previous}`); previous = error;
    evidence.push({ refinement, error, maximumTubeIntervals: result.maximumTubeIntervals });
  }
  assert.ok(previous * data.nt < .05); t.diagnostic(JSON.stringify(evidence));
});

test('curved references reject changed observation nodes, sampled folds, and invalid resource or mass controls', () => {
  const data = conformalPolynomialFixture({ nx: 1, nt: 1 });
  assert.throws(() => createCoonsGridGeometry({ ...data, boundary: {
    bottom: k => { const value = data.boundary.bottom(k); value.point.y += .01; return value; },
  } }), /does not match its prescribed grid nodes/);
  assert.throws(() => solveCurvedHarmonicReference({ ...data, geometry: {
    at: (...args) => { const value = data.geometry.at(...args); value.point.x += .01; return value; },
  } }), /does not match its observation nodes/);
  const square = { nodes: [[{ x: 0, y: 0 }, { x: 0, y: 1 }], [{ x: 1, y: 0 }, { x: 1, y: 1 }]], massFlows: [1] };
  const folded = createCoonsGridGeometry({ ...square, boundary: {
    bottom: s => ({ point: { x: s, y: 5 * s * (1 - s) }, derivative: { x: 1, y: 5 - 10 * s } }),
  } });
  // Positive corner nodes alone cannot make this bowed cell admissible.
  assert.throws(() => solveCurvedHarmonicReference({ ...square, geometry: folded }), /Nonpositive.*Jacobian/);
  assert.throws(() => solveCurvedHarmonicReference({ ...data, massFlows: [0] }), /Invalid/);
  assert.throws(() => solveCurvedHarmonicReference(data, { refinement: 9 }), /Invalid/);
  assert.throws(() => solveCurvedHarmonicReference(data, { refinement: 4, maxUnknowns: 1 }), /exceeding/);
  const boundaryOnly = solveCurvedHarmonicReference(data);
  assert.equal(boundaryOnly.unknowns, 0); assert.equal(boundaryOnly.maximumTubeIntervals, 0);
  assert.deepEqual(boundaryOnly.values, [[0, 1], [0, 1]]);
  assert.equal(boundaryOnly.physicalAcceptance, false);
  assert.equal(boundaryOnly.geometryCheck.globallyCertified, false);
});

test('curved-reference mass coordinates are invariant under rigid motion and uniform geometry or mass scaling', () => {
  const data = conformalPolynomialFixture({ nx: 3, nt: 3 });
  const geometry = createCoonsGridGeometry(data), factor = .002, theta = .71;
  const rotate = p => ({ x: factor * (Math.cos(theta) * p.x - Math.sin(theta) * p.y),
    y: factor * (Math.sin(theta) * p.x + Math.cos(theta) * p.y) });
  const point = p => { const q = rotate(p); return { x: q.x + .017, y: q.y - .025 }; };
  const transformed = { nodes: data.nodes.map(row => row.map(point)), massFlows: data.massFlows.map(m => 7 * m),
    geometry: { at: (...args) => { const q = geometry.at(...args); return { point: point(q.point), ds: rotate(q.ds), dt: rotate(q.dt) }; } } };
  const original = solveCurvedHarmonicReference({ ...data, geometry }, { refinement: 2 });
  const moved = solveCurvedHarmonicReference(transformed, { refinement: 2 });
  original.values.forEach((row, i) => row.forEach((value, j) => assert.ok(Math.abs(value - moved.values[i][j]) < 2e-13)));
  assert.ok(original.maximumTubeIntervals > 1e-6, 'This must exercise a nonzero field correction.');
});
