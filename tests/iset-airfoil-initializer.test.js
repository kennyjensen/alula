import test from 'node:test';
import assert from 'node:assert/strict';
import { createIsetAirfoilInitialGrid, locateIsetLeadingEdge } from '../src/geometry/tests/iset-airfoil-initializer.js';
import { naca4 } from '../src/geometry/airfoil.js';
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { solveLinear } from '../src/numerics/linear.js';

const close = (a, b, tolerance = 3e-12) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);
const fixture = () => ({ points: naca4('0012', 40), inletStations: [0, .3, .7, .9, .97, 1],
  outletStations: [0, .03, .1, .3, .7, 1],
  surfaceStations: { upper: Array.from({ length: 9 }, (_, i) => .5 * (1 - Math.cos(Math.PI * i / 8))),
    lower: Array.from({ length: 9 }, (_, i) => .5 * (1 - Math.cos(Math.PI * i / 8))) },
  transverseStations: [0, .03, .15, .5, .5, .85, .97, 1] });

test('ISET supplied physical stations share LE/TE endpoints and use the same geometric cuts on both surfaces', () => {
  const input = fixture(), before = structuredClone(input), r = createIsetAirfoilInitialGrid(input);
  assert.deepEqual(input, before);
  assert.equal(r.leadingIndex, 5); assert.equal(r.trailingIndex, 13); assert.equal(r.nodes[0].length, 19);
  for (let i = 0; i < r.nodes[0].length; i++) {
    if (i <= r.leadingIndex || i >= r.trailingIndex) assert.deepEqual(r.outline.upper[i], r.outline.lower[i]);
    if (i < r.leadingIndex) close(r.outline.upper[i].x, r.leading.point.x + 4 * (input.inletStations[i] - 1));
    if (i > r.trailingIndex) close(r.outline.upper[i].x, 1 + 4 * input.outletStations[i - r.trailingIndex]);
  }
  assert.deepEqual(r.outline.upper[r.leadingIndex], r.leading.point);
  assert.equal(r.diagnostics.originalFortranExecuted, false); assert.equal(r.diagnostics.flowSolved, false);
  assert.deepEqual(r.diagnostics.stages, ['NORMIT', 'SPLNIT', 'OUTLIN', 'NORLIN', 'RESPLI']);
});

test('NORLIN outer x is uniform in index and independent of supplied inlet/body station spacing', () => {
  const input = fixture(), a = createIsetAirfoilInitialGrid(input);
  input.inletStations = [0, .2, .4, .6, .8, 1];
  input.surfaceStations.upper = input.surfaceStations.upper.map(f => f ** 1.2);
  const b = createIsetAirfoilInitialGrid(input);
  assert.deepEqual(a.farfield, b.farfield);
  assert.notDeepEqual(a.outline, b.outline); assert.notDeepEqual(a.xi, b.xi);
  a.farfield.upper.forEach((p, i) => {
    close(p.x, -4 + i / 2); close(p.y, 2);
    close(a.farfield.lower[i].y, -2);
  });
});

test('NORLIN mass-coordinate interpolation, nonconstant farfield y and region orientation follow the topology break', () => {
  const r = createIsetAirfoilInitialGrid({ ...fixture(), inletSlope: .12, outletSlope: -.03,
    transverseStations: [10, 11, 14, 20, 20, 27, 30] });
  assert.deepEqual(r.distributions.transverse, [0, .05, .2, .5, .5, .85, 1]);
  assert.deepEqual(r.nodes.map(group => group[0].length), [3, 4]);
  r.nodes.forEach((group, g) => {
    const weights = r.massFlows[g], total = weights.reduce((s, w) => s + w, 0);
    group.forEach((row, i) => {
      let eta = 0;
      row.forEach((p, j) => {
        const a = g ? r.outline.upper[i] : r.farfield.lower[i], b = g ? r.farfield.upper[i] : r.outline.lower[i];
        close(p.x, (1 - eta) * a.x + eta * b.x); close(p.y, (1 - eta) * a.y + eta * b.y);
        eta += (weights[j] ?? 0) / total;
      });
      close(r.farfield.upper[i].y - r.farfield.lower[i].y, 4);
    });
  });
  assert.notEqual(r.farfield.upper[0].y, r.farfield.upper.at(-1).y);
  // No forced vertical inlet: outer endpoints use normalized chord limits,
  // whereas the cut follows the actual slope-dependent nose point.
  assert.ok(Math.abs(r.outline.upper[0].x - r.farfield.upper[0].x) > 1e-4);
});

test('RESPLI resamples the actual wall nodes before ELLIP and preserves two natural, nonperiodic TE endpoints', () => {
  const input = fixture(); input.surfaceStations.lower = input.surfaceStations.lower.map(f => f ** 1.3);
  const r = createIsetAirfoilInitialGrid(input), { points, curve, leadingParameter, surfaceFractions } = r.resampled;
  assert.equal(points.length, 17); assert.equal(leadingParameter, curve.knots[8]);
  for (const side of ['upper', 'lower']) for (let k = 0; k <= 8; k++) {
    const p = curve.branch(side, surfaceFractions[side][k], leadingParameter).point;
    const expected = r.outline[side][r.leadingIndex + k];
    close(p.x, expected.x); close(p.y, expected.y);
  }
  for (const s of [0, curve.length]) {
    close(curve.evaluate(s).secondDerivative.x, 0); close(curve.evaluate(s).secondDerivative.y, 0);
  }
  assert.notEqual(curve.evaluate(0).derivative.y, curve.evaluate(curve.length).derivative.y);
  assert.notEqual(r.curve.length, curve.length, 'RESPLI must recompute polygon lengths on the new samples');
});

test('the natural contour agrees with the thesis first-derivative spline equations via an independent dense solve', () => {
  const curve = createIsetAirfoilInitialGrid(fixture()).resampled.curve, n = curve.knots.length;
  for (const key of ['x', 'y']) {
    const s = curve.knots, f = s.map(v => curve.evaluate(v).point[key]);
    const a = new Float64Array(n * n), rhs = new Float64Array(n);
    a[0] = 2; a[1] = 1; rhs[0] = 3 * (f[1] - f[0]) / (s[1] - s[0]);
    a[(n - 1) * n + n - 2] = 1; a[n * n - 1] = 2;
    rhs[n - 1] = 3 * (f[n - 1] - f[n - 2]) / (s[n - 1] - s[n - 2]);
    for (let i = 1; i < n - 1; i++) {
      const left = s[i] - s[i - 1], right = s[i + 1] - s[i];
      a[i * n + i - 1] = right; a[i * n + i] = 2 * (left + right); a[i * n + i + 1] = left;
      rhs[i] = 3 * (left * (f[i + 1] - f[i]) / right + right * (f[i] - f[i - 1]) / left);
    }
    const derivative = solveLinear(a, rhs);
    derivative.forEach((v, i) => close(v, curve.evaluate(s[i]).derivative[key]));
  }
});

test('normalization preserves translation/scale and reports the source surface-prefix truncation', () => {
  const input = fixture(), original = createIsetAirfoilInitialGrid(input);
  input.points = input.points.map(p => ({ x: 3 * p.x + 7, y: 3 * p.y - 4 }));
  for (const name of ['inletStations', 'outletStations', 'transverseStations']) input[name] = input[name].map(v => 5 * v + 2);
  input.surfaceStations.upper.push(1.5, 2);
  const r = createIsetAirfoilInitialGrid(input);
  assert.deepEqual(r.diagnostics.surfacePrefixTruncation, { upper: 2, lower: 0 });
  assert.deepEqual(r.diagnostics.geometryNormalization, { origin: { x: 7, y: -4 }, chord: 3 });
  r.nodes.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => {
    close(p.x, original.nodes[g][i][j].x); close(p.y, original.nodes[g][i][j].y);
  })));
});

test('SPLNIT retains its sampled tangent interpolation and rejects unresolved source inputs', () => {
  // Projection crosses zero at a known knot; no nonlinear/panel refinement.
  const curve = { knots: [0, 1, 2], evaluate: s => ({ point: { x: (s - 1) ** 2, y: s }, derivative: { x: 2 * (s - 1), y: 1 } }) };
  close(locateIsetLeadingEdge(curve, .5).parameter, .75);
  assert.throws(() => locateIsetLeadingEdge(curve, 4), /did not find/);
  assert.throws(() => createIsetAirfoilInitialGrid({ ...fixture(), transverseStations: [0, .2, .5, .8, 1] }), /exactly one/);
  assert.throws(() => createIsetAirfoilInitialGrid({ ...fixture(), transverseStations: [0, .2, .2, .5, .5, .8, 1] }), /exactly one/);
  assert.throws(() => createIsetAirfoilInitialGrid({ ...fixture(), inletStations: [0, .5, .4, 1] }), /ordered/);
  assert.throws(() => createIsetAirfoilInitialGrid({ ...fixture(), leadingEdge: 'sharp' }), /rounded/);
  const open = fixture(); open.points.at(-1).y += 1e-13;
  assert.throws(() => createIsetAirfoilInitialGrid(open), /closed trailing edge/);
  assert.throws(() => createIsetAirfoilInitialGrid({ ...fixture(), points: fixture().points.toReversed() }), /orientation/);
  assert.throws(() => createIsetAirfoilInitialGrid({ ...fixture(),
    transverseStations: [-1e16, -5e15, -.125, .125, 5e15, 1e16] }), /collapse/);
});

test('the small isolated ISET seed enters the source indexed-y ELLIP equations without a panel or Euler solve', () => {
  const r = createIsetAirfoilInitialGrid(fixture()), before = structuredClone(r.nodes);
  const results = r.nodes.map((nodes, g) => smoothEllipticStreamtubeGrid(createEllipticStreamtubeGrid({ nodes,
    massFlows: r.massFlows[g], streamwiseCoordinates: r.xi, discretization: 'giles-1985',
    boundaryConditions: { [g ? 'upper' : 'lower']: 'giles-indexed-y' } }), { maxSweeps: 160 }));
  assert.deepEqual(r.nodes, before);
  for (const [g, q] of results.entries()) {
    assert.equal(q.converged, true, q.reason); assert.equal(q.quality.valid, true);
    q.nodes.forEach((row, i) => row.forEach((p, j) => {
      if (!i || i === q.nodes.length - 1 || !j || j === row.length - 1) {
        if (!i || i === q.nodes.length - 1 || j === (g ? 0 : row.length - 1)) assert.deepEqual(p, before[g][i][j]);
        if (!j || j === row.length - 1) assert.equal(p.y, before[g][i][j].y);
      }
    }));
  }
});
