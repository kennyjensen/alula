import test from 'node:test';
import assert from 'node:assert/strict';
import { createPolynomialGridGeometry } from '../src/geometry/polynomial-grid-geometry.js';
import { conformalPolynomialGridFixture } from './fixtures/polynomial-grid.js';

const close = (a, b, tolerance = 2e-12) => assert.ok(Math.hypot(a.x - b.x, a.y - b.y) < tolerance,
  `${JSON.stringify(a)} != ${JSON.stringify(b)}`);
const sharedEdges = (geometry, nx, nt) => {
  for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) for (const q of [0, .17, .51, .89, 1]) {
    if (i + 1 < nx) close(geometry.at(i, j, 1, q).point, geometry.at(i + 1, j, 0, q).point);
    if (j + 1 < nt) close(geometry.at(i, j, q, 1).point, geometry.at(i, j + 1, q, 0).point);
  }
};

test('canonical tensor-Q2 coefficients reproduce the conformal map and share complete exact faces', () => {
  const data = conformalPolynomialGridFixture(), model = createPolynomialGridGeometry(data);
  assert.deepEqual(model.degree, { s: 2, t: 2 });
  for (let i = 0; i < data.nx; i++) for (let j = 0; j < data.nt; j++) {
    if (i) assert.deepEqual(data.controlPoints[i][j][0], data.controlPoints[i - 1][j][2]);
    if (j) data.controlPoints[i][j].forEach((row, a) => assert.deepEqual(row[0], data.controlPoints[i][j - 1][a][2]));
    for (const s of [0, .27, .6, 1]) for (const t of [0, .31, .73, 1]) {
      const exact = data.analyticAt(i, j, s, t), actual = model.at(i, j, s, t);
      close(actual.point, exact.point); close(actual.ds, exact.ds); close(actual.dt, exact.dt);
    }
  }
  sharedEdges(model, data.nx, data.nt);
  data.nodes.forEach((row, i) => row.forEach((p, j) => assert.ok(Math.abs(data.psiAt(p) - j / data.nt) < 1e-14)));
});

test('polynomial correction vanishes at observation nodes and has consistent analytic derivatives', () => {
  const data = conformalPolynomialGridFixture(), model = createPolynomialGridGeometry(data), h = 1e-6;
  for (let i = 0; i < data.nx; i++) for (let j = 0; j < data.nt; j++) {
    for (const s of [0, 1]) for (const t of [0, 1]) assert.deepEqual(model.correction(i, j, s, t).point, { x: 0, y: 0 });
    const s = .27, t = .43, q = model.correction(i, j, s, t);
    const a = model.correction(i, j, s + h, t).point, b = model.correction(i, j, s - h, t).point;
    const c = model.correction(i, j, s, t + h).point, d = model.correction(i, j, s, t - h).point;
    close(q.ds, { x: (a.x - b.x) / (2 * h), y: (a.y - b.y) / (2 * h) }, 1e-10);
    close(q.dt, { x: (c.x - d.x) / (2 * h), y: (c.y - d.y) / (2 * h) }, 1e-10);
  }
});

test('Q1 nodal perturbations preserve shared polynomial faces and onGrid snapshots its input', () => {
  const data = conformalPolynomialGridFixture(), model = createPolynomialGridGeometry(data);
  const moved = data.nodes.map((row, i) => row.map((p, j) => ({ ...p,
    y: p.y + (!i || !j || i === data.nx || j === data.nt ? 0 : .008 * Math.sin(i + j)),
  })));
  const geometry = model.onGrid(moved); sharedEdges(geometry, data.nx, data.nt);
  for (let i = 0; i < data.nx; i++) for (let j = 0; j < data.nt; j++) {
    const s = .31, t = .47, weights = [(1 - s) * (1 - t), s * (1 - t), s * t, (1 - s) * t];
    const ij = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]], original = model.at(i, j, s, t).point;
    const dy = ij.reduce((sum, [a, b], k) => sum + weights[k] * (moved[a][b].y - data.nodes[a][b].y), 0);
    close(geometry.at(i, j, s, t).point, { x: original.x, y: original.y + dy });
  }
  const before = geometry.at(0, 0, .4, .7); moved[1][1].y += 10;
  assert.deepEqual(geometry.at(0, 0, .4, .7), before);
});

test('nonconforming polynomial control corners and complete faces are rejected', () => {
  const data = conformalPolynomialGridFixture();
  const corners = structuredClone(data.controlPoints); corners[0][0][0][0].x += .001;
  assert.throws(() => createPolynomialGridGeometry({ ...data, controlPoints: corners }), /corners/);
  const streamwise = structuredClone(data.controlPoints); streamwise[1][0][0][1].y += .001;
  assert.throws(() => createPolynomialGridGeometry({ ...data, controlPoints: streamwise }), /streamwise face/);
  const transverse = structuredClone(data.controlPoints); transverse[0][1][1][0].x += .001;
  assert.throws(() => createPolynomialGridGeometry({ ...data, controlPoints: transverse }), /transverse face/);
  const invalid = structuredClone(data.controlPoints); invalid[0][0][1][1].x = NaN;
  assert.throws(() => createPolynomialGridGeometry({ ...data, controlPoints: invalid }), /degrees/);
});

test('whole-cell certificates accept conformal geometry and reject reversed guides or a folded map', () => {
  const data = conformalPolynomialGridFixture({ nx: 3, nt: 3 }), model = createPolynomialGridGeometry(data);
  const good = model.quality(data.nodes, data.directions);
  assert.equal(good.valid, true, JSON.stringify(good)); assert.equal(good.positiveCells, true); assert.equal(good.transverse, true);
  assert.ok(good.minimumJacobian > 0); assert.ok(good.minimumTransversality > 0);
  assert.equal(good.cellCertificates.length, data.nx * data.nt);
  const reversed = data.directions.map(row => row.map(() => ({ x: 0, y: -1 })));
  const badGuide = model.quality(data.nodes, reversed);
  assert.equal(badGuide.valid, false); assert.equal(badGuide.transverse, false); assert.ok(badGuide.invalidCells.length > 0);
  const folded = data.nodes.map(row => row.map(p => ({ x: -p.x, y: p.y })));
  const badMap = model.quality(folded, data.directions);
  assert.equal(badMap.valid, false); assert.equal(badMap.positiveCells, false); assert.ok(badMap.invalidCells.length > 0);
});
