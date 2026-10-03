import test from 'node:test';
import assert from 'node:assert/strict';
import { quadLaplaceMatrix } from '../src/numerics/quad-laplace.js';
import { solveAlternatingScalarLines } from '../src/numerics/alternating-scalar-lines.js';
import { assembleTangentialCoordinate, redistributeStreamtubeTangentially } from '../src/geometry/streamtube-tangential-redistribution.js';
import { sparseDense, sparseProduct, sparseMatrix, sparseAdd } from '../src/numerics/sparse.js';
import { solveLinear } from '../src/numerics/linear.js';

const near = (a, b, tolerance = 2e-12) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);
const grid = (nx = 8, nt = 6, amplitude = .08) => Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => {
  const u = i / nx, v = j / nt, interior = i > 0 && i < nx && j > 0 && j < nt;
  return { x: u + (interior ? amplitude * Math.sin(Math.PI * u) * Math.sin(Math.PI * v) : 0),
    y: v + (interior ? .3 * amplitude * Math.sin(2 * Math.PI * u) * Math.sin(Math.PI * v) : 0) };
}));

test('2x2 STIFF reconstruction recovers the closed-form square matrix and exact linear-field energy', () => {
  const square = [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: 1, y: 1 }, { x: 0, y: 1 }];
  const expected = [4, -1, -2, -1, -1, 4, -1, -2, -2, -1, 4, -1, -1, -2, -1, 4].map(v => v / 6);
  const a = quadLaplaceMatrix(square, { quadratureOrder: 2 }); a.forEach((v, i) => near(v, expected[i]));
  for (const vertices of [square, [{ x: -.2, y: .1 }, { x: 1.1, y: -.3 }, { x: 1.4, y: 1.2 }, { x: .1, y: 1.4 }]]) {
    const matrix = quadLaplaceMatrix(vertices, { quadratureOrder: 2 }), u = vertices.map(p => p.x + 2 * p.y);
    let energy = 0, area = 0;
    for (let r = 0; r < 4; r++) {
      near(matrix.subarray(4 * r, 4 * r + 4).reduce((a, b) => a + b, 0), 0);
      const next = vertices[(r + 1) % 4]; area += .5 * (vertices[r].x * next.y - vertices[r].y * next.x);
      for (let c = 0; c < 4; c++) { near(matrix[4 * r + c], matrix[4 * c + r]); energy += u[r] * matrix[4 * r + c] * u[c]; }
    }
    near(energy, 5 * area);
  }
  assert.deepEqual(quadLaplaceMatrix(square), quadLaplaceMatrix(square, { quadratureOrder: 3 }));
  assert.throws(() => quadLaplaceMatrix(square, { quadratureOrder: 1 }), /quadrature/);
});

test('assembled coordinate equation reproduces an exact physical harmonic field on a distorted grid', () => {
  const nodes = grid(), a = assembleTangentialCoordinate(nodes), dense = sparseDense(a.matrix);
  const expected = [];
  for (let i = 1; i < a.nx; i++) for (let j = 1; j < a.nt; j++) expected.push(nodes[i][j].x * a.nx * (1 / a.nx) ** .25 - a.baseline[i]);
  const rhs = sparseProduct(a.matrix, expected);
  rhs.forEach((v, i) => near(v, a.rhs[i], 5e-14));
  solveLinear(dense, a.rhs).forEach((v, i) => near(v, expected[i], 5e-14));
  for (let r = 0; r < a.matrix.n; r++) for (let c = 0; c < a.matrix.n; c++) near(dense[r * a.matrix.n + c], dense[c * a.matrix.n + r]);
});

test('fixed alternating tridiagonal pairs agree with independent dense line solves and converge to a full dense solve', () => {
  const a = assembleTangentialCoordinate(grid()), dense = sparseDense(a.matrix), n = a.matrix.n;
  const id = (i, j) => (i - 1) * (a.nt - 1) + j - 1;
  const lineIds = [...Array.from({ length: a.nx - 1 }, (_, i) => Array.from({ length: a.nt - 1 }, (_, j) => id(i + 1, j + 1))),
    ...Array.from({ length: a.nt - 1 }, (_, j) => Array.from({ length: a.nx - 1 }, (_, i) => id(i + 1, j + 1)))];
  for (const pairs of [1, 5]) {
    const expected = new Float64Array(n);
    for (let pair = 0; pair < pairs; pair++) for (const ids of lineIds) {
      const b = ids.map(r => a.rhs[r] - expected.reduce((s, v, c) => s + dense[r * n + c] * v, 0));
      const block = ids.flatMap(r => ids.map(c => dense[r * n + c]));
      const correction = solveLinear(block, b);
      ids.forEach((r, k) => { expected[r] += correction[k]; });
    }
    const actual = solveAlternatingScalarLines(a.matrix, a.rhs, { ...a, pairs });
    actual.x.forEach((v, i) => near(v, expected[i], 2e-14));
    assert.equal(actual.history.length, pairs + 1); assert.equal(actual.pairs, pairs);
  }
  const exact = solveLinear(dense, a.rhs), converged = solveAlternatingScalarLines(a.matrix, a.rhs, { ...a, pairs: 100 });
  assert.equal(converged.converged, true); converged.x.forEach((v, i) => near(v, exact[i]));
  const none = solveAlternatingScalarLines(a.matrix, a.rhs, { ...a, pairs: 0 });
  assert.equal(none.converged, false); assert.equal(none.relativeResidual, 1);
});

test('SMOVE formula preserves boundaries and the uniform grid, and every displacement follows the original centered tangent', () => {
  const uniform = grid(8, 6, 0), unchanged = redistributeStreamtubeTangentially(uniform);
  assert.ok(unchanged.maxDisplacement < 2e-15);
  const nodes = grid(), before = structuredClone(nodes), result = redistributeStreamtubeTangentially(nodes);
  assert.deepEqual(nodes, before); assert.equal(result.solution.pairs, 5); assert.ok(result.maxDisplacement > .01);
  for (let i = 0; i < nodes.length; i++) for (let j = 0; j < nodes[0].length; j++) {
    const p = nodes[i][j], q = result.nodes[i][j];
    if (i === 0 || j === 0 || i === nodes.length - 1 || j === nodes[0].length - 1) assert.deepEqual(q, p);
    else {
      const a = nodes[i - 1][j], b = nodes[i + 1][j];
      near((q.x - p.x) * (b.y - a.y) - (q.y - p.y) * (b.x - a.x), 0, 5e-17);
    }
  }
});

test('redistribution and its fourth-root coordinate increments respect rotation, translation, scale and bank choice', () => {
  const nodes = grid(7, 5), angle = .7, c = Math.cos(angle), s = Math.sin(angle), scale = 3;
  const transform = p => ({ x: scale * (c * p.x - s * p.y) + 4, y: scale * (s * p.x + c * p.y) - 2 });
  const movedInput = nodes.map(row => row.map(transform));
  for (const referenceBank of [0, 5]) {
    const a = redistributeStreamtubeTangentially(nodes, { referenceBank }), b = redistributeStreamtubeTangentially(movedInput, { referenceBank });
    a.nodes.forEach((row, i) => row.forEach((p, j) => { const q = transform(p); near(b.nodes[i][j].x, q.x); near(b.nodes[i][j].y, q.y); }));
    const ca = assembleTangentialCoordinate(nodes, { referenceBank }), cb = assembleTangentialCoordinate(movedInput, { referenceBank });
    ca.increments.forEach((v, i) => near(cb.increments[i], scale ** .25 * v));
  }
});

test('invalid grids and nonlocal or unresolved line systems fail without modifying input', () => {
  const bad = grid(3, 3), before = structuredClone(bad); bad[1][1] = { x: 4, y: 4 };
  const invalid = structuredClone(bad);
  assert.throws(() => redistributeStreamtubeTangentially(bad), /Invalid quadrilateral/); assert.deepEqual(bad, invalid);
  assert.throws(() => assembleTangentialCoordinate(before, { referenceBank: 1 }), /reference bank/);
  const a = assembleTangentialCoordinate(grid(3, 3));
  assert.throws(() => solveAlternatingScalarLines(a.matrix, a.rhs, { ...a, pairs: -1 }), /controls/);
  const singular = structuredClone(a.matrix); singular.values.fill(0);
  assert.throws(() => solveAlternatingScalarLines(singular, a.rhs, a), /pivot/);
  const remote = sparseMatrix([[8], [], [], [], [], [], [], [], []]); sparseAdd(remote, 0, 8, 1);
  assert.throws(() => solveAlternatingScalarLines(remote, new Float64Array(9), { nx: 4, nt: 4 }), /nine-point/);
});

test('free farfield bank rows reproduce a harmonic field with zero natural normal flux', () => {
  const nx = 8, nt = 6;
  for (const fixedBanks of [[false, true], [true, false], [false, false]]) {
    const nodes = grid(nx, nt, 0);
    for (let i = 1; i < nx; i++) for (let j = 0; j <= nt; j++) {
      if (j === 0 && fixedBanks[0] || j === nt && fixedBanks[1]) continue;
      nodes[i][j].x += .04 * Math.sin(Math.PI * i / nx) * (1 + .2 * j / nt);
    }
    // xi is proportional to physical x. The free horizontal banks have
    // zero normal flux; the fixed bank and inlet/outlet values are exact.
    // Supply an unperturbed reference bank even when its coordinates move:
    // for the fully free test use no perturbation on bank zero.
    const referenceBank = fixedBanks[0] ? 0 : fixedBanks[1] ? nt : 0;
    if (!fixedBanks.some(Boolean)) for (let i = 0; i <= nx; i++) nodes[i][0].x = i / nx;
    const a = assembleTangentialCoordinate(nodes, { fixedBanks, referenceBank }), expected = [];
    for (let i = 1; i < nx; i++) for (let j = a.firstStreamline; j <= a.lastStreamline; j++)
      expected.push(nodes[i][j].x * nx * (1 / nx) ** .25 - a.baseline[i]);
    sparseProduct(a.matrix, expected).forEach((v, i) => near(v, a.rhs[i], 1e-13));
    const solved = solveAlternatingScalarLines(a.matrix, a.rhs, { ...a.lineDimensions, pairs: 100 });
    assert.equal(solved.converged, true); solved.x.forEach((v, i) => near(v, expected[i]));
    const r = redistributeStreamtubeTangentially(nodes, { fixedBanks, referenceBank });
    let freeMotion = 0;
    for (let i = 0; i <= nx; i++) for (const j of [0, nt]) {
      if (i === 0 || i === nx || fixedBanks[j === 0 ? 0 : 1]) assert.deepEqual(r.nodes[i][j], nodes[i][j]);
      else { near(r.nodes[i][j].y, nodes[i][j].y); freeMotion = Math.max(freeMotion, Math.abs(r.nodes[i][j].x - nodes[i][j].x)); }
    }
    assert.ok(freeMotion > .01);
  }
});
