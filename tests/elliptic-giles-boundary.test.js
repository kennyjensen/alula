// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { solveLinear } from '../src/numerics/linear.js';

const clone = value => structuredClone(value);
const close = (actual, expected, tolerance = 3e-11) => assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} != ${expected}`);
const masses = eta => eta.slice(1).map((v, j) => v - eta[j]);
const xi = [0, 0.09, 0.27, 0.55, 0.79, 1];
const eta = [0, 0.12, 0.39, 0.7, 1];
const makeSystem = (nodes, boundaryConditions = {}, extra = {}) => createEllipticStreamtubeGrid({
  nodes, massFlows: masses(eta), streamwiseCoordinates: xi,
  discretization: 'giles-1985', boundaryConditions, ...extra,
});

function patch(boundaryConditions, perturb = true) {
  return xi.map((u, i) => eta.map((v, j) => {
    const fixedColumn = i === 0 || i === xi.length - 1;
    const copyMode = mode => ['giles-vertical', 'giles-indexed-y'].includes(mode);
    const freeX = !fixedColumn && (j === 0 ? copyMode(boundaryConditions.lower)
      : j === eta.length - 1 ? copyMode(boundaryConditions.upper) : true);
    const interior = !fixedColumn && j > 0 && j < eta.length - 1;
    return { x: 2 * u + (perturb && freeX ? 0.04 * Math.sin(Math.PI * u) * Math.cos(0.7 * Math.PI * v) : 0),
      y: 0.8 + 1.2 * v + (perturb && interior ? 0.016 * Math.sin(Math.PI * u) * Math.sin(Math.PI * v) : 0) };
  }));
}

// Independent transcription of ELLIP's centered secants and edge-slope
// second differences. No production metric/operator/line helpers are used.
function frozenOperator(grid, i, j) {
  const left = xi[i] - xi[i - 1], right = xi[i + 1] - xi[i];
  const bottom = eta[j] - eta[j - 1], top = eta[j + 1] - eta[j];
  const hx = (left + right) / 2, hy = (bottom + top) / 2;
  const sx = (grid[i + 1][j].x - grid[i - 1][j].x) / (2 * hx);
  const sy = (grid[i + 1][j].y - grid[i - 1][j].y) / (2 * hx);
  const tx = (grid[i][j + 1].x - grid[i][j - 1].x) / (2 * hy);
  const ty = (grid[i][j + 1].y - grid[i][j - 1].y) / (2 * hy);
  const alpha = tx * tx + ty * ty, beta = sx * tx + sy * ty, gamma = sx * sx + sy * sy;
  return (state, key) => {
    const middle = state[i][j][key];
    const xx = ((state[i + 1][j][key] - middle) / right - (middle - state[i - 1][j][key]) / left) / hx;
    const yy = ((state[i][j + 1][key] - middle) / top - (middle - state[i][j - 1][key]) / bottom) / hy;
    const xy = (state[i + 1][j + 1][key] - state[i - 1][j + 1][key]
      - state[i + 1][j - 1][key] + state[i - 1][j - 1][key]) / (4 * hx * hy);
    return alpha * xx - 2 * beta * xy + gamma * yy;
  };
}

function denseSweep(initial, omega) {
  const expected = clone(initial), n = xi.length - 2, nt = eta.length - 1;
  for (let i = 1; i <= n; i++) expected[i][0].x = expected[i][1].x;
  for (let j = 1; j < nt; j++) {
    const operators = Array.from({ length: n }, (_, row) => frozenOperator(expected, row + 1, j));
    const matrix = new Float64Array(n * n);
    // Since the frozen operator is linear, unit-basis columns give its
    // exact dense restriction to this line without finite-difference error.
    for (let col = 0; col < n; col++) {
      const basis = expected.map(row => row.map(() => ({ x: 0, y: 0 })));
      basis[col + 1][j].x = 1;
      for (let row = 0; row < n; row++) matrix[row * n + col] = operators[row](basis, 'x');
    }
    for (const key of ['x', 'y']) {
      const correction = solveLinear(matrix, operators.map(operator => -operator(expected, key)));
      for (let i = 1; i <= n; i++) expected[i][j][key] += omega * correction[i - 1];
    }
  }
  for (let i = 1; i <= n; i++) expected[i][nt].x = expected[i][nt - 1].x;
  return expected;
}

test('Giles horizontal boundaries follow the dense sweep oracle: lower lagged, upper fresh, no omega on copies', () => {
  const boundaryConditions = { lower: 'giles-vertical', upper: 'giles-vertical' };
  const nodes = patch(boundaryConditions), before = clone(nodes), system = makeSystem(nodes, boundaryConditions);
  const omega = 1.27, expected = denseSweep(nodes, omega), actual = system.sweep(nodes, omega);
  let largestChange = 0;
  actual.nodes.forEach((row, i) => row.forEach((p, j) => {
    close(p.x, expected[i][j].x); close(p.y, expected[i][j].y);
    largestChange = Math.max(largestChange, Math.abs(p.x - nodes[i][j].x), Math.abs(p.y - nodes[i][j].y));
    if (i === 0 || i === xi.length - 1) assert.deepEqual(p, nodes[i][j]);
    if (j === 0 || j === eta.length - 1) assert.equal(p.y, nodes[i][j].y);
  }));
  const nt = eta.length - 1;
  for (let i = 1; i < xi.length - 1; i++) {
    assert.equal(actual.nodes[i][0].x, nodes[i][1].x);
    assert.equal(actual.nodes[i][nt].x, actual.nodes[i][nt - 1].x);
  }
  assert.ok(Math.abs(actual.nodes[2][0].x - actual.nodes[2][1].x) > 1e-6, 'fixture distinguishes a lagged lower copy from a fresh one');
  assert.ok(Math.abs(actual.nodes[2][nt].x - nodes[2][nt - 1].x) > 1e-6, 'fixture distinguishes a fresh upper copy from a lagged one');
  close(actual.maxUpdate, largestChange / system.lengthScale);
  assert.deepEqual(nodes, before);
  assert.deepEqual(boundaryConditions, { lower: 'giles-vertical', upper: 'giles-vertical' });
  assert.deepEqual(system.coordinateEquations.boundaryConditions, boundaryConditions);
});

test('Giles indexed-y boundaries follow the dense oracle with nonconstant outer y and ordered x copies', () => {
  const boundaryConditions = { lower: 'giles-indexed-y', upper: 'giles-indexed-y' };
  const nodes = patch(boundaryConditions).map((row, i) => row.map((p, j) => ({
    x: p.x, y: p.y + .17 * xi[i] + .035 * Math.sin(2 * Math.PI * xi[i]) * (.4 + .3 * eta[j]),
  })));
  const before = clone(nodes), system = makeSystem(nodes, boundaryConditions), omega = 1.27;
  const expected = denseSweep(nodes, omega), actual = system.sweep(nodes, omega), nt = eta.length - 1;
  for (const j of [0, nt]) assert.ok(nodes.some(row => row[j].y !== nodes[0][j].y));
  let largestChange = 0;
  actual.nodes.forEach((row, i) => row.forEach((p, j) => {
    close(p.x, expected[i][j].x); close(p.y, expected[i][j].y);
    largestChange = Math.max(largestChange, Math.abs(p.x - nodes[i][j].x), Math.abs(p.y - nodes[i][j].y));
    if (!i || i === xi.length - 1) assert.deepEqual(p, nodes[i][j]);
    if (!j || j === nt) assert.equal(p.y, nodes[i][j].y);
  }));
  for (let i = 1; i < xi.length - 1; i++) {
    assert.equal(actual.nodes[i][0].x, nodes[i][1].x);
    assert.equal(actual.nodes[i][nt].x, actual.nodes[i][nt - 1].x);
  }
  assert.ok(Math.abs(actual.nodes[2][0].x - actual.nodes[2][1].x) > 1e-6,
    'Nonconstant-y fixture distinguishes the lagged lower copy.');
  assert.ok(Math.abs(actual.nodes[2][nt].x - nodes[2][nt - 1].x) > 1e-6,
    'Nonconstant-y fixture distinguishes the fresh upper copy.');
  close(actual.maxUpdate, largestChange / system.lengthScale);
  const measured = system.residuals(actual.nodes);
  let expectedResidual = 0;
  for (const { i, j, x } of measured.boundaryRows) {
    const adjacent = j === 0 ? 1 : nt - 1;
    const expected = (actual.nodes[i][j].x - actual.nodes[i][adjacent].x)
      / Math.abs(eta[j] - eta[adjacent]) / system.lengthScale;
    close(x, expected, 1e-14);
    expectedResidual = Math.max(expectedResidual, Math.abs(expected));
  }
  assert.equal(measured.boundaryRows.length, 2 * (xi.length - 2));
  close(measured.boundaryResidual, expectedResidual, 1e-14);
  assert.ok(system.quality(actual.nodes).valid);
  assert.deepEqual(nodes, before);
  assert.deepEqual(system.coordinateEquations.boundaryConditions, boundaryConditions);
});

test('indexed-y affine sloped boundaries measure x_eta, not the physical normal condition', () => {
  const boundaryConditions = { lower: 'giles-indexed-y', upper: 'giles-indexed-y' };
  for (const shear of [0, .3]) {
    const nodes = xi.map(u => eta.map(v => ({ x: 2 * u + shear * v, y: .8 + 1.2 * v + .2 * u })));
    const system = makeSystem(nodes, boundaryConditions), r = system.residuals(nodes);
    assert.ok(r.interiorResidual < 2e-13);
    close(r.boundaryResidual, shear / system.lengthScale, 2e-15);
    for (const row of r.boundaryRows)
      close(row.x, (row.j === 0 ? -shear : shear) / system.lengthScale, 2e-15);
    const result = smoothEllipticStreamtubeGrid(system, { maxSweeps: 0, tolerance: 1e-11 });
    assert.equal(result.converged, shear === 0);
    assert.equal(result.history.length, 1);
    assert.match(result.formulation, /each original indexed boundary y/);
    assert.match(result.formulation, /without a fixed physical curve constraint/);
    if (!shear) {
      // The boundary tangent is (2,.2); its dot product with this vertical
      // crossline is nonzero even though the literal ELLIP x_eta row is zero.
      const p = nodes[2][0], q = nodes[2][1];
      assert.ok(Math.abs(2 * (q.x - p.x) + .2 * (q.y - p.y)) > .01);
      assert.deepEqual(result.nodes, nodes);
    }
  }
});

test('indexed-y mode allows only selected interior boundary x to move and rejects any indexed y change', () => {
  const nodes = xi.map(u => eta.map(v => ({ x: 2 * u, y: .8 + 1.2 * v + .2 * u })));
  for (const side of ['lower', 'upper']) {
    const j = side === 'lower' ? 0 : eta.length - 1, other = j === 0 ? eta.length - 1 : 0;
    const system = makeSystem(nodes, { [side]: 'giles-indexed-y' });
    const allowed = clone(nodes); allowed[2][j].x += .01;
    assert.doesNotThrow(() => system.residuals(allowed));
    assert.doesNotThrow(() => system.sweep(allowed));
    for (const [i, row, key] of [[2, j, 'y'], [2, other, 'x'], [2, other, 'y'],
      [0, j, 'x'], [0, j, 'y'], [xi.length - 1, j, 'x'], [xi.length - 1, j, 'y']]) {
      const changed = clone(nodes); changed[i][row][key] += 1e-12;
      assert.throws(() => system.residuals(changed), /boundary.*fixed/i);
      assert.throws(() => system.sweep(changed), /boundary.*fixed/i);
    }
    assert.throws(() => makeSystem(nodes, { [side]: 'giles-vertical' }), /horizontal/i);
  }
});

test('affine and perturbed rectangles recover with lower, upper, or both horizontal boundaries free', () => {
  for (const boundaryConditions of [{ lower: 'giles-vertical' }, { upper: 'giles-vertical' },
    { lower: 'giles-vertical', upper: 'giles-vertical' }]) {
    const exact = patch(boundaryConditions, false);
    const affine = smoothEllipticStreamtubeGrid(makeSystem(exact, boundaryConditions), { maxSweeps: 0, tolerance: 1e-11 });
    assert.equal(affine.converged, true, affine.reason);
    assert.equal(affine.history.length, 1);
    const nodes = patch(boundaryConditions), before = clone(nodes), system = makeSystem(nodes, boundaryConditions);
    const result = smoothEllipticStreamtubeGrid(system, { maxSweeps: 500, tolerance: 1e-11, omega: 1.2 });
    assert.equal(result.converged, true, `${JSON.stringify(boundaryConditions)}: ${result.reason}`);
    assert.equal(result.quality.valid, true);
    assert.ok(result.history.length > 1);
    for (const entry of result.history) {
      assert.ok(Number.isFinite(entry.interiorResidual) && Number.isFinite(entry.boundaryResidual));
      assert.equal(entry.residual, Math.max(entry.interiorResidual, entry.boundaryResidual));
    }
    result.nodes.forEach((row, i) => row.forEach((p, j) => {
      close(p.x, exact[i][j].x, 3e-10); close(p.y, exact[i][j].y, 3e-10);
      if (i === 0 || i === xi.length - 1) assert.deepEqual(p, nodes[i][j]);
      if (j === 0 || j === eta.length - 1) {
        assert.equal(p.y, nodes[i][j].y);
        const side = j === 0 ? 'lower' : 'upper';
        if (!boundaryConditions[side]) assert.deepEqual(p, nodes[i][j]);
      }
    }));
    assert.deepEqual(nodes, before);
    assert.deepEqual(system.coordinateEquations.boundaryConditions, { lower: 'fixed', upper: 'fixed', ...boundaryConditions });
  }
});

test('an exact sheared affine interior cannot converge at sweep zero with nonzero farfield x_eta', () => {
  const nodes = xi.map(u => eta.map(v => ({ x: 2 * u + 0.3 * v, y: v })));
  const boundaryConditions = { lower: 'giles-vertical', upper: 'giles-vertical' };
  const system = makeSystem(nodes, boundaryConditions), residuals = system.residuals(nodes);
  const scale = Math.hypot(2.3, 1), expected = 0.3 / scale;
  close(system.lengthScale, scale, 1e-15);
  assert.ok(residuals.interiorResidual < 2e-13);
  close(residuals.boundaryResidual, expected, 2e-15);
  assert.equal(residuals.residual, residuals.boundaryResidual);
  assert.equal(residuals.rows.length, (xi.length - 2) * (eta.length - 2));
  assert.equal(residuals.boundaryRows.length, 2 * (xi.length - 2));
  for (const row of residuals.boundaryRows) {
    assert.ok(row.i > 0 && row.i < xi.length - 1);
    assert.ok(row.j === 0 || row.j === eta.length - 1);
    close(row.x, row.j === 0 ? -expected : expected, 2e-15);
  }
  const stopped = smoothEllipticStreamtubeGrid(system, { maxSweeps: 0, tolerance: 1e-11 });
  assert.equal(stopped.converged, false);
  assert.equal(stopped.history.length, 1);
  assert.equal(stopped.history[0].boundaryResidual, residuals.boundaryResidual);
  const scaled = nodes.map(row => row.map(p => ({ x: 7 * p.x, y: 7 * p.y })));
  const scaledResidual = makeSystem(scaled, boundaryConditions, { massFlows: masses(eta).map(m => m * 13) }).residuals(scaled);
  close(scaledResidual.boundaryResidual, expected, 2e-15);
});

test('only selected interior farfield x coordinates may change; endpoint, y, and other boundary data remain fixed', () => {
  const nodes = patch({ lower: 'giles-vertical' }), system = makeSystem(nodes, { lower: 'giles-vertical' });
  const allowed = clone(nodes); allowed[2][0].x += 0.01;
  assert.doesNotThrow(() => system.residuals(allowed));
  for (const [i, j, key] of [[2, 0, 'y'], [2, eta.length - 1, 'x'], [0, 1, 'x'], [xi.length - 1, 0, 'x']]) {
    const invalid = clone(nodes); invalid[i][j][key] += 0.001;
    assert.throws(() => system.residuals(invalid), /boundary|fixed/i);
  }
  const fixed = makeSystem(nodes);
  assert.deepEqual(fixed.coordinateEquations.boundaryConditions, { lower: 'fixed', upper: 'fixed' });
  assert.throws(() => fixed.residuals(allowed), /boundary|fixed/i);
});

test('Giles farfield options reject curved or sloped selected rows and unknown conditions', () => {
  const flat = patch({}, false);
  for (const side of ['lower', 'upper']) {
    const j = side === 'lower' ? 0 : eta.length - 1;
    const curved = clone(flat); curved[2][j].y += 1e-15;
    assert.throws(() => makeSystem(curved, { [side]: 'giles-vertical' }), /horizontal/i);
    const sloped = clone(flat); sloped.forEach(row => { row[j].y += 0.1 * row[j].x; });
    assert.throws(() => makeSystem(sloped, { [side]: 'giles-vertical' }), /horizontal/i);
    assert.doesNotThrow(() => makeSystem(sloped, { [side]: 'fixed' }));
  }
  for (const boundaryConditions of [null, [], 'giles-vertical', { lower: 'unknown' }, { upper: true }, { left: 'fixed' }])
    assert.throws(() => makeSystem(flat, boundaryConditions), /boundary|condition|option/i);
});
