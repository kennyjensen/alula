import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';

const fixture = () => Array.from({ length: 3 }, (_, i) => Array.from({ length: 3 }, (_, j) => ({
  x: i === 0 ? 0 : i === 2 ? 1 : j === 1 ? .9 : .1, y: j / 2,
})));

test('a convex SLOR step bounds a real overshoot without changing its harmonic equations or fixed boundaries', () => {
  const nodes = fixture(), before = structuredClone(nodes);
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: [1, 1], discretization: 'giles-1985' });
  assert.ok(system.quality(nodes).valid);
  const full = system.sweep(nodes, 1.9);
  assert.equal(system.quality(full.nodes).valid, false);
  const snapshots = [];
  const result = smoothEllipticStreamtubeGrid(system, { omega: 1.9, maxSweeps: 1, requireConvex: true,
    onSweep: (history, state) => snapshots.push(structuredClone({ history, state })) });
  assert.equal(result.history[0].fraction, .5);
  assert.ok(result.quality.valid);
  assert.ok(snapshots.every(s => system.quality(s.state).valid));
  assert.equal(result.converged, false);
  assert.equal(result.reason, 'sweep limit');
  assert.equal(result.coordinateEquations.xi, 'Laplace');
  assert.equal(result.coordinateEquations.eta, 'Laplace');
  result.nodes.forEach((row, i) => row.forEach((p, j) => {
    if (i !== 1 || j !== 1) assert.deepEqual(p, before[i][j]);
  }));
  assert.deepEqual(nodes, before); assert.deepEqual(system.initial, before);
});

test('bounded SLOR recovers the independently known harmonic root and reports its actual residual', () => {
  const system = createEllipticStreamtubeGrid({ nodes: fixture(), massFlows: [1, 1], discretization: 'giles-1985' });
  const result = smoothEllipticStreamtubeGrid(system, { omega: 1.9, maxSweeps: 500, tolerance: 1e-11, requireConvex: true });
  assert.ok(result.converged, result.reason);
  // With the centered secants, both metrics are one and the single x row
  // is proportional to (0 + 1 - 2*x) + (.1 + .1 - 2*x), giving x=.3.
  assert.ok(Math.abs(result.nodes[1][1].x - .3) < 1e-11);
  assert.equal(result.nodes[1][1].y, .5);
  assert.ok(system.residuals(result.nodes).residual <= 1e-11);
  assert.ok(result.history.every(h => h.invalidCells === 0));
});

test('an already folded initial grid cannot be published as a convex SLOR iterate', () => {
  const nodes = fixture(); nodes[1][1].x = -1;
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: [1, 1] });
  let publications = 0;
  const result = smoothEllipticStreamtubeGrid(system, { requireConvex: true, onSweep: () => publications++ });
  assert.equal(result.converged, false); assert.equal(publications, 0);
  assert.match(result.reason, /Invalid starting state/);
  assert.deepEqual(result.nodes, nodes);
  assert.throws(() => smoothEllipticStreamtubeGrid(system, { requireConvex: 1 }), /controls/);
});

test('the exact default seven-tube slat passage reaches a healthy Giles harmonic grid with unchanged physical boundaries', () => {
  const { nodes, massFlows, xi } = JSON.parse(fs.readFileSync(new URL('./fixtures/three-slat-harmonic-grid.json', import.meta.url)));
  const before = structuredClone(nodes), system = createEllipticStreamtubeGrid({ nodes, massFlows,
    streamwiseCoordinates: xi, discretization: 'giles-1985' });
  assert.equal(system.nt, 7); assert.equal(system.nx, 229);
  const result = smoothEllipticStreamtubeGrid(system, { omega: 1, maxSweeps: 100, tolerance: 1e-9, requireConvex: true });
  assert.ok(result.converged, result.reason); assert.ok(result.quality.minCornerSine > .2);
  assert.ok(result.history.every(h => h.invalidCells === 0));
  assert.ok(system.residuals(result.nodes).residual <= 1e-9);
  assert.deepEqual(system.xi, xi);
  result.nodes.forEach((row, i) => row.forEach((p, j) => {
    if (!i || i === system.nx || !j || j === system.nt) assert.deepEqual(p, before[i][j]);
  }));
  assert.deepEqual(nodes, before);
});

// An affine map solves the continuous inverse-Laplace equations exactly.
// Thin mass-coordinate bands amplify rounding of its stored x coordinates.
const thinAffineGrid = () => {
  const xi = [0, .1, .3, .6, .8, 1], eta = [0, .5, .50001, .50002, 1];
  return { nodes: xi.map(u => eta.map(v => ({ x: 4 + 4 * u + .1234567 * v, y: .1 * u + v }))),
    massFlows: eta.slice(1).map((v, j) => v - eta[j]), streamwiseCoordinates: xi, discretization: 'giles-1985' };
};

test('coordinate precision is reported separately from a tolerance root in a known affine grid', () => {
  const input = thinAffineGrid(), system = createEllipticStreamtubeGrid(input);
  assert.ok(system.residuals(input.nodes).residual > 1e-8);
  for (const backend of ['wasm', 'javascript']) {
    const result = smoothEllipticStreamtubeGrid(system, { omega: 1, maxSweeps: 20, tolerance: 1e-10,
      requireConvex: true, detectCoordinateRoundoff: true, backend });
    assert.equal(result.converged, false);
    assert.equal(result.reason, 'coordinate precision limit');
    assert.ok(result.coordinatePrecision.limited);
    assert.ok(result.coordinatePrecision.undampedUpdate <= Number.EPSILON);
    assert.equal(result.coordinatePrecision.residual, system.residuals(result.nodes).residual);
    assert.deepEqual(result.nodes, input.nodes, 'the exact affine geometry needs no correction');
    assert.ok(result.quality.valid);
  }
});

test('roundoff detection cannot accept a disturbed grid or an artificially damped update', () => {
  const input = thinAffineGrid(); input.nodes[2][2].x += .001;
  const system = createEllipticStreamtubeGrid(input);
  assert.equal(system.coordinateRoundoff(input.nodes, 1e-10).limited, false);
  for (const omega of [1, 1e-16]) {
    const result = smoothEllipticStreamtubeGrid(system, { omega, maxSweeps: 1, tolerance: 1e-10,
      detectCoordinateRoundoff: true });
    assert.equal(result.converged, false);
    assert.equal(result.coordinatePrecision, undefined);
  }
  const controlled = createEllipticStreamtubeGrid({ ...thinAffineGrid(),
    streamwiseSource: thinAffineGrid().nodes.map(row => row.map(() => .1)) });
  assert.equal(controlled.coordinateRoundoff(controlled.initial, 1e-10), null);
});
