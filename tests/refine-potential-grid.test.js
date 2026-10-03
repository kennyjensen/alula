import test from 'node:test';
import assert from 'node:assert/strict';
import { refinePotentialGrid } from '../src/geometry/tests/refine-potential-grid.js';
import { naca4 } from '../src/geometry/airfoil.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { createInitialStreamtubeTopology } from '../src/geometry/streamtube-topology.js';
import { createPanelStreamtubeGrid } from '../src/euler/streamtube-body-initializer.js';

const map = (u, v) => ({ x: u, y: v + .12 * Math.sin(Math.PI * u) * (1 - v) });
const coordinate = (u, v) => ({ x: u + .3 * u * u, y: v });
function fixture() {
  const nx = 6, nt = 4, levels = Array.from({ length: nt + 1 }, (_, j) => Math.expm1(j / nt) / Math.expm1(1));
  const at = (i, j) => {
    const k = Math.min(nt - 1, Math.floor(j)), v = levels[k] + (j - k) * (levels[k + 1] - levels[k]);
    return { point: map(i / nx, v), coordinate: coordinate(i / nx, v) };
  };
  const nodes = Array.from({ length: nx + 1 }, (_, i) => levels.map(v => map(i / nx, v)));
  const coordinates = Array.from({ length: nx + 1 }, (_, i) => levels.map(v => coordinate(i / nx, v)));
  return { nx, nt, levels, at, nodes, coordinates };
}

test('physical boundary refinement exactly reproduces a curved Coons map and preserves every coarse node and label', () => {
  const f = fixture(), before = structuredClone({ nodes: f.nodes, coordinates: f.coordinates }); let calls = 0;
  const r = refinePotentialGrid(f, { boundaryAt: (i, j) => {
    assert.ok(!Number.isInteger(i) || !Number.isInteger(j), 'Existing boundary nodes must be copied without resampling'); calls++; return f.at(i, j);
  } });
  assert.ok(calls > 0); assert.equal(r.quality.valid, true); assert.equal(r.computationalQuality.valid, true);
  r.nodes.forEach((row, i) => row.forEach((p, j) => {
    const exact = f.at(i / 2, j / 2), w = r.coordinates[i][j];
    assert.ok(Math.hypot(p.x - exact.point.x, p.y - exact.point.y) < 5e-15);
    assert.ok(Math.hypot(w.x - exact.coordinate.x, w.y - exact.coordinate.y) < 5e-15);
    if (i % 2 === 0 && j % 2 === 0) { assert.deepEqual(p, f.nodes[i / 2][j / 2]); assert.deepEqual(w, f.coordinates[i / 2][j / 2]); }
  }));
  assert.deepEqual({ nodes: f.nodes, coordinates: f.coordinates }, before);
});

test('refinement refuses inconsistent transverse labels, infeasible parents and excessive work', () => {
  const f = fixture();
  assert.throws(() => refinePotentialGrid(f, { maxNodes: 1, boundaryAt: () => assert.fail('Budget check must precede boundary work') }), /budget/);
  assert.throws(() => refinePotentialGrid(f, { boundaryAt: (i, j) => { const p = f.at(i, j); p.coordinate.y += .01; return p; } }), /streamfunction labels/);
  const folded = structuredClone(f.nodes); folded[2][2].x += 5;
  assert.throws(() => refinePotentialGrid({ nodes: folded, coordinates: f.coordinates }, { boundaryAt: f.at }), /positive parent/);
  assert.throws(() => refinePotentialGrid(f, { factor: 0, boundaryAt: f.at }), /controls/);
});

test('panel boundary sampling refines the true surface parameter and evaluates its potential instead of splitting potential uniformly', () => {
  const input = { elements: [{ points: naca4('0012', 80) }], alpha: 4 };
  const r = createPanelStreamtubeGrid(createInitialStreamtubeTopology(input, { surfaceIntervals: 16, tubes: 7 }), { recordPotentialCoordinates: true });
  const before = structuredClone({ nodes: r.nodes, coordinates: r.potentialCoordinates }), body = r.input.bodies[0], i = body.leadingIndex, j = r.nodes[0][0].length - 1;
  const midpoint = r.potentialBoundaryAt(0, i + .5, j), fractions = r.system.fractions[0].lower;
  const curve = createContourCurve(body.points), exact = curve.branch('lower', .5 * (fractions[0] + fractions[1]), body.stagnationParameter);
  assert.ok(Math.hypot(midpoint.point.x - exact.point.x, midpoint.point.y - exact.point.y) < 1e-13);
  const a = r.potentialCoordinates[0][i][j], b = r.potentialCoordinates[0][i + 1][j];
  assert.ok(midpoint.coordinate.x > a.x && midpoint.coordinate.x < b.x);
  assert.ok(Math.abs((midpoint.coordinate.x - a.x) / (b.x - a.x) - .5) > .05, 'Stagnation potential cannot stand in for physical curve distance');
  assert.equal(midpoint.coordinate.y, a.y);
  // These eighth-interval samples lie nearer the edge than the original
  // one-fifth-interval trace seed. They must be integrated, not extrapolated.
  for (const station of [body.leadingIndex - .125, body.trailingIndex + .125]) {
    const value = r.potentialBoundaryAt(0, station, j), k = Math.floor(station), t = station - k;
    const x = r.nodes[0][k][j].x + t * (r.nodes[0][k + 1][j].x - r.nodes[0][k][j].x);
    assert.ok(Math.abs(value.point.x - x) < 3e-9);
    assert.ok(Math.abs(r.guideField.streamfunctionAt(value.point) - value.coordinate.y) < 5e-8);
    assert.equal(r.guideField.admissibleNode(value.point), true);
  }
  const anchor = r.potentialBoundaryAt(0, i, j); assert.deepEqual(anchor.point, r.nodes[0][i][j]); anchor.point.x += 1;
  assert.deepEqual({ nodes: r.nodes, coordinates: r.potentialCoordinates }, before);
  assert.throws(() => r.potentialBoundaryAt(0, i, 1), /boundary index/);
  assert.throws(() => r.potentialBoundaryAt(-1, 0, 0), /boundary index/);
  assert.throws(() => r.potentialBoundaryAt(0, i, j + 1), /boundary index/);
});
