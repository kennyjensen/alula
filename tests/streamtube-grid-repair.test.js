import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { untangleQuadrilaterals } from '../src/geometry/untangle-quadrilaterals.js';
import { repairStreamtubeGrid, smoothStreamtubeGrid } from '../src/euler/streamtube-grid-repair.js';
import { streamtubeMeshConnectivity } from '../src/geometry/streamtube-mesh-connectivity.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const fixture = JSON.parse(readFileSync(new URL('./fixtures/folded-main-flap-patch.json', import.meta.url)));
const mesh = () => ({ ...fixture, fixed: new Set(fixture.fixed) });
const distance = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);
const badCells = (vertices, cells) => cells.filter(ids => {
  const p = ids.map(i => vertices[i]);
  // Independent corner triangle orientation; no production quality helper.
  return p.some((a, i) => {
    const b = p[(i + 1) % 4], c = p[(i + 2) % 4];
    return !((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) > 0);
  });
}).length;

test('single free corner repair equals the independently derived closest feasible point', () => {
  const vertices = [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: -.2, y: 1 }, { x: 0, y: 1 }];
  const r = untangleQuadrilaterals({ vertices, cells: [[0, 1, 2, 3]], fixed: new Set([0, 1, 3]) });
  // The active inequalities are x >= .001*.2 and
  // x+y-1 >= .001*.2*sqrt(1.2^2+1), using initial edge products.
  const expected = { x: .0002, y: 1 + .0002 * Math.hypot(1.2, 1) - .0002 };
  assert.equal(r.converged, true, r.reason); assert.ok(distance(r.vertices[2], expected) < 2e-12);
  [0, 1, 3].forEach(i => assert.deepEqual(r.vertices[i], vertices[i]));
  assert.equal(badCells(r.vertices, [[0, 1, 2, 3]]), 0);
});

test('fixed and movement-limited infeasible grids cannot be reported as repaired; valid meshes are unchanged', () => {
  const vertices = [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: -.2, y: 1 }, { x: 0, y: 1 }], cells = [[0, 1, 2, 3]];
  const fixed = untangleQuadrilaterals({ vertices, cells, fixed: new Set([0, 1, 2, 3]) });
  assert.equal(fixed.converged, false); assert.equal(fixed.quality.valid, false); assert.deepEqual(fixed.vertices, vertices);
  const limited = untangleQuadrilaterals({ vertices, cells, fixed: new Set([0, 1, 3]) }, { displacementLimit: .01, maxSweeps: 5 });
  assert.equal(limited.converged, false); assert.equal(limited.quality.valid, false);
  const valid = vertices.map((p, i) => i === 2 ? { x: 1, y: 1 } : p);
  const identity = untangleQuadrilaterals({ vertices: valid, cells });
  assert.equal(identity.converged, true); assert.deepEqual(identity.vertices, valid); assert.equal(identity.history.length, 1);
  assert.throws(() => untangleQuadrilaterals({ vertices, cells: [[0, 1, 1, 3]] }), /Invalid/);
  assert.throws(() => untangleQuadrilaterals({ vertices, cells }, { minimumCorner: 0 }), /Invalid/);
});

test('frozen 28-cell main/flap patch removes nine invalid cells with 11 free vertices and unchanged boundaries', () => {
  const input = mesh(), before = structuredClone(input), r = untangleQuadrilaterals(input);
  assert.equal(input.vertices.length, 48); assert.equal(input.cells.length, 28); assert.equal(input.vertices.length - input.fixed.size, 11);
  assert.equal(badCells(input.vertices, input.cells), 9);
  assert.equal(r.converged, true, r.reason); assert.equal(badCells(r.vertices, input.cells), 0);
  assert.deepEqual(input, before); input.fixed.forEach(i => assert.deepEqual(r.vertices[i], input.vertices[i]));
  assert.ok(r.maxDisplacement < .0035); assert.ok(r.quality.minCornerSine > 4e-4);
  // The minimax update should not recreate a larger normalized defect,
  // including when a feasible polygon collapses to a segment or point.
  r.history.slice(1).forEach((h, i) => assert.ok(h.maximumDefect <= r.history[i].maximumDefect + 1e-10));
  assert.match(r.status, /flow equations have not been solved/);
});

test('repair is invariant to rotation, translation and units, including short-edge roundoff', () => {
  const input = mesh(), a = untangleQuadrilaterals(input), angle = .4, length = 3;
  const move = p => ({ x: 7 + length * (p.x * Math.cos(angle) - p.y * Math.sin(angle)),
    y: -2 + length * (p.x * Math.sin(angle) + p.y * Math.cos(angle)) });
  const b = untangleQuadrilaterals({ ...input, vertices: input.vertices.map(move) });
  assert.equal(a.converged, true); assert.equal(b.converged, true, b.reason);
  a.vertices.forEach((p, i) => assert.ok(distance(move(p), b.vertices[i]) < 2e-10));
  assert.ok(Math.abs(length * a.maxDisplacement - b.maxDisplacement) < 2e-10);
});

test('streamtube wrapper preserves walls, outer boundaries and shared cuts, and rejects disconnected input', () => {
  const system = createStreamtubeBodySystem(intrinsicBodyFixture({ elements: 2 })), nodes = system.decode(system.initial).nodes;
  const before = structuredClone(nodes), chart = system.geometryChart(), r = repairStreamtubeGrid({ system, nodes });
  assert.equal(r.converged, true); assert.deepEqual(nodes, before); assert.deepEqual(system.geometryChart(), chart);
  assert.deepEqual(r.nodes, nodes); assert.equal(streamtubeMeshSnapshot({ system, nodes: r.nodes }).quality.valid, true);
  for (let b = 0; b < 2; b++) for (let i = 0; i <= system.layout.nx; i++) if (!system.layout.active(b, i))
    assert.deepEqual(r.nodes[b][i].at(-1), r.nodes[b + 1][i][0]);
  const bad = structuredClone(nodes); bad[1][0][0] = { ...bad[1][0][0], y: bad[1][0][0].y + .01 };
  assert.throws(() => streamtubeMeshConnectivity(system.layout, bad), /Disconnected/);
  const smoothed = smoothStreamtubeGrid({ system, nodes }, { maxSweeps: 2 });
  assert.deepEqual(system.geometryChart(), chart); assert.deepEqual(nodes, before);
  const originalMesh = streamtubeMeshConnectivity(system.layout, nodes), resultMesh = streamtubeMeshConnectivity(system.layout, smoothed.nodes);
  assert.deepEqual(resultMesh.cells, originalMesh.cells);
  originalMesh.fixed.forEach(i => assert.deepEqual(resultMesh.vertices[i], originalMesh.vertices[i]));
  assert.equal(streamtubeMeshSnapshot({ system, nodes: smoothed.nodes }).quality.valid, true);
});

test('frozen first Newton direction reproduces the gap corner collapse without a flow solve', () => {
  const report = JSON.parse(readFileSync(new URL('../docs/streamtube-grid-repair-diagnosis.json', import.meta.url)));
  const c = report.firstDirection.corner;
  assert.equal(report.repair.converged, true); assert.equal(report.flow.converged, false);
  assert.ok(report.firstDirection.linearResidual < 1e-10);
  const determinantAt = step => {
    const p = c.points.map((v, i) => ({ x: v.x + step * c.motion[i].x, y: v.y + step * c.motion[i].y }));
    const [a, b, d] = [0, 1, 2].map(k => p[(c.corner + k) % 4]);
    return (b.x - a.x) * (d.y - a.y) - (b.y - a.y) * (d.x - a.x);
  };
  assert.ok(Math.abs(determinantAt(0) - c.det) < 1e-19);
  // Symmetric differences are exact for the quadratic determinant along
  // these frozen physical tangents (apart from floating-point roundoff).
  const h = 1e-3, derivative = (determinantAt(h) - determinantAt(-h)) / (2 * h);
  assert.ok(Math.abs(derivative - c.derivative) < 1e-14);
  assert.ok(determinantAt(.0006) > 0); assert.ok(determinantAt(.0008) < 0);
  let interior = 0, global = 0;
  for (const term of c.contributions) {
    if (term.col === term.node.column) interior += term.derivative * term.step;
    else global += term.derivative * term.step;
  }
  assert.ok(Math.abs(interior + global - c.derivative) < 1e-17);
  assert.ok(interior < 0); assert.ok(global > 0);
  // This diagnoses an accurate Newton direction leaving the convex set;
  // it does not certify a finite update of the nonlinear geometry chart.
});
