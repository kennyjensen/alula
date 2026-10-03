import test from 'node:test';
import assert from 'node:assert/strict';
import { naca4, transform } from '../src/geometry/airfoil.js';
import { solveStreamtubeAssembly } from '../src/euler/streamtube-result.js';
import { streamtubeAuditRegions } from '../src/geometry/streamtube-grid-audit.js';

test('the explicit potential-block cosine main/flap SLOR mesh and wall subdivision preserve the two holes without a flow solve', () => {
  const input = { elements: [{ points: naca4('2412', 160) }, { points: transform(naca4('0012', 160), { chord: .3, x: .94, y: -.08, angle: -15 }) }],
    alpha: 4, gridEllipticSmoothing: true, gridCrosslinePlacement: 'potential',
    gridInletIntervals: 32, gridOutletIntervals: 32 }, before = structuredClone(input);
  let initialMesh;
  const r = solveStreamtubeAssembly(input, { meshOnly: true,
    onMesh: (mesh, stage) => { if (stage === 'initial') initialMesh = mesh; },
    onIteration: () => assert.fail('Mesh preview must not enter Euler') });
  assert.deepEqual(input, before); assert.equal(r.status, 'mesh-ready');
  assert.equal(r.mesh.initialization.flowSolved, false); assert.equal(r.flow, undefined);
  const { mesh } = r, init = mesh.initialization;
  assert.equal(init.gridSmoothing.converged, true);
  assert.deepEqual(init.tubes, [10, 13, 10]);
  assert.equal(init.gridSpacing.inlet.intervals, 32); assert.equal(init.gridSpacing.outlet.intervals, 32);
  assert.ok(init.potentialCrosslines.turnResolution.some(r => r.addedPoints > 0));
  // Check all farfield intervals, including the formerly discontinuous
  // LE/TE joins. This is a spacing screen, not a flow-accuracy assertion.
  const regions = streamtubeAuditRegions(initialMesh);
  for (const [g, side] of [[0, 0], [2, init.tubes[2]]]) {
    const points = regions[g].nodes.map(row => row[side]);
    const length = points.slice(1).map((p, i) => Math.hypot(p.x - points[i].x, p.y - points[i].y));
    for (let i = 1; i < length.length; i++) assert.ok(Math.max(length[i] / length[i - 1], length[i - 1] / length[i]) <= 1.5);
  }
  const edges = new Map();
  for (const cell of mesh.cells) {
    assert.equal(new Set(cell).size, 4);
    for (let k = 0; k < 4; k++) {
      const [a, b, c] = [0, 1, 2].map(j => mesh.vertices[cell[(k + j) % 4]]);
      assert.ok((b.x - a.x) * (c.y - b.y) > (b.y - a.y) * (c.x - b.x));
      const u = cell[k], v = cell[(k + 1) % 4], key = u < v ? `${u}:${v}` : `${v}:${u}`;
      edges.set(key, (edges.get(key) ?? 0) + 1);
    }
  }
  assert.ok([...edges.values()].every(n => n === 1 || n === 2));
  assert.equal(mesh.vertices.length - edges.size + mesh.cells.length, -1);
});
