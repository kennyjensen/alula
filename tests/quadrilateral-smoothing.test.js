import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { quadrilateralCornerShape } from '../src/geometry/quadrilateral-shape.js';
import { smoothQuadrilaterals } from '../src/geometry/smooth-quadrilaterals.js';
import { untangleQuadrilaterals } from '../src/geometry/untangle-quadrilaterals.js';

const fixture = JSON.parse(readFileSync(new URL('./fixtures/folded-main-flap-patch.json', import.meta.url)));
const distance = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);
const rectangle = () => {
  const vertices = Array.from({ length: 9 }, (_, i) => ({ x: 2 * (i % 3), y: Math.floor(i / 3) }));
  vertices[4] = { x: 2.35, y: .8 };
  const cells = [[0, 1, 4, 3], [1, 2, 5, 4], [3, 4, 7, 6], [4, 5, 8, 7]];
  return { vertices, cells, fixed: new Set([0, 1, 2, 3, 5, 6, 7, 8]) };
};

test('corner shape gradient and Hessian agree with independent differences in all three vertex positions', () => {
  const points = [{ x: -.2, y: .1 }, { x: .8, y: .3 }, { x: 1, y: 1.2 }], aspect = 2.3, h = 2e-6;
  const value = p => {
    const [a, b, c] = p, u = [b.x - a.x, b.y - a.y], v = [c.x - b.x, c.y - b.y];
    return (Math.hypot(...u) ** 2 / aspect + aspect * Math.hypot(...v) ** 2) / (2 * (u[0] * v[1] - u[1] * v[0]));
  };
  for (let slot = 0; slot < 3; slot++) {
    const r = quadrilateralCornerShape(points, aspect, slot);
    for (const [k, axis] of ['x', 'y'].entries()) {
      const a = structuredClone(points), b = structuredClone(points); a[slot][axis] += h; b[slot][axis] -= h;
      assert.ok(Math.abs((value(a) - value(b)) / (2 * h) - r.gradient[k]) < 2e-9);
      const ga = quadrilateralCornerShape(a, aspect, slot).gradient, gb = quadrilateralCornerShape(b, aspect, slot).gradient;
      for (let row = 0; row < 2; row++) assert.ok(Math.abs((ga[row] - gb[row]) / (2 * h) - r.hessian[2 * row + k]) < 2e-8);
    }
    assert.ok(r.hessian[0] > 0 && r.hessian[0] * r.hessian[3] - r.hessian[1] ** 2 > 0);
  }
  const orthogonal = [{ x: 0, y: 0 }, { x: 2, y: 0 }, { x: 2, y: 1 }];
  assert.equal(quadrilateralCornerShape(orthogonal, 2, 0).value, 1);
  assert.throws(() => quadrilateralCornerShape([...orthogonal].reverse(), 2, 0), /positive/);
});

test('fixed rectangular boundary recovers its independently known orthogonal grid minimum', () => {
  const input = rectangle(), before = structuredClone(input), aspects = input.cells.map(() => [2, .5, 2, .5]);
  const r = smoothQuadrilaterals(input, { aspects });
  assert.equal(r.converged, true, r.reason); assert.ok(distance(r.vertices[4], { x: 2, y: 1 }) < 2e-8);
  assert.ok(Math.abs(r.quality.energy - 16) < 1e-12); assert.deepEqual(input, before);
  input.fixed.forEach(id => assert.deepEqual(r.vertices[id], input.vertices[id]));
  assert.match(r.status, /flow equations have not been solved/);
  assert.throws(() => smoothQuadrilaterals(input, { aspects: [[1]] }), /Invalid/);
});

test('repaired main/flap patch decreases shape energy while keeping every corner positive and boundaries fixed', () => {
  const fixed = new Set(fixture.fixed), repaired = untangleQuadrilaterals({ ...fixture, fixed });
  const input = { ...fixture, fixed, vertices: repaired.vertices }, before = structuredClone(input);
  const r = smoothQuadrilaterals(input);
  assert.ok(r.quality.energy < .5 * r.history[0].energy);
  for (const [i, h] of r.history.entries()) {
    assert.ok(h.minCornerSine > 0); if (i) assert.ok(h.energy <= r.history[i - 1].energy * (1 + 1e-13));
  }
  for (const cell of input.cells) {
    const p = cell.map(i => r.vertices[i]);
    p.forEach((a, i) => { const b = p[(i + 1) % 4], c = p[(i + 2) % 4];
      assert.ok((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) > 0); });
  }
  assert.deepEqual(input, before); fixed.forEach(id => assert.deepEqual(r.vertices[id], input.vertices[id]));
  assert.throws(() => smoothQuadrilaterals({ ...fixture, fixed }), /positive/);
});

test('shape smoothing preserves rotation, translation and physical units', () => {
  const input = rectangle(), angle = .4, scale = 3;
  const move = p => ({ x: 7 + scale * (p.x * Math.cos(angle) - p.y * Math.sin(angle)),
    y: -2 + scale * (p.x * Math.sin(angle) + p.y * Math.cos(angle)) });
  const a = smoothQuadrilaterals(input), b = smoothQuadrilaterals({ ...input, vertices: input.vertices.map(move) });
  assert.equal(a.converged, true, a.reason); assert.equal(b.converged, true, b.reason);
  a.vertices.forEach((p, i) => assert.ok(distance(move(p), b.vertices[i]) < 2e-8));
  assert.ok(Math.abs(a.quality.energy - b.quality.energy) < 1e-12);
});
