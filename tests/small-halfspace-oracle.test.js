import test from 'node:test';
import assert from 'node:assert/strict';
import { projectSmallHalfspaces } from './oracles/small-halfspace-projection.js';

const face = (normal, lower) => ({ normal: Float64Array.from(normal), lower });
const close = (a, b, t = 2e-12) => assert.ok(Math.abs(a - b) < t, `${a} != ${b}`);

test('the independent QR face oracle resolves intersecting, redundant and replaced active constraints', () => {
  for (const [faces, expected] of [
    [[face([1, 0], 1), face([0, 1], 2)], [1, 2]],
    [[face([1, 0], 1), face([2, 0], 4), face([1, 1], 1)], [2, 0]],
    [[face([1, 1], 2), face([1, 0], -.1)], [1, 1]],
    [[face([1, 0], -1), face([0, 1], -1)], [0, 0]]
  ]) {
    const before = structuredClone(faces), r = projectSmallHalfspaces(faces);
    assert.equal(r.converged, true, r.reason); r.point.forEach((v, i) => close(v, expected[i]));
    assert.deepEqual(faces, before); assert.ok(r.maximumPrimalToleranceRatio <= 1);
  }
});

test('row-span QR resolves near-parallel normals without squaring their conditioning', () => {
  // These faces meet at (1,1). The row Gram matrix rounds away the 1e-16
  // angular distinction, whereas QR retains the 1e-8 transverse direction.
  const e = 1e-8, faces = [face([1, e], 1 + e), face([1, -e], 1 - e)];
  // Positive multipliers for both faces are impossible at (1,1): the true
  // projection activates only the stronger face, x=(1+e)*(1,e)/(1+e²).
  const r = projectSmallHalfspaces(faces);
  assert.equal(r.converged, true, r.reason);
  close(r.point[0], (1 + e) / (1 + e * e)); close(r.point[1], e * (1 + e) / (1 + e * e));
  assert.equal(r.rank, 2); assert.equal(r.active, 1);
});

test('inconsistent faces and an insufficient norm bound cannot certify a projection', () => {
  assert.equal(projectSmallHalfspaces([face([1, 0], 1), face([-1, 0], 1)]).converged, false);
  assert.equal(projectSmallHalfspaces([face([1, 1], 2)], { maximumNorm: 1 }).converged, false);
  assert.throws(() => projectSmallHalfspaces(Array.from({ length: 13 }, () => face([1], 1))), /1–12/);
});
