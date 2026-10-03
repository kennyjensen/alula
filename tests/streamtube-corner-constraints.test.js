import test from 'node:test';
import assert from 'node:assert/strict';
import { streamtubeCornerConstraints } from '../src/euler/streamtube-corner-constraints.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

test('all signed corner areas and Cartesian derivatives match independent quadrilateral differences', () => {
  const x = Float64Array.of(0, 0, 2, .1, 1.6, 1, -.2, .8);
  const ids = [[0, 3], [1, 2]], euler = { layout: { n: 8 },
    decode: z => ({ nodes: [ids.map(row => row.map(k => ({ x: z[2 * k], y: z[2 * k + 1] })))] }),
    geometryDerivatives: () => [ids.map(row => row.map(k => new Map([[2 * k, { x: 1, y: 0 }], [2 * k + 1, { x: 0, y: 1 }]])))] };
  const independent = z => Array.from({ length: 4 }, (_, k) => {
    const a = k, b = (k + 1) % 4, c = (k + 2) % 4;
    // Twice signed triangle area, in a different algebraic form.
    return z[2 * a] * (z[2 * b + 1] - z[2 * c + 1]) + z[2 * b] * (z[2 * c + 1] - z[2 * a + 1])
      + z[2 * c] * (z[2 * a + 1] - z[2 * b + 1]);
  });
  const corners = streamtubeCornerConstraints(euler, x);
  for (let col = 0; col < x.length; col++) {
    const plus = x.slice(), minus = x.slice(), h = 1e-4; plus[col] += h; minus[col] -= h;
    const a = independent(plus), b = independent(minus);
    corners.forEach((corner, k) => {
      assert.ok(Math.abs(corner.value - independent(x)[k]) < 1e-14);
      assert.ok(Math.abs((a[k] - b[k]) / (2 * h) - (corner.gradient.get(col) ?? 0)) < 1e-10);
    });
  }
  // A deliberately sizeable angular margin makes its norm derivative
  // measurable; the operational 1e-12 margin uses the identical formula.
  const margin = .15, bounded = streamtubeCornerConstraints(euler, x, { minimumSine: margin });
  const boundedAreas = z => independent(z).map((area, k) => {
    const a = k, b = (k + 1) % 4, c = (k + 2) % 4;
    return area - margin * Math.hypot(z[2 * b] - z[2 * a], z[2 * b + 1] - z[2 * a + 1])
      * Math.hypot(z[2 * c] - z[2 * b], z[2 * c + 1] - z[2 * b + 1]);
  });
  for (let col = 0; col < x.length; col++) {
    const plus = x.slice(), minus = x.slice(), h = 1e-5; plus[col] += h; minus[col] -= h;
    const a = boundedAreas(plus), b = boundedAreas(minus);
    bounded.forEach((corner, k) => assert.ok(Math.abs((a[k] - b[k]) / (2 * h) - (corner.gradient.get(col) ?? 0)) < 1e-9));
  }
  x[4] = -.5; x[5] = .2;
  assert.ok(streamtubeCornerConstraints(euler, x).some(c => c.value < 0), 'folded geometry must remain signed');
});

test('coupled corner derivatives include normal motion, stagnation and all surface/wake displacement chains', () => {
  const system = createCoupledStreamtubeBody(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }));
  const { euler, bl, ne } = system, x = system.initial, baseline = system.residual(x), before = x.slice();
  const constraints = streamtubeCornerConstraints(euler, x.subarray(0, ne), { displacementMap: bl.thicknessMap });
  const directions = [
    new Map(euler.layout.positions.map(({ column }, i) => [column, .02 * Math.sin(i + 1)])),
    new Map(euler.layout.globals.stagnation.filter(c => c !== null).map((c, i) => [c, .02 * (i + 1)])),
    new Map(bl.stations.filter(s => s.kind === 'surface').map(s => [ne + 4 * s.id + 2, .2 * Math.sin(s.id + 1)])),
    new Map(bl.stations.filter(s => s.kind !== 'surface').map(s => [ne + 4 * s.id + 2, .2 * Math.cos(s.id + 1)])),
  ];
  const areas = z => {
    euler.setDisplacement(bl.thicknesses(z.subarray(ne)));
    const { nodes } = euler.decode(z.subarray(0, ne));
    return constraints.map(({ g, i, j, k }) => {
      const quad = [nodes[g][i][j], nodes[g][i + 1][j], nodes[g][i + 1][j + 1], nodes[g][i][j + 1]];
      const [a, b, c] = [quad[k], quad[(k + 1) % 4], quad[(k + 2) % 4]];
      return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    });
  };
  for (const direction of directions) {
    assert.ok(direction.size > 0);
    const h = 1e-4, samples = [-2, -1, 1, 2].map(m => areas(x.map((v, col) => v + m * h * (direction.get(col) ?? 0))));
    let maximum = 0;
    constraints.forEach((c, i) => {
      const exact = [...c.gradient].reduce((s, [col, d]) => s + d * (direction.get(col) ?? 0), 0);
      const fd = (samples[0][i] - 8 * samples[1][i] + 8 * samples[2][i] - samples[3][i]) / (12 * h);
      maximum = Math.max(maximum, Math.abs(exact));
      assert.ok(Math.abs(fd - exact) < 2e-9 * Math.max(1, Math.abs(exact)), `${c.g}/${c.i}/${c.j}/${c.k}: ${fd} != ${exact}`);
    });
    assert.ok(maximum > 0, 'the checked direction must move at least one corner');
  }
  assert.deepEqual(x, before); assert.deepEqual(system.residual(x), baseline);
});
