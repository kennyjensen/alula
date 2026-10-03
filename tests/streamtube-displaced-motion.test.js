import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { streamtubeMotionDirections } from '../src/euler/streamtube-geometry.js';

test('a restored coupled chart follows each displaced wake centerline, including the shifted TE', () => {
  const f = JSON.parse(fs.readFileSync(new URL('fixtures/default-coupled-ises-wake-fold.json', import.meta.url)));
  const system = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: Float64Array.from(f.initialBL) });
  const x = system.initial, v = system.evaluate(x), { layout, conditions } = system.euler;
  for (const key in f.families) assert.ok(Math.abs(v.families[key] - f.families[key]) < 1e-10);
  const check = state => {
    const flow = system.evaluate(state).outer, chart = new Map(system.euler.geometryChart().map(p => [p.column, p.normal]));
    for (const node of layout.positions.filter(p => p.kind === 'cut')) {
      const left = Math.max(0, node.i - 1), right = Math.min(layout.nx, node.i + 1), b = node.body;
      // Independently form the centerline chord as the sum of physical
      // upper/lower bank chords; any normal must have zero tangent product.
      const lower = flow.nodes[b], upper = flow.nodes[b + 1];
      const tx = lower[right].at(-1).x - lower[left].at(-1).x + upper[right][0].x - upper[left][0].x;
      const ty = lower[right].at(-1).y - lower[left].at(-1).y + upper[right][0].y - upper[left][0].y;
      const n = chart.get(node.column), length = Math.hypot(tx, ty);
      assert.ok(Math.abs((tx * n.x + ty * n.y) / length) < 2e-12);
      assert.ok(tx * n.y - ty * n.x > 0);
    }
  };
  check(x);
  const body = layout.bodies[0], i = body.trailingIndex + 1, col = layout.nodes[0][i].at(-1).column;
  const raw = v.outer.undisplacedNodes[0], a = raw[i - 1].at(-1), b = raw[i + 1].at(-1);
  const n = system.euler.geometryChart().find(p => p.column === col).normal;
  assert.ok(Math.abs(((b.x - a.x) * n.x + (b.y - a.y) * n.y) / Math.hypot(b.x - a.x, b.y - a.y)) > .6,
    'The frozen case must distinguish displaced and bare-TE directions.');
  const h = 1e-7, plus = x.slice(), minus = x.slice(); plus[col] += h; minus[col] -= h;
  const p = system.euler.decode(plus.subarray(0, system.ne)), m = system.euler.decode(minus.subarray(0, system.ne));
  for (const key of ['x', 'y']) {
    const derivative = (p.nodes[0][i].at(-1)[key] + p.nodes[1][i][0][key] - m.nodes[0][i].at(-1)[key] - m.nodes[1][i][0][key]) / (4 * h);
    assert.ok(Math.abs(derivative - conditions.lengthScale * n[key]) < 2e-9);
  }
  const rebased = system.rebase(x); check(rebased);
  const after = system.evaluate(rebased);
  assert.ok(Math.max(...v.residual.map((r, k) => Math.abs(r - after.residual[k]))) < 1e-10);
  v.outer.nodes.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => {
    const q = after.outer.nodes[g][i][j]; assert.ok(Math.hypot(p.x - q.x, p.y - q.y) < 1e-12);
  })));
  const angle = .73, c = Math.cos(angle), s = Math.sin(angle), move = p => ({ x: 4 + 3.2 * (c * p.x - s * p.y), y: -2 + 3.2 * (s * p.x + c * p.y) });
  const transformed = streamtubeMotionDirections(layout, v.outer.nodes.map(g => g.map(row => row.map(move))), 'body-stations');
  for (const { column, normal } of system.euler.geometryChart()) {
    const rotated = transformed.get(column);
    assert.ok(Math.abs(rotated.x - c * normal.x + s * normal.y) < 2e-11);
    assert.ok(Math.abs(rotated.y - s * normal.x - c * normal.y) < 2e-11);
  }
});
