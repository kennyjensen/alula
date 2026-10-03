import test from 'node:test';
import assert from 'node:assert/strict';
import { streamtubeWakeGap } from '../src/euler/streamtube-wake-geometry.js';
import { createStreamtubeBodyLayout } from '../src/euler/streamtube-body-layout.js';
const close = (a, b, tol = 2e-10) => assert.ok(Math.abs(a - b) < tol, `${a} != ${b}`);

test('normal wake gap permits tangential offset, rotates and scales, and is regular at zero thickness', () => {
  for (const angle of [0, .6, 2]) for (const scale of [.03, 1, 40]) for (const gap of [0, 1e-12, .2]) {
    const t = { x: Math.cos(angle), y: Math.sin(angle) }, n = { x: -t.y, y: t.x };
    const point = (s, d) => ({ x: 2 + scale * (s * t.x + d * n.x), y: -3 + scale * (s * t.y + d * n.y) });
    const lower = [-1, 0, 1].map(s => point(s, 0)), upper = [-1, 0, 1].map(s => point(s + .17, gap));
    const r = streamtubeWakeGap(lower, upper); close(r.gap, scale * gap); close(r.tangentialOffset, scale * .17);
    close(r.apply(lower.map(() => ({ x: 0, y: 0 })), upper.map(() => n)), 1);
    close(streamtubeWakeGap(upper, lower).gap, -scale * gap);
    assert.ok(Math.hypot(upper[1].x - lower[1].x, upper[1].y - lower[1].y) > r.gap);
  }
});

test('gap derivative includes all bank coordinates and the rotating mean secant', () => {
  const banks = [[{ x: -.7, y: -.2 }, { x: .1, y: .04 }, { x: 1.3, y: .5 }],
    [{ x: -.6, y: .04 }, { x: .4, y: .3 }, { x: 1.6, y: .65 }]];
  const r = streamtubeWakeGap(...banks), h = 1e-6;
  for (let b = 0; b < 2; b++) for (let i = 0; i < 3; i++) for (const key of ['x', 'y']) {
    const p = structuredClone(banks), m = structuredClone(banks), d = banks.map(row => row.map(() => ({ x: 0, y: 0 })));
    p[b][i][key] += h; m[b][i][key] -= h; d[b][i][key] = 1;
    close(r.apply(...d), (streamtubeWakeGap(...p).gap - streamtubeWakeGap(...m).gap) / (2 * h));
  }
  assert.throws(() => streamtubeWakeGap([], []), /stencil/);
  assert.throws(() => streamtubeWakeGap(Array(3).fill({ x: 0, y: 0 }), Array(3).fill({ x: 0, y: 1 })), /Degenerate/);
});

test('independent downstream banks add exactly one unknown and one equation per station including the outlet', () => {
  for (const densityUnknowns of [false, true]) {
    const options = { segments: 12, tubes: [2, 3, 2], bodies: [{ leadingIndex: 3, trailingIndex: 8 }, { leadingIndex: 5, trailingIndex: 10 }], densityUnknowns };
    const a = createStreamtubeBodyLayout(options), b = createStreamtubeBodyLayout({ ...options, independentWakeBanks: true });
    assert.equal(b.n, a.n + 6); assert.equal(b.rows.length, b.n); assert.equal(b.rowCounts.wakeGap, 6);
    assert.equal(b.rowCounts.endTangency, a.rowCounts.endTangency);
    for (let body = 0; body < 2; body++) {
      assert.equal(b.nodes[body][0].at(-1).column, b.nodes[body + 1][0][0].column);
      assert.notEqual(b.nodes[body][12].at(-1).column, b.nodes[body + 1][12][0].column);
      assert.equal(b.rows.filter(r => r.kind === 'endTangency' && r.i === 12 && r.node.body === body).length, 1);
    }
  }
});
