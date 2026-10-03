import test from 'node:test';
import assert from 'node:assert/strict';
import { captureStreamtubeInletFractions, adjustStreamtubeInlets, dekinkStreamtubeInteriors } from '../src/geometry/streamtube-grid-maintenance.js';
import { streamtubeGridConvexity } from '../src/geometry/streamtube-convex-step.js';

const grid = () => Array.from({ length: 3 }, (_, g) => Array.from({ length: 6 }, (_, i) => Array.from({ length: 4 }, (_, j) => ({ x: i, y: g + j / 3 }))));
const bodies = [{ leadingIndex: 3 }, { leadingIndex: 4 }];
const near = (a, b) => assert.ok(Math.abs(a - b) < 2e-13, `${a} != ${b}`);

test('inlet adjustment restores known straight-cut spacing on every body and preserves both shared banks', () => {
  const initial = grid(), fractions = captureStreamtubeInletFractions(initial, bodies), moved = structuredClone(initial);
  moved[0][1][3].x += .2; moved[1][1][0].x += .2;
  moved[1][2][3].x -= .3; moved[2][2][0].x -= .3;
  const before = structuredClone(moved), r = adjustStreamtubeInlets(moved, bodies, fractions);
  assert.deepEqual(moved, before); near(r.maxArcErrorBefore, .3); assert.ok(r.maxArcErrorAfter < 1e-13);
  r.nodes.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => {
    near(p.x, initial[g][i][j].x); near(p.y, initial[g][i][j].y);
  })));
});

test('DEKINK uses original neighbors in one pass, preserving cuts, walls and domain endpoints', () => {
  const nodes = grid(); nodes[0][2][1].x = .75; nodes[0][3][1].x = 2;
  nodes[0][2][0].x = .75; // A bank kink is deliberately left alone.
  const before = structuredClone(nodes), r = dekinkStreamtubeInteriors(nodes);
  assert.deepEqual(nodes, before); assert.equal(r.repairs.length, 2);
  near(r.nodes[0][1][1].x, .375); near(r.nodes[0][2][1].x, 1.5);
  assert.deepEqual(r.nodes[0][2][0], nodes[0][2][0]);
  for (let g = 0; g < nodes.length; g++) for (let i = 0; i < nodes[g].length; i++) for (let j = 0; j < 4; j++)
    if (!r.repairs.some(p => p.group === g && p.i === i && p.j === j)) assert.deepEqual(r.nodes[g][i][j], nodes[g][i][j]);
});

test('maintenance preserves rotation/scale covariance and rejects invalid cut data', () => {
  const nodes = grid(), angle = .6, c = Math.cos(angle), s = Math.sin(angle);
  const transform = p => ({ x: 3 * (c * p.x - s * p.y) + 4, y: 3 * (s * p.x + c * p.y) - 2 });
  const fractions = captureStreamtubeInletFractions(nodes, bodies);
  nodes[0][1][3].x += .2; nodes[1][1][0].x += .2; nodes[0][2][1].x = .75;
  const a = dekinkStreamtubeInteriors(adjustStreamtubeInlets(nodes, bodies, fractions).nodes);
  const rotated = nodes.map(group => group.map(row => row.map(transform)));
  const b = dekinkStreamtubeInteriors(adjustStreamtubeInlets(rotated, bodies, fractions).nodes);
  a.nodes.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => { const q = transform(p); near(b.nodes[g][i][j].x, q.x); near(b.nodes[g][i][j].y, q.y); })));
  const bad = structuredClone(nodes); bad[0][1][3].y += .01;
  assert.throws(() => captureStreamtubeInletFractions(bad, bodies), /Disconnected/);
  assert.throws(() => adjustStreamtubeInlets(nodes, bodies, [[0, .5, .4, 1], fractions[1]]), /fractions/);
});

test('inlet redistribution follows the curved cut without extrapolating through its neighbor', () => {
  const cut = [{ x: 0, y: 1 }, { x: 1, y: 1 }, { x: 2, y: 2 }];
  const nodes = [cut.map((p, i) => [{ x: i, y: 0 }, { ...p }]),
    cut.map((p, i) => [{ ...p }, { x: i, y: 3 }])];
  const before = structuredClone(nodes), body = [{ leadingIndex: 2 }], fractions = [[0, .9, 1]];
  assert.equal(streamtubeGridConvexity(nodes).valid, true);
  const moved = adjustStreamtubeInlets(nodes, body, fractions, { preserveConvexity: true });
  const fraction = (.9 * (1 + Math.SQRT2) - 1) / Math.SQRT2;
  near(moved.nodes[0][1][1].x, 1 + fraction);
  near(moved.nodes[0][1][1].y, 1 + fraction);
  assert.equal(streamtubeGridConvexity(moved.nodes).valid, true);
  assert.equal(moved.convexity, undefined, 'the curve-following movement needs no clipping here');
  assert.equal(moved.reparameterization.method, 'polyline-arclength');
  assert.deepEqual(moved.nodes[0][1][1], moved.nodes[1][1][0]);
  assert.deepEqual(nodes, before);
});

test('inlet redistribution retains clipping when the curve-following alternative would also fold a cell', () => {
  const cut = [1, .25, .25, 1].map((y, x) => ({ x, y }));
  const nodes = [cut.map((p, i) => [{ x: i, y: 0 }, { ...p }]),
    cut.map((p, i) => [{ ...p }, { x: i, y: 2 }])];
  const before = structuredClone(nodes), body = [{ leadingIndex: 3 }], fractions = [[0, .01, .1, 1]];
  assert.equal(streamtubeGridConvexity(nodes).valid, true);
  const full = adjustStreamtubeInlets(nodes, body, fractions);
  assert.equal(streamtubeGridConvexity(full.nodes).valid, false);
  const limited = adjustStreamtubeInlets(nodes, body, fractions, { preserveConvexity: true });
  assert.equal(streamtubeGridConvexity(limited.nodes).valid, true);
  assert.ok(limited.convexity.scale > 0 && limited.convexity.scale < 1);
  assert.equal(limited.reparameterization, undefined);
  assert.ok(limited.maxArcErrorAfter < limited.maxArcErrorBefore);
  for (let i = 0; i < 4; i++) assert.deepEqual(limited.nodes[0][i][1], limited.nodes[1][i][0]);
  for (let g = 0; g < 2; g++) for (const i of [0, 3]) assert.deepEqual(limited.nodes[g][i], nodes[g][i]);
  assert.deepEqual(nodes, before);
});

test('DEKINK cannot move an interior streamline through its neighboring bank in a convex curved passage', () => {
  const nodes = [Array.from({ length: 3 }, (_, i) => [1, 2, 3].map(r => ({
    x: r * Math.cos(-i * 2 * Math.PI / 3), y: r * Math.sin(-i * 2 * Math.PI / 3),
  })))];
  const before = structuredClone(nodes);
  assert.equal(streamtubeGridConvexity(nodes).valid, true);
  const full = dekinkStreamtubeInteriors(nodes);
  assert.equal(full.repairs.length, 1);
  assert.equal(streamtubeGridConvexity(full.nodes).valid, false);
  const limited = dekinkStreamtubeInteriors(nodes, { preserveConvexity: true });
  assert.equal(streamtubeGridConvexity(limited.nodes).valid, true);
  assert.ok(limited.convexity.scale > 0 && limited.convexity.scale < 1);
  assert.deepEqual(limited.repairs, full.repairs);
  for (let i = 0; i < 3; i++) for (const j of [0, 2]) assert.deepEqual(limited.nodes[0][i][j], nodes[0][i][j]);
  for (const i of [0, 2]) assert.deepEqual(limited.nodes[0][i], nodes[0][i]);
  assert.deepEqual(nodes, before);
});
