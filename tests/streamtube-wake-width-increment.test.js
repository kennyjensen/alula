// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { incrementIndependentWakeWidths } from '../src/euler/streamtube-wake-geometry.js';
import { streamtubeWakeGap } from '../src/euler/streamtube-wake-geometry.js';
import { extendWarmBoundaryIncrements } from '../src/euler/streamtube-displacement.js';

const close = (a, b, tolerance = 3e-12) => assert.ok(Math.abs(a - b) <= tolerance, `${a} != ${b}`);
const center = (nodes, body, i) => ({ x: .5 * (nodes[body][i].at(-1).x + nodes[body + 1][i][0].x),
  y: .5 * (nodes[body][i].at(-1).y + nodes[body + 1][i][0].y) });
const gap = (nodes, b, i, nx) => streamtubeWakeGap(...[nodes[b].map(row => row.at(-1)), nodes[b + 1].map(row => row[0])]
  .map(bank => [bank[i - 1], bank[i], bank[Math.min(nx, i + 1)]]));

function fixture() {
  const nx = 7, tubes = [2, 3, 2], bodies = [{ trailingIndex: 2 }, { trailingIndex: 3 }];
  const banks = bodies.map((body, b) => Array.from({ length: nx + 1 }, (_, i) => {
    const c = { x: .19 * i + .009 * Math.sin(i), y: b + .013 * i * i },
      length = Math.hypot(1, .026 * i), t = { x: 1 / length, y: .026 * i / length }, n = { x: -t.y, y: t.x },
      width = .025 + .001 * i, offset = .022 - .003 * i;
    return { lower: { x: c.x - .5 * (width * n.x + offset * t.x), y: c.y - .5 * (width * n.y + offset * t.y) },
      upper: { x: c.x + .5 * (width * n.x + offset * t.x), y: c.y + .5 * (width * n.y + offset * t.y) } };
  }));
  const nodes = tubes.map((nt, g) => Array.from({ length: nx + 1 }, (_, i) => {
    const lo = g ? banks[g - 1][i].upper : { x: .19 * i, y: -1 }, hi = g < 2 ? banks[g][i].lower : { x: .19 * i, y: 3 };
    return Array.from({ length: nt + 1 }, (_, j) => ({ x: (1 - j / nt) * lo.x + j / nt * hi.x,
      y: (1 - j / nt) * lo.y + j / nt * hi.y, tag: `${g}/${i}/${j}` }));
  }));
  const beforeWidths = bodies.map((b, k) => Array.from({ length: nx - b.trailingIndex }, (_, i) => .025 + .001 * i + .003 * k));
  const afterWidths = beforeWidths.map((row, b) => row.map((d, i) => d + (i % 2 ? -.002 : .003) * (b + 1)));
  return { layout: { nx, tubes, bodies, independentWakeBanks: true }, nodes, beforeWidths, afterWidths };
}

test('normal width increments preserve centers, tangent offsets and existing gap residual on both curved wakes', () => {
  const input = fixture(), old = structuredClone(input), result = incrementIndependentWakeWidths(input), { nx, tubes, bodies } = input.layout;
  assert.deepEqual(input, old); assert.equal(result.diagnostics.changedNodes, 18);
  for (const [body, spec] of bodies.entries()) for (let i = spec.trailingIndex + 1; i <= nx; i++) {
    const a = gap(input.nodes, body, i, nx), b = gap(result.nodes, body, i, nx), k = i - spec.trailingIndex - 1;
    close(a.gap - input.beforeWidths[body][k], b.gap - input.afterWidths[body][k]);
    close(a.tangentialOffset, b.tangentialOffset);
    for (const axis of ['x', 'y']) close(center(input.nodes, body, i)[axis], center(result.nodes, body, i)[axis]);
  }
  for (let g = 0; g < tubes.length; g++) for (let i = 0; i <= nx; i++) for (let j = 0; j <= tubes[g]; j++) {
    const b = j === 0 && g > 0 ? g - 1 : j === tubes[g] && g < bodies.length ? g : null;
    if (b === null || i <= bodies[b].trailingIndex) assert.deepEqual(result.nodes[g][i][j], input.nodes[g][i][j]);
  }
  assert.equal(result.diagnostics.zeroesGapResidual, false); assert.equal(result.diagnostics.geometryAdmissibilityChecked, false);
});

test('zero changes return the original coordinates exactly without reapplying absolute wake thickness', () => {
  const input = fixture(); input.afterWidths = structuredClone(input.beforeWidths);
  input.nodes[0][0][0].x = -0;
  const old = structuredClone(input), result = incrementIndependentWakeWidths(input);
  assert.equal(result.nodes, input.nodes); assert.deepEqual(input, old);
  assert.equal(result.diagnostics.active, false); assert.equal(result.diagnostics.changedNodes, 0);
  assert.ok(Object.is(result.nodes[0][0][0].x, -0));
});

test('physical total-width updates rotate and scale and do not add the finite-base width twice', () => {
  const original = fixture(), base = incrementIndependentWakeWidths(original);
  for (const angle of [.7, 2.1]) for (const scale of [.03, 6]) {
    const rotate = p => ({ ...p, x: 2 + scale * (Math.cos(angle) * p.x - Math.sin(angle) * p.y),
      y: -3 + scale * (Math.sin(angle) * p.x + Math.cos(angle) * p.y) });
    const input = { ...original, nodes: original.nodes.map(g => g.map(r => r.map(rotate))),
      beforeWidths: original.beforeWidths.map(r => r.map(d => scale * d)), afterWidths: original.afterWidths.map(r => r.map(d => scale * d)) };
    const result = incrementIndependentWakeWidths(input);
    result.nodes.forEach((g, a) => g.forEach((r, i) => r.forEach((p, j) => {
      const expected = rotate(base.nodes[a][i][j]); close(p.x, expected.x); close(p.y, expected.y);
    })));
  }
  // Adding an unchanged dead-air component to both TOTAL widths cancels in
  // the increment. The caller owns any changed arc-dependent base-gap model.
  const shifted = { ...original, beforeWidths: original.beforeWidths.map(r => r.map(d => d + .004)),
    afterWidths: original.afterWidths.map(r => r.map(d => d + .004)) };
  const result = incrementIndependentWakeWidths(shifted);
  result.nodes.forEach((g, a) => g.forEach((r, i) => r.forEach((p, j) => {
    close(p.x, base.nodes[a][i][j].x); close(p.y, base.nodes[a][i][j].y);
  })));
});

test('boundary increments extend through a shared passage using physical masses and retain interior detail', () => {
  const input = fixture(), target = incrementIndependentWakeWidths(input).nodes, masses = [[1, 9], [1, 2, 7], [3, 1]];
  const nodes = extendWarmBoundaryIncrements({ sourceNodes: input.nodes, targetNodes: target, masses });
  for (let g = 0; g < nodes.length; g++) for (let i = 0; i <= input.layout.nx; i++) {
    let eta = 0;
    for (let j = 0; j <= masses[g].length; j++) {
      if (j) eta += masses[g][j - 1] / masses[g].reduce((a, b) => a + b, 0);
      for (const axis of ['x', 'y']) {
        const deltaLo = target[g][i][0][axis] - input.nodes[g][i][0][axis],
          deltaHi = target[g][i].at(-1)[axis] - input.nodes[g][i].at(-1)[axis];
        close(nodes[g][i][j][axis], input.nodes[g][i][j][axis] + (1 - eta) * deltaLo + eta * deltaHi);
      }
      assert.equal(nodes[g][i][j].tag, input.nodes[g][i][j].tag);
    }
  }
});

test('reject wrong layouts, nonfinite or negative widths, bad dimensions and degenerate wake normals', () => {
  for (const mutate of [a => { a.layout.independentWakeBanks = false; }, a => { a.afterWidths[0][0] = NaN; },
    a => { a.afterWidths[0][0] = -.001; }, a => { a.beforeWidths[0].pop(); }, a => { a.nodes[0][0].pop(); },
    a => { for (const i of [2, 4]) { a.nodes[0][i][2] = { x: 1, y: 0 }; a.nodes[1][i][0] = { x: 1, y: 1 }; } }]) {
    const input = fixture(); mutate(input); const old = structuredClone(input);
    assert.throws(() => incrementIndependentWakeWidths(input)); assert.deepEqual(input, old);
  }
});
