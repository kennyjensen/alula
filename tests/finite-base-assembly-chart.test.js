// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { buildReliabilityCase } from '../scripts/validation/solver-reliability-cases.js';
import { solveInviscid, velocityAt } from '../src/inviscid/linear-vortex.js';
import { createDisplacementOperator } from '../src/inviscid/displacement.js';
import { streamfunctionAt, streamfunctionBranchIncrement } from '../src/inviscid/streamfunction.js';
import { createContourTopology } from '../src/geometry/contour-topology.js';
import { assemblyFiniteBaseSourceChart } from '../src/inviscid/finite-base-source-chart.js';
import { finiteBaseSourceChart, xfoilSurfaceDerivatives } from '../src/inviscid/finite-base-influence.js';

const input = buildReliabilityCase({ preset: 'three', mode: 'streamtube-bl', changes: {} }).caseData;
const panel = solveInviscid({ ...input, mach: 0, boundaryCondition: 'streamfunction' });

test('finite-base source chart avoids downstream solids consistently in both panel operators', () => {
  const operator = createDisplacementOperator({ elements: input.elements, alpha: input.alpha, wakeCount: 24, wakeInitialization: 'inviscid' });
  const velocity = operator.velocityField(new Float64Array(operator.total));
  for (const p of [{ x: -.1, y: .1 }, { x: .3, y: .3 }, { x: 1.4, y: -.3 }]) {
    const v = velocityAt(p, panel.field), w = velocity(p);
    assert.ok(Math.hypot(v.u - w.u, v.v - w.v) < 1e-10);
  }
  for (const element of input.elements) {
    const topology = createContourTopology(element.points, element);
    const psi = topology.surface.points.map(p => streamfunctionAt(p, panel.field));
    assert.ok(Math.max(...psi) - Math.min(...psi) < 1e-11);
    // Mid-panel values are approximate, but there must be no source-sized
    // inter-node jump where the original slat ray crossed the main nose.
    for (let i = 1; i < topology.surface.points.length; i++) {
      const a = topology.surface.points[i - 1], b = topology.surface.points[i];
      const mid = streamfunctionAt({ x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 }, panel.field);
      assert.ok(Math.abs(mid - psi[i]) < 1e-4);
    }
  }
});

test('source chart selection preserves unobstructed charts and commutes with rigid transformations', () => {
  const t = createContourTopology(input.elements[0].points, input.elements[0]);
  const d = xfoilSurfaceDerivatives(t.surface.points).derivatives;
  const tx = -d[0].x + d.at(-1).x, ty = -d[0].y + d.at(-1).y, length = Math.hypot(tx, ty);
  const preferred = { x: tx / length, y: ty / length }, others = input.elements.slice(1).map(e => e.points);
  assert.deepEqual(assemblyFiniteBaseSourceChart(t, preferred, []), finiteBaseSourceChart(t, preferred));
  const chart = assemblyFiniteBaseSourceChart(t, preferred, others);
  const angle = .73, c = Math.cos(angle), s = Math.sin(angle);
  const rotate = p => ({ x: c * p.x - s * p.y, y: s * p.x + c * p.y });
  const transform = p => { const r = rotate(p); return { x: 3 * r.x + 2, y: 3 * r.y - 4 }; };
  const transformed = createContourTopology(input.elements[0].points.map(transform), input.elements[0]);
  const result = assemblyFiniteBaseSourceChart(transformed, rotate(preferred), others.map(p => p.map(transform)));
  const expected = rotate(chart.direction), origin = transform(chart.origin);
  assert.ok(Math.hypot(result.direction.x - expected.x, result.direction.y - expected.y) < 1e-12);
  assert.ok(Math.hypot(result.origin.x - origin.x, result.origin.y - origin.y) < 1e-12);
});

test('source-sheet transport uses exact signed flux and is independent of path subdivision', () => {
  const field = { basePanels: [{ cutOrigin: { x: 0, y: 0 }, cutDirection: { x: 1, y: 0 }, sourceStrength: 2, length: .3 }] };
  const path = [{ x: 1, y: -1 }, { x: 1, y: 1 }];
  assert.equal(streamfunctionBranchIncrement(path, field), -.6);
  assert.equal(streamfunctionBranchIncrement(path.toReversed(), field), .6);
  assert.equal(streamfunctionBranchIncrement([path[0], { x: 1, y: 0 }, path[1]], field), -.6);
  assert.equal(streamfunctionBranchIncrement([{ x: -1, y: -1 }, { x: -1, y: 1 }], field), 0);
  assert.equal(streamfunctionBranchIncrement([...path, path[0]], field), 0);
});
