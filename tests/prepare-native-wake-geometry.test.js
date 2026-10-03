// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { prepareNativeWakeGeometry } from '../scripts/validation/prepare-native-wake-geometry.js';

const mean = (a, b) => ({ x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 });
const center = (n, b, i) => mean(n[b][i].at(-1), n[b + 1][i][0]);
const distance = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);
const close = (a, b, tolerance = 3e-14) => assert.ok(Math.abs(a - b) <= tolerance, `${a} differs from ${b}`);
const dot = (a, b) => a.x * b.x + a.y * b.y;
const separation = (n, b, i) => ({ x: n[b + 1][i][0].x - n[b][i].at(-1).x,
  y: n[b + 1][i][0].y - n[b][i].at(-1).y });
const frame = (n, b, i) => {
  const a = center(n, b, i - 1), z = center(n, b, Math.min(n[0].length - 1, i + 1));
  const length = distance(a, z), t = { x: (z.x - a.x) / length, y: (z.y - a.y) / length };
  return { t, n: { x: -t.y, y: t.x } };
};
const scene = () => {
  const nx = 7, masses = [[1, 3], [2, 1, 4], [1, 5]], trailingIndices = [2, 4];
  const pair = (b, i) => {
    const slope = .04 + .02 * b, den = Math.hypot(1, slope), t = { x: 1 / den, y: slope / den };
    const c = { x: i, y: 2 * (b + 1) + slope * i }, tau = b === 0 ? .12 : -.06, gap = .2 + .02 * b + .001 * i;
    const d = { x: tau * t.x - gap * t.y, y: tau * t.y + gap * t.x };
    return [{ x: c.x - d.x / 2, y: c.y - d.y / 2 }, { x: c.x + d.x / 2, y: c.y + d.y / 2 }];
  };
  const rawNodes = masses.map((mass, g) => Array.from({ length: nx + 1 }, (_, i) => {
    const lower = g === 0 ? { x: i, y: -1 } : pair(g - 1, i)[1];
    const upper = g === 2 ? { x: i, y: 8 } : pair(g, i)[0];
    return Array.from({ length: mass.length + 1 }, (_, j) => ({ x: lower.x + j / mass.length * (upper.x - lower.x),
      y: lower.y + j / mass.length * (upper.y - lower.y) }));
  }));
  const targetNodes = structuredClone(rawNodes);
  trailingIndices.forEach((te, b) => {
    for (let i = 1; i <= te; i++) {
      const weight = i / te, lower = targetNodes[b][i].at(-1), upper = targetNodes[b + 1][i][0];
      lower.x += (.03 + .01 * b) * weight; lower.y += (.01 - .005 * b) * weight;
      upper.x += (.01 - .005 * b) * weight; upper.y -= (.07 - .015 * b) * weight;
    }
  });
  const solidTrailingCenters = trailingIndices.map((te, b) => {
    const c = center(rawNodes, b, te); return { x: c.x - .1, y: c.y - .03 };
  });
  const wakeGaps = trailingIndices.map((te, b) => Array.from({ length: nx - te }, (_, k) => .3 + .025 * b + .015 * k));
  return { rawNodes, targetNodes, trailingIndices, wakeGaps, masses, solidTrailingCenters, lengthScale: 2.5 };
};

test('two staggered wakes preserve tangential offsets while matching gaps and translated centers', () => {
  const input = scene(), before = structuredClone(input), { nodes, diagnostics } = prepareNativeWakeGeometry(input);
  const { rawNodes, targetNodes, trailingIndices, wakeGaps } = input;
  trailingIndices.forEach((te, b) => {
    const oldTe = center(rawNodes, b, te), newTe = center(targetNodes, b, te);
    for (let i = 0; i <= te; i++) {
      assert.deepEqual(nodes[b][i].at(-1), targetNodes[b][i].at(-1));
      assert.deepEqual(nodes[b + 1][i][0], targetNodes[b + 1][i][0]);
    }
    for (let i = te + 1; i < nodes[0].length; i++) {
      const old = center(rawNodes, b, i), next = center(nodes, b, i), f = frame(nodes, b, i), oldF = frame(rawNodes, b, i);
      close(next.x - old.x, newTe.x - oldTe.x); close(next.y - old.y, newTe.y - oldTe.y);
      close(dot(separation(nodes, b, i), f.n), wakeGaps[b][i - te - 1]);
      close(dot(separation(nodes, b, i), f.t), dot(separation(rawNodes, b, i), oldF.t));
      assert.ok(Math.abs(dot(separation(nodes, b, i), f.t)) > .05, 'Retained tangential offset was not averaged away.');
    }
  });
  for (let i = 0; i < nodes[0].length; i++) {
    assert.deepEqual(nodes[0][i][0], rawNodes[0][i][0]); assert.deepEqual(nodes.at(-1)[i].at(-1), rawNodes.at(-1)[i].at(-1));
  }
  assert.deepEqual(input, before);
  assert.equal(diagnostics.operations.boundaryExtensions, 1);
  for (const d of diagnostics.bodies) {
    close(d.maxGapError, 0); close(d.maxTangentialOffsetError, 0); close(d.maxCenterIntervalLengthError, 0);
  }
});

test('shared passage receives both wake-bank increments once in cumulative physical mass fractions', () => {
  const input = scene(), { nodes } = prepareNativeWakeGeometry(input), i = 6, g = 1;
  const old = input.rawNodes[g][i], next = nodes[g][i];
  for (const [j, eta] of [[1, 2 / 7], [2, 3 / 7]]) for (const k of ['x', 'y']) {
    close(next[j][k], old[j][k] + (1 - eta) * (next[0][k] - old[0][k]) + eta * (next.at(-1)[k] - old.at(-1)[k]));
  }
});

test('bare-TE BL arc shift and changed outlet direction remain explicit and unqualified', () => {
  const input = scene(), { nodes, diagnostics } = prepareNativeWakeGeometry(input);
  diagnostics.bodies.forEach((d, b) => {
    const i = input.trailingIndices[b] + 1, bare = input.solidTrailingCenters[b];
    const shift = distance(center(nodes, b, i), bare) - distance(center(input.rawNodes, b, i), bare);
    close(d.bareTeArc.physicalShift, shift); close(d.bareTeArc.kernelShift, shift / input.lengthScale);
    assert.ok(Math.abs(shift) > .001, 'Manufactured translation must exercise the fixed bare-TE anchor.');
    assert.equal(d.outlet.qualified, false);
    assert.ok(d.outlet.banks.some(v => Math.abs(v.directionChangeRadians) > .001));
  });
  assert.equal(diagnostics.converged, false); assert.equal(diagnostics.geometryQualified, false);
  assert.equal(diagnostics.gasQualified, false); assert.equal(diagnostics.physicalAcceptance, false);
});

test('unchanged walls and existing gaps are an exact identity with fresh output storage', () => {
  const input = scene(); input.targetNodes = input.rawNodes;
  input.wakeGaps = input.trailingIndices.map((te, b) => Array.from({ length: input.rawNodes[0].length - te - 1 }, (_, k) => {
    const i = te + 1 + k; return dot(separation(input.rawNodes, b, i), frame(input.rawNodes, b, i).n);
  }));
  const before = structuredClone(input), { nodes } = prepareNativeWakeGeometry(input);
  assert.deepEqual(nodes, input.rawNodes); assert.deepEqual(input, before);
  assert.notEqual(nodes[0][0][0], input.rawNodes[0][0][0]);
});

test('rotation, translation and physical scaling preserve the construction and BL-unit arc shift', () => {
  const original = scene(), result = prepareNativeWakeGeometry(original), angle = .61, scale = 3.7;
  const transform = p => ({ x: 11 + scale * (Math.cos(angle) * p.x - Math.sin(angle) * p.y),
    y: -7 + scale * (Math.sin(angle) * p.x + Math.cos(angle) * p.y) });
  const transformed = { ...original, rawNodes: original.rawNodes.map(g => g.map(r => r.map(transform))),
    targetNodes: original.targetNodes.map(g => g.map(r => r.map(transform))),
    solidTrailingCenters: original.solidTrailingCenters.map(transform),
    wakeGaps: original.wakeGaps.map(r => r.map(g => scale * g)), masses: original.masses.map(r => r.map(m => scale * m)),
    lengthScale: scale * original.lengthScale };
  const next = prepareNativeWakeGeometry(transformed);
  result.nodes.forEach((g, gi) => g.forEach((r, i) => r.forEach((p, j) => close(distance(transform(p), next.nodes[gi][i][j]), 0, 1e-13))));
  next.diagnostics.bodies.forEach((d, b) => close(d.bareTeArc.kernelShift, result.diagnostics.bodies[b].bareTeArc.kernelShift, 1e-14));
});

test('tiny exterior target chart roundoff is explicit and never displaces the raw exterior', () => {
  const input = scene(); input.targetNodes[0][2][0].x += 2 * Number.EPSILON;
  input.targetNodes[2][4].at(-1).y += 8 * Number.EPSILON;
  const before = structuredClone(input), result = prepareNativeWakeGeometry(input);
  assert.ok(result.diagnostics.exteriorBanks.targetDeparture > 0);
  assert.ok(result.diagnostics.exteriorBanks.targetDeparture <= result.diagnostics.exteriorBanks.coordinateRoundoffTolerance);
  for (let i = 0; i < input.rawNodes[0].length; i++) {
    assert.deepEqual(result.nodes[0][i][0], input.rawNodes[0][i][0]);
    assert.deepEqual(result.nodes[2][i].at(-1), input.rawNodes[2][i].at(-1));
  }
  assert.deepEqual(input, before);
});

test('rejects invalid gaps, mass, indices, changed exterior banks and degenerate wake centers without mutation', () => {
  const changes = [
    i => { i.wakeGaps[0][0] = -1; }, i => { i.wakeGaps[0][0] = NaN; },
    i => { i.masses[1][0] = 0; }, i => { i.trailingIndices[0] = 7; },
    i => { i.targetNodes[0][3][0].x += .1; }, i => { i.targetNodes[1][0][1].y = Infinity; },
    i => { i.solidTrailingCenters.pop(); },
    i => { const c = center(i.rawNodes, 0, 2), d = separation(i.rawNodes, 0, 4);
      i.rawNodes[0][4][2] = { x: c.x - d.x / 2, y: c.y - d.y / 2 };
      i.rawNodes[1][4][0] = { x: c.x + d.x / 2, y: c.y + d.y / 2 }; },
  ];
  for (const change of changes) {
    const input = scene(); change(input); const before = structuredClone(input);
    assert.throws(() => prepareNativeWakeGeometry(input)); assert.deepEqual(input, before);
  }
});
