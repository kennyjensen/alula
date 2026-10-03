// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem as create } from '../src/euler/streamtube-body.js';
import { finiteBaseBodyFixture } from './fixtures/finite-base-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
const fixture = () => {
  const input = finiteBaseBodyFixture({ bodySegments: 4, tubes: 2 });
  input.bodies[0].points[0].x += .0002; input.bodies[0].points.at(-1).x = input.bodies[0].points[0].x;
  input.bodies[0].wakeTangentialReference = 'material-te'; return input;
};

test('material-TE geometry and the complete Euler Jacobian include the tangential projection derivative', () => {
  const system = create(fixture()), state = system.initial, h = 3e-7;
  const d = Float64Array.from(state, (_, i) => Math.sin(.41 * (i + 1)));
  const p = state.map((v, i) => v + h * d[i]), m = state.map((v, i) => v - h * d[i]);
  const a = system.decode(p).nodes, z = system.decode(m).nodes, maps = system.geometryDerivatives(state);
  let geometryError = 0;
  for (let g = 0; g < a.length; g++) for (let i = 0; i < a[g].length; i++) for (let j = 0; j < a[g][i].length; j++) for (const k of ['x', 'y']) {
    const exact = [...maps[g][i][j]].reduce((sum, [col, v]) => sum + v[k] * d[col], 0);
    geometryError = Math.max(geometryError, Math.abs(exact - (a[g][i][j][k] - z[g][i][j][k]) / (2 * h)));
  }
  assert.ok(geometryError < 2e-8, `${geometryError}`);
  const J = system.jacobian(state), rp = system.residual(p), rm = system.residual(m);
  let rowError = 0;
  for (let row = 0; row < state.length; row++) {
    let exact = 0; for (let col = 0; col < state.length; col++) exact += J[row * state.length + col] * d[col];
    rowError = Math.max(rowError, Math.abs(exact - (rp[row] - rm[row]) / (2 * h)) / Math.max(1, Math.abs(exact)));
  }
  assert.ok(rowError < 4e-6, `${rowError}`);
});

test('the declared centerline chart rejects omitted markers and transfers displaced TE banks exactly into independent mode', () => {
  const input = fixture(), base = create(input); input.displacement = structuredClone(base.displacement);
  input.displacement.surfaces[0].upper.fill(.0003); input.displacement.surfaces[0].lower.fill(.0001);
  input.displacement.surfaces[0].upper[0] = input.displacement.surfaces[0].lower[0] = 0;
  input.displacement.wakes[0] = input.displacement.wakes[0].map(v => v + .0004);
  const source = create(input), nodes = source.decode(source.initial).nodes;
  const omitted = structuredClone(input); delete omitted.bodies[0].wakeTangentialReference;
  const plain = create(omitted);
  assert.throws(() => plain.adoptGeometry(plain.initial, nodes), /Restart wall or wake banks/);
  assert.throws(() => source.adoptGeometry(source.initial, plain.decode(plain.initial).nodes), /Restart wall or wake banks/);
  omitted.bodies[0].wakeTangentialReference = 'unknown'; assert.throws(() => create(omitted), /Unknown tangential wake reference/);
  const independentInput = { ...JSON.parse(JSON.stringify(input)), wakeGeometry: 'independent-banks', wakeDisplacementMotion: 'te-center', wakeOutlet: 'banks' };
  const target = create(independentInput), state = target.adoptGeometry(target.initial, nodes);
  assert.deepEqual(target.decode(state).nodes, nodes);
  const restored = create(JSON.parse(JSON.stringify(independentInput)));
  assert.deepEqual(restored.decode(restored.adoptGeometry(restored.initial, nodes)).nodes, nodes);
  const maps = target.geometryDerivatives(state, { includeDisplacement: true }), h = 2e-7;
  for (const side of ['upper', 'lower']) {
    const param = target.displacementParameters.findIndex(p => p.kind === 'surface' && p.side === side
      && p.index === input.bodies[0].trailingIndex - input.bodies[0].leadingIndex);
    const varied = [1, -1].map(sign => {
      const d = structuredClone(input.displacement); d.surfaces[0][side][d.surfaces[0][side].length - 1] += sign * h;
      target.setDisplacement(d); return target.decode(state).nodes;
    });
    let max = 0;
    for (let g = 0; g < nodes.length; g++) for (let i = 0; i < nodes[g].length; i++) for (let j = 0; j < nodes[g][i].length; j++) for (const k of ['x', 'y'])
      max = Math.max(max, Math.abs((maps[g][i][j].get(target.layout.n + param)?.[k] ?? 0)
        - (varied[0][g][i][j][k] - varied[1][g][i][j][k]) / (2 * h)));
    assert.ok(max < 2e-8, `${side}: ${max}`);
  }
  target.setDisplacement(input.displacement); assert.deepEqual(target.decode(state).nodes, nodes);
});

test('a zero solid TE vector reduces exactly to the original displacement geometry', () => {
  const input = intrinsicBodyFixture({ bodySegments: 4, tubes: 2 });
  input.displacement = { surfaces: input.bodies.map(b => ({ upper: Array(b.trailingIndex - b.leadingIndex + 1).fill(0),
    lower: Array(b.trailingIndex - b.leadingIndex + 1).fill(0) })), wakes: input.bodies.map(b => Array(input.outerLower.length - 1 - b.trailingIndex).fill(0)) };
  const before = create(input); input.bodies[0].wakeTangentialReference = 'material-te'; const after = create(input);
  assert.deepEqual(after.decode(after.initial).nodes, before.decode(before.initial).nodes);
  assert.deepEqual(after.geometryDerivatives(after.initial), before.geometryDerivatives(before.initial));
});
