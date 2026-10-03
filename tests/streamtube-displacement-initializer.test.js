import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { extendStreamtubeDisplacement } from '../src/euler/streamtube-displacement.js';
import { evaluateStreamtubeCell } from '../src/euler/streamtube-cell.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const fixture = amplitude => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  input.displacement = { surfaces: input.bodies.map(b => ({
    upper: Array.from({ length: b.trailingIndex - b.leadingIndex + 1 }, (_, i) => amplitude * (1 + .2 * i)),
    lower: Array.from({ length: b.trailingIndex - b.leadingIndex + 1 }, (_, i) => amplitude * (1 + .1 * i)) })),
    wakes: input.bodies.map(b => Array(input.outerLower.length - 1 - b.trailingIndex).fill(3 * amplitude)) };
  return input;
};

test('displacement-grid extension preserves zero offsets, prescribed walls/wake gaps and fixed farfield/inlet boundaries', () => {
  for (const amplitude of [0, .0003]) {
    const input = fixture(amplitude), system = createStreamtubeBodySystem(input), x = system.initial;
    const before = system.decode(x), chart = system.geometryChart(), candidate = extendStreamtubeDisplacement(system, x);
    assert.deepEqual(system.decode(x), before); assert.deepEqual(system.geometryChart(), chart);
    if (!amplitude) assert.deepEqual(candidate, before.nodes);
    for (let g = 0; g < candidate.length; g++) assert.deepEqual(candidate[g][0], before.nodes[g][0]);
    for (let i = 0; i <= system.layout.nx; i++) {
      assert.deepEqual(candidate[0][i][0], before.nodes[0][i][0]);
      assert.deepEqual(candidate.at(-1)[i].at(-1), before.nodes.at(-1)[i].at(-1));
    }
    // Adoption independently validates the prescribed wall normals and
    // every wake bank against the reconstructed centerline.
    const check = createStreamtubeBodySystem(input), adopted = check.adoptGeometry(x, candidate), after = check.decode(adopted);
    after.nodes.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => {
      assert.ok(Math.hypot(p.x - candidate[g][i][j].x, p.y - candidate[g][i][j].y) < 1e-13);
    })));
    system.layout.bodies.forEach((body, b) => {
      for (const side of ['upper', 'lower']) for (let i = body.leadingIndex; i <= body.trailingIndex; i++) {
        const g = side === 'lower' ? b : b + 1, j = side === 'lower' ? system.layout.tubes[b] : 0;
        assert.deepEqual(candidate[g][i][j], before.nodes[g][i][j]);
      }
      for (let i = body.trailingIndex + 1; i <= system.layout.nx; i++) {
        const a = candidate[b][i].at(-1), c = candidate[b + 1][i][0];
        assert.ok(Math.abs(Math.hypot(c.x - a.x, c.y - a.y) - input.displacement.wakes[b][i - body.trailingIndex - 1]) < 1e-13);
      }
    });
  }
});

test('the saved wide outer tube pressure failure responds to local geometry curvature, without changing gas state or pressure guards', () => {
  const { cell } = JSON.parse(readFileSync(new URL('./fixtures/displaced-outer-te-cell.json', import.meta.url)));
  assert.throws(() => evaluateStreamtubeCell(cell), /interface pressure/);
  const straightened = structuredClone(cell), [a, p, b] = straightened.upper;
  const t = (p.x - a.x) / (b.x - a.x); p.y = (1 - t) * a.y + t * b.y;
  const checked = evaluateStreamtubeCell(straightened);
  assert.ok(checked.interfacePressure.lower > 0 && checked.interfacePressure.upper > 0);
  assert.ok(checked.states.every(s => s.machSquared < 1));
  assert.deepEqual(straightened.densities, cell.densities); assert.equal(straightened.massFlow, cell.massFlow);
  assert.ok(Math.abs(checked.streamwiseResidual) > 1e-3, 'a locally admissible cell is not a converged flow');
});
