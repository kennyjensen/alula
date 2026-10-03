import test from 'node:test';
import assert from 'node:assert/strict';
import { initializePanelStreamtubeBody, createPanelStreamtubeGrid } from '../src/euler/streamtube-body-initializer.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

test('panel/potential initialization builds admissible one/two-body grids and preserves the supplied topology', () => {
  for (const elements of [1, 2]) {
    const input = intrinsicBodyFixture({ elements, alpha: 2, bodySegments: 16, tubes: 5, tubeGrowth: 3, surfaceSpacing: 'cosine' });
    const original = structuredClone(input), r = initializePanelStreamtubeBody(input);
    assert.deepEqual(input, original, 'initializer mutated its caller');
    assert.equal(r.system.layout.elements, elements);
    assert.ok(r.diagnostics.maxStreamfunctionDrift < 2e-8);
    assert.ok(r.diagnostics.maxSurfaceStreamfunctionDefect < 1e-4);
    assert.ok(r.system.admissible(r.initial));
    assert.ok(r.diagnostics.initialEuler.maxMach < .6);
    assert.ok(r.diagnostics.initialEuler.residualByFamily.inletDensity < 3e-15);
    const { nodes, system } = r;
    for (let b = 0; b < elements; b++) for (let i = 0; i <= system.layout.nx; i++) {
      if (!system.layout.active(b, i)) assert.deepEqual(nodes[b][i].at(-1), nodes[b + 1][i][0]);
      else for (const side of ['lower', 'upper']) {
        const range = r.input.bodies[b], expected = system.curves[b].branch(side, system.fractions[b][side][i - range.leadingIndex], range.stagnationParameter).point;
        assert.deepEqual(side === 'lower' ? nodes[b][i].at(-1) : nodes[b + 1][i][0], expected);
      }
    }
  }
});

test('streamline integration refinement and reference scaling preserve the initialized physical grid', () => {
  const input = intrinsicBodyFixture({ elements: 2, alpha: 2, bodySegments: 16, tubes: 5, tubeGrowth: 3, surfaceSpacing: 'cosine' });
  const a = initializePanelStreamtubeBody(input), b = initializePanelStreamtubeBody(input, { relativeTolerance: 2e-10 });
  let change = 0;
  for (let g = 0; g < a.nodes.length; g++) for (let i = 0; i < a.nodes[g].length; i++) for (let j = 0; j < a.nodes[g][i].length; j++) {
    const p = a.nodes[g][i][j], q = b.nodes[g][i][j]; change = Math.max(change, Math.hypot(p.x - q.x, p.y - q.y));
  }
  assert.ok(change < 2e-7, `Grid changes by ${change} under trace refinement`);
  assert.ok(b.diagnostics.maxStreamfunctionDrift < a.diagnostics.maxStreamfunctionDrift);
  const scale = points => points.map(p => ({ x: 2 * p.x, y: 2 * p.y }));
  const doubled = { ...input, bodies: input.bodies.map(body => ({ ...body, points: scale(body.points), stagnationParameter: 2 * body.stagnationParameter })),
    outerLower: scale(input.outerLower), outerUpper: scale(input.outerUpper), cutPaths: input.cutPaths.map(scale) };
  const c = initializePanelStreamtubeBody(doubled);
  for (let g = 0; g < a.nodes.length; g++) for (let i = 0; i < a.nodes[g].length; i++) for (let j = 0; j < a.nodes[g][i].length; j++) {
    const p = a.nodes[g][i][j], q = c.nodes[g][i][j]; assert.ok(Math.hypot(2 * p.x - q.x, 2 * p.y - q.y) < 1e-7);
  }
});

test('a surface-refined but normally coarse grid is explicitly rejected instead of clipping sonic flux', () => {
  const input = intrinsicBodyFixture({ elements: 2, alpha: 2, bodySegments: 32, tubes: 5, tubeGrowth: 3, surfaceSpacing: 'cosine' });
  assert.throws(() => initializePanelStreamtubeBody(input), /Initial section.*sonic/);
  assert.throws(() => initializePanelStreamtubeBody(input, { seedDistance: 0 }), /controls/);
});

test('dividing and wake traces cover refined cut intervals and are insensitive to a smaller seed distance', () => {
  const input = intrinsicBodyFixture({ elements: 1, alpha: 2, bodySegments: 32, tubes: 7, tubeGrowth: 3,
    surfaceSpacing: 'cosine', cutSpacing: 'surface-matched' }), original = structuredClone(input);
  const a = createPanelStreamtubeGrid(input), b = createPanelStreamtubeGrid(input, { seedDistance: 1e-4 });
  assert.deepEqual(input, original);
  assert.ok(a.diagnostics.profiles[0].seedDistance < .002);
  assert.ok(b.diagnostics.profiles[0].seedDistance < a.diagnostics.profiles[0].seedDistance);
  for (const r of [a, b]) assert.ok(r.diagnostics.maxStreamfunctionDrift < 2e-8);
  let difference = 0;
  for (let g = 0; g < a.nodes.length; g++) for (let i = 0; i < a.nodes[g].length; i++) for (let j = 0; j < a.nodes[g][i].length; j++)
    difference = Math.max(difference, Math.hypot(a.nodes[g][i][j].x - b.nodes[g][i][j].x, a.nodes[g][i][j].y - b.nodes[g][i][j].y));
  assert.ok(difference < 2e-7, `Changing the trace seed moves the grid by ${difference}`);
});
