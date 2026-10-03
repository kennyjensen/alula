import test from 'node:test';
import assert from 'node:assert/strict';
import { initializePanelStreamtubeBody, createPanelStreamtubeGrid } from '../src/euler/streamtube-body-initializer.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const fixture = () => ({ ...intrinsicBodyFixture({ elements: 2, alpha: 2, bodySegments: 16, tubes: 5, tubeGrowth: 3,
  surfaceSpacing: 'cosine' }), streamwiseMode: 'isentropic' });
const controls = { crosslinePlacement: 'potential', outerCrosslineSpread: .15 };
const distance = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);

test('potential spacing uses independent wall fractions while preserving exact body endpoints and connected cuts', () => {
  const input = fixture(), before = structuredClone(input), r = initializePanelStreamtubeBody(input, controls);
  assert.deepEqual(input, before); assert.equal(r.system.admissible(r.initial), true);
  assert.ok(r.diagnostics.maxStreamfunctionDrift < 2e-8);
  const mesh = streamtubeMeshSnapshot(r); assert.equal(mesh.quality.valid, true);
  const edges = new Set();
  for (const cell of mesh.cells) for (let i = 0; i < 4; i++) {
    const [a, b] = [cell[i], cell[(i + 1) % 4]].sort((a, b) => a - b); edges.add(`${a}:${b}`);
  }
  assert.equal(mesh.vertices.length - edges.size + mesh.cells.length, -1);
  for (let b = 0; b < 2; b++) {
    const body = r.input.bodies[b], curve = r.system.curves[b], fractions = body.surfaceFractions;
    assert.ok(fractions.upper.some((f, i) => Math.abs(f - fractions.lower[i]) > 1e-3));
    for (let i = 0; i <= r.system.layout.nx; i++) {
      const lower = r.nodes[b][i].at(-1), upper = r.nodes[b + 1][i][0];
      if (!r.system.layout.active(b, i)) assert.deepEqual(lower, upper);
      else {
        assert.ok(distance(lower, curve.branch('lower', fractions.lower[i - body.leadingIndex], body.stagnationParameter).point) < 1e-14);
        assert.ok(distance(upper, curve.branch('upper', fractions.upper[i - body.leadingIndex], body.stagnationParameter).point) < 1e-14);
      }
    }
    assert.deepEqual(r.nodes[b][body.trailingIndex].at(-1), input.bodies[b].points[0]);
    assert.deepEqual(r.nodes[b][body.leadingIndex].at(-1), r.nodes[b + 1][body.leadingIndex][0]);
  }
});

test('matched potential grids preserve units and converge under independent tracing and seed refinement', () => {
  const input = fixture(), a = createPanelStreamtubeGrid(input, controls);
  const b = createPanelStreamtubeGrid(input, { ...controls, relativeTolerance: 2e-10, seedDistance: 1e-4 });
  const scale = points => points.map(p => ({ x: 2 * p.x, y: 2 * p.y }));
  const doubled = { ...input, bodies: input.bodies.map(body => ({ ...body, points: scale(body.points), stagnationParameter: 2 * body.stagnationParameter })),
    outerLower: scale(input.outerLower), outerUpper: scale(input.outerUpper), cutPaths: input.cutPaths.map(scale) };
  const c = createPanelStreamtubeGrid(doubled, controls);
  assert.equal(a.system.layout.nx, b.system.layout.nx); assert.equal(a.system.layout.nx, c.system.layout.nx);
  let difference = 0;
  for (let g = 0; g < a.nodes.length; g++) for (let i = 0; i < a.nodes[g].length; i++) for (let j = 0; j < a.nodes[g][i].length; j++) {
    const p = a.nodes[g][i][j]; difference = Math.max(difference, distance(p, b.nodes[g][i][j]));
    assert.ok(distance({ x: 2 * p.x, y: 2 * p.y }, c.nodes[g][i][j]) < 1e-7);
  }
  assert.ok(difference < 2e-7, `Trace/seed refinement moves the grid by ${difference}`);
  assert.throws(() => createPanelStreamtubeGrid(input, { crosslinePlacement: 'unknown' }), /placement/);
  assert.throws(() => createPanelStreamtubeGrid(input, { ...controls, outerCrosslineSpread: 1.1 }), /spreading/);
  assert.throws(() => createPanelStreamtubeGrid(input, { ...controls, maxCrosslineIntervals: 4 }), /budget/);
});

test('refined potential cuts extend actual streamlines toward edges without extrapolating the grid', () => {
  const input = { ...intrinsicBodyFixture({ elements: 1, alpha: 2, bodySegments: 64, tubes: 9, tubeGrowth: 3,
    surfaceSpacing: 'cosine' }), streamwiseMode: 'isentropic' };
  const a = createPanelStreamtubeGrid(input, controls), b = createPanelStreamtubeGrid(input, { ...controls, seedDistance: 1e-5 });
  assert.ok(a.diagnostics.cutTraceExtensions > 0);
  assert.equal(a.system.layout.nx, b.system.layout.nx);
  assert.equal(streamtubeMeshSnapshot(a).quality.valid, true);
  let difference = 0;
  for (let g = 0; g < a.nodes.length; g++) for (let i = 0; i < a.nodes[g].length; i++) for (let j = 0; j < a.nodes[g][i].length; j++)
    difference = Math.max(difference, distance(a.nodes[g][i][j], b.nodes[g][i][j]));
  assert.ok(difference < 2e-7, `Extended cuts depend on the trace seed by ${difference}`);
});
