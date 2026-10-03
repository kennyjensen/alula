// SPDX-License-Identifier: GPL-2.0-or-later
// Tiny intrinsic charts only. No panel solve, SLOR, gas evaluation or Newton.
import test from 'node:test';
import assert from 'node:assert/strict';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { refineStreamtubeBody } from '../src/euler/streamtube-refinement.js';

function fixture(elements = 1) {
  const input = intrinsicBodyFixture({ bodySegments: 4, tubes: 2, contourPanels: 40, elements });
  input.bodies.forEach((body, b) => { body.element = elements - b - 1; });
  const system = createStreamtubeBodySystem(input);
  system.evaluate = () => { throw new Error('Metadata tests must not evaluate gas/flow.'); };
  return { input, system, nodes: system.decode(system.initial).nodes };
}

test('actual grid count metadata overrides stale requested or parent diagnostics', () => {
  const { system, nodes } = fixture();
  const metadata = streamtubeMeshSnapshot({ system, nodes, diagnostics: {
    groups: 99, streamwiseSegments: 1, tubes: [99], surfaceIntervals: [{ element: 4, intervals: 999 }],
    gridSpacing: { surfaceIntervals: 128 },
  } });
  const intervalCount = system.layout.bodies[0].trailingIndex - system.layout.bodies[0].leadingIndex;
  assert.equal(metadata.initialization.groups, nodes.length);
  assert.equal(metadata.initialization.streamwiseSegments, nodes[0].length - 1);
  assert.deepEqual(metadata.initialization.tubes, nodes.map(g => g[0].length - 1));
  assert.deepEqual(metadata.initialization.surfaceIntervals, [{ element: 0, intervals: intervalCount }]);
  assert.equal(metadata.initialization.gridSpacing.surfaceIntervals, 128, 'requested/parent diagnostic stays separate');
  assert.equal(metadata.cells.length, (nodes[0].length - 1) * nodes.reduce((n, g) => n + g[0].length - 1, 0));
  assert.equal(metadata.initialization.flowSolved, false);
});

test('multiple element counts keep the physical solver order and explicit element identities', () => {
  const { system, nodes } = fixture(2), snapshot = streamtubeMeshSnapshot({ system, nodes });
  assert.deepEqual(snapshot.initialization.elementOrder, [1, 0]);
  assert.deepEqual(snapshot.initialization.surfaceIntervals.map(b => b.element), [1, 0]);
  for (const [b, body] of system.layout.bodies.entries()) {
    const counts = ['lower', 'upper'].map(side => nodes[side === 'lower' ? b : b + 1]
      .slice(body.leadingIndex, body.trailingIndex + 1).length);
    assert.deepEqual(counts, [snapshot.initialization.surfaceIntervals[b].intervals + 1,
      snapshot.initialization.surfaceIntervals[b].intervals + 1]);
  }
  assert.deepEqual(snapshot.initialization.tubes.map(n => n + 1), nodes.map(g => g[0].length));
});

test('refined metadata reports the child chart without a flow evaluation', () => {
  const { input, system } = fixture();
  const refined = refineStreamtubeBody(input, system, { initial: system.initial, streamwiseFactor: 2,
    normalSubdivisions: system.layout.tubes.map(n => Array(n).fill(2)), initializeFlow: false });
  const snapshot = streamtubeMeshSnapshot({ ...refined, nodes: refined.initialEuler.nodes,
    diagnostics: { tubes: system.layout.tubes, surfaceIntervals: [{ element: 0, intervals: 4 }] } });
  assert.equal(snapshot.initialization.streamwiseSegments, 2 * system.layout.nx);
  assert.deepEqual(snapshot.initialization.tubes, system.layout.tubes.map(n => 2 * n));
  assert.deepEqual(snapshot.initialization.surfaceIntervals, system.layout.bodies.map(body => ({
    element: body.element, intervals: 2 * (body.trailingIndex - body.leadingIndex),
  })));
});
