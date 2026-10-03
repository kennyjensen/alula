import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { streamtubeMeshQuality, streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { streamtubeMeshConnectivity } from '../src/geometry/streamtube-mesh-connectivity.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { twoActiveFiniteBaseWakes } from './fixtures/two-active-finite-base-wakes.js';

// Retain the former published-cell scan as an independent equivalence oracle.
function publishedQuality(layout, nodes) {
  const { vertices, cells } = streamtubeMeshConnectivity(layout, nodes);
  let minArea = Infinity, minCornerSine = Infinity; const invalidCells = [];
  cells.forEach((cell, index) => {
    const p = cell.map(id => vertices[id]); let twiceArea = 0, valid = true;
    for (let j = 0; j < 4; j++) {
      const a = p[j], b = p[(j + 1) % 4], c = p[(j + 2) % 4];
      twiceArea += (a.x - p[0].x) * (b.y - p[0].y) - (a.y - p[0].y) * (b.x - p[0].x);
      const ux = b.x - a.x, uy = b.y - a.y, vx = c.x - b.x, vy = c.y - b.y;
      const sine = (ux * vy - uy * vx) / (Math.hypot(ux, uy) * Math.hypot(vx, vy));
      minCornerSine = Math.min(minCornerSine, sine); if (!(sine > 1e-12)) valid = false;
    }
    minArea = Math.min(minArea, .5 * twiceArea); if (!valid || !(twiceArea > 0)) invalidCells.push(index);
  });
  return { valid: !invalidCells.length, invalidCells, minArea, minCornerSine };
}

function fixtures() {
  const result = [1, 2].map(elements => {
    const system = createStreamtubeBodySystem(intrinsicBodyFixture({ elements, bodySegments: 4, contourPanels: 40 }));
    return { system, nodes: system.decode(system.initial).nodes };
  });
  const { system, x } = twoActiveFiniteBaseWakes();
  result.push({ system: system.euler, nodes: system.evaluate(x).outer.nodes });
  return result;
}

test('logical quality equals the published-cell oracle for valid, flat and folded grids', () => {
  for (const { system, nodes: original } of fixtures()) for (const shape of ['original', 'flat', 'folded', 'near-flat']) {
    const nodes = structuredClone(original);
    if (shape !== 'original') {
      const a = nodes[0][1][0], b = nodes[0][1][1];
      const fraction = shape === 'flat' ? 0 : shape === 'folded' ? -1 : 1e-14;
      nodes[0][1][1] = { x: a.x + fraction * (b.x - a.x), y: a.y + fraction * (b.y - a.y) };
    }
    const expected = publishedQuality(system.layout, nodes);
    assert.deepEqual(streamtubeMeshQuality({ system, nodes }), expected);
    assert.deepEqual(streamtubeMeshSnapshot({ system, nodes }).quality, expected);
    if (shape === 'flat' || shape === 'folded') assert.equal(expected.valid, false);
  }
});

test('quality-only acceptance still rejects disconnected shared cuts and endpoints', () => {
  for (const { system, nodes } of fixtures()) {
    const { indices } = streamtubeMeshConnectivity(system.layout, nodes), seen = new Set();
    let tested = 0;
    for (let g = 0; g < indices.length; g++) for (let i = 0; i < indices[g].length; i++) for (let j = 0; j < indices[g][i].length; j++) {
      const id = indices[g][i][j];
      if (seen.has(id)) {
        const changed = structuredClone(nodes);
        changed[g][i][j] = { ...changed[g][i][j], y: changed[g][i][j].y + 1e-8 };
        assert.throws(() => publishedQuality(system.layout, changed), /Disconnected/);
        assert.throws(() => streamtubeMeshQuality({ system, nodes: changed }), /Disconnected/);
        tested++;
      }
      seen.add(id);
    }
    assert.ok(tested > 0);
  }
});
