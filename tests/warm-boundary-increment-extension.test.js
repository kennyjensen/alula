// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { extendWarmBoundaryIncrements } from '../src/euler/streamtube-displacement.js';

const fixture = () => ({
  sourceNodes: [
    [[{ x: -0, y: 0 }, { x: .1, y: 1 }, { x: .2, y: 2 }], [{ x: 1, y: 0 }, { x: 1.1, y: 1 }, { x: 1.2, y: 2 }]],
    [[{ x: .2, y: 2 }, { x: .3, y: 3 }, { x: .4, y: 4 }, { x: .5, y: 5 }], [{ x: 1.2, y: 2 }, { x: 1.3, y: 3 }, { x: 1.4, y: 4 }, { x: 1.5, y: 5 }]],
    [[{ x: .5, y: 5 }, { x: .6, y: 6 }, { x: .7, y: 7 }], [{ x: 1.5, y: 5 }, { x: 1.6, y: 6 }, { x: 1.7, y: 7 }]],
  ], masses: [[1, 3], [1, 2, 1], [3, 1]],
});
const near = (a, b) => assert.ok(Math.abs(a - b) < 2e-14, `${a} differs from ${b}`);

test('zero boundary increments preserve every physical coordinate exactly and detach points', () => {
  const f = fixture(), targetNodes = structuredClone(f.sourceNodes), before = structuredClone(f);
  const result = extendWarmBoundaryIncrements({ ...f, targetNodes });
  assert.deepEqual(result, f.sourceNodes); assert.deepEqual(f, before);
  for (let g = 0; g < result.length; g++) for (let i = 0; i < result[g].length; i++)
    for (let j = 0; j < result[g][i].length; j++) assert.notEqual(result[g][i][j], f.sourceNodes[g][i][j]);
  assert.ok(Object.is(result[0][0][0].x, -0));
});

test('two body boundaries blend once into their common passage with unequal masses', () => {
  const f = fixture(), targetNodes = structuredClone(f.sourceNodes), before = structuredClone(f);
  for (let i = 0; i < 2; i++) {
    targetNodes[0][i][2].x += .04; targetNodes[0][i][2].y += .2;
    targetNodes[1][i][0] = { ...targetNodes[0][i][2] };
    targetNodes[1][i][3].x -= .08; targetNodes[1][i][3].y -= .1;
    targetNodes[2][i][0] = { ...targetNodes[1][i][3] };
    // Target interior coordinates must not be treated as another displacement.
    targetNodes[1][i][1].y = 90;
  }
  const targetBefore = structuredClone(targetNodes), result = extendWarmBoundaryIncrements({ ...f, targetNodes });
  for (let i = 0; i < 2; i++) {
    near(result[0][i][1].y - f.sourceNodes[0][i][1].y, .05);
    near(result[1][i][1].y - f.sourceNodes[1][i][1].y, .125);
    near(result[1][i][2].y - f.sourceNodes[1][i][2].y, -.025);
    near(result[2][i][1].y - f.sourceNodes[2][i][1].y, -.025);
    near(result[1][i][1].x - f.sourceNodes[1][i][1].x, .01);
    assert.deepEqual(result[0][i][2], result[1][i][0]);
    assert.deepEqual(result[1][i][3], result[2][i][0]);
  }
  assert.deepEqual(f, before); assert.deepEqual(targetNodes, targetBefore);
});

test('independent noncoincident wake banks and outer endpoints remain exact', () => {
  const f = fixture(), targetNodes = structuredClone(f.sourceNodes);
  targetNodes[1][1][0].y += .03; // retained distinct bank, not silently averaged
  const result = extendWarmBoundaryIncrements({ ...f, targetNodes });
  for (let g = 0; g < result.length; g++) for (let i = 0; i < result[g].length; i++) {
    assert.deepEqual(result[g][i][0], targetNodes[g][i][0]);
    assert.deepEqual(result[g][i].at(-1), targetNodes[g][i].at(-1));
  }
  assert.notDeepEqual(result[0][1].at(-1), result[1][1][0]);
});

test('rigid translation, rotation, and length/mass scaling commute with extension', () => {
  const f = fixture(), targetNodes = structuredClone(f.sourceNodes);
  targetNodes[0][1][2].y += .12;
  const result = extendWarmBoundaryIncrements({ ...f, targetNodes });
  const transform = p => ({ x: 3 - 2 * p.y, y: 4 + 2 * p.x });
  const map = nodes => nodes.map(g => g.map(row => row.map(transform)));
  const transformed = extendWarmBoundaryIncrements({ sourceNodes: map(f.sourceNodes), targetNodes: map(targetNodes), masses: f.masses.map(g => g.map(m => 7 * m)) });
  const expected = map(result);
  for (let g = 0; g < result.length; g++) for (let i = 0; i < result[g].length; i++) for (let j = 0; j < result[g][i].length; j++) {
    near(transformed[g][i][j].x, expected[g][i][j].x); near(transformed[g][i][j].y, expected[g][i][j].y);
  }
});

test('invalid mass and topology are rejected', () => {
  const f = fixture(), targetNodes = structuredClone(f.sourceNodes);
  assert.throws(() => extendWarmBoundaryIncrements({ ...f, targetNodes, masses: [[0, 1], ...f.masses.slice(1)] }), /positive/);
  assert.throws(() => extendWarmBoundaryIncrements({ ...f, targetNodes: targetNodes.slice(1) }), /matching/);
  targetNodes[1][0].pop(); assert.throws(() => extendWarmBoundaryIncrements({ ...f, targetNodes }), /dimensions/);
});
