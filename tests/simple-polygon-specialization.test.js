// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { requirePositiveSimplePolygon as check } from '../src/geometry/simple-polygon.js';
import { requirePositiveSimplePolygon as reference } from './oracles/simple-polygon-reference.js';
function outcome(fn, p) { try { return { area: fn(p) }; } catch (e) { return { error: e.message }; } }
let seed = 15485863;
const random = () => ((seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0) / 2 ** 32);

test('fixed polygon kernels preserve signed area and domain rejection on ordinary and extreme cells', () => {
  for (const n of [4, 6]) for (let k = 0; k < 12000; k++) {
    const scale = [1, 1e-150, 1e150][k % 3], offset = k % 7 === 0 ? 1e12 : 0;
    const p = Array.from({ length: n }, (_, i) => {
      const a = 2 * Math.PI * i / n, r = .1 + random();
      return { x: offset + scale * r * Math.cos(a), y: scale * r * Math.sin(a) };
    });
    if (k % 4 === 0) [p[1], p[2]] = [p[2], p[1]];
    if (k % 11 === 0) p[2] = { ...p[1] };
    if (k % 13 === 0) p.reverse();
    assert.deepEqual(outcome(check, p), outcome(reference, p), JSON.stringify(p));
  }
  for (const p of [
    [{x:0,y:0},{x:1,y:0},{x:.5,y:0},{x:0,y:1}],
    [{x:0,y:0},{x:1,y:0},{x:1,y:1},{x:NaN,y:1}],
    [{x:-1e308,y:0},{x:1e308,y:0},{x:1e308,y:1},{x:-1e308,y:1}],
    [{x:0,y:0},{x:1,y:0},{x:1,y:1},{x:0,y:1},{x:0,y:.5},{x:.5,y:0}],
  ]) assert.deepEqual(outcome(check,p),outcome(reference,p));
});
