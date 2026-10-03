import test from 'node:test';
import assert from 'node:assert/strict';
import { machContours } from '../src/ui/mach-contours.js';
const line = (tube, y, values, group = 0) => ({group, tube, points: [0, 1, 2].map(x => ({ x, y })), machNumbers: values});
test('linear Mach field gives exact contour locations without changing samples', () => {
 const flow = { lines: [line(0, 0, [.5, 1.5]), line(1, 1, [.5, 1.5])] };
 const before = structuredClone(flow), sonic = machContours(flow).find(c => c.level === 1);
 assert.ok(sonic.segments.length);
 for (const segment of sonic.segments) for (const p of segment) assert.ok(Math.abs(p.x - 1) < 1e-14);
 assert.deepEqual(flow, before);
});
test('contours never bridge distinct passages or missing Mach samples', () => {
 assert.deepEqual(machContours({lines: [line(0, 0, [.5, 1.5]), line(1, 1, [.5, 1.5], 1)]}), []);
 assert.deepEqual(machContours({lines: [line(0, 0, [null, null]), line(1, 1, [.5, 1.5])]}), []);
 assert.deepEqual(machContours({lines: [line(0, 0, [.5, .5]), line(1, 1, [.5, .5])]}), []);
});
