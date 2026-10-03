import test from 'node:test';
import assert from 'node:assert/strict';
import { solveInviscid, velocityAt } from '../src/inviscid/linear-vortex.js';
import { traceStreamlines } from '../src/inviscid/streamlines.js';
import { naca4, segmentsTouch } from '../src/geometry/airfoil.js';
test('display streamlines follow computed velocity and never cross body panels', () => {
  const r = solveInviscid({ elements: [{ points: naca4('0012', 60) }], alpha: 4 });
  const lines = traceStreamlines(r, { xMin: -0.5, xMax: 1.5, yMin: -0.5, yMax: 0.5 }, 9);
  assert.equal(lines.length, 9); assert.ok(lines.every(line => line.length > 10));
  for (const line of lines) for (let i = 1; i < line.length; i++) {
    const a = line[i - 1]; const b = line[i];
    for (const panel of r.field.panels) assert.equal(segmentsTouch(a, b, panel.a, panel.b, 1e-12), false);
    const mid = { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 };
    const v = velocityAt(mid, r.field);
    const error = Math.abs((b.x - a.x) * v.v - (b.y - a.y) * v.u) / (Math.hypot(b.x - a.x, b.y - a.y) * Math.hypot(v.u, v.v));
    assert.ok(error < 0.02);
  }
});
