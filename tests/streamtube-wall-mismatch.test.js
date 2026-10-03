import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { naca4, transform, pointInside } from '../src/geometry/airfoil.js';
import { solveInviscid, velocityAt } from '../src/inviscid/linear-vortex.js';
import { streamfunctionAt } from '../src/inviscid/streamfunction.js';
import { isentropicSectionDensity } from '../src/euler/streamtube-initial-state.js';

test('frozen near-wall sections isolate mass/geometry mismatch without tracing or solving an airfoil grid', () => {
  const { cases } = JSON.parse(readFileSync(new URL('./fixtures/streamtube-wall-mismatch.json', import.meta.url)));
  for (const c of cases) {
    const mainFlap = c.name === 'main-flap', n = c.conditions.contourPanels;
    const elements = mainFlap ? [{ points: naca4('2412', n) },
      { points: transform(naca4('0012', n), { chord: .3, x: .94, y: -.08, angle: -15 }) }]
      : [{ points: naca4('0012', n) }, { points: transform(naca4('0012', n), { chord: .3, x: -.4, y: .2 }) }];
    const { field } = solveInviscid({ elements, alpha: c.conditions.alpha, boundaryCondition: c.conditions.panelBoundaryCondition });
    for (const s of c.sections) {
      const [a, b, d, e] = s.points, psi = s.points.map(p => streamfunctionAt(p, field));
      psi.forEach((value, i) => assert.ok(Math.abs(value - s.recordedStreamfunctions[i]) < 1e-11));
      const tx = (b.x + e.x - a.x - d.x) / 2, ty = (b.y + e.y - a.y - d.y) / 2;
      const ax = (d.x + e.x - a.x - b.x) / 2, ay = (d.y + e.y - a.y - b.y) / 2;
      const area = (tx * ay - ty * ax) / Math.hypot(tx, ty);
      // Cross-line endpoint differences give actual panel flux at each
      // station. Their mean diagnoses the error relative to allocated mass;
      // it is not a replacement mass assignment for the Euler solver.
      const actualMass = .5 * (psi[2] - psi[0] + psi[3] - psi[1]);
      assert.ok(area > 0 && actualMass > 0);
      assert.ok(s.massFlow / actualMass > (mainFlap ? 1.07 : 1.8));
      // The tracer controls local position error, not a global psi bound.
      // Here the fluid trace drift is over 1000 times smaller than the
      // wall/mass mismatch; tighter tracing alone cannot remove the latter.
      assert.ok(Math.abs(psi[2] - psi[3]) < .001 * (s.massFlow - actualMass));
      if (!mainFlap) assert.throws(() => isentropicSectionDensity({ massFlux: s.massFlow / area, mach: .2 }), /sonic/);
    }
    if (!mainFlap) {
      const s = c.sections[0], center = s.points.reduce((p, q) => ({ x: p.x + q.x / 4, y: p.y + q.y / 4 }), { x: 0, y: 0 });
      assert.ok(elements.some(e => pointInside(center, e.points)));
      // Mid-cell panel velocity is an interior-field query here. The actual
      // fluid nodes have substantial speed, so do not interpret it as stall.
      const speed = p => { const q = velocityAt(p, field); return Math.hypot(q.u, q.v); };
      assert.ok(speed(center) < .02);
      s.points.slice(2).forEach(p => assert.ok(speed(p) > 1.7 && !elements.some(e => pointInside(p, e.points))));
    }
  }
});
