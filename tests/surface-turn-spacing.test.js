import test from 'node:test';
import assert from 'node:assert/strict';
import { resolveSurfaceTurns } from '../src/geometry/surface-turn-spacing.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { naca4, transform } from '../src/geometry/airfoil.js';

test('turn refinement resolves an analytic parabola and preserves supplied fractions', () => {
  const curve = { length: 2, knots: [0, 1, 2], evaluate: s => ({ derivative: { x: 1, y: 4 * (s - 1) }, secondDerivative: { x: 0, y: 4 } }) };
  for (const side of ['upper', 'lower']) {
    const r = resolveSurfaceTurns({ curve, side, stagnation: 1, fractions: [0, .4, 1] });
    assert.ok(r.addedPoints > 0); assert.ok(r.fractions.includes(.4));
    for (let i = 1; i < r.fractions.length; i++) {
      // Exact tangent angle for (x,y)=(s,2(s-1)^2).
      const angle = f => Math.atan(4 * f * (side === 'upper' ? -1 : 1));
      assert.ok(Math.abs(angle(r.fractions[i]) - angle(r.fractions[i - 1])) <= Math.PI / 12 + 1e-14);
    }
  }
  assert.throws(() => resolveSurfaceTurns({ curve, side: 'upper', stagnation: 1, fractions: [0, 1], maxPoints: 2 }), /budget/);
});

test('spline tangent-cone refinement is invariant to units and resolves the full cambered nose', () => {
  const points = naca4('2412', 80), curve = createContourCurve(points), moved = createContourCurve(transform(points, { chord: 3, angle: 37, x: 2, y: -1 }));
  const fractions = Array.from({ length: 9 }, (_, i) => .5 * (1 - Math.cos(Math.PI * i / 8)));
  for (const side of ['upper', 'lower']) {
    const r = resolveSurfaceTurns({ curve, side, stagnation: .53 * curve.length, fractions });
    const m = resolveSurfaceTurns({ curve: moved, side, stagnation: .53 * moved.length, fractions });
    assert.deepEqual(m.fractions, r.fractions);
    for (let i = 1; i < r.fractions.length; i++) {
      const directions = Array.from({ length: 129 }, (_, k) => {
        const f = r.fractions[i - 1] + k / 128 * (r.fractions[i] - r.fractions[i - 1]);
        return curve.branch(side, f, .53 * curve.length).derivative;
      });
      const a = directions[0], angles = directions.map(b => Math.atan2(a.x * b.y - a.y * b.x, a.x * b.x + a.y * b.y));
      assert.ok(Math.max(...angles) - Math.min(...angles) <= Math.PI / 12 + 1e-12);
    }
  }
});
