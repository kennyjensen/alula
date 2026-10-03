import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createStagnationGuideConnector } from '../src/inviscid/stagnation-guide-connector.js';
import { pointInside, segmentsTouch } from '../src/geometry/airfoil.js';
import { potentialDifference } from '../src/inviscid/streamfunction.js';

const field = (u = 1, v = 0) => ({ u, v, panels: [], gamma: [] });
const straight = () => ({ join: { x: -1, y: 0, potential: 3 }, anchor: { x: 0, y: 0 },
  anchorDirection: { x: 1, y: 0 }, stagnationPotential: 4, field: field(), tolerance: 1e-12 });
const close = (a, b, tolerance = 2e-12) => assert.ok(Math.abs(a - b) <= tolerance, `Absolute error ${Math.abs(a - b)}`);

test('straight uniform field has exact geometry, common potential and inverse at both endpoints and interiors', () => {
  const input = straight(), before = structuredClone(input), c = createStagnationGuideConnector(input);
  for (const x of [-1, -.95, -.6, -.1, 0]) {
    const p = c.atX(x); close(p.x, x); close(p.y, 0); close(p.potential, 4 + x); close(p.derivative, 1);
    close(c.atPotential(4 + x).x, x);
  }
  assert.deepEqual(input, before);
  assert.deepEqual(c.anchor, { x: 0, y: 0, potential: 4 });
});

test('curved geometric path uses physical potential instead of treating its curve parameter as potential', () => {
  const input = { ...straight(), join: { x: -1, y: -.3, potential: 2.85 },
    anchorDirection: { x: 1, y: 1.2 }, field: field(1, .5) };
  const c = createStagnationGuideConnector(input);
  for (const x of [-1, -.8, -.5, -.2, 0]) {
    const p = c.atX(x); close(p.potential, 4 + p.x + .5 * p.y);
    close(c.atPotential(p.potential).x, x);
  }
  const dx = 1e-6, slope = (c.atX(-dx).y - c.anchor.y) / -dx;
  close(slope, 1.2, 3e-6);
  assert.match(c.diagnostics.interpretation, /not an exact constant-streamfunction/);
});

test('length scaling and translation preserve a geometric curve and its physical-potential differences', () => {
  const a = { ...straight(), join: { x: -1, y: -.3, potential: 2.85 },
    anchorDirection: { x: 1, y: 1.2 }, field: field(1, .5) }, original = createStagnationGuideConnector(a);
  const scale = 7, shift = { x: 3, y: -4 }, transform = p => ({ x: scale * p.x + shift.x, y: scale * p.y + shift.y });
  const b = { ...a, join: { ...transform(a.join), potential: scale * a.join.potential + 9 },
    anchor: transform(a.anchor), stagnationPotential: scale * a.stagnationPotential + 9, tolerance: scale * a.tolerance };
  const transformed = createStagnationGuideConnector(b);
  for (const x of [-1, -.75, -.5, -.25, 0]) {
    const p = original.atX(x), q = transformed.atX(scale * x + shift.x), expected = transform(p);
    close(q.x, expected.x); close(q.y, expected.y); close(q.potential, scale * p.potential + 9);
  }
});

function square(x0, x1, y0, y1) {
  const points = [{ x: x0, y: y0 }, { x: x1, y: y0 }, { x: x1, y: y1 }, { x: x0, y: y1 }];
  points.push({ ...points[0] });
  const panels = points.slice(1).map((b, k) => {
    const a = points[k], length = Math.hypot(b.x - a.x, b.y - a.y);
    return { a, b, length, tx: (b.x - a.x) / length, ty: (b.y - a.y) / length, node: k, element: 0 };
  });
  return { points, panels };
}

test('Bezier hull certificate rejects an obstacle hit by the curved path even when its endpoint chord is clear', () => {
  const obstacle = square(.45, .55, .2, .3);
  const a = { x: 0, y: 0 }, b = { x: 1, y: 0 };
  assert.equal(obstacle.panels.some(p => segmentsTouch(a, b, p.a, p.b)), false);
  assert.throws(() => createStagnationGuideConnector({ join: { ...a, potential: 0 }, anchor: b,
    anchorDirection: { x: 1, y: -1 }, stagnationPotential: 1,
    field: { ...field(1, 1), panels: obstacle.panels, gamma: Array(5).fill(0) },
    admissible: p => !pointInside(p, obstacle.points) }), /cannot be certified clear/);
});

test('reverse flow, nonmonotone physical potential, unrelated gauges and out-of-range queries remain rejected', () => {
  assert.throws(() => createStagnationGuideConnector({ ...straight(), field: field(-1) }), /forward/);
  assert.throws(() => createStagnationGuideConnector({ ...straight(), field: field(1, 1),
    anchorDirection: { x: 1, y: -2 } }), /nonincreasing panel potential/);
  assert.throws(() => createStagnationGuideConnector({ ...straight(), stagnationPotential: 6 }), /common potential gauge/);
  const c = createStagnationGuideConnector(straight());
  assert.throws(() => c.atX(.1), /outside/);
  assert.throws(() => c.atPotential(4.1), /outside/);
});

test('frozen RAE connector fills all missing near-LE requests on an exterior monotone common potential branch', () => {
  const input = JSON.parse(fs.readFileSync(new URL('./fixtures/rae-stagnation-guide.json', import.meta.url)));
  const before = structuredClone(input), fluid = p => !input.bodies.some(b => pointInside(p, b));
  const segment = (a, b) => !input.field.panels.some(p => segmentsTouch(a, b, p.a, p.b, 1e-12));
  const c = createStagnationGuideConnector({ ...input, tolerance: 2e-9, admissible: fluid, admissibleSegment: segment });
  assert.ok(c.diagnostics.certificateLeaves > 0);
  assert.ok(c.diagnostics.minimumPotentialDerivative > .003);
  assert.ok(c.diagnostics.maximumStreamfunctionDefect > 1e-5, 'Keep the initializer approximation visible.');
  assert.ok(Math.abs(c.diagnostics.joinPotentialError) < 3e-11);
  let previous = c.join;
  for (const x of input.missingRequests) {
    const p = c.atX(x); close(p.x, x, 1e-15); assert.equal(fluid(p), true);
    assert.ok(p.potential > previous.potential && p.potential < c.anchor.potential);
    assert.equal(segment(previous, p), true);
    close(c.anchor.potential - p.potential, potentialDifference(p, c.anchor, input.field), 3e-14);
    close(c.atPotential(p.potential).x, p.x, 2e-9);
    previous = p;
  }
  assert.deepEqual(c.anchor, { ...input.anchor, potential: input.stagnationPotential });
  assert.deepEqual(input, before);
});
