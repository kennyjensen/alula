import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { untangleQuadrilaterals } from '../src/geometry/untangle-quadrilaterals.js';
import { isentropicSectionDensity } from '../src/euler/streamtube-initial-state.js';
import { repairSubsonicStreamtubeGrid } from '../src/euler/streamtube-grid-repair.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const vertices = [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: 1, y: .2 }, { x: 0, y: .2 }];
const cells = [[0, 1, 2, 3]], fixed = new Set([0, 1, 3]);
const width = p => {
  // Independent shoelace area divided by streamwise midpoint distance.
  const area = .5 * p.reduce((a, q, i) => a + q.x * p[(i + 1) % 4].y - q.y * p[(i + 1) % 4].x, 0);
  return 2 * area / Math.hypot(p[1].x + p[2].x - p[0].x - p[3].x, p[1].y + p[2].y - p[0].y - p[3].y);
};

test('one free vertex recovers the independently minimized section-width constraint', () => {
  // Solve the boundary conic explicitly for y(x), then minimize distance in
  // one dimension. This oracle does not use half-plane cuts or gradients.
  const y = x => {
    const a = .91, b = .4 * x + .036, c = .04 * x * x - .09 * ((1 + x) ** 2 + .04);
    return (-b + Math.sqrt(b * b - 4 * a * c)) / (2 * a);
  };
  const objective = x => (x - 1) ** 2 + (y(x) - .2) ** 2;
  let lo = .5, hi = 1.5;
  for (let i = 0; i < 100; i++) {
    const a = (2 * lo + hi) / 3, b = (lo + 2 * hi) / 3;
    if (objective(a) < objective(b)) hi = b; else lo = a;
  }
  const x = .5 * (lo + hi), result = untangleQuadrilaterals({ vertices, cells, fixed }, { minimumNormalAreas: [.3] });
  assert.equal(result.converged, true, result.reason);
  assert.ok(Math.hypot(result.vertices[2].x - x, result.vertices[2].y - y(x)) < 2e-7);
  assert.ok(Math.abs(width(result.vertices) - .3) < 1e-10);
  fixed.forEach(i => assert.deepEqual(result.vertices[i], vertices[i]));
});

test('normal-width repair respects fixed/movement limits and transforms with physical units', () => {
  const original = structuredClone(vertices), controls = { minimumNormalAreas: [.3] };
  for (const options of [{ fixed: new Set([0, 1, 2, 3]) }, { fixed, displacementLimit: .01 }]) {
    const r = untangleQuadrilaterals({ vertices, cells, fixed: options.fixed }, { ...controls, ...options, maxSweeps: 5 });
    assert.equal(r.converged, false); assert.ok(width(r.vertices) < .3);
    assert.deepEqual(r.quality.invalidNormalAreaCells, [0]);
  }
  const a = untangleQuadrilaterals({ vertices, cells, fixed }, controls);
  const move = p => ({ x: 3 + 2 * (p.x * Math.cos(.4) - p.y * Math.sin(.4)), y: -1 + 2 * (p.x * Math.sin(.4) + p.y * Math.cos(.4)) });
  const b = untangleQuadrilaterals({ vertices: vertices.map(move), cells, fixed }, { minimumNormalAreas: [.6] });
  assert.equal(b.converged, true); assert.ok(Math.abs(width(b.vertices) - .6) < 2e-10);
  a.vertices.forEach((p, i) => assert.ok(Math.hypot(move(p).x - b.vertices[i].x, move(p).y - b.vertices[i].y) < 2e-7));
  assert.deepEqual(vertices, original);
  for (const minimumNormalAreas of [[], [-1], [NaN]])
    assert.throws(() => untangleQuadrilaterals({ vertices, cells }, { minimumNormalAreas }), /Invalid/);
});

test('frozen 24-cell main/flap patch obtains positive corners and admissible subsonic section states with fixed boundaries', () => {
  const f = JSON.parse(readFileSync(new URL('./fixtures/coarse-normal-area-patch.json', import.meta.url)));
  const { mach, gamma } = f.conditions;
  // Compute the sonic bound independently from local Mach = 1.
  const temperature = (1 + .5 * (gamma - 1) * mach * mach) / (1 + .5 * (gamma - 1));
  const fluxLimit = f.massFluxFraction * temperature ** (1 / (gamma - 1) + .5) / mach;
  const minimumNormalAreas = f.massFlows.map(m => m / fluxLimit);
  const badCorners = p => p.filter((a, i) => {
    const b = p[(i + 1) % 4], c = p[(i + 2) % 4];
    return (b.x - a.x) * (c.y - b.y) <= (b.y - a.y) * (c.x - b.x);
  }).length;
  assert.ok(f.cells.some(c => badCorners(c.map(i => f.vertices[i])) > 0));
  const r = untangleQuadrilaterals({ ...f, fixed: new Set(f.fixed) }, { minimumCorner: .2, minimumNormalAreas });
  assert.equal(r.converged, true, r.reason); assert.ok(r.maxDisplacement < .001);
  f.fixed.forEach(i => assert.deepEqual(r.vertices[i], f.vertices[i]));
  f.cells.forEach((c, i) => {
    const p = c.map(id => r.vertices[id]); assert.equal(badCorners(p), 0);
    const area = width(p), massFlux = f.massFlows[i] / area;
    assert.ok(massFlux < fluxLimit);
    const rho = isentropicSectionDensity({ ...f.conditions, massFlux }), q = massFlux / rho;
    const t = 1 + .5 * (gamma - 1) * mach * mach * (1 - q * q);
    assert.ok(mach * mach * q * q / t < 1);
    assert.ok(Math.abs(rho ** (gamma - 1) / t - 1) < 3e-14);
  });
});

test('subsonic body wrapper assigns every cell its own tube mass limit and preserves the input chart', () => {
  const system = createStreamtubeBodySystem(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }));
  const { nodes, allocation } = system.decode(system.initial), chart = system.geometryChart(), before = structuredClone(nodes);
  const r = repairSubsonicStreamtubeGrid({ system, nodes });
  const { mach, gamma } = system.conditions;
  const t = (1 + .5 * (gamma - 1) * mach * mach) / (1 + .5 * (gamma - 1));
  const limit = .99 * t ** (1 / (gamma - 1) + .5) / mach;
  let index = 0;
  for (let g = 0; g < nodes.length; g++) for (let i = 0; i < system.layout.nx; i++) for (let j = 0; j < system.layout.tubes[g]; j++) {
    assert.ok(Math.abs(r.minimumNormalAreas[index++] * limit / allocation.groups[g][j].massFlow - 1) < 1e-14);
  }
  assert.equal(index, r.minimumNormalAreas.length); assert.deepEqual(nodes, before); assert.deepEqual(system.geometryChart(), chart);
  assert.equal(r.converged, true); assert.match(r.status, /flow equations have not been solved/);
  for (const massFluxFraction of [0, 1, NaN])
    assert.throws(() => repairSubsonicStreamtubeGrid({ system, nodes }, { massFluxFraction }), /conditions/);
});
