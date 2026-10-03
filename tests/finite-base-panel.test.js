import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
import { solveInviscid, velocityAt } from '../src/inviscid/linear-vortex.js';
import { streamfunctionAt, potentialDifference } from '../src/inviscid/streamfunction.js';
import { xfoilSurfaceDerivatives, basePanelStreamfunctionBasis, finiteBaseSourceChart } from '../src/inviscid/finite-base-influence.js';
import { makePanel, sourceVelocity } from '../src/inviscid/panel.js';
import { naca4, transform } from '../src/geometry/airfoil.js';
import { getBenchmarkAirfoil } from '../src/geometry/benchmark-airfoils.js';
import { createContourTopology } from '../src/geometry/contour-topology.js';

const native = JSON.parse(fs.readFileSync('tests/fixtures/fortran/finite-base-panel.json'));
const close = (a, b, tolerance = 2e-10) => assert.ok(Math.abs(a - b) <= tolerance, `${a} != ${b}; error ${Math.abs(a - b)}`);
function element(surface, segments = 1, bulge = 0) {
  const points = structuredClone(surface), lowerIndex = points.length - 1, a = points.at(-1), b = points[0];
  for (let k = 1; k <= segments; k++) {
    const t = k / segments;
    points.push(k === segments ? { ...b } : { x: a.x + t * (b.x - a.x) + bulge * Math.sin(Math.PI * t), y: a.y + t * (b.y - a.y) });
  }
  return { points, trailingEdge: { kind: 'finite-base', upperIndex: 0, lowerIndex } };
}
const baseline = () => element(native.cases[0].points);
const fieldChecks = (r, points) => {
  const h = 2e-6;
  const derivative = (point, axis, fn) => [-2, -1, 1, 2].map(k => fn({ ...point, [axis]: point[axis] + k * h }))
    .reduce((sum, v, i) => sum + [1, -8, 8, -1][i] * v, 0) / (12 * h);
  for (const p of points) {
    const q = velocityAt(p, r.field);
    close(derivative(p, 'x', z => streamfunctionAt(z, r.field)), -q.v, 1e-8);
    close(derivative(p, 'y', z => streamfunctionAt(z, r.field)), q.u, 1e-8);
    close(derivative(p, 'x', z => potentialDifference(p, z, r.field)), q.u, 1e-8);
    close(derivative(p, 'y', z => potentialDifference(p, z, r.field)), q.v, 1e-8);
  }
};

test('finite-base nodal strengths, base projections, force and off-body velocity match executed original Fortran', () => {
  assert.equal(native.passed, true);
  for (const c of native.cases) {
    const e = element(c.points), before = structuredClone(e), r = solveInviscid({ elements: [e], alpha: c.alpha });
    assert.deepEqual(e, before); assert.equal(r.status, 'solved'); assert.equal(r.boundaryCondition, 'streamfunction');
    const d = xfoilSurfaceDerivatives(c.points).derivatives;
    c.expected.nodes.forEach((q, i) => {
      close(r.field.gamma[i], -q.gamma); close(d[i].x, q.dxds); close(d[i].y, q.dyds);
    });
    close(r.field.basePanels[0].sourceStrength, c.expected.sourceStrength);
    close(r.field.basePanels[0].vortexStrength, -c.expected.vortexStrength);
    for (const key of ['cl', 'cm']) close(r[key], c.expected[key]);
    close(r.diagnostics.pressureDrag, c.expected.pressureDrag);
    c.expected.fields.forEach(p => { const q = velocityAt(p, r.field); close(q.u, p.u); close(q.v, p.v); });
    assert.equal(r.diagnostics.unknowns, c.points.length + 1);
    assert.equal(r.diagnostics.surfaceRows, c.points.length); assert.equal(r.diagnostics.kuttaRows, 1);
    for (const p of c.points) close(streamfunctionAt(p, r.field), r.diagnostics.surfaceStreamfunctions[0]);
  }
});

test('every collinear base segment is retained and subdivision preserves the complete potential solution', () => {
  const e1 = baseline(), e32 = element(native.cases[0].points, 32), a = solveInviscid({ elements: [e1], alpha: 4 }), b = solveInviscid({ elements: [e32], alpha: 4 });
  assert.deepEqual(b.elements[0].points, e32.points);
  assert.equal(b.field.panels.length, 40); assert.equal(b.field.basePanels.length, 32);
  assert.equal(b.panelCount, 72); assert.equal(b.diagnostics.unknowns, 42);
  assert.equal(b.elements[0].cp.filter(p => p.base).length, 32);
  a.field.gamma.forEach((g, i) => close(g, b.field.gamma[i], 2e-11));
  for (const key of ['cl', 'cm']) close(a[key], b[key], 2e-11);
  for (const key of ['baseSourceFlux', 'baseCirculation']) close(a.diagnostics[key], b.diagnostics[key], 2e-12);
  for (const p of native.cases[0].queries) { close(streamfunctionAt(p, a.field), streamfunctionAt(p, b.field), 2e-12); }
  // Inside the old geometric strip, both sides of the single common cut are
  // fluid chart points even at the former per-panel midpoint cuts.
  const insideStrip = [{ x: 1.2, y: .0005 }, { x: 1.2, y: -.0005 }];
  for (const p of insideStrip) close(streamfunctionAt(p, a.field), streamfunctionAt(p, b.field), 2e-12);
  fieldChecks(b, insideStrip);
});

test('noncollinear base geometry is used as supplied while integrated TECALC projections telescope', () => {
  const e = element(native.cases[1].points, 8, .0008), r = solveInviscid({ elements: [e], alpha: 7 });
  assert.deepEqual(r.elements[0].points, e.points);
  const t = r.elements[0].finiteBase.teDerivative.vector, a = e.points[e.trailingEdge.lowerIndex], b = e.points[0];
  const difference = r.field.gamma[0] - r.field.gamma[e.trailingEdge.lowerIndex], dx = b.x - a.x, dy = b.y - a.y;
  close(r.diagnostics.baseSourceFlux, -.5 * difference * (t.x * dy - t.y * dx), 2e-15);
  close(r.diagnostics.baseCirculation, -.5 * difference * (t.x * dx + t.y * dy), 2e-15);
  const straight = solveInviscid({ elements: [element(native.cases[1].points, 8)], alpha: 7 });
  // The combined projected source/vortex differential is analytic and its
  // exterior integral depends on endpoints. This model retains the supplied
  // cap but does not impose no-penetration on each base segment.
  r.field.gamma.forEach((g, i) => close(g, straight.field.gamma[i], 2e-11));
  assert.ok(Math.max(...r.field.basePanels.map(p => p.sourceStrength)) - Math.min(...r.field.basePanels.map(p => p.sourceStrength)) > .01);
  for (const p of native.cases[1].queries) {
    let u = 0, v = 0;
    for (const b of r.field.basePanels) {
      const n = 128;
      for (let i = 0; i <= n; i++) {
        const s = i / n, dx = p.x - b.a.x - s * (b.b.x - b.a.x), dy = p.y - b.a.y - s * (b.b.y - b.a.y);
        const w = (i === 0 || i === n ? 1 : i % 2 ? 4 : 2) * b.length / (3 * n * 2 * Math.PI * (dx * dx + dy * dy));
        u += w * (b.sourceStrength * dx - b.vortexStrength * dy);
        v += w * (b.sourceStrength * dy + b.vortexStrength * dx);
      }
    }
    const all = velocityAt(p, r.field), surface = velocityAt(p, { ...r.field, basePanels: [] });
    close(all.u - surface.u, u, 2e-12); close(all.v - surface.v, v, 2e-12);
  }
  fieldChecks(r, native.cases[1].queries);
});

test('complete finite-base field has consistent velocity/psi/potential and the intended nonzero flux and source-cut jump', () => {
  const r = solveInviscid({ elements: [baseline()], alpha: 4 });
  fieldChecks(r, native.cases[0].queries);
  let flux = 0, circulation = 0;
  const radius = 2, count = 192, da = 2 * Math.PI / count;
  for (let i = 0; i < count; i++) {
    const a = da * (i + .5), nx = Math.cos(a), ny = Math.sin(a), q = velocityAt({ x: .5 + radius * nx, y: radius * ny }, r.field);
    flux += (q.u * nx + q.v * ny) * radius * da;
    circulation += (-q.u * ny + q.v * nx) * radius * da;
  }
  assert.ok(r.diagnostics.baseSourceFlux > .001);
  close(flux, r.diagnostics.baseSourceFlux, 3e-12);
  close(circulation, r.elements[0].circulation, 3e-12);
  assert.throws(() => streamfunctionAt({ x: 1.2, y: 0 }, r.field), /source-cut ray/);
  const bottom = { x: 1.2, y: -.02 }, top = { x: 1.2, y: .02 }, n = 128, dy = (top.y - bottom.y) / n;
  let integral = 0;
  for (let i = 0; i <= n; i++) integral += (i === 0 || i === n ? 1 : i % 2 ? 4 : 2) * velocityAt({ x: 1.2, y: bottom.y + i * dy }, r.field).u * dy / 3;
  close(streamfunctionAt(top, r.field) - streamfunctionAt(bottom, r.field) - integral, -r.diagnostics.baseSourceFlux, 2e-12);
});

test('individual source charts differentiate correctly inside the geometric base strip away from their actual ray', () => {
  const b = { ...makePanel({ x: 0, y: -.5 }, { x: 0, y: .5 }), cutDirection: { x: 1, y: 0 } };
  for (const p of [{ x: 2, y: .1 }, { x: 2, y: -.1 }, { x: .1, y: .1 }, { x: .1, y: -.1 }]) {
    const h = 1e-6, dx = (basePanelStreamfunctionBasis({ ...p, x: p.x + h }, b).source - basePanelStreamfunctionBasis({ ...p, x: p.x - h }, b).source) / (2 * h);
    const dy = (basePanelStreamfunctionBasis({ ...p, y: p.y + h }, b).source - basePanelStreamfunctionBasis({ ...p, y: p.y - h }, b).source) / (2 * h);
    const q = sourceVelocity(p, b); assert.equal(basePanelStreamfunctionBasis(p, b).onSourceCut, false);
    close(dx, -q.v, 2e-10); close(dy, q.u, 2e-10);
  }
});

test('source charts preserve the native gauge and handle retained NLR base coordinate jitter without a flow solve', () => {
  const inputs = [baseline(), ...getBenchmarkAirfoil('nlr7301').elements], before = structuredClone(inputs);
  inputs.forEach((e, index) => {
    const topology = createContourTopology(e.points, { trailingEdge: e.trailingEdge }), d = xfoilSurfaceDerivatives(topology.surface.points).derivatives;
    const t = { x: .5 * (-d[0].x + d.at(-1).x), y: .5 * (-d[0].y + d.at(-1).y) }, length = Math.hypot(t.x, t.y);
    t.x /= length; t.y /= length;
    const chart = finiteBaseSourceChart(topology, t);
    if (index === 0) { assert.equal(chart.upstreamShift, 0); assert.deepEqual(chart.origin, topology.trailingEdge.center); }
    else assert.ok(chart.upstreamShift > 0 && chart.upstreamShift < .05 * topology.trailingEdge.gapLength);
    close((chart.origin.x - chart.originalCenter.x) * t.y - (chart.origin.y - chart.originalCenter.y) * t.x, 0, 2e-16);
    // Every original segment is still used, and each source gradient in
    // nearby fluid agrees with its independent analytic velocity.
    for (const panel of topology.base.panels) {
      const p = { ...makePanel(panel.start, panel.end), cutDirection: t, cutOrigin: chart.origin };
      const q = { x: chart.originalCenter.x + .1 * t.x - .0002 * t.y, y: chart.originalCenter.y + .1 * t.y + .0002 * t.x }, h = 1e-6;
      const dx = (basePanelStreamfunctionBasis({ ...q, x: q.x + h }, p).source - basePanelStreamfunctionBasis({ ...q, x: q.x - h }, p).source) / (2 * h);
      const dy = (basePanelStreamfunctionBasis({ ...q, y: q.y + h }, p).source - basePanelStreamfunctionBasis({ ...q, y: q.y - h }, p).source) / (2 * h);
      const velocity = sourceVelocity(q, p);
      close(dx, -velocity.v, 2e-10); close(dy, velocity.u, 2e-10);
    }
  });
  assert.deepEqual(inputs, before);
});

test('finite and sharp elements have an explicit square system with independent constants and all cross-influences', () => {
  const elements = [baseline(), { points: transform(naca4('0012', 40), { chord: .25, x: 1.05, y: -.1 }) }];
  const r = solveInviscid({ elements, alpha: 4 }), reverse = solveInviscid({ elements: elements.slice().reverse(), alpha: 4 });
  assert.equal(r.status, 'solved'); assert.equal(r.diagnostics.unknowns, 84);
  assert.equal(r.diagnostics.surfaceRows + r.diagnostics.teProbes.length + r.diagnostics.kuttaRows, 84);
  assert.equal(r.diagnostics.teProbes.length, 1);
  const p = r.diagnostics.teProbes[0], q = velocityAt(p.point, r.field);
  close(q.u * p.tangent.x + q.v * p.tangent.y, 0);
  close(r.cl, reverse.cl); close(r.cm, reverse.cm);
  assert.ok(Math.abs(r.diagnostics.surfaceStreamfunctions[0] - r.diagnostics.surfaceStreamfunctions[1]) > .01);
});

test('sharp-only complete results remain exact to the archived implementation for both original modes', async () => {
  const source = fs.readFileSync('docs/nlr-finite-base/before/panel/linear-vortex.js.txt', 'utf8')
    .replace(/from '(\.[^']+)'/g, (_, p) => `from '${pathToFileURL(resolve('src/inviscid', p)).href}'`);
  const old = await import('data:text/javascript;base64,' + Buffer.from(source).toString('base64'));
  for (const boundaryCondition of ['normal-velocity', 'streamfunction']) {
    const input = { elements: [{ points: naca4('2412', 40) }, { points: transform(naca4('0012', 40), { chord: .3, x: .94, y: -.08, angle: -15 }) }], alpha: 4, boundaryCondition };
    assert.deepEqual(solveInviscid(input), old.solveInviscid(input));
  }
});
