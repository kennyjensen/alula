import test from 'node:test';
import assert from 'node:assert/strict';
import { solveInviscid, velocityAt, vortexBasis } from '../src/inviscid/linear-vortex.js';
import { solveInviscid as solveSources, sourceVelocity, makePanel } from '../src/inviscid/panel.js';
import { naca4, transform } from '../src/geometry/airfoil.js';
import { joukowski } from './fixtures/joukowski.js';

test('source and linear-vortex kernels match direct off-body numerical integration', () => {
  const panel = makePanel({ x: 0.2, y: -0.1 }, { x: 1.1, y: 0.4 });
  const point = { x: -0.3, y: 0.9 };
  const q = [{ u: 0, v: 0 }, { u: 0, v: 0 }]; const source = { u: 0, v: 0 };
  const samples = 10000;
  for (let i = 0; i < samples; i++) {
    const t = (i + 0.5) / samples;
    const dx = point.x - (panel.a.x + t * (panel.b.x - panel.a.x));
    const dy = point.y - (panel.a.y + t * (panel.b.y - panel.a.y));
    const weight = panel.length / samples / (2 * Math.PI * (dx * dx + dy * dy));
    source.u += dx * weight; source.v += dy * weight;
    [1 - t, t].forEach((shape, k) => { q[k].u -= dy * shape * weight; q[k].v += dx * shape * weight; });
  }
  const actual = vortexBasis(point, panel);
  for (let k = 0; k < 2; k++) for (const component of ['u', 'v']) assert.ok(Math.abs(q[k][component] - actual[k][component]) < 1e-9);
  const actualSource = sourceVelocity(point, panel);
  assert.ok(Math.abs(source.u - actualSource.u) < 1e-9); assert.ok(Math.abs(source.v - actualSource.v) < 1e-9);
});
test('nonlifting cylinder matches analytic Cp and zero net force', () => {
  const n = 80;
  const points = Array.from({ length: n + 1 }, (_, i) => ({ x: Math.cos(2 * Math.PI * i / n), y: Math.sin(2 * Math.PI * i / n) }));
  const r = solveSources({ elements: [{ points, lifting: false }] });
  for (const p of r.elements[0].cp) assert.ok(Math.abs(p.cp - (1 - 4 * Math.sin(Math.atan2(p.y, p.x)) ** 2)) < 1e-11);
  assert.ok(Math.abs(r.cl) < 1e-12); assert.ok(Math.abs(r.diagnostics.pressureDrag) < 1e-12);
});
test('lifting Joukowski airfoil converges toward exact Cp, lift and zero drag', () => {
  let previous = { lift: Infinity, cp: Infinity, drag: Infinity };
  for (const n of [80, 160, 320]) {
    const exact = joukowski(n);
    const r = solveInviscid({ elements: [exact], alpha: 4 });
    let squared = 0; let weight = 0;
    r.elements[0].cp.forEach((p, i) => { squared += (p.cp - exact.cpAtAngle(2 * Math.PI * (i + 0.5) / n)) ** 2 * p.length; weight += p.length; });
    const errors = { lift: Math.abs(r.cl - exact.cl), cp: Math.sqrt(squared / weight), drag: Math.abs(r.diagnostics.pressureDrag) };
    for (const key of ['lift', 'cp', 'drag']) assert.ok(errors[key] < previous[key], `${key} must improve at ${n} panels`);
    assert.ok(r.diagnostics.normalVelocityResidual < 1e-10); assert.ok(r.diagnostics.kuttaResidual < 1e-10);
    previous = errors;
  }
  assert.ok(previous.lift < 0.002); assert.ok(previous.cp < 0.01); assert.ok(previous.drag < 0.0011);
});
test('NACA symmetry, incidence signs and reference transformations are consistent', () => {
  const points = naca4('0012', 120);
  const solve = (alpha, p = points, rest = {}) => solveInviscid({ elements: [{ points: p }], alpha, ...rest });
  const zero = solve(0); const positive = solve(4); const negative = solve(-4);
  assert.ok(Math.abs(zero.cl) < 1e-10);
  assert.ok(positive.cl > 0.4 && positive.cl < 0.6);
  assert.ok(Math.abs(positive.cl + negative.cl) < 1e-10);
  const moved = solve(4, transform(points, { chord: 3, x: 5, y: -2 }), { referenceChord: 3, momentReference: { x: 5.75, y: -2 } });
  assert.ok(Math.abs(moved.cl - positive.cl) < 1e-9); assert.ok(Math.abs(moved.cm - positive.cm) < 1e-9);
  const rotated = solve(17, transform(points, { angle: 13 }), { momentReference: transform([{ x: 0.25, y: 0 }], { angle: 13 })[0] });
  assert.ok(Math.abs(rotated.cl - positive.cl) < 1e-9); assert.ok(Math.abs(rotated.cm - positive.cm) < 1e-9);
});
test('multielement coupling is simultaneous, order-independent, and vanishes at large distance', () => {
  const main = { points: naca4('2412', 120) };
  const flap = { points: transform(naca4('0012', 80), { chord: 0.3, angle: -15, x: 0.94, y: -0.08 }) };
  const both = solveInviscid({ elements: [main, flap], alpha: 4 });
  const reverse = solveInviscid({ elements: [flap, main], alpha: 4 });
  const single = solveInviscid({ elements: [main], alpha: 4 });
  const isolatedFlap = solveInviscid({ elements: [flap], alpha: 4 });
  assert.ok(Math.abs(both.cl - reverse.cl) < 1e-9);
  assert.ok(Math.abs(both.cl - single.cl - isolatedFlap.cl) > 0.05);
  assert.ok(both.diagnostics.normalVelocityResidual < 1e-10); assert.ok(both.diagnostics.kuttaResidual < 1e-10);
  const far = solveInviscid({ elements: [main, { points: transform(flap.points, { y: 100 }) }], alpha: 4 });
  assert.ok(Math.abs(far.cl - single.cl - isolatedFlap.cl) < 0.002);
  const speed = velocityAt({ x: 1000, y: 1000 }, both.field);
  assert.ok(Math.abs(speed.u - Math.cos(4 * Math.PI / 180)) < 0.001);
  assert.ok(Math.abs(speed.v - Math.sin(4 * Math.PI / 180)) < 0.001);
});
test('unsupported physics and invalid inputs cannot silently produce results', () => {
  const elements = [{ points: naca4() }];
  assert.throws(() => solveInviscid({ elements, mach: 0.3 }), /Only incompressible/);
  assert.throws(() => solveInviscid({ elements, viscous: true }), /Only incompressible/);
  assert.throws(() => solveInviscid({ elements, alpha: NaN }), /Invalid/);
  assert.throws(() => solveInviscid({ elements: [...elements, ...elements] }), /intersect/);
  assert.equal(solveInviscid({ elements }).cd, null);
});
