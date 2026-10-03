import test from 'node:test';
import assert from 'node:assert/strict';
import { solveInviscid, velocityAt } from '../src/inviscid/linear-vortex.js';
import { streamfunctionAt } from '../src/inviscid/streamfunction.js';
import { naca4, transform } from '../src/geometry/airfoil.js';
import { joukowski } from './fixtures/joukowski.js';

const solve = options => solveInviscid({ ...options, boundaryCondition: 'streamfunction' });
const close = (a, b, tolerance = 2e-10) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);

test('nodal streamfunction equations impose independent constants and sharp-TE conditions on every element', () => {
  const elements = [{ points: naca4('2412', 80) },
    { points: transform(naca4('0012', 60), { chord: .3, x: .94, y: -.08, angle: -15 }) }];
  const r = solve({ elements, alpha: 4 }); assert.equal(r.status, 'solved');
  assert.equal(r.field.gamma.length, 142, 'streamfunction constants must not become vortex strengths');
  assert.ok(Math.abs(r.diagnostics.surfaceStreamfunctions[0] - r.diagnostics.surfaceStreamfunctions[1]) > .01);
  r.elements.forEach((element, e) => {
    for (const p of element.points) close(streamfunctionAt(p, r.field), r.diagnostics.surfaceStreamfunctions[e], 2e-12);
    const { point, tangent } = r.diagnostics.teProbes[e], q = velocityAt(point, r.field);
    close(q.u * tangent.x + q.v * tangent.y, 0, 2e-12);
    const panels = r.field.panels.filter(p => p.element === e);
    close(r.field.gamma[panels[0].node] + r.field.gamma[panels.at(-1).node + 1], 0, 2e-12);
  });
  assert.equal(r.diagnostics.normalVelocityResidual, undefined, 'mid-panel no-penetration is not the imposed condition');
});

test('streamfunction panel solution refines toward independent Joukowski pressure, circulation and forces', () => {
  let previous = { lift: Infinity, pressure: Infinity, drag: Infinity, circulation: Infinity };
  for (const n of [80, 160, 320]) {
    const exact = joukowski(n), r = solve({ elements: [exact], alpha: 4 });
    let square = 0, length = 0;
    r.elements[0].cp.forEach((p, i) => { square += p.length * (p.cp - exact.cpAtAngle(2 * Math.PI * (i + .5) / n)) ** 2; length += p.length; });
    const error = { lift: Math.abs(r.cl - exact.cl), pressure: Math.sqrt(square / length), drag: Math.abs(r.diagnostics.pressureDrag),
      circulation: Math.abs(r.diagnostics.circulationLift - exact.cl) };
    for (const key of Object.keys(error)) assert.ok(error[key] < .6 * previous[key], `${key} failed refinement at ${n}: ${error[key]}`);
    previous = error;
  }
  assert.ok(previous.lift < .002 && previous.pressure < .01 && previous.drag < .0011 && previous.circulation < 1e-5);
});

test('streamfunction solution preserves symmetry, units, rotation and multielement permutation', () => {
  const points = naca4('0012', 80), positive = solve({ elements: [{ points }], alpha: 4 });
  close(solve({ elements: [{ points }], alpha: -4 }).cl, -positive.cl);
  close(solve({ elements: [{ points }], alpha: 0 }).cl, 0);
  const movement = { chord: 3, x: 5, y: -2, angle: 13 };
  const moved = solve({ elements: [{ points: transform(points, movement) }], alpha: 17, referenceChord: 3,
    momentReference: transform([{ x: .25, y: 0 }], movement)[0] });
  close(moved.cl, positive.cl); close(moved.cm, positive.cm);
  const elements = [{ points }, { points: transform(naca4('0012', 60), { chord: .3, x: -.4, y: .2 }) }];
  const both = solve({ elements, alpha: 2 }), reverse = solve({ elements: elements.toReversed(), alpha: 2 });
  close(both.cl, reverse.cl); close(both.cm, reverse.cm);
  both.elements.forEach((e, i) => close(e.circulation, reverse.elements[1 - i].circulation));
  const isolated = elements.map(element => solve({ elements: [element], alpha: 2 }));
  assert.ok(Math.abs(both.cl - isolated.reduce((sum, r) => sum + r.cl, 0)) > .01);
  const far = solve({ elements: [elements[0], { points: transform(elements[1].points, { y: 100 }) }], alpha: 2 });
  close(far.cl, isolated.reduce((sum, r) => sum + r.cl, 0), .002);
  assert.throws(() => solveInviscid({ elements, boundaryCondition: 'unknown' }), /boundary condition/);
});
