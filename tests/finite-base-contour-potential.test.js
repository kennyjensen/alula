import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { createContourPotential } from '../src/inviscid/contour-potential.js';
import { refineSurfaceStagnation } from '../src/inviscid/surface-stagnation.js';
import { createSurfaceContourCurve } from '../src/geometry/contour-topology.js';
import { potentialDifference } from '../src/inviscid/streamfunction.js';
import { velocityAt } from '../src/inviscid/linear-vortex.js';
import { makePanel } from '../src/inviscid/panel.js';

const saved = JSON.parse(fs.readFileSync(new URL('../docs/nlr-finite-base/nlr7301-panel-final.json', import.meta.url)));
const diagnostic = JSON.parse(fs.readFileSync(new URL('../docs/nlr-finite-base/contour-field/diagnostic.json', import.meta.url)));
const close = (a, b, tolerance = 2e-12) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);
const exterior = { length: 1, knots: [0, .3, .7, 1],
  evaluate: s => ({ point: { x: -.2 + 1.6 * s, y: .5 + .02 * s }, derivative: { x: 1.6, y: .02 } }) };

test('retained NLR full-field contour phase and directional velocity agree, including both finite bases', t => {
  const field = structuredClone(saved.result.field), before = structuredClone(field), s0 = .43;
  const branch = createContourPotential({ curve: exterior, field, stagnation: s0, stagnationPotential: 2 });
  let potentialError = 0, velocityError = 0;
  for (const s of [.13, .51, .77, .93]) {
    const p = exterior.evaluate(s), direct = potentialDifference(exterior.evaluate(s0).point, p.point, field);
    potentialError = Math.max(potentialError, Math.abs(branch.phase(s) - 2 - direct));
    const h = 1e-5, derivative = (-branch.phase(s + 2 * h) + 8 * branch.phase(s + h)
      - 8 * branch.phase(s - h) + branch.phase(s - 2 * h)) / (12 * h), q = velocityAt(p.point, field);
    velocityError = Math.max(velocityError, Math.abs(derivative - q.u * p.derivative.x - q.v * p.derivative.y));
  }
  assert.ok(potentialError < 3e-12); assert.ok(velocityError < 1e-9);
  assert.deepEqual(field, before);
  const expected = branch.phase(.5291);
  field.gamma.fill(8); field.panels[0].a.x = -9;
  for (const p of field.basePanels) {
    p.a.x -= .7; p.b.y += .3; p.sourceStrength *= 7; p.vortexStrength *= 11;
    p.cutDirection.x = -99; p.cutOrigin.y = 99;
  }
  assert.equal(branch.phase(.5291), expected, 'uncached queries must retain the same complete copied field as cached prefixes');
  t.diagnostic(JSON.stringify({ potentialError, velocityError, retainedBasePanels: before.basePanels.length }));
});

test('finite-base closed-body membership rejects interior paths and both actual incompatible NLR C2 contours', () => {
  const field = saved.result.field;
  // This small path lies wholly inside the main element and does not cross
  // a wetted sheet. Omitting the base incorrectly made that body an open sheet.
  const inside = { length: 1, knots: [0, .5, 1], evaluate: s => ({ point: { x: .4 + .01 * s, y: .01 }, derivative: { x: .01, y: 0 } }) };
  assert.throws(() => createContourPotential({ curve: inside, field, stagnation: .4, stagnationPotential: 0 }),
    e => e.code === 'CONTOUR_POTENTIAL_PATH_INCOMPATIBLE' && /inside a panel body/.test(e.message));
  saved.input.elements.forEach((body, b) => {
    const curve = createSurfaceContourCurve(body.points, body), proof = diagnostic.profiles[b];
    const root = refineSurfaceStagnation({ curve, field, ...proof.roots[0], chord: 1 });
    assert.equal(root.parameter, proof.stagnation.parameter);
    assert.throws(() => createContourPotential({ curve, field, stagnation: root.parameter, stagnationPotential: 0 }), e => {
      assert.equal(e.code, 'CONTOUR_POTENTIAL_PATH_INCOMPATIBLE'); assert.match(e.message, /inside a panel body/);
      assert.deepEqual(e.interval, proof.current.interval); return true;
    });
  });
});

test('base sheets have the same explicit boundary limits and complete-chain validation as wetted sheets', () => {
  const field = structuredClone(saved.result.field), p = field.basePanels[0];
  const center = { x: .5 * (p.a.x + p.b.x), y: .5 * (p.a.y + p.b.y) };
  const curve = { length: 1, knots: [0, 1], evaluate: () => ({ point: center, derivative: { x: 0, y: 1 } }) };
  assert.throws(() => createContourPotential({ curve, field, stagnation: .5, stagnationPotential: 0 }),
    e => e.code === 'CONTOUR_POTENTIAL_PATH_INCOMPATIBLE' && /one-sided limit/.test(e.message));
  const missing = structuredClone(field); missing.basePanels.pop();
  assert.throws(() => createContourPotential({ curve: exterior, field: missing, stagnation: .4, stagnationPotential: 0 }), /complete ordered closed body/);
  const invalid = structuredClone(field); invalid.basePanels[0].sourceStrength = NaN;
  assert.throws(() => createContourPotential({ curve: exterior, field: invalid, stagnation: .4, stagnationPotential: 0 }), /Invalid.*finite-base field/);
});

test('sharp-field phase values and diagnostics retain exact archived behavior', async () => {
  const source = fs.readFileSync(new URL('../docs/nlr-finite-base/contour-field/contour-potential.js.before.txt', import.meta.url), 'utf8');
  const archived = await import('data:text/javascript;base64,' + Buffer.from(source.replace(/from\s+(['"])(\.[^'"]+)\1/g,
    (_, quote, p) => `from ${quote}${pathToFileURL(path.resolve('src/inviscid', p)).href}${quote}`)).toString('base64'));
  const points = [{ x: -1, y: -1 }, { x: 1, y: -1 }, { x: 1, y: 1 }, { x: -1, y: 1 }, { x: -1, y: -1 }];
  const field = { u: .4, v: .2, gamma: Float64Array.of(1, 1, 1, 1, 1),
    panels: points.slice(0, -1).map((p, i) => ({ ...makePanel(p, points[i + 1], 0), node: i })) };
  const curve = { length: 2 * Math.PI, knots: [0, Math.PI / 2, Math.PI, 3 * Math.PI / 2, 2 * Math.PI],
    evaluate: s => ({ point: { x: 3 * Math.cos(s), y: 3 * Math.sin(s) }, derivative: { x: -3 * Math.sin(s), y: 3 * Math.cos(s) } }) };
  const input = { curve, field, stagnation: Math.PI, stagnationPotential: 2 };
  const a = createContourPotential(input), b = archived.createContourPotential(input);
  assert.deepEqual(a.diagnostics, b.diagnostics);
  for (const s of [0, .121, .79, 1.111, Math.PI, 4.333, curve.length]) assert.equal(a.phase(s), b.phase(s));
});
