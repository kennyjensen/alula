import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createContourPotential } from '../src/inviscid/contour-potential.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { solveInviscid } from '../src/inviscid/linear-vortex.js';
import { makePanel } from '../src/inviscid/panel.js';
import { potentialDifference } from '../src/inviscid/streamfunction.js';

const close = (a, b, tolerance = 2e-12) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b} (limit ${tolerance})`);
const uniform = { u: 2, v: -.3, panels: [], gamma: [] };
const parabola = { length: 1, knots: [0, .2, .7, 1],
  evaluate: s => ({ point: { x: s, y: s * s }, derivative: { x: 1, y: 2 * s } }) };
const squareField = () => {
  const points = [{ x: -1, y: -1 }, { x: 1, y: -1 }, { x: 1, y: 1 }, { x: -1, y: 1 }, { x: -1, y: -1 }];
  return { u: .4, v: .2, gamma: Float64Array.of(1, 1, 1, 1, 1),
    panels: points.slice(0, -1).map((p, i) => ({ ...makePanel(p, points[i + 1], 0), node: i })) };
};

test('cached and arbitrary contour phases reproduce a uniform analytic field and preserve the supplied origin', () => {
  const stagnation = .413, origin = 7.125, before = structuredClone(uniform);
  const result = createContourPotential({ curve: parabola, stagnation, stagnationPotential: origin, field: uniform });
  assert.equal(result.phase(stagnation), origin);
  for (const s of [0, .015, .2, .412999, .413001, .7, .987, 1])
    close(result.phase(s), origin + 2 * (s - stagnation) - .3 * (s * s - stagnation * stagnation));
  assert.deepEqual(uniform, before);
  assert.equal(result.diagnostics.sampleCount, 9);
  assert.equal(result.diagnostics.adaptiveSubdivisions, 0);
  assert.equal(result.diagnostics.sampledMonotone, false, 'A supplied point need not be a stagnation point of the analytic field.');
  assert.throws(() => result.phase(-.1), /outside/);
  assert.throws(() => result.phase(NaN), /outside/);
});

test('analytic contour branches retain circulation instead of forcing a single-valued TE potential', () => {
  const field = squareField(), curve = { length: 2 * Math.PI, knots: [0, Math.PI / 2, Math.PI, 3 * Math.PI / 2, 2 * Math.PI],
    evaluate: s => ({ point: { x: 3 * Math.cos(s), y: 3 * Math.sin(s) }, derivative: { x: -3 * Math.sin(s), y: 3 * Math.cos(s) } }) };
  const args = { curve, field, stagnation: Math.PI, stagnationPotential: 2 };
  const two = createContourPotential({ ...args, subdivisions: 2 }), four = createContourPotential({ ...args, subdivisions: 4 });
  close(two.phase(curve.length) - two.phase(0), 8);
  close(two.diagnostics.contourPotentialIncrement, 8);
  for (const s of [0, .121, 1.111, Math.PI, 4.333, curve.length]) close(two.phase(s), four.phase(s));
  const expected = two.phase(.79);
  field.gamma.fill(100); field.u = 99; field.panels[0].a.x = -9;
  close(two.phase(.79), expected, 1e-14, 'Cached and uncached values use the same copied panel field.');
});

test('exact oblique panel endpoints have continuous exterior potential limits under rotation', () => {
  for (const angle of [0, .41, 1.17, 2.36, 3.8, 5.1]) {
    const rotate = p => ({ x: p.x * Math.cos(angle) - p.y * Math.sin(angle),
      y: p.x * Math.sin(angle) + p.y * Math.cos(angle) });
    const a = rotate({ x: 1, y: 0 }), b = rotate({ x: 0, y: 1 });
    const panel = { ...makePanel(a, b, 0), node: 0 }, field = { panels: [panel], gamma: [1, .7], u: 0, v: 0 };
    for (const endpoint of [a, b]) {
      let previous = Infinity;
      for (const h of [1e-3, 1e-5, 1e-7]) {
        const exterior = { x: endpoint.x + h * panel.nx, y: endpoint.y + h * panel.ny };
        const delta = potentialDifference(exterior, endpoint, field);
        assert.ok(Math.abs(delta) < previous / 20, `Endpoint limit jumps at angle ${angle}: ${delta}`);
        close(delta, -potentialDifference(endpoint, exterior, field));
        previous = Math.abs(delta);
      }
      assert.ok(previous < 5e-7);
    }
  }
});

test('short-path subdivision resolves chord visibility and rejects a curve inside the solid', () => {
  const field = squareField(), curve = { length: 1, knots: [0, 1],
    evaluate: s => ({ point: { x: -2 + 4 * s, y: 8 * s * (1 - s) }, derivative: { x: 4, y: 8 - 16 * s } }) };
  const result = createContourPotential({ curve, field, stagnation: .125, stagnationPotential: 0, subdivisions: 1 });
  assert.ok(result.diagnostics.adaptiveSubdivisions > 0);
  const a = curve.evaluate(0).point, b = curve.evaluate(1).point, top = { x: 0, y: 3 };
  close(result.phase(1) - result.phase(0), potentialDifference(a, top, field) + potentialDifference(top, b, field));
  const incompatible = { length: 1, knots: [0, .5, 1],
    evaluate: s => ({ point: { x: -2 + 4 * s, y: 0 }, derivative: { x: 4, y: 0 } }) };
  assert.throws(() => createContourPotential({ curve: incompatible, field, stagnation: .1, stagnationPotential: 0 }),
    error => error.code === 'CONTOUR_POTENTIAL_PATH_INCOMPATIBLE' && /inside/.test(error.message));
  const open = { u: 0, v: 0, gamma: [1, 1],
    panels: [{ ...makePanel({ x: 0, y: -1 }, { x: 0, y: 1 }, 0), node: 0 }] };
  assert.throws(() => createContourPotential({ curve: incompatible, field: open, stagnation: .1, stagnationPotential: 0 }),
    error => error.code === 'CONTOUR_POTENTIAL_PATH_INCOMPATIBLE' && /one-sided limit/.test(error.message));
});

test('contour potential rejects malformed inputs without mutating the curve', () => {
  const args = { curve: parabola, field: uniform, stagnation: .4, stagnationPotential: 0 };
  for (const subdivisions of [0, 1.5, 33, Infinity]) assert.throws(() => createContourPotential({ ...args, subdivisions }), /Invalid/);
  for (const stagnation of [0, 1, NaN]) assert.throws(() => createContourPotential({ ...args, stagnation }), /Invalid/);
  assert.throws(() => createContourPotential({ ...args, curve: { ...parabola, knots: [0, .7, .2, 1] } }), /knots/);
  assert.throws(() => createContourPotential({ ...args, field: { ...uniform, u: NaN } }), /field/);
  const knots = [...parabola.knots];
  createContourPotential(args);
  assert.deepEqual(parabola.knots, knots);
});

test('all four default surface branches agree under subdivision and recover both TE circulation jumps', t => {
  // One panel solve only: no streamline tracing, grid preparation, or SLOR.
  const fixture = JSON.parse(fs.readFileSync(new URL('../docs/mset-block-skeleton-fixture.json', import.meta.url)));
  const panel = solveInviscid({ elements: fixture.input.bodies.map(b => ({ points: b.points })), alpha: fixture.input.alpha,
    boundaryCondition: fixture.diagnostics.panelBoundaryCondition });
  const reports = [], results = [];
  fixture.input.bodies.forEach((body, b) => {
    const curve = createContourCurve(body.points), origin = fixture.diagnostics.profiles[b].stagnationPotential;
    const args = { curve, field: panel.field, stagnation: body.stagnationParameter, stagnationPotential: origin };
    const two = createContourPotential({ ...args, subdivisions: 2 }), four = createContourPotential({ ...args, subdivisions: 4 });
    assert.equal(two.phase(body.stagnationParameter), origin);
    let maximumDifference = 0;
    for (const side of ['upper', 'lower']) for (const f of [0, .0001, .07, .19, .53, .93, 1]) {
      const s = curve.branch(side, f, body.stagnationParameter).parameter;
      maximumDifference = Math.max(maximumDifference, Math.abs(two.phase(s) - four.phase(s)));
    }
    assert.ok(maximumDifference < 2e-10);
    const loop = two.phase(curve.length) - two.phase(0);
    close(loop, panel.elements[b].circulation, 2e-10);
    const monotonicity = two.diagnostics.nearStagnation;
    assert.ok(Number.isFinite(monotonicity.derivative));
    assert.equal(monotonicity.sampledMonotone, true,
      'The prepared stagnation point must be the zero of the analytic C2-surface tangential derivative.');
    assert.ok(Math.abs(monotonicity.derivative) < 1e-9);
    assert.equal(two.diagnostics.sampledMonotone, true);
    assert.equal(two.diagnostics.branches.upper.negativeIncrements, 0);
    assert.equal(two.diagnostics.branches.lower.negativeIncrements, 0);
    reports.push({ body: b, maximumDifference, loop, circulation: panel.elements[b].circulation,
      nearStagnation: monotonicity, branches: two.diagnostics.branches });
    results.push({ curve, potential: two });
  });
  // Independent direct passage path. The saved scalar endpoint reference
  // predates the exact-endpoint branch correction, so recompute from the
  // current analytic field instead of preserving that known defect.
  const gap = results[1].potential.phase(results[1].curve.length) - results[0].potential.phase(0);
  const direct = potentialDifference(results[0].curve.evaluate(0).point,
    results[1].curve.evaluate(results[1].curve.length).point, panel.field);
  t.diagnostic(JSON.stringify({ gap, direct, reports }));
  close(gap, direct, 2e-10);
});
