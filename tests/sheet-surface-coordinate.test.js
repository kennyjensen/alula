// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import fs from 'node:fs';
import { createSheetSurfaceCoordinate } from '../src/inviscid/sheet-surface-coordinate.js';
import { naca4, pointInside, transform } from '../src/geometry/airfoil.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { getBenchmarkAirfoil } from '../src/geometry/benchmark-airfoils.js';
import { makePanel } from '../src/inviscid/panel.js';
import { createContourPotential } from '../src/inviscid/contour-potential.js';

const close = (actual, expected, tolerance = 128 * Number.EPSILON * Math.max(1, Math.abs(expected))) =>
  assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} != ${expected}; tolerance ${tolerance}`);
function fixture({ points = naca4('0012', 20), strength = (s, length) => 3.25 * (s - .46 * length), element = 0, node = 0 } = {}) {
  const curve = createContourCurve(points), knots = curve.knots;
  const field = { u: 1, v: 0, gamma: [...Array(node).fill(99), ...knots.map(s => strength(s, curve.length))],
    panels: points.slice(0, -1).map((p, i) => ({ ...makePanel(p, points[i + 1], element), node: node + i })) };
  return { points, knots, field, element };
}

test('linear manufactured strength has the exact zero, quadratic primitive and independent TE jump', () => {
  const args = fixture(), before = structuredClone(args), length = args.knots.at(-1), root = .46 * length, gauge = 7.125;
  const result = createSheetSurfaceCoordinate({ ...args, stagnationPotential: gauge });
  close(result.stagnationParameter, root); assert.equal(result.phase(result.stagnationParameter), gauge);
  assert.equal(result.derivative(result.stagnationParameter), 0);
  for (const s of [0, .001, .1, .33, .46, .461, .8, 1].map(f => f * length).concat(args.knots)) {
    close(result.derivative(s), 3.25 * (s - root));
    close(result.phase(s), gauge + 3.25 / 2 * (s - root) ** 2);
  }
  // Closed-form circulation of this linear strength, not the implementation's panel sum.
  const circulation = 3.25 / 2 * ((length - root) ** 2 - root ** 2);
  close(result.phase(length) - result.phase(0), circulation);
  close(result.diagnostics.circulation, circulation);
  close(result.diagnostics.stagnationSlope, 3.25);
  assert.equal(result.diagnostics.sampledMonotone, true);
  assert.deepEqual(args, before);
});

test('a nodal zero keeps two exact one-sided quadratics and a declared mean strain slope', () => {
  const points = naca4('0012', 20), knots = createContourCurve(points).knots, root = knots[8];
  const args = fixture({ points, strength: s => (s <= root ? 4 : 2) * (s - root) });
  const result = createSheetSurfaceCoordinate(args);
  assert.equal(result.stagnationParameter, root);
  assert.equal(result.diagnostics.stagnation.kind, 'nodal-zero');
  close(result.diagnostics.stagnationSlopes.upper, 4); close(result.diagnostics.stagnationSlopes.lower, 2);
  close(result.diagnostics.stagnationSlope, 3);
  for (const delta of [-.4, -.000001, -1e-10, 0, 1e-10, .000001, .4]) {
    const s = root + delta, factor = delta <= 0 ? 4 : 2;
    close(result.phase(s), factor / 2 * (s - root) ** 2, Math.max(1e-29, 64 * Number.EPSILON * (s - root) ** 2));
    close(result.derivative(s), factor * (s - root));
  }
});

test('coordinate obeys gauge, geometry-unit, rigid-transform and strength sensitivity identities', () => {
  const args = fixture(), base = createSheetSurfaceCoordinate(args), length = base.length;
  const scaled = fixture({ points: transform(args.points, { chord: 7, angle: 31, x: 2, y: -.3 }),
    strength: (s, l) => 3.25 * length * (s / l - .46) });
  const moved = createSheetSurfaceCoordinate(scaled);
  for (const f of [0, .17, .46, .69, 1]) {
    close(moved.derivative(f * moved.length), base.derivative(f * length));
    close(moved.phase(f * moved.length), 7 * base.phase(f * length), 2e-13);
  }
  close(moved.stagnationParameter, 7 * base.stagnationParameter);
  const perturbation = .13, slope = 3.25;
  const changed = createSheetSurfaceCoordinate(fixture({ strength: (s, l) => slope * (s - .46 * l) + perturbation }));
  close(changed.stagnationParameter, base.stagnationParameter - perturbation / slope);
  const changedRoot = .46 * length - perturbation / slope;
  close(changed.phase(.8 * length), slope / 2 * (.8 * length - changedRoot) ** 2);
  const doubled = createSheetSurfaceCoordinate(fixture({ strength: (s, l) => 2 * slope * (s - .46 * l) }));
  close(doubled.phase(.2 * length), 2 * base.phase(.2 * length));
  const gauge = createSheetSurfaceCoordinate({ ...args, stagnationPotential: -12 });
  close(gauge.phase(.7 * length), base.phase(.7 * length) - 12);
});

test('selected element, independent TE strengths and returned closures remain isolated from caller mutation', () => {
  const args = fixture({ element: 3, node: 7 }), result = createSheetSurfaceCoordinate(args);
  const s = .31 * result.length, oldPhase = result.phase(s), oldDerivative = result.derivative(s);
  assert.notEqual(args.field.gamma[7], args.field.gamma.at(-1), 'the repeated geometric TE has independent sheet strengths');
  args.points[1].x += 99; args.knots[1] = 99; args.field.gamma.fill(99); args.field.panels[0].node = 99;
  result.diagnostics.stagnation.parameter = 99;
  assert.equal(result.phase(s), oldPhase); assert.equal(result.derivative(s), oldDerivative);
});

test('malformed geometry, panel correspondence, controls and unresolved or nonmonotone strength roots reject', () => {
  const failures = [
    a => { a.points = a.points.slice().reverse(); },
    a => { a.points.at(-1).x += 1e-12; },
    a => { a.knots[0] = .01; }, a => { a.knots[3] += 1e-8; }, a => { a.knots[3] = a.knots[2]; },
    a => { a.field.panels[2].node++; }, a => { a.field.panels.reverse(); },
    a => { a.field.panels[2].tx *= -1; }, a => { a.field.panels[2].length *= 1.01; },
    a => { a.field.panels[2].b = { x: 0, y: 0 }; }, a => { a.field.gamma[2] = NaN; },
    a => { a.field.gamma.pop(); }, a => { a.element = 1; }, a => { a.stagnationPotential = Infinity; },
    a => { a.field.gamma.fill(1); }, a => { a.field.gamma.fill(-1); },
    a => { a.field.gamma[7] = a.field.gamma[8] = 0; },
    a => { a.field.gamma[2] = 1; a.field.gamma[3] = -1; }, // second incoming zero
    a => { a.field.gamma[0] = 1; }, // one incoming zero but an invalid upper branch
  ];
  for (const change of failures) { const args = fixture(); change(args); assert.throws(() => createSheetSurfaceCoordinate(args)); }
  const good = createSheetSurfaceCoordinate(fixture());
  for (const s of [-1, good.length + 1, NaN, Infinity]) {
    assert.throws(() => good.phase(s), /outside/); assert.throws(() => good.derivative(s), /outside/);
  }
});

test('RAE reference fields reject reversed normal-velocity TE strengths and accept the source-consistent streamfunction BIE', t => {
  const preset = getBenchmarkAirfoil('rae2822'), points = preset.elements[0].points, before = structuredClone(points);
  const sourceHash = createHash('sha256').update(fs.readFileSync(preset.provenance.sourcePath)).digest('hex');
  assert.equal(sourceHash, preset.provenance.sourceHash);
  // Retain the actual failed normal-velocity field as a regression. The
  // prescribed two-field producer solved each unchanged formulation once;
  // this test reuses those full fields instead of repeating their LU work.
  const audit = JSON.parse(fs.readFileSync('docs/rae2822/panel-coordinate-audit.json', 'utf8'));
  for (const path of ['src/inviscid/linear-vortex.js', 'src/inviscid/panel.js', 'src/inviscid/streamfunction.js', 'src/geometry/airfoil.js'])
    assert.equal(createHash('sha256').update(fs.readFileSync(path)).digest('hex'), audit.sourceHashes[path]);
  assert.deepEqual(audit.input.points, points);
  const curve = createContourCurve(points);
  const rejected = audit.cases.find(c => c.boundaryCondition === 'normal-velocity');
  assert.equal(rejected.coordinate.accepted, false);
  assert.ok(rejected.endpointGamma[0] > 0 && rejected.endpointGamma[1] < 0);
  assert.throws(() => createSheetSurfaceCoordinate({ points, knots: curve.knots, field: rejected.panel.field }), /not monotone/);
  const panel = audit.cases.find(c => c.boundaryCondition === 'streamfunction').panel;
  assert.equal(panel.status, 'solved');
  const result = createSheetSurfaceCoordinate({ points, knots: curve.knots, field: panel.field, element: 0 });
  let insideCount = 0, maximumDifference = 0, sumSquares = 0;
  const count = panel.panelCount;
  panel.elements[0].cp.forEach((cp, i) => {
    const s = .5 * (curve.knots[i] + curve.knots[i + 1]);
    const difference = result.derivative(s) - cp.qt;
    maximumDifference = Math.max(maximumDifference, Math.abs(difference)); sumSquares += difference * difference;
    if (pointInside(curve.evaluate(s).point, points)) insideCount++;
    assert.ok(Number.isFinite(result.phase(s)));
    assert.ok(s < result.stagnationParameter ? result.derivative(s) <= 0 : result.derivative(s) >= 0);
  });
  assert.ok(insideCount > 0, 'the RAE spline/polygon incompatibility must remain visible');
  assert.throws(() => createContourPotential({ curve, field: panel.field, stagnation: result.stagnationParameter, stagnationPotential: 0 }),
    e => e.code === 'CONTOUR_POTENTIAL_PATH_INCOMPATIBLE' && /inside/.test(e.message));
  close(result.phase(result.length) - result.phase(0), panel.elements[0].circulation, 2e-14);
  assert.deepEqual(points, before);
  t.diagnostic(JSON.stringify({ preset: 'rae2822', sourceHash, panelCount: count, alpha: 4,
    boundaryCondition: 'streamfunction', panelResidual: panel.diagnostics.linearResidual,
    stagnationParameter: result.stagnationParameter, stagnationPoint: curve.evaluate(result.stagnationParameter).point,
    circulation: panel.elements[0].circulation, coordinate: result.diagnostics,
    c2MidpointsInsidePanelPolygon: insideCount,
    gammaVsExteriorCollocationQt: { maximumAbsoluteDifference: maximumDifference, rmsDifference: Math.sqrt(sumSquares / count),
      interpretation: 'Measured approximation discrepancy, not an acceptance tolerance or an exact-field equality claim.' } }));
});
