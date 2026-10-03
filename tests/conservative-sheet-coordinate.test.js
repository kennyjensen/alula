import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { createSurfaceContourCurve } from '../src/geometry/contour-topology.js';
import { createSheetSurfaceCoordinate } from '../src/inviscid/sheet-surface-coordinate.js';

const saved = JSON.parse(fs.readFileSync('docs/solver-reliability/gui-defaults-current/30p30n/recovery-panel-source.json'));
const points = saved.topology.bodies[0].points, curve = createContourCurve(points);
const input = { points, knots: curve.knots, field: { ...saved.panel.field, gamma: Float64Array.from(Object.values(saved.panel.field.gamma)) }, element: 0 };

test('conservative policy leaves already monotone exact sheet coordinates unchanged', () => {
  for (const element of [1, 2]) {
    const body = saved.topology.bodies[element];
    const c = body.trailingEdge?.kind === 'finite-base' ? createSurfaceContourCurve(body.points, body) : createContourCurve(body.points);
    const args = { points: body.points, knots: c.knots, field: input.field, element, trailingEdge: body.trailingEdge };
    const old = createSheetSurfaceCoordinate(args), value = createSheetSurfaceCoordinate({ ...args, reconstruction: 'conservative-panel' });
    assert.deepEqual(value.diagnostics, old.diagnostics);
    assert.equal(value.stagnationParameter, old.stagnationParameter);
    for (const s of [...c.knots, ...c.knots.slice(1).map((z, i) => .5 * (c.knots[i] + z))]) {
      assert.equal(value.phase(s), old.phase(s));
      assert.equal(value.derivative(s), old.derivative(s));
    }
  }
});

test('slat cove coordinate preserves every original panel integral and the true incoming stagnation quadratic', () => {
  const before = JSON.stringify(input);
  assert.throws(() => createSheetSurfaceCoordinate(input), error => error.code === 'SHEET_SURFACE_NONMONOTONE');
  const value = createSheetSurfaceCoordinate({ ...input, reconstruction: 'conservative-panel' });
  assert.deepEqual(value.diagnostics.reconstructedPanels.map(p => p.index), [125, 126]);
  const gamma = input.field.gamma;
  for (let i = 0; i < points.length - 1; i++) {
    const a = curve.knots[i], b = curve.knots[i + 1], expected = (b - a) * (.5 * gamma[i] + .5 * gamma[i + 1]);
    assert.ok(Math.abs(value.phase(b) - value.phase(a) - expected) < 4e-15, `panel ${i}`);
    for (let j = 0; j <= 8; j++) {
      const s = a + j / 8 * (b - a), d = value.derivative(s);
      assert.ok(s < value.stagnationParameter ? d <= 0 : d >= 0, `sign ${i}/${j}`);
    }
  }
  const root = value.stagnationParameter, slope = value.diagnostics.stagnationSlope;
  for (const delta of [-1e-7, -1e-8, 1e-8, 1e-7]) {
    const s = root + delta, ds = s - root;
    assert.equal(value.phase(s), .5 * slope * ds * ds);
    assert.equal(value.derivative(s), slope * ds);
  }
  const shifted = createSheetSurfaceCoordinate({ ...input, reconstruction: 'conservative-panel', stagnationPotential: 2 });
  for (const s of curve.knots) assert.ok(Math.abs(shifted.phase(s) - value.phase(s) - 2) < 1e-15);
  assert.ok(Math.abs(value.diagnostics.contourPotentialIncrement - value.diagnostics.circulation) < 1e-14);
  assert.equal(JSON.stringify(input), before);
});

test('conservative reconstruction rejects resolved reverse flux, invalid TE signs and false controls', () => {
  const clone = () => ({ ...input, field: { ...input.field, gamma: input.field.gamma.slice() } });
  for (const mutate of [a => { a.field.gamma[125] = -1; a.field.gamma[126] = -1; },
    a => { a.field.gamma[0] = 1; }, a => { a.field.gamma[125] = 0; a.field.gamma[126] = 0; }]) {
    const a = clone(); mutate(a); const before = JSON.stringify(a);
    assert.throws(() => createSheetSurfaceCoordinate({ ...a, reconstruction: 'conservative-panel' }));
    assert.equal(JSON.stringify(a), before);
  }
  assert.throws(() => createSheetSurfaceCoordinate({ ...input, reconstruction: 'clip' }), /Unknown/);
});
