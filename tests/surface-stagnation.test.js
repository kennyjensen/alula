import test from 'node:test';
import assert from 'node:assert/strict';
import { makePanel } from '../src/inviscid/panel.js';
import { velocityAt } from '../src/inviscid/linear-vortex.js';
import { refineSurfaceStagnation } from '../src/inviscid/surface-stagnation.js';
import { createContourPotential } from '../src/inviscid/contour-potential.js';

function junction({ scale = 1, angle = 0, middleStrength = 0 } = {}) {
  const rotate = p => ({ x: scale * (p.x * Math.cos(angle) - p.y * Math.sin(angle)),
    y: scale * (p.x * Math.sin(angle) + p.y * Math.cos(angle)) });
  const points = [{ x: -.5, y: .5 }, { x: 0, y: 0 }, { x: .5, y: .5 }].map(rotate);
  // Linear strengths vanish at the central vertex when middleStrength=0.
  // The exact weighted endpoint velocity is finite and the uniform term
  // cancels it, but the existing unit-strength velocity basis is undefined.
  const v = Math.SQRT2 / (2 * Math.PI);
  const field = { u: -v * Math.sin(angle), v: v * Math.cos(angle), gamma: [-1, middleStrength, 1],
    panels: points.slice(0, -1).map((p, i) => ({ ...makePanel(p, points[i + 1], 0), node: i })) };
  const curve = { length: 1, knots: [0, .5, 1], evaluate: s => ({
    point: rotate({ x: s - .5, y: 0 }), derivative: rotate({ x: 1, y: 0 }) }) };
  return { curve, field, left: .2, right: .8, chord: scale, relativeTolerance: 1e-9 };
}

test('an exact panel vertex is returned as a physically narrow potential-minimum bracket without inventing its velocity', () => {
  for (const controls of [{}, { scale: 7, angle: .43 }]) {
    const args = junction(controls), before = structuredClone(args.field);
    assert.throws(() => velocityAt(args.curve.evaluate(.5).point, args.field), /singular/);
    const result = refineSurfaceStagnation(args);
    assert.equal(result.kind, 'potential-minimum-bracket');
    assert.equal(result.parameter, .5); assert.equal(result.derivative, null);
    assert.equal(result.derivativeStatus, 'undefined at exact panel vertex');
    assert.ok(result.bracket.left < .5 && result.bracket.right > .5);
    assert.ok(result.bracket.leftDerivative < 0 && result.bracket.rightDerivative > 0);
    assert.ok(result.bracket.physicalWidth <= args.relativeTolerance * args.chord);
    assert.ok(result.bracket.physicalWidth >= args.chord * (result.bracket.right - result.bracket.left) * (1 - 1e-14));
    assert.ok(result.comparison.vertexBelowEndpoints);
    assert.ok(result.comparison.leftPotential > 0 && result.comparison.rightPotential > 0);
    const branch = createContourPotential({ ...args, stagnation: result.parameter, stagnationPotential: 0 });
    assert.equal(branch.diagnostics.nearStagnation.derivative, null);
    assert.equal(branch.diagnostics.nearStagnation.sampledMonotone, true);
    assert.equal(branch.diagnostics.sampledMonotone, true);
    assert.deepEqual(args.field, before);
  }
});

test('a nonzero corner density can move the minimum away from the vertex and is not snapped to symmetry', () => {
  const args = junction({ middleStrength: .1 }), result = refineSurfaceStagnation(args);
  assert.equal(result.kind, 'velocity-zero');
  assert.ok(Math.abs(result.parameter - .5) > .001);
  assert.ok(Number.isFinite(result.derivative));
  assert.ok(Math.abs(result.derivative) <= result.velocityTolerance);
});

test('a curved bracket is narrowed using an arc bound even when its endpoint distance already meets tolerance', () => {
  const args = junction();
  args.field.v = -1;
  args.curve = { length: 1, knots: [0, .5, 1], evaluate: s => ({
    point: { x: s - .5, y: -((s - .5) ** 2) }, derivative: { x: 1, y: -2 * (s - .5) } }) };
  // Original endpoint distance is .6, while exact parabola arc length is
  // .3*sqrt(1.36)+asinh(.6)/2 > .61.
  const exactOriginalArc = .3 * Math.sqrt(1.36) + Math.asinh(.6) / 2;
  assert.ok(exactOriginalArc > .61);
  const result = refineSurfaceStagnation({ ...args, relativeTolerance: .61 });
  assert.equal(result.kind, 'potential-minimum-bracket');
  assert.ok(result.iterations > 0);
  assert.ok(result.bracket.physicalWidth <= .61);
  const primitive = s => { const x = s - .5; return .5 * x * Math.sqrt(1 + 4 * x * x) + .25 * Math.asinh(2 * x); };
  assert.ok(result.bracket.physicalWidth >= primitive(result.bracket.right) - primitive(result.bracket.left));
});

test('ordinary off-vertex roots retain a finite tangential-velocity residual check', () => {
  const curve = { length: 1, knots: [0, 1], evaluate: s => ({ point: { x: s, y: (s - .37) ** 2 },
    derivative: { x: 1, y: 2 * (s - .37) } }) };
  const result = refineSurfaceStagnation({ curve, field: { u: 0, v: 1, panels: [], gamma: [] },
    left: .1, right: .8, chord: 1 });
  assert.equal(result.kind, 'velocity-zero'); assert.equal(result.derivativeStatus, 'finite');
  assert.ok(Math.abs(result.parameter - .37) < 1e-15);
  assert.ok(Math.abs(result.derivative) <= result.velocityTolerance);
  assert.equal(result.comparison, null);
});

test('neither an iteration limit nor an absent measured sign bracket is accepted as a corner minimum', () => {
  assert.throws(() => refineSurfaceStagnation({ ...junction(), maxIterations: 1 }), /physical width/);
  const args = junction(); args.field = { ...args.field, u: 10 };
  assert.throws(() => refineSurfaceStagnation(args), /finite signed/);
});
