// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { limitStreamtubeGridStep, interpolateStreamtubeGridNodes, streamtubeGridConvexity,
  assertConvexStreamtubeGrid, StreamtubeGridStepError, STREAMTUBE_CONVEX_MINIMUM_CORNER_SINE } from '../src/geometry/streamtube-convex-step.js';
import { STREAMTUBE_MINIMUM_CORNER_SINE } from '../src/euler/streamtube-mesh-preview.js';
import { potentialGridQuality } from '../src/geometry/potential-plane-grid.js';

const square = () => [[[{ x: 0, y: 0 }, { x: 0, y: 1 }], [{ x: 1, y: 0 }, { x: 1, y: 1 }]]];
const map = (nodes, fn) => nodes.map(group => group.map(row => row.map(p => fn(p))));
const near = (actual, expected, tolerance = 2e-14) => assert.ok(Math.abs(actual - expected) <= tolerance * Math.max(1, Math.abs(expected)), `${actual} != ${expected}`);
function independentPathCheck(from, to, step) {
  for (let k = 0; k <= 100; k++) for (const group of interpolateStreamtubeGridNodes(from, to, step * k / 100))
    assert.equal(potentialGridQuality(group).valid, true, `Invalid grid at sample ${k}.`);
}

test('safe translation and shear retain the exact full update and do not mutate either grid', () => {
  const from = square(), to = map(from, p => ({ x: p.x + .3 * p.y + 7, y: 1.2 * p.y - 2, marker: 'proposal' }));
  const beforeFrom = structuredClone(from), beforeTo = structuredClone(to);
  const result = limitStreamtubeGridStep(from, to);
  assert.equal(result.step, 1); assert.equal(result.limited, false); assert.equal(result.limiter, null);
  assert.equal(result.cellsChecked, 1); assert.equal(result.quality.valid, true);
  assert.deepEqual(interpolateStreamtubeGridNodes(from, to, result.step), to);
  assert.notEqual(interpolateStreamtubeGridNodes(from, to, result.step)[0][0][0], to[0][0][0]);
  assert.deepEqual(from, beforeFrom); assert.deepEqual(to, beforeTo);
  independentPathCheck(from, to, result.step);
});

test('the first folding corner gives an analytic bound before the concave proposal', () => {
  const from = square(), to = square(); to[0][1][1] = { x: .2, y: .2 };
  assert.equal(streamtubeGridConvexity(to).valid, false);
  const result = limitStreamtubeGridStep(from, to);
  near(result.limiter.boundaryStep, .625); near(result.step, .95 * .625);
  assert.deepEqual({ ...result.limiter, boundaryStep: .625 }, { group: 0, i: 0, tube: 0, corner: 1, boundaryStep: .625, reason: 'corner-orientation' });
  independentPathCheck(from, to, result.step);
  assert.equal(streamtubeGridConvexity(interpolateStreamtubeGridNodes(from, to, .625 + 1e-7)).valid, false);
});

test('positive endpoint cells cannot conceal an intervening inverted interval', () => {
  const from = square(), to = map(from, p => ({ x: -2 * p.x, y: -.5 * p.y }));
  // Each orientation is (1 - 3t)(1 - 1.5t): positive endpoints,
  // negative between t=1/3 and t=2/3. Endpoint-only guards miss this.
  assert.equal(streamtubeGridConvexity(to).valid, true);
  assert.equal(streamtubeGridConvexity(interpolateStreamtubeGridNodes(from, to, .5)).valid, false);
  const result = limitStreamtubeGridStep(from, to);
  near(result.limiter.boundaryStep, 1 / 3); near(result.step, .95 / 3);
  independentPathCheck(from, to, result.step);
});

test('a double orientation root blocks passage through a collapsed grid even without a sign change', () => {
  const from = square(), to = map(from, p => ({ x: -p.x, y: -p.y }));
  // (1 - 2t)^2 is never negative; its double zero still closes the domain.
  const result = limitStreamtubeGridStep(from, to);
  near(result.limiter.boundaryStep, .5); near(result.step, .475);
  assert.equal(streamtubeGridConvexity(to).valid, true);
  assert.equal(streamtubeGridConvexity(interpolateStreamtubeGridNodes(from, to, .5)).valid, false);
  independentPathCheck(from, to, result.step);
});

test('near-double-root discrimination retains strictly convex near-180-degree rotations', () => {
  for (const delta of [1e-6, 1e-8, 1e-10]) {
    const angle = Math.PI - delta, from = square();
    const to = map(from, p => ({ x: Math.cos(angle) * p.x - Math.sin(angle) * p.y,
      y: Math.sin(angle) * p.x + Math.cos(angle) * p.y }));
    const result = limitStreamtubeGridStep(from, to);
    assert.equal(result.step, 1); assert.equal(result.limited, false);
    independentPathCheck(from, to, 1);
  }
});

test('near-double roots still detect a narrow inverted interval with positive endpoints', () => {
  for (const gap of [1e-6, 1e-8, 1e-10]) {
    const from = square(), to = map(from, p => ({ x: -p.x, y: -(1 + gap) * p.y }));
    const result = limitStreamtubeGridStep(from, to);
    near(result.limiter.boundaryStep, 1 / (2 + gap), 2e-12);
    independentPathCheck(from, to, result.step);
  }
});

test('near-linear orientation polynomials retain their finite root without cancellation', () => {
  for (const epsilon of [0, 1e-15, 1e-10, 1e-4]) {
    const from = square(), to = map(from, p => ({ x: -p.x, y: (1 + epsilon) * p.y }));
    const result = limitStreamtubeGridStep(from, to);
    near(result.limiter.boundaryStep, .5); near(result.step, .475);
    independentPathCheck(from, to, result.step);
  }
});

test('corner bounds are invariant under physical units and rotations', () => {
  const a = square(), b = square(); b[0][1][1] = { x: .2, y: .2 };
  for (const [sx, sy, angle] of [[1e-120, 1e-120, .3], [1e120, 1e120, .3]]) {
    const transform = p => ({ x: sx * p.x * Math.cos(angle) - sy * p.y * Math.sin(angle),
      y: sx * p.x * Math.sin(angle) + sy * p.y * Math.cos(angle) });
    const from = map(a, transform), to = map(b, transform), result = limitStreamtubeGridStep(from, to);
    near(result.limiter.boundaryStep, .625); near(result.step, .59375);
    assert.equal(result.quality.valid, true);
    assert.equal(assertConvexStreamtubeGrid(interpolateStreamtubeGridNodes(from, to, result.step)).valid, true);
  }
});

test('independent edge normalization handles unequal edge scales without determinant underflow', () => {
  for (const [sx, sy] of [[1e-200, 1e-200], [1e200, 1e200], [1e-120, 1e120], [1e120, 1e-120]]) {
    const from = map(square(), p => ({ x: sx * p.x, y: sy * p.y }));
    const to = map(from, p => ({ x: -p.x, y: p.y })), result = limitStreamtubeGridStep(from, to);
    near(result.limiter.boundaryStep, .5); near(result.step, .475);
    assert.equal(assertConvexStreamtubeGrid(interpolateStreamtubeGridNodes(from, to, result.step)).valid, true);
  }
});

test('maximumStep is absolute, safety is applied only when limited, and zero steps retain the start', () => {
  const from = square(), to = map(from, p => ({ x: -p.x, y: p.y }));
  const short = limitStreamtubeGridStep(from, to, { maximumStep: .4 });
  assert.equal(short.step, .4); assert.equal(short.limited, false);
  const limited = limitStreamtubeGridStep(from, to, { maximumStep: .8, safetyFraction: .9 });
  near(limited.step, .45); near(limited.limiter.boundaryStep, .5);
  const zero = limitStreamtubeGridStep(from, to, { maximumStep: 0 });
  assert.equal(zero.step, 0); assert.equal(zero.limited, false);
  assert.deepEqual(interpolateStreamtubeGridNodes(from, to, 0), from);
});

test('the earliest cell among all passages controls one common step', () => {
  const from = [square()[0], ...map(square(), p => ({ x: p.x + 3, y: p.y }))];
  const to = structuredClone(from); to[1][1][1] = { x: 2, y: -1 };
  const result = limitStreamtubeGridStep(from, to);
  assert.equal(result.cellsChecked, 2); assert.equal(result.limiter.group, 1);
  near(result.limiter.boundaryStep, .25); near(result.step, .2375);
  independentPathCheck(from, to, result.step);
});

test('actual-coordinate checks retain the existing sine threshold and reject degenerate starts with typed diagnostics', () => {
  assert.equal(STREAMTUBE_CONVEX_MINIMUM_CORNER_SINE, STREAMTUBE_MINIMUM_CORNER_SINE);
  const from = square(), to = square(); to[0][1][1] = { x: 1.5, y: 1e-14 };
  // Positive orientations alone do not satisfy the existing numerical
  // domain. The final-coordinate check must still enforce its threshold.
  const result = limitStreamtubeGridStep(from, to);
  assert.equal(result.limited, true); assert.equal(result.limiter.reason, 'corner-sine');
  assert.ok(result.validationBacktracks > 0);
  assert.ok(result.quality.minCornerSine > STREAMTUBE_MINIMUM_CORNER_SINE);
  assert.throws(() => limitStreamtubeGridStep(to, from), error => error instanceof StreamtubeGridStepError
    && error.code === 'streamtube-grid-invalid-start' && error.details.group === 0 && error.details.sine <= 1e-12);
  assert.throws(() => assertConvexStreamtubeGrid(to), error => error instanceof StreamtubeGridStepError
    && error.code === 'streamtube-grid-nonconvex' && /i=0, group=0, tube=0/.test(error.message)
    && error.diagnostics.cell.i === 0 && error.diagnostics.cell.group === 0 && error.diagnostics.cell.tube === 0);
});

test('malformed, mismatched and nonfinite coordinates fail explicitly without mutating inputs', () => {
  const from = square(), before = structuredClone(from);
  for (const to of [[], [[]], [[from[0][0]]], [[from[0][0], [{ x: NaN, y: 0 }, { x: 1, y: 1 }]]]])
    assert.throws(() => limitStreamtubeGridStep(from, to), StreamtubeGridStepError);
  for (const options of [{ maximumStep: -1 }, { maximumStep: 2 }, { maximumStep: NaN }, { safetyFraction: 1 }, { safetyFraction: 0 }])
    assert.throws(() => limitStreamtubeGridStep(from, from, options), error => error.code === 'streamtube-grid-step-control');
  assert.deepEqual(from, before);
});
