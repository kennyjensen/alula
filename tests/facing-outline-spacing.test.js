import test from 'node:test';
import assert from 'node:assert/strict';
import { createMonotoneCubicMap } from '../src/numerics/monotone-cubic.js';
import { reconcileFacingArcSlopes } from '../src/geometry/facing-arc-spacing.js';
import { reconcileFacingOutlineSlopes } from '../src/geometry/facing-outline-spacing.js';

const close = (a, b, tol = 2e-13) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const map = (knots, values, slopes) => ({ knots, values,
  slopes: slopes ?? createMonotoneCubicMap(knots, values).slopes });
const interpolants = (maps, fit) => maps.map((m, i) => createMonotoneCubicMap(m.knots, m.values,
  { derivatives: 'prescribed', slopes: fit.slopes[i] }));
const normalized = (m, x, a, b) => (m.value(x) - m.value(a)) / (m.value(b) - m.value(a));

test('forest fit agrees with the existing two-wall solution, including active bounds and disjoint maps', () => {
  for (const maps of [
    [map([0, 1, 2], [0, .8, 1]), map([1, 2, 3], [0, .5, 1])],
    [map([0, 1, 2], [0, .01, 1.01]), map([1, 2, 3], [0, 1, 1.01])],
    [map([0, 1], [0, 2]), map([2, 3], [0, 1])],
  ]) {
    const before = structuredClone(maps), old = reconcileFacingArcSlopes(...maps);
    const fit = reconcileFacingOutlineSlopes({ maps, pairs: [[0, 1]] });
    fit.slopes.forEach((row, m) => row.forEach((v, i) => close(v, [old.first, old.second][m][i])));
    fit.knots.forEach((k, i) => {
      close(k.pairs[0].targetLogRatio, old.knots[i].targetLogRatio);
      close(k.components[0].commonLogAdjustment, old.knots[i].commonLogAdjustment);
    });
    assert.equal(fit.exactModernMsetLaw, false); assert.deepEqual(maps, before);
    interpolants(maps, fit);
  }
});

test('three bodies merge duplicate cut variables and minimize changes of all six occurrences without order dependence', () => {
  const secants = [2, 2, 4, 4, 8, 8], original = [1.2, 2.2, 3.6, 4.4, 6.4, 8.8];
  const maps = secants.map((s, i) => map([0, 1, 2], [10 * i, 10 * i + s, 10 * i + 2 * s], Array(3).fill(original[i])));
  const pairs = [[1, 2], [3, 4]], shared = [0, 2].flatMap(rank => [[0, 1], [2, 3], [4, 5]].map(maps => ({ rank, maps })));
  const before = structuredClone({ maps, pairs, shared });
  const fit = reconcileFacingOutlineSlopes({ maps, pairs, shared });
  const factor = Math.exp(original.reduce((sum, v, i) => sum + Math.log(v / secants[i]), 0) / 6);
  for (const i of [0, 2]) {
    for (let m = 0; m < 6; m++) close(fit.slopes[m][i], secants[m] * factor);
    for (const [a, b] of [[0, 1], [2, 3], [4, 5]]) assert.equal(fit.slopes[a][i], fit.slopes[b][i]);
    const component = fit.knots.find(k => k.rank === i).components[0];
    assert.equal(component.occurrences, 6); assert.equal(component.variables.length, 3);
    assert.equal(component.commonLogAdjustment, 0);
  }
  const order = [4, 1, 5, 2, 0, 3], inverse = index => order.indexOf(index);
  const reordered = reconcileFacingOutlineSlopes({ maps: order.map(i => maps[i]),
    pairs: pairs.toReversed().map(p => p.toReversed().map(inverse)),
    shared: shared.toReversed().map(s => ({ rank: s.rank, maps: s.maps.toReversed().map(inverse) })) });
  fit.slopes.forEach((row, m) => row.forEach((v, i) => close(v, reordered.slopes[inverse(m)][i])));
  assert.deepEqual({ maps, pairs, shared }, before);
  interpolants(maps, fit);
});

test('compatible block ratios give the same normalized cubic at every facing block with physical hits and C1 joins fixed', () => {
  const knots = [0, .4, 1.3, 2], values = [0, .2, 1, 1.4];
  const maps = [1, 3, .4].map((scale, m) => {
    const result = map(knots, values.map(v => 7 * m + scale * v));
    result.slopes = result.slopes.map((s, i) => s * (.4 + .1 * ((i + m) % 3)));
    return result;
  });
  const fit = reconcileFacingOutlineSlopes({ maps, pairs: [[0, 1], [1, 2]] }), curves = interpolants(maps, fit);
  for (let k = 0; k < knots.length - 1; k++) for (const t of [0, .01, .1, .3, .7, .95, 1]) {
    const a = knots[k], b = knots[k + 1], x = a + t * (b - a);
    for (const c of curves.slice(1)) close(normalized(curves[0], x, a, b), normalized(c, x, a, b), 2e-13);
  }
  for (const [m, curve] of curves.entries()) for (let i = 0; i < knots.length; i++) {
    assert.equal(curve.value(knots[i]), maps[m].values[i]);
    close(curve.evaluate(knots[i]).derivative, fit.slopes[m][i]);
    if (i && i < knots.length - 1) close(curve.evaluate(knots[i] - 1e-9).derivative,
      curve.evaluate(knots[i] + 1e-9).derivative, 2e-7);
  }
  assert.ok(fit.knots.every(k => k.pairs.every(p => p.rmsIncidentLogMismatch < 2e-14)));
});

test('incompatible incident block ratios record the log compromise while keeping positive C1 interpolants', () => {
  const maps = [map([0, 1, 2], [0, 1, 5]), map([0, 1, 2], [0, 4, 5]), map([0, 1, 2], [0, 8, 10])];
  const fit = reconcileFacingOutlineSlopes({ maps, pairs: [[0, 1], [1, 2]] });
  const middle = fit.knots.find(k => k.rank === 1), pair = middle.pairs[0];
  assert.equal(pair.targetLogRatio, 0); close(fit.slopes[0][1], fit.slopes[1][1]);
  close(pair.rmsIncidentLogMismatch, Math.log(4));
  const objective = ratio => pair.incidentLogRatios.reduce((s, target) => s + (ratio - target) ** 2, 0);
  assert.ok(objective(0) < objective(.1) && objective(0) < objective(-.1));
  const curves = interpolants(maps, fit);
  assert.ok(Math.abs(normalized(curves[0], .4, 0, 1) - normalized(curves[1], .4, 0, 1)) > .01);
  for (const curve of curves) {
    for (let k = 0; k <= 40; k++) assert.ok(curve.evaluate(k / 20).derivative > 0);
    close(curve.evaluate(1 - 1e-9).derivative, curve.evaluate(1 + 1e-9).derivative, 2e-7);
  }
});

test('unit changes commute with the fit; merged cut copies use the same metric unit conversion', () => {
  const maps = [map([0, 1, 2], [0, .01, 1.01]), map([1, 2, 3], [0, 1, 1.01])];
  const fit = reconcileFacingOutlineSlopes({ maps, pairs: [[0, 1]] }), factors = [7, 4];
  const convert = (m, factor, i) => ({ knots: m.knots.map(v => 3 + 2 * v),
    values: m.values.map(v => 9 * i + factor * v), slopes: m.slopes.map(v => factor * v / 2) });
  const scaled = reconcileFacingOutlineSlopes({ maps: maps.map((m, i) => convert(m, factors[i], i)), pairs: [[0, 1]] });
  fit.slopes.forEach((row, m) => row.forEach((v, i) => close(scaled.slopes[m][i], factors[m] * v / 2)));
  const sharedMaps = [map([0, 1], [0, 2], [.7, 1]), map([0, 1], [5, 7], [1.4, 1.6])];
  const shared = [{ rank: 0, maps: [0, 1] }], sharedFit = reconcileFacingOutlineSlopes({ maps: sharedMaps, pairs: [], shared });
  const converted = reconcileFacingOutlineSlopes({ maps: sharedMaps.map((m, i) => convert(m, 5, i)), pairs: [],
    shared: [{ rank: 3, maps: [0, 1] }] });
  sharedFit.slopes.forEach((row, m) => row.forEach((v, i) => close(converted.slopes[m][i], 2.5 * v)));
  assert.equal(converted.slopes[0][0], converted.slopes[1][0]);
});

test('invalid shared ranks, nonpositive cones, self demands and cycles fail explicitly', () => {
  const a = map([0, 1, 2], [0, 1, 2]), b = map([0, 1, 2], [0, 2, 4]);
  for (const shared of [null, [{ rank: .5, maps: [0, 1] }], [{ rank: 1, maps: [0, 2] }],
    [{ rank: 1, maps: [0, 0] }], [{ rank: Infinity, maps: [0, 1] }], [{ rank: 1, maps: [0] }]])
    assert.throws(() => reconcileFacingOutlineSlopes({ maps: [a, b], pairs: [], shared }), /shared|Shared/);
  for (const pairs of [null, [[0, 0]], [[0, 2]], [[0, 1], [1, 0]], [[0, 1, 0]]])
    assert.throws(() => reconcileFacingOutlineSlopes({ maps: [a, b], pairs }), /pair/);
  assert.throws(() => reconcileFacingOutlineSlopes({ maps: [a, { ...b, slopes: [1, 0, 1] }], pairs: [[0, 1]] }), /Invalid/);
  assert.throws(() => reconcileFacingOutlineSlopes({ maps: [map([0, 1], [0, 1e308], [1, 1])], pairs: [] }), /cone/);
  assert.throws(() => reconcileFacingOutlineSlopes({ maps: [map([0, 1], [0, 1], [3, 1])], pairs: [] }), /cone/);
  assert.throws(() => reconcileFacingOutlineSlopes({ maps: [a, map([.5, 1.5], [0, 1])], pairs: [[0, 1]] }), /shared physical hit/);
  assert.throws(() => reconcileFacingOutlineSlopes({ maps: [a, b], pairs: [[0, 1]],
    shared: [{ rank: 1, maps: [0, 1] }] }), /self-edge/);
  const redundant = reconcileFacingOutlineSlopes({ maps: [a, structuredClone(a)], pairs: [[0, 1]],
    shared: [{ rank: 1, maps: [0, 1] }] });
  assert.equal(redundant.knots[1].pairs[0].redundantSharedConstraint, true);
  assert.throws(() => reconcileFacingOutlineSlopes({ maps: [a, b, a], pairs: [[0, 1], [1, 2], [2, 0]] }), /cycles.*unsupported/);
});
