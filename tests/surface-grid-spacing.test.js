import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { transform } from '../src/geometry/airfoil.js';
import { measureSurfaceGridSpacing, streamtubeSurfaceNodes } from '../src/geometry/tests/surface-grid-spacing.js';

const chord = { leadingEdge: { x: 0, y: 0 }, trailingEdge: { x: 1, y: 0 } };
const uniform = Array.from({ length: 51 }, (_, i) => ({ x: i / 50, y: 0 }));

test('uniform mid-chord spacing has unit length ratio and zero variation', () => {
  const r = measureSurfaceGridSpacing({ ...chord, points: uniform });
  assert.equal(r.passed, true); assert.equal(r.intervals.length, 30);
  assert.ok(Math.abs(r.lengthRatio - 1) < 1e-13);
  assert.ok(r.coefficientVariation < 1e-13);
  assert.ok(Math.abs(r.mean - .02) < 1e-15);
});

test('cosine surface distributions remain nearly uniform over the middle chord at several counts', () => {
  for (const n of [16, 32, 64, 128]) {
    const points = Array.from({ length: n + 1 }, (_, i) => {
      const x = .5 * (1 - Math.cos(Math.PI * i / n));
      return { x, y: .03 * Math.sin(Math.PI * x) };
    });
    const r = measureSurfaceGridSpacing({ ...chord, points });
    assert.equal(r.passed, true, `cosine ${n}: ${JSON.stringify(r)}`);
    assert.ok(r.lengthRatio < 1.3);
  }
});

test('gradual pile-up fails even when every neighboring interval ratio passes 1.5', () => {
  const points = Array.from({ length: 101 }, (_, i) => ({
    x: Math.expm1(i * Math.log(1.05)) / Math.expm1(100 * Math.log(1.05)), y: 0,
  }));
  const r = measureSurfaceGridSpacing({ ...chord, points });
  assert.ok(r.maximumAdjacentRatio < 1.051);
  assert.ok(r.lengthRatio > 3);
  assert.equal(r.status, 'excessive-clustering'); assert.equal(r.passed, false);
});

test('length ratios and chord windows are invariant under rotation, translation and scale', () => {
  const controls = { chord: 2.4, angle: 37, x: -3, y: .7 };
  const points = uniform.map(p => ({ x: p.x, y: .03 * Math.sin(Math.PI * p.x) }));
  const a = measureSurfaceGridSpacing({ ...chord, points });
  const [leadingEdge, trailingEdge] = transform([chord.leadingEdge, chord.trailingEdge], controls);
  const b = measureSurfaceGridSpacing({ leadingEdge, trailingEdge, points: transform(points, controls) });
  assert.deepEqual(b.intervals.map(p => p.from), a.intervals.map(p => p.from));
  for (const key of ['lengthRatio', 'coefficientVariation', 'mean']) assert.ok(Math.abs(a[key] - b[key]) < 1e-12);
});

test('sparse, incomplete, reversed and degenerate windows cannot silently pass', () => {
  assert.equal(measureSurfaceGridSpacing({ ...chord, points: uniform.filter((_, i) => i % 10 === 0) }).status, 'insufficient-resolution');
  assert.equal(measureSurfaceGridSpacing({ ...chord, points: uniform.slice(12, 30) }).status, 'insufficient-coverage');
  assert.throws(() => measureSurfaceGridSpacing({ ...chord, points: [...uniform].reverse() }), /reverse/);
  assert.throws(() => measureSurfaceGridSpacing({ ...chord, points: [uniform[0], uniform[0]] }), /distinct/);
  assert.throws(() => measureSurfaceGridSpacing({ ...chord, points: uniform, trailingEdge: chord.leadingEdge }), /nonzero/);
  for (const extra of [{ range: [.8, .2] }, { maximumLengthRatio: .5 }, { maximumCoefficientVariation: NaN }])
    assert.throws(() => measureSurfaceGridSpacing({ ...chord, points: uniform, ...extra }), /inputs/);
});

test('wall extraction selects both sides and all LE-through-TE nodes from unequal passage connectivity', () => {
  const nx = 5, tubes = [2, 3], vertices = [], cells = [];
  for (let g = 0; g < tubes.length; g++) {
    const offset = vertices.length, nt = tubes[g], id = (i, j) => offset + i * (nt + 1) + j;
    for (let i = 0; i <= nx; i++) for (let j = 0; j <= nt; j++) vertices.push({ x: i, y: 10 * g + j });
    for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) cells.push([id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)]);
  }
  const mesh = { vertices, cells, initialization: { tubes, streamwiseSegments: nx,
    potentialCrosslines: { bodies: [{ leadingIndex: 1, trailingIndex: 4 }] } } };
  for (const [side, y] of [['upper', 10], ['lower', 2]])
    assert.deepEqual(streamtubeSurfaceNodes(mesh, { body: 0, side }), [1, 2, 3, 4].map(x => ({ x, y })));
  assert.throws(() => streamtubeSurfaceNodes(mesh, { body: 1, side: 'upper' }), /connectivity/);
});

test('retained actual default exposes main-upper and flap clustering missed by the earlier screen', () => {
  const fixture = JSON.parse(readFileSync(new URL('./fixtures/default-grid-surface-spacing.json', import.meta.url)));
  const results = fixture.surfaces.map(surface => ({ ...surface, audit: measureSurfaceGridSpacing(surface) }));
  assert.deepEqual(results.filter(s => !s.audit.passed).map(s => [s.element, s.side]), [[1, 'upper'], [1, 'lower'], [0, 'upper']]);
  const main = results.find(s => s.element === 0 && s.side === 'upper');
  assert.ok(main.audit.lengthRatio > 17 && main.audit.lengthRatio < 18);
  assert.ok(main.audit.smallest.start > .7 && main.audit.smallest.end < .8);
  // This tests detection of a retained failure. Acceptance is a separate
  // command and must keep failing until the current generated grid improves.
});
