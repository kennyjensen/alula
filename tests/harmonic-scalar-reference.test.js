import test from 'node:test';
import assert from 'node:assert/strict';
import { solveHarmonicScalarReference } from '../src/geometry/tests/harmonic-scalar-reference.js';
import { solveHarmonicGridReference } from '../src/geometry/harmonic-grid-audit.js';

const rectangular = n => Array.from({ length: n + 1 }, (_, i) => Array.from({ length: n + 1 }, (_, j) => {
  const x = i / n, y = j / n, d = .035 * Math.sin(Math.PI * x) * Math.sin(Math.PI * y);
  return { x: x + d, y: y - .3 * d };
}));

test('arbitrary Dirichlet scalar labels recover an affine physical field regardless of interior seed', () => {
  const nodes = rectangular(5), exact = p => 3 + 2 * p.x - .4 * p.y;
  const labels = nodes.map((row, i) => row.map((p, j) => exact(p) + (i && i < 5 && j && j < 5 ? .2 * Math.cos(i + 2 * j) : 0)));
  const saved = structuredClone({ nodes, labels });
  for (const refinement of [1, 2, 4]) {
    const result = solveHarmonicScalarReference({ nodes, labels }, { refinement });
    result.values.forEach((row, i) => row.forEach((v, j) => assert.ok(Math.abs(v - exact(nodes[i][j])) < 3e-13)));
    assert.ok(result.maximumError > .1);
  }
  assert.deepEqual({ nodes, labels }, saved);
  const constant = solveHarmonicScalarReference({ nodes, labels: nodes.map(row => row.map(() => 12)) });
  assert.equal(constant.maximumError, 0);
});

test('scalar reference agrees with the existing mass-coordinate reference and preserves scaling', () => {
  const nodes = rectangular(5), labels = nodes.map(row => row.map((_, j) => j / 5));
  const a = solveHarmonicScalarReference({ nodes, labels }, { refinement: 2 });
  const b = solveHarmonicGridReference({ nodes, massFlows: Array(5).fill(.2) }, { refinement: 2 });
  const transformed = nodes.map(row => row.map(p => ({ x: -2 + 3 * (.8 * p.x - .6 * p.y), y: 7 + 3 * (.6 * p.x + .8 * p.y) })));
  const c = solveHarmonicScalarReference({ nodes: transformed, labels: labels.map(row => row.map(v => 2 * v + 3)) }, { refinement: 2 });
  a.values.forEach((row, i) => row.forEach((v, j) => {
    assert.ok(Math.abs(v - b.values[i][j].eta) < 2e-13);
    assert.ok(Math.abs(c.values[i][j] - (2 * v + 3)) < 3e-13);
  }));
});

test('nonlinear analytic harmonic field converges under physical-grid refinement', t => {
  const exact = p => Math.sin(Math.PI * p.x) * Math.sinh(Math.PI * p.y) / Math.sinh(Math.PI);
  let previous = Infinity; const evidence = [];
  for (const n of [8, 16, 32, 64]) {
    const nodes = rectangular(n), labels = nodes.map((row, i) => row.map((p, j) => !i || i === n || !j || j === n ? exact(p) : 0));
    const result = solveHarmonicScalarReference({ nodes, labels }, { refinement: 2 });
    let error = 0;
    result.values.forEach((row, i) => row.forEach((v, j) => { error = Math.max(error, Math.abs(v - exact(nodes[i][j]))); }));
    assert.ok(error < .35 * previous); previous = error; evidence.push({ n, error });
  }
  t.diagnostic(JSON.stringify(evidence)); assert.ok(previous < .0005);
});

test('scalar reference rejects malformed data, inverted geometry and excessive budgets', () => {
  const nodes = rectangular(4), labels = nodes.map(row => row.map(p => p.x));
  for (const invalid of [labels.slice(1), labels.map(row => row.slice(1)), labels.map(row => row.map(() => NaN)), Array(5)])
    assert.throws(() => solveHarmonicScalarReference({ nodes, labels: invalid }), /Invalid/);
  assert.throws(() => solveHarmonicScalarReference({ nodes, labels }, { refinement: 4, maxUnknowns: 10 }), /budget/);
  assert.throws(() => solveHarmonicScalarReference({ nodes: nodes.toReversed(), labels }), /Invalid/);
});
