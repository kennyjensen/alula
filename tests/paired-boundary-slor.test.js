import test from 'node:test';
import { solveBlockTridiagonal } from '../src/numerics/block-tridiagonal.js';
import assert from 'node:assert/strict';
import { createEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { smoothPairedBoundaryGrid } from '../src/geometry/paired-boundary-slor.js';

const uniformFlow = (sourceForm, streamwiseSourceDiscretization) => {
  const xi = Array.from({ length: 9 }, (_, i) => i / 8), eta = [0, .04, .16, .4, .7, .91, 1];
  const nodes = xi.map(u => eta.map(e => ({ x: u + .04 * Math.sin(Math.PI * u) * e * (1 - e),
    y: e + .025 * Math.sin(Math.PI * u) * Math.sin(Math.PI * e) })));
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: eta.slice(1).map((e, j) => e - eta[j]), streamwiseCoordinates: xi,
    discretization: 'giles-1985', streamwiseSourceDiscretization, lineLinearization: 'full-metrics', lineGrouping: 'boundary-pairs', lineSearch: 'armijo',
    orthogonalBoundaryControl: { sourceForm, background: nodes.map(row => row.map(() => 0)), decay: { lower: 4, upper: 4 } } });
  return { system, nodes, xi, eta };
};
for (const sourceForm of ['poisson', 'metric-stretch']) for (const scheme of ['centered', 'grape-1980'])
  test(`paired ${sourceForm}/${scheme} relaxation recovers exact uniform-flow coordinates with positive cells`, () => {
    const { system, nodes, xi, eta } = uniformFlow(sourceForm, scheme), original = structuredClone(nodes);
    // Residual-only metrics skip derivative construction, not physics checks.
    assert.deepEqual(system.residuals(nodes), system.residuals(nodes, system.metrics(nodes)));
    const result = smoothPairedBoundaryGrid(system, { maxSweeps: 150, tolerance: 1e-9 });
    assert.equal(result.converged, true, result.reason); assert.ok(result.residual <= 1e-9);
    assert.ok(result.history.every(h => !h.invalidCells)); assert.equal(result.quality.valid, true);
    result.history.forEach(h => {
      for (const step of h.rowSteps ?? []) {
        if (step.skipped) { assert.ok(step.residual <= step.tolerance); continue; }
        assert.ok(step.merit <= (1 - 2e-4 * step.fraction) * step.baseMerit);
        assert.ok(step.minCornerSine > 0);
        assert.equal(step.trials.filter(t => t.accepted).length, 1);
      }
    });
    result.nodes.forEach((row, i) => row.forEach((p, j) => {
      // Uniform physical potential/streamfunction: phi=x, psi=y. These
      // errors are independent of the solver's inverse-grid residual.
      assert.ok(Math.abs(p.x - xi[i]) < 2e-8); assert.ok(Math.abs(p.y - eta[j]) < 2e-8);
      if (!i || i === xi.length - 1 || !j || j === eta.length - 1) assert.deepEqual(p, nodes[i][j]);
    }));
    assert.deepEqual(nodes, original);
  });

test('paired relaxation leaves an invalid starting grid unchanged and rejects unsupported settings', () => {
  const { system, nodes } = uniformFlow('poisson', 'centered'), invalid = structuredClone(nodes);
  invalid[3][2].y = -.4;
  const result = smoothPairedBoundaryGrid(system, { initial: invalid });
  assert.equal(result.converged, false); assert.match(result.reason, /Invalid starting state: folded grid/);
  assert.equal(result.history.length, 1); assert.deepEqual(result.nodes, invalid);
  assert.throws(() => smoothPairedBoundaryGrid(system, { maxSweeps: -1 }), /settings/);
  assert.throws(() => smoothPairedBoundaryGrid(system, { tolerance: 0 }), /settings/);
});

// Independent full-grid replay: every trial recomputes all metrics and owns all
// points, so it detects stale metric caches, workspace reuse and trial aliasing.
function fullGridSweep(system, grid, groups, omega) {
  let next = structuredClone(grid);
  const accepted = [];
  for (const rows of groups) {
    const matrix = system.linearizeRows(next, rows), correction = solveBlockTridiagonal(matrix);
    const baseMerit = .5 * matrix.rhs.reduce((sum, row, i) => sum + row.reduce((s, r, k) => s + (r / matrix.scales[i][k]) ** 2, 0), 0);
    let found = false;
    for (let halving = 0; halving <= 20; halving++) {
      const fraction = omega * 2 ** -halving, trial = structuredClone(next);
      for (let i = 1; i < grid.length - 1; i++) for (const [r, j] of rows.entries()) for (const [k, key] of ['x', 'y'].entries())
        trial[i][j][key] += fraction * correction[i - 1][2 * r + k];
      if (!system.quality(trial).valid) continue;
      let residual;
      try { residual = system.residuals(trial, system.metrics(trial)); } catch { continue; }
      // Undo residual normalization, then apply the fixed line scaling.
      const metrics = system.metrics(trial);
      let merit = 0;
      for (let i = 1; i < grid.length - 1; i++) for (const [r, j] of rows.entries()) for (const [k, key] of ['x', 'y'].entries()) {
        const row = residual.rows.find(row => row.i === i && row.j === j);
        const scale = metrics[i][j].alpha + metrics[i][j].gamma;
        // Recover the unnormalized operator using the public physical scale.
        merit += .5 * (row[key] * scale * system.lengthScale / matrix.scales[i - 1][2 * r + k]) ** 2;
      }
      if (merit <= (1 - 2e-4 * fraction) * baseMerit) {
        next = trial; accepted.push(fraction); found = true; break;
      }
    }
    assert.ok(found, 'full-grid reference found an admissible decreasing step');
  }
  return { nodes: next, fractions: accepted };
}

for (const sourceForm of ['poisson', 'metric-stretch']) for (const sides of [['lower'], ['lower', 'upper']])
  test(`incremental paired trials match full-grid replay: ${sourceForm}, ${sides}`, () => {
    const nx = 12, nt = 11;
    const xi = Array.from({ length: nx + 1 }, (_, i) => i / nx);
    const nodes = xi.map(u => Array.from({ length: nt + 1 }, (_, j) => {
      const e = j / nt;
      return { x: u + .06 * Math.sin(Math.PI * u) * e * (1 - e),
        y: e + .04 * Math.sin(Math.PI * u) * Math.sin(Math.PI * e) };
    }));
    const system = createEllipticStreamtubeGrid({ nodes, massFlows: Array(nt).fill(1 / nt), streamwiseCoordinates: xi,
      discretization: 'giles-1985', lineLinearization: 'full-metrics', lineGrouping: 'boundary-pairs', lineSearch: 'armijo',
      orthogonalBoundaryControl: { sourceForm, sides, background: nodes.map(row => row.map(() => 0)),
        corners: { lower: [4], ...(sides.length === 2 ? { upper: [7] } : {}) } } });
    const groups = [[1, 2], [3], [4], [5], [6], [7], [8], ...(sides.length === 2 ? [[9, 10]] : [[9], [10]])];
    let grid = structuredClone(nodes), backtracked = false;
    const snapshot = structuredClone(grid);
    // An overrelaxed step exercises rejected trials as well as accepted ones.
    for (let iteration = 0; iteration < 3; iteration++) {
      const original = structuredClone(grid), expected = fullGridSweep(system, grid, groups, 1.9999);
      const actual = system.sweep(grid, 1.9999);
      assert.deepEqual(grid, original, 'sweep must not mutate caller coordinates');
      assert.deepEqual(actual.nodes, expected.nodes);
      assert.deepEqual(actual.rowSteps.map(row => row.fraction), expected.fractions);
      backtracked ||= actual.rowSteps.some(row => row.trials.length > 1);
      grid = actual.nodes;
    }
    assert.deepEqual(nodes, snapshot);
    assert.ok(backtracked, 'exercise rejected trial restoration');
  });
