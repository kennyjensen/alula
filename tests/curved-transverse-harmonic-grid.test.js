import test from 'node:test';
import assert from 'node:assert/strict';
import { createPolynomialGridGeometry } from '../src/geometry/polynomial-grid-geometry.js';
import { createTransverseHarmonicGrid, smoothTransverseHarmonicGrid } from '../src/geometry/transverse-harmonic-grid.js';
import { solveCurvedHarmonicReference } from '../src/geometry/tests/curved-harmonic-reference.js';
import { conformalPolynomialGridFixture } from './fixtures/polynomial-grid.js';

const boundary = (i, j, data) => !i || !j || i === data.nx || j === data.nt;
const perturb = data => data.nodes.map((row, i) => row.map((p, j) => ({ ...p,
  y: p.y + (boundary(i, j, data) ? 0 : .025 * Math.sin(Math.PI * i / data.nx) * Math.sin(Math.PI * j / data.nt)),
})));

test('certified curved Newton recovers the exact conformal streamfunction with either fixed quadrature rule', t => {
  const data = conformalPolynomialGridFixture({ nx: 4, nt: 3 }), curvedGeometry = createPolynomialGridGeometry(data);
  const nodes = perturb(data), before = structuredClone(nodes), results = [];
  for (const quadratureOrder of [3, 5]) {
    const system = createTransverseHarmonicGrid({ ...data, nodes, curvedGeometry, quadratureOrder });
    const result = smoothTransverseHarmonicGrid(system, { tolerance: 1e-11 });
    assert.equal(result.converged, true, JSON.stringify({ reason: result.reason, history: result.history }));
    assert.equal(result.quality.valid, true); assert.equal(result.physicsValidated, false);
    assert.match(result.geometryModel, /curved correction/);
    assert.equal(result.linearBackend, 'scalar streamwise SLOR');
    let error = 0, displacementError = 0;
    result.nodes.forEach((row, i) => row.forEach((p, j) => {
      assert.equal(p.x, nodes[i][j].x);
      if (boundary(i, j, data)) assert.deepEqual(p, nodes[i][j]);
      error = Math.max(error, Math.abs(data.psiAt(p) - j / data.nt));
      displacementError = Math.max(displacementError, Math.hypot(p.x - data.nodes[i][j].x, p.y - data.nodes[i][j].y));
    }));
    assert.ok(error < 1e-10, `Physical psi error ${error}`);
    assert.ok(displacementError < 2e-10, `Exact map discrepancy ${displacementError}`);
    const independent = solveCurvedHarmonicReference({ ...data, nodes: result.nodes, geometry: curvedGeometry.onGrid(result.nodes) },
      { refinement: 2, quadratureOrder: 5 });
    assert.ok(independent.maximumTubeIntervals < 1e-9, JSON.stringify(independent.linear));
    results.push({ quadratureOrder, iterations: result.history.length - 1, error, displacementError,
      independentTubeIntervals: independent.maximumTubeIntervals, minimumJacobian: result.quality.minimumJacobian,
      minimumTransversality: result.quality.minimumTransversality });
  }
  assert.deepEqual(nodes, before); t.diagnostic(JSON.stringify(results));
});

test('curved smoothing requires the fixed certified geometry and exact prescribed boundary positions', () => {
  const data = conformalPolynomialGridFixture({ nx: 3, nt: 3 }), geometry = createPolynomialGridGeometry(data);
  assert.throws(() => createTransverseHarmonicGrid({ ...data, geometryCorrection: geometry.correction }), /certified polynomial/);
  assert.throws(() => createTransverseHarmonicGrid({ ...data, curvedGeometry: { ...geometry } }), /whole-cell certificates/);
  const shiftedBoundary = structuredClone(data.nodes); shiftedBoundary[0][1].y += .001;
  assert.throws(() => createTransverseHarmonicGrid({ ...data, nodes: shiftedBoundary, curvedGeometry: geometry }), /boundaries must match/);
  const reversed = data.directions.map(row => row.map(() => ({ x: 0, y: -1 })));
  assert.throws(() => createTransverseHarmonicGrid({ ...data, directions: reversed, curvedGeometry: geometry }), /certified positive/);
  const system = createTransverseHarmonicGrid({ ...data, nodes: perturb(data), curvedGeometry: geometry });
  assert.throws(() => { system.initial[0][1].y += .001; }, TypeError);
  assert.throws(() => { system.directions[1][1].x = 1; }, TypeError);
  assert.throws(() => { system.initial[0] = []; }, TypeError);
  const offGuide = structuredClone(system.initial); offGuide[1][1].x += .001;
  assert.throws(() => system.evaluate(offGuide), /guide line/);
  const boundaryMove = structuredClone(system.initial); boundaryMove[0][1].y += .001;
  assert.throws(() => system.evaluate(boundaryMove), /boundaries/);
  const failed = smoothTransverseHarmonicGrid(system, { linearSolve: () => ({ converged: false }) });
  assert.equal(failed.converged, false); assert.equal(failed.reason, 'linear solve failed');
  assert.deepEqual(failed.nodes, system.initial);
});

test('positive straight-sided nodes cannot hide a folded polynomial bow from the solver', () => {
  const data = conformalPolynomialGridFixture({ nx: 3, nt: 3 }), controlPoints = structuredClone(data.controlPoints);
  // An interior Bernstein coefficient leaves every face and observation
  // node intact while folding the curved map inside its first cell.
  controlPoints[0][0][1][1].y += 1;
  const curvedGeometry = createPolynomialGridGeometry({ ...data, controlPoints });
  assert.equal(curvedGeometry.quality(data.nodes, data.directions).valid, false);
  assert.throws(() => createTransverseHarmonicGrid({ ...data, curvedGeometry }), /certified positive/);
});
