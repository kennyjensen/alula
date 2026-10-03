import test from 'node:test';
import assert from 'node:assert/strict';
import { createTransverseHarmonicGrid, smoothTransverseHarmonicGrid } from '../src/geometry/transverse-harmonic-grid.js';
import { solveHarmonicGridReference } from '../src/geometry/harmonic-grid-audit.js';

const b = .12, clusteredX = u => u + .6 * u * (1 - u);
const psi = p => p.y + b * (p.x * p.x - p.y * p.y);
const inverseY = (x, eta) => 2 * (eta - b * x * x) / (1 + Math.sqrt(1 - 4 * b * (eta - b * x * x)));
const boundary = (i, j, nx, nt) => !i || !j || i === nx || j === nt;
const masses = eta => eta.slice(1).map((v, j) => v - eta[j]);
const vertical = nodes => nodes.map(row => row.map(() => ({ x: 0, y: 1 })));
const fixed = (actual, initial) => actual.forEach((row, i) => row.forEach((p, j) => {
  assert.equal(p.x, initial[i][j].x);
  if (boundary(i, j, actual.length - 1, row.length - 1)) assert.deepEqual(p, initial[i][j]);
}));
function fixture(nx, nt = nx) {
  const eta = Array.from({ length: nt + 1 }, (_, j) => j / nt);
  const exact = Array.from({ length: nx + 1 }, (_, i) => eta.map(v => {
    const x = clusteredX(i / nx); return { x, y: inverseY(x, v) };
  }));
  const nodes = exact.map((row, i) => row.map((p, j) => boundary(i, j, nx, nt) ? { ...p } : {
    x: p.x, y: (1 - eta[j]) * row[0].y + eta[j] * row[nt].y
      + .01 * Math.sin(Math.PI * i / nx) * Math.sin(Math.PI * eta[j]),
  }));
  return { nodes, directions: vertical(nodes), massFlows: masses(eta), eta };
}

test('harmonic streamfunction Newton preserves clustered straight flow without imposing harmonic xi', () => {
  const eta = [0, .03, .11, .29, .58, .84, 1];
  const nodes = Array.from({ length: 13 }, (_, i) => eta.map(y => ({ x: clusteredX(i / 12), y })));
  const data = { nodes, massFlows: masses(eta), directions: vertical(nodes) }, before = structuredClone(nodes);
  const system = createTransverseHarmonicGrid(data), result = smoothTransverseHarmonicGrid(system);
  assert.equal(result.converged, true, result.reason); assert.equal(system.quality(result.nodes).valid, true);
  assert.deepEqual(result.nodes, nodes); assert.deepEqual(nodes, before);
  const reference = solveHarmonicGridReference({ nodes: result.nodes, massFlows: data.massFlows }, { refinement: 1 });
  assert.ok(reference.maximum.tubeIntervals < 1e-11);
  assert.ok(reference.maximum.crosslineIntervals > .1, 'The independently prescribed clustered xi is intentionally not harmonic.');
});

test('direct harmonic residual converges to the exact physical quadratic streamfunction under refinement', t => {
  let previous = Infinity; const evidence = [];
  for (const n of [4, 8, 16]) {
    const data = fixture(n), before = structuredClone(data.nodes), system = createTransverseHarmonicGrid(data);
    const result = smoothTransverseHarmonicGrid(system, { tolerance: 1e-10 });
    assert.equal(result.converged, true, JSON.stringify({ n, reason: result.reason, last: result.history.at(-1) }));
    assert.equal(system.quality(result.nodes).valid, true); fixed(result.nodes, data.nodes); assert.deepEqual(data.nodes, before);
    let error = 0, initialError = 0;
    result.nodes.forEach((row, i) => row.forEach((p, j) => {
      error = Math.max(error, Math.abs(psi(p) - data.eta[j]));
      initialError = Math.max(initialError, Math.abs(psi(data.nodes[i][j]) - data.eta[j]));
    }));
    assert.ok(error < initialError / 10);
    assert.ok(error < .35 * previous, `Physical psi error ${error} after ${previous}`); previous = error;
    evidence.push({ n, error, initialError, iterations: result.history.length - 1,
      minCornerSine: system.quality(result.nodes).minCornerSine });
  }
  assert.ok(previous < 2e-4); t.diagnostic(JSON.stringify(evidence));
});

test('a separate fixed-geometry finite-element solve reproduces converged eta labels', t => {
  const data = fixture(8), system = createTransverseHarmonicGrid(data);
  const result = smoothTransverseHarmonicGrid(system, { tolerance: 1e-10 });
  assert.equal(result.converged, true, result.reason);
  // This assembles and solves the physical Q1 Laplace field on the final
  // fixed geometry. Only eta is relevant; xi remains an independent guide.
  const reference = solveHarmonicGridReference({ nodes: result.nodes, massFlows: data.massFlows }, { refinement: 1 });
  let error = 0;
  reference.values.forEach(row => row.forEach((p, j) => { error = Math.max(error, Math.abs(p.eta - data.eta[j])); }));
  assert.ok(error < 1e-9, `Independent discrete eta discrepancy ${error}`);
  assert.ok(reference.maximum.tubeIntervals < 1e-8);
  t.diagnostic(JSON.stringify({ etaError: error, tubeIntervals: reference.maximum.tubeIntervals,
    etaLinearRelativeResidual: reference.linear.eta.relativeResidual }));
});

test('harmonic streamfunction limits and failed or falsely certified linear solves cannot pass', () => {
  const data = fixture(5), system = createTransverseHarmonicGrid(data);
  const stopped = smoothTransverseHarmonicGrid(system, { maxIterations: 0 });
  assert.equal(stopped.converged, false); assert.equal(stopped.reason, 'iteration limit');
  assert.deepEqual(stopped.nodes, system.initial);
  const failed = smoothTransverseHarmonicGrid(system, { linearSolve: () => ({ converged: false }) });
  assert.equal(failed.converged, false); assert.match(failed.reason, /linear/); assert.deepEqual(failed.nodes, system.initial);
  const uncertified = smoothTransverseHarmonicGrid(system, {
    linearSolve: () => ({ x: new Float64Array(system.n), converged: true, relativeResidual: 0 }),
  });
  assert.equal(uncertified.converged, false); assert.match(uncertified.reason, /linear/); assert.deepEqual(uncertified.nodes, system.initial);
});

test('direct harmonic grid moves preserve boundaries, guide lines and original inputs', () => {
  const data = fixture(5, 4), before = structuredClone(data.nodes), system = createTransverseHarmonicGrid(data);
  const delta = Float64Array.from({ length: system.n }, (_, k) => .002 * Math.cos(k));
  const moved = system.move(system.initial, delta, .5);
  fixed(moved, data.nodes); assert.deepEqual(system.initial, before); assert.deepEqual(data.nodes, before);
  assert.notDeepEqual(moved, before);
  const offGuide = structuredClone(system.initial); offGuide[2][2].x += .001;
  assert.throws(() => system.evaluate(offGuide), /guide line/);
  const movedBoundary = structuredClone(system.initial); movedBoundary[0][2].y += .001;
  assert.throws(() => system.evaluate(movedBoundary), /boundaries/);
});
