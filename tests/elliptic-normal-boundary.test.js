import test from 'node:test';
import assert from 'node:assert/strict';
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { createNormalGraphBoundary } from '../src/geometry/normal-graph-boundary.js';

const makeCurve = (points, slopes) => ({ points, slopes });

test('normal-curve SLOR reduces to the literal horizontal Giles boundary update', () => {
  const nx = 8, nt = 4;
  const nodes = Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => ({
    x: i / nx + (i && i !== nx ? .02 * Math.sin(Math.PI * i / nx) * (1 + j / nt) : 0), y: j / nt,
  })));
  const base = { nodes, massFlows: Array(nt).fill(1 / nt), discretization: 'giles-1985' };
  const horizontal = createEllipticStreamtubeGrid({ ...base, boundaryConditions: { lower: 'giles-vertical', upper: 'giles-vertical' } });
  const curve = createEllipticStreamtubeGrid({ ...base, boundaryConditions: { lower: 'normal-curve', upper: 'normal-curve' },
    boundaryCurves: { lower: makeCurve(nodes.map(row => row[0]), Array(nx + 1).fill(0)),
      upper: makeCurve(nodes.map(row => row[nt]), Array(nx + 1).fill(0)) } });
  for (const omega of [.8, 1.3]) {
    assert.deepEqual(curve.sweep(nodes, omega).nodes, horizontal.sweep(nodes, omega).nodes);
    assert.deepEqual(curve.residuals(nodes), horizontal.residuals(nodes));
  }
});

test('normal farfields recover a rotated affine harmonic map with fixed ends', () => {
  const nx = 8, nt = 4, slope = .2, xi = [0, .04, .12, .25, .4, .57, .73, .88, 1];
  const eta = [0, .1, .3, .65, 1];
  const exact = xi.map(u => eta.map(v => ({ x: u - slope * v, y: slope * u + v })));
  const boundaryCurves = Object.fromEntries([['lower', 0], ['upper', nt]].map(([side, j]) =>
    [side, makeCurve(exact.map(row => row[j]), Array(nx + 1).fill(slope))]));
  const nodes = exact.map((row, i) => row.map((p, j) => {
    const d = i && i !== nx ? .015 * Math.sin(Math.PI * xi[i]) : 0;
    return { x: p.x + d, y: p.y + slope * d };
  }));
  const saved = structuredClone(nodes), savedCurves = structuredClone(boundaryCurves);
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: eta.slice(1).map((v, j) => v - eta[j]),
    streamwiseCoordinates: xi, discretization: 'giles-1985',
    boundaryConditions: { lower: 'normal-curve', upper: 'normal-curve' }, boundaryCurves });
  const result = smoothEllipticStreamtubeGrid(system, { maxSweeps: 300, tolerance: 1e-10 });
  assert.equal(result.converged, true, result.reason);
  assert.ok(result.history.at(-1).boundaryResidual < 1e-10);
  result.nodes.forEach((row, i) => row.forEach((p, j) => {
    assert.ok(Math.hypot(p.x - exact[i][j].x, p.y - exact[i][j].y) < 1e-9);
    if (!i || i === nx) assert.deepEqual(p, nodes[i][j]);
  }));
  assert.deepEqual(nodes, saved); assert.deepEqual(boundaryCurves, savedCurves);
});

test('curved normal-farfield SLOR refines toward an independently inverted conformal flow', t => {
  // z=w+a*w^2, derivative real part >=1 on the unit square, hence injective.
  // The top streamline is curved, with dy/dx=2a/(1+2au).
  const a = .1, map = (u, v) => ({ x: u + a * (u * u - v * v), y: v + 2 * a * u * v });
  const knots = Array.from({ length: 65 }, (_, i) => i / 64);
  const descriptor = makeCurve(knots.map(u => map(u, 1)), knots.map(u => 2 * a / (1 + 2 * a * u)));
  const boundary = createNormalGraphBoundary(descriptor), cases = [];
  for (const [nx, nt] of [[8, 4], [16, 8], [32, 16]]) {
    const nodes = Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => map(i / nx, j / nt)));
    const system = createEllipticStreamtubeGrid({ nodes, massFlows: Array(nt).fill(1 / nt), discretization: 'giles-1985',
      boundaryConditions: { upper: 'normal-curve' }, boundaryCurves: { upper: descriptor } });
    const result = smoothEllipticStreamtubeGrid(system, { maxSweeps: 800, tolerance: 1e-10 });
    assert.equal(result.converged, true, result.reason); assert.equal(result.quality.valid, true);
    let coordinate = 0, xi = 0, eta = 0, boundaryError = 0;
    result.nodes.forEach((row, i) => row.forEach((p, j) => {
      if (!i || i === nx || !j) assert.deepEqual(p, nodes[i][j]);
      if (j === nt) boundaryError = Math.max(boundaryError, Math.abs(p.y - boundary.evaluate(p.x).point.y));
      // Closed-form complex square root: w=(sqrt(1+4az)-1)/(2a).
      const re = 1 + 4 * a * p.x, im = 4 * a * p.y, realRoot = Math.sqrt((Math.hypot(re, im) + re) / 2);
      const u = (realRoot - 1) / (2 * a), v = im / (4 * a * realRoot);
      const reconstructed = map(u, v);
      assert.ok(Math.hypot(reconstructed.x - p.x, reconstructed.y - p.y) < 4e-15);
      xi = Math.max(xi, Math.abs(u - i / nx)); eta = Math.max(eta, Math.abs(v - j / nt));
      coordinate = Math.max(coordinate, Math.hypot(p.x - nodes[i][j].x, p.y - nodes[i][j].y));
    }));
    assert.ok(boundaryError < 1e-14);
    const record = { nx, nt, sweeps: result.history.length - 1, coordinate, xi, eta,
      boundaryError, residual: result.history.at(-1).residual };
    cases.push(record); t.diagnostic(JSON.stringify(record));
  }
  for (let k = 1; k < cases.length; k++) for (const key of ['coordinate', 'xi', 'eta'])
    assert.ok(cases[k][key] < .7 * cases[k - 1][key], `${key} did not refine: ${JSON.stringify(cases)}`);
});
