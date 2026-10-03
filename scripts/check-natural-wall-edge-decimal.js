// SPDX-License-Identifier: GPL-2.0-or-later
// Independent arithmetic and chart checks for the failed thin-wake row.
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { numericalSourceHashes, changedSources, sha256 } from './validation/provenance.js';
const source = 'docs/current-natural-graded-wall-checked-seed.json';
const output = 'docs/current-natural-wall-edge-decimal.json';
assert.ok(!fs.existsSync(output));
const sourceHash = sha256(source), saved = JSON.parse(fs.readFileSync(source)), f = saved.seed;
assert.deepEqual(changedSources(saved.sourceHashes), []);
const sourceHashes = numericalSourceHashes(['scripts/check-natural-wall-edge-decimal.js', 'tests/oracles/decimal-streamtube-edge.py']);
const start = performance.now(), s = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
const x = s.initial.slice(), phases = s.bl.snapshotActive(), v = s.evaluate(x);
assert.deepEqual(v.families, saved.initial.families);
const d = x.map((v, i) => (i < s.ne ? .001 : Math.max(.01, Math.abs(v))) * Math.sin(.43 * i + .2));
const station = s.bl.stations.find(p => p.kind === 'wake' && p.i === 321);
assert.ok(station && station.i < s.euler.layout.nx);
const row = s.ne + 4 * station.id + 3, matrix = s.jacobian(x), analytic = sparseProduct(matrix, d)[row];
const maps = s.euler.geometryDerivatives(x.subarray(0, s.ne));
const tangent = (g, i, j) => {
  const out = { x: 0, y: 0 };
  for (const [col, p] of maps[g][i][j]) { assert.ok(col < s.ne); out.x += p.x * d[col]; out.y += p.y * d[col]; }
  return out;
};
const point = p => [p.x, p.y];
let chartError = 0;
const cells = ['upper', 'lower'].map(side => {
  const g = side === 'lower' ? station.body : station.body + 1;
  const tube = side === 'lower' ? s.euler.layout.tubes[g] - 1 : 0;
  const indices = [station.i - 1, station.i, station.i + 1];
  const nodes = [tube, tube + 1].map(j => indices.map(i => point(v.outer.nodes[g][i][j])));
  const nodeTangents = [tube, tube + 1].map(j => indices.map(i => point(tangent(g, i, j))));
  const baseMass = v.outer.allocation.groups[g][tube].massFlow;
  for (const h of [1e-3, 5e-4]) {
    const plus = s.euler.decode(x.subarray(0, s.ne).map((v, i) => v + h * d[i]));
    const minus = s.euler.decode(x.subarray(0, s.ne).map((v, i) => v - h * d[i]));
    [tube, tube + 1].forEach((j, bank) => indices.forEach((i, k) => {
      for (const [coord, c] of [['x', 0], ['y', 1]]) chartError = Math.max(chartError,
        Math.abs((plus.nodes[g][i][j][coord] - minus.nodes[g][i][j][coord]) / (2 * h) - nodeTangents[bank][k][c]));
    }));
    assert.equal(plus.allocation.groups[g][tube].massFlow, baseMass);
    assert.equal(minus.allocation.groups[g][tube].massFlow, baseMass);
  }
  const c = v.outer.cells[station.i - 1][g][tube];
  return { nodes, nodeTangents, mass: baseMass, massTangent: 0, densities: c.states.map(q => q.rho),
    logDensityTangents: [station.i - 1, station.i].map(i => d[s.euler.layout.densityIndex(i, g, tube)]) };
});
assert.ok(chartError < 1e-11);
const input = { cells, gamma: s.euler.conditions.gamma, h0: s.euler.conditions.h0,
  ue: x[row], ueTangent: d[row] };
const oracle = JSON.parse(execFileSync('python3', ['tests/oracles/decimal-streamtube-edge.py'], { input: JSON.stringify(input), encoding: 'utf8' }));
const baselineError = Math.abs(Number(oracle.residual) - v.residual[row]);
assert.ok(baselineError < 5e-10);
const checks = oracle.checks.map(q => ({ ...q, error: Math.abs(Number(q.derivative) - analytic) / Math.max(1, Math.abs(analytic)) }));
assert.ok(checks.every(q => q.error < 1e-8));
assert.deepEqual(s.initial, x); assert.deepEqual(s.bl.snapshotActive(), phases); assert.deepEqual(s.evaluate(x).families, v.families);
assert.deepEqual(changedSources(sourceHashes), []); assert.equal(sha256(source), sourceHash);
const report = { date: new Date().toISOString(), physicalAcceptance: false, passed: true, sourceHashes,
  source: { path: source, sha256: sourceHash },
  scope: 'One failed mixed-direction wake edge row. Validate its affine coordinate tangents independently at large steps, then evaluate the signed-bank-bend stencil in 80-digit Decimal arithmetic. Complete coupled matrix supplies the derivative being tested; no nonlinear solve or physical acceptance.',
  station, row, analytic, chartError, baselineError, input, oracle, checks, parentUnchanged: true, seconds: (performance.now() - start) / 1000 };
fs.writeFileSync(output, JSON.stringify(report) + '\n');
console.log(JSON.stringify({ output, passed: true, seconds: report.seconds, row, analytic, chartError, baselineError, checks }));
