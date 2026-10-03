// SPDX-License-Identifier: GPL-2.0-or-later
// Bounded simultaneous solve of the saved admissible default two-element
// initializer. Every Euler, BL and wake unknown participates in each step.
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { numericalSourceHashes, changedSources, sha256 } from './validation/provenance.js';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { directChannelConservation } from '../tests/oracles/streamtube.js';

const file = 'docs/accepted-coupled-initialization.json', saved = JSON.parse(fs.readFileSync(file));
const started = performance.now(), report = { date: new Date().toISOString(), input: { path: file, sha256: sha256(file) },
  sourceHashes: numericalSourceHashes(['scripts/check-accepted-coupled-trust-region.js', 'tests/oracles/streamtube.js']),
  controls: { stepMethod: 'dogleg', initialTrustRadius: 1, maxIterations: 24, tolerance: 1e-8 },
  physicalAcceptance: false, scope: 'Full 7079-unknown moving Euler/grid/BL/wake solve from a saved admissible guess, with no reinitialization or changed residual equations.' };
try {
  const seed = saved.restart, system = createCoupledStreamtubeBody(seed.input, {
    ...seed.options, initialBL: Float64Array.from(seed.initialBL), initialEuler: seed.initialEuler });
  report.initial = { n: system.n, ne: system.ne, stations: system.bl.stations.length,
    surfaces: system.bl.surfaces.length, wakes: system.bl.wakes.length, families: system.evaluate(system.initial).families };
  report.result = solveCoupledStreamtubeBody(system, { ...report.controls, onIteration: h => console.log(JSON.stringify(h)) });
  const r = report.result, state = system.evaluate(r.x);
  report.restart = { input: seed.input, options: seed.options, initialBL: Array.from(r.x.slice(system.ne)),
    initialEuler: { x: Array.from(r.x.slice(0, system.ne)), nodes: state.outer.nodes, undisplacedNodes: state.outer.undisplacedNodes } };
  report.conservation = r.flow.nodes.map((nodes, g) => directChannelConservation({ nodes,
    sections: r.flow.sections.map(row => row[g]), cells: r.flow.cells.map(row => row[g]) }, system.euler.conditions.gamma));
  report.passed = r.converged && report.conservation.every(c => ['maxLocal', 'total', 'internalCancellation'].every(k => c[k].every(v => Math.abs(v) < 2e-9)));
} catch (error) { report.error = error.stack ?? error.message; report.passed = false; }
assert.deepEqual(changedSources(report.sourceHashes), []); assert.equal(sha256(file), report.input.sha256);
report.seconds = (performance.now() - started) / 1000;
fs.writeFileSync('docs/accepted-coupled-trust-region.json', JSON.stringify(report, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v, 2) + '\n');
console.log(JSON.stringify({ seconds: report.seconds, initial: report.initial, passed: report.passed,
  reason: report.result?.reason, families: report.result?.families, quality: report.result?.mesh.quality, error: report.error }));
if (!report.passed) process.exitCode = 1;
