// SPDX-License-Identifier: GPL-2.0-or-later
// Continue the accepted two-element state without rebuilding its initializer.
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { numericalSourceHashes, sha256, changedSources } from './validation/provenance.js';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { directChannelConservation } from '../tests/oracles/streamtube.js';

const file = 'docs/accepted-grid-euler-handoff.json', saved = JSON.parse(fs.readFileSync(file)), started = performance.now();
const report = { date: new Date().toISOString(), input: { path: file, sha256: sha256(file) },
  sourceHashes: numericalSourceHashes(['scripts/continue-accepted-grid-euler.js', 'tests/oracles/streamtube.js']),
  controls: { maxIterations: 24, tolerance: 1e-10, stepMethod: 'dogleg',
    initialTrustRadius: saved.final.history.at(-1).trustRadius }, physicalAcceptance: false,
  scope: 'Bounded continuation of the actual accepted two-element Euler chart. No panel tracing, smoothing, or viscosity. Conservative polygon fluxes are independently measured.' };
try {
  const { input, initialEuler } = saved.restart, system = createStreamtubeBodySystem(input);
  const initial = system.adoptGeometry(Float64Array.from(initialEuler.x), initialEuler.nodes);
  report.initialDiagnostics = system.evaluate(initial).diagnostics;
  assert.ok(Math.abs(report.initialDiagnostics.residual - saved.final.diagnostics.residual) < 1e-11);
  const flow = solveStreamtubeBody(system, { initial, ...report.controls,
    onIteration: h => console.log(JSON.stringify(h)) });
  report.final = { converged: flow.converged, reason: flow.reason, history: flow.history,
    diagnostics: flow.diagnostics, linearDiagnostics: flow.linearDiagnostics,
    rejectedTrials: flow.rejectedTrials, quality: streamtubeMeshSnapshot({ system, nodes: flow.nodes }).quality };
  report.conservation = flow.nodes.map((nodes, g) => directChannelConservation({ nodes,
    sections: flow.sections.map(row => row[g]), cells: flow.cells.map(row => row[g]) }, system.conditions.gamma));
  report.restart = { input, initialEuler: { x: Array.from(flow.x), nodes: flow.nodes } };
  report.surfaces = flow.surfaces;
  report.passed = flow.converged && report.final.quality.valid && report.conservation.every(c =>
    ['maxLocal', 'total', 'internalCancellation'].every(key => c[key].every(v => Math.abs(v) < 2e-9)));
} catch (error) { report.error = error.message; report.passed = false; }
assert.deepEqual(changedSources(report.sourceHashes), []); assert.equal(sha256(file), report.input.sha256);
report.seconds = (performance.now() - started) / 1000;
fs.writeFileSync('docs/accepted-grid-euler-continuation.json', JSON.stringify(report, null, 2) + '\n');
console.log(JSON.stringify({ passed: report.passed, error: report.error, final: report.final, seconds: report.seconds }));
if (!report.passed) process.exitCode = 1;
