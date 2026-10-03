// SPDX-License-Identifier: GPL-2.0-or-later
// Recover the accepted two-element chart with the feasible Newton-ray step,
// then initialize all BL/wake displacements if Euler reaches its root.
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { numericalSourceHashes, sha256, changedSources } from './validation/provenance.js';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { extendStreamtubeDisplacement } from '../src/euler/streamtube-displacement.js';
import { initializeStreamtubeDensities } from '../src/euler/streamtube-initial-state.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { directChannelConservation } from '../tests/oracles/streamtube.js';

const file = 'docs/accepted-grid-euler-handoff.json', saved = JSON.parse(fs.readFileSync(file)), started = performance.now();
const report = { date: new Date().toISOString(), input: { path: file, sha256: sha256(file) },
  sourceHashes: numericalSourceHashes(['scripts/check-accepted-flow-recovery.js', 'tests/oracles/streamtube.js']),
  eulerControls: { maxIterations: 40, tolerance: 1e-10, stepMethod: 'dogleg', initialTrustRadius: saved.final.history.at(-1).trustRadius },
  blOptions: { reynolds: 1e6, ncrit: 9, tripFractions: [[.05, .05], [.05, .05]] },
  coupledControls: { maxIterations: 16, tolerance: 1e-8 }, physicalAcceptance: false,
  scope: 'Same accepted two-element chart, restored three-update state and unchanged conservation equations. Feasible Newton ray remains inside the trust radius. Full four-surface/two-wake initialization/solve attempted only after Euler convergence; fixed material trips, no shocks or confluence.' };
const save = () => fs.writeFileSync('docs/accepted-flow-recovery.json', JSON.stringify(report, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v, 2) + '\n');
const conservation = (flow, gamma) => flow.nodes.map((nodes, g) => directChannelConservation({ nodes,
  sections: flow.sections.map(row => row[g]), cells: flow.cells.map(row => row[g]) }, gamma));
try {
  const { input, initialEuler } = saved.restart, system = createStreamtubeBodySystem(input);
  const initial = system.adoptGeometry(Float64Array.from(initialEuler.x), initialEuler.nodes);
  assert.ok(Math.abs(system.evaluate(initial).diagnostics.residual - saved.final.diagnostics.residual) < 1e-11);
  const flow = solveStreamtubeBody(system, { initial, ...report.eulerControls,
    onIteration: h => console.log(JSON.stringify({ stage: 'euler', ...h })) });
  report.euler = { converged: flow.converged, reason: flow.reason, history: flow.history,
    diagnostics: flow.diagnostics, linearDiagnostics: flow.linearDiagnostics, rejectedTrials: flow.rejectedTrials,
    quality: streamtubeMeshSnapshot({ system, nodes: flow.nodes }).quality,
    conservation: conservation(flow, system.conditions.gamma) };
  report.eulerRestart = { input, initialEuler: { x: flow.x, nodes: flow.nodes } };
  report.eulerSeconds = (performance.now() - started) / 1000; save();
  if (!flow.converged) throw new Error(`Euler remains unconverged: ${flow.reason}`);
  assert.ok(report.euler.conservation.every(c => ['maxLocal', 'total', 'internalCancellation'].every(key => c[key].every(v => Math.abs(v) < 2e-9))));
  const seeded = createCoupledStreamtubeBody(input, { ...report.blOptions, initialEuler: flow });
  const bl = seeded.initial.slice(seeded.ne), state = seeded.initial.slice(0, seeded.ne), euler = seeded.euler;
  report.blSeed = { x: bl, initialization: seeded.initialization, surfaces: seeded.bl.surfaces, wakes: seeded.bl.wakes };
  euler.setDisplacement(seeded.bl.thicknesses(bl));
  const extended = extendStreamtubeDisplacement(euler, state);
  report.displacementMesh = streamtubeMeshSnapshot({ system: euler, nodes: extended }); save();
  assert.ok(report.displacementMesh.quality.valid, 'Initial BL displacement folds the extended grid');
  const gas = initializeStreamtubeDensities(euler, euler.adoptGeometry(state, extended)), value = euler.evaluate(gas);
  const coupled = createCoupledStreamtubeBody(input, { ...report.blOptions, initialBL: bl, initialEuler: { ...value, x: gas } });
  report.coupledInitial = { x: coupled.initial, ne: coupled.ne, families: coupled.evaluate(coupled.initial).families };
  report.coupledRestart = { input, options: report.blOptions, initialBL: bl, initialEuler: { ...value, x: gas } }; save();
  const result = solveCoupledStreamtubeBody(coupled, { ...report.coupledControls,
    onIteration: h => console.log(JSON.stringify({ stage: 'coupled', ...h })) });
  report.coupledResult = result; report.coupledConservation = conservation(result.flow, euler.conditions.gamma);
  report.passed = result.converged && report.coupledConservation.every(c =>
    ['maxLocal', 'total', 'internalCancellation'].every(key => c[key].every(v => Math.abs(v) < 2e-9)));
} catch (error) { report.error = error.message; report.passed = false; }
assert.deepEqual(changedSources(report.sourceHashes), []); assert.equal(sha256(file), report.input.sha256);
report.seconds = (performance.now() - started) / 1000; save();
console.log(JSON.stringify({ passed: report.passed, error: report.error, euler: report.euler?.diagnostics,
  coupled: report.coupledResult?.families, seconds: report.seconds }));
if (!report.passed) process.exitCode = 1;
