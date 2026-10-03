// SPDX-License-Identifier: GPL-2.0-or-later
// Change only the edge matching of a retained coupled root, with a bounded
// simultaneous restart and independently integrated conservation diagnostics.
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { directChannelConservation } from '../tests/oracles/streamtube.js';
import { numericalSourceHashes, changedSources, sha256 } from './validation/provenance.js';

const [source = 'docs/coupled-isentropic-radius-recovery.json', output = 'docs/coupled-section-velocity-root.json', iterations = '8'] = process.argv.slice(2);
assert.notEqual(source, output);
const old = JSON.parse(fs.readFileSync(source)), started = performance.now();
const report = { date: new Date().toISOString(), source: { path: source, sha256: sha256(source) },
  sourceHashes: numericalSourceHashes(['scripts/validate-streamtube-edge-matching.js', 'tests/oracles/streamtube.js']),
  scope: 'Section-speed edge matching on the same grid, trips, gas and Reynolds number as a retained root. Full simultaneous Euler/BL/displacement/wake restart; physical profile/force acceptance remains independent.',
  physicalAcceptance: false, failures: [], input: old.input, options: { ...old.options, edgeMatching: 'section-velocity' },
  controls: { maxIterations: Number(iterations), tolerance: 1e-10, stepMethod: 'dogleg', initialTrustRadius: 1 } };
assert.ok(Number.isInteger(report.controls.maxIterations) && report.controls.maxIterations >= 0 && report.controls.maxIterations <= 24);
const save = () => fs.writeFileSync(output, JSON.stringify(report, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v instanceof Map ? [...v] : v, 2) + '\n');
try {
  assert.ok(old.result?.flow && old.result?.x, 'A retained complete iterate is required.');
  const continuing = old.options?.edgeMatching === 'section-velocity';
  if (!continuing) assert.equal(old.result.converged, true);
  const parent = createCoupledStreamtubeBody(old.input, { ...old.options, initialEuler: old.result.flow, initialBL: old.result.x.slice(old.initial.ne) });
  assert.ok(parent.admissible(parent.initial));
  if (!continuing) assert.ok(Math.max(...Array.from(parent.residual(parent.initial), Math.abs)) < 1e-10);
  report.continuingSectionMatching = continuing;
  const system = createCoupledStreamtubeBody(report.input, { ...report.options, initialEuler: old.result.flow, initialBL: old.result.x.slice(old.initial.ne) });
  assert.ok(system.admissible(system.initial));
  report.initial = { ne: system.ne, x: system.initial, families: system.evaluate(system.initial).families }; save();
  report.result = solveCoupledStreamtubeBody(system, { ...report.controls, onIteration: h => console.log(JSON.stringify(h)) });
  report.conservation = report.result.flow.nodes.map((nodes, g) => directChannelConservation({ nodes,
    sections: report.result.flow.sections.map(row => row[g]), cells: report.result.flow.cells.map(row => row[g]) }, system.euler.conditions.gamma));
  if (!report.result.converged) report.failures.push(`Coupled equations: ${report.result.reason}`);
  for (const c of report.conservation) for (const key of ['maxLocal', 'total', 'internalCancellation'])
    for (const i of system.euler.conditions.streamwiseMode === 'momentum' || key === 'internalCancellation' ? [0, 1, 2, 3] : [0, 3])
      if (Math.abs(c[key][i]) > 2e-9) report.failures.push(`Independent conservation: ${key}[${i}]`);
  if (system.euler.conditions.streamwiseMode === 'isentropic') {
    const { gamma, mach, pInf } = system.euler.conditions;
    const p0 = pInf * (1 + .5 * (gamma - 1) * mach ** 2) ** (gamma / (gamma - 1));
    report.totalPressureError = Math.max(...report.result.flow.sections.flat(2).map(s =>
      Math.abs(s.p * (1 + .5 * (gamma - 1) * s.q ** 2 / (gamma * s.p / s.rho)) ** (gamma / (gamma - 1)) / p0 - 1)));
    if (report.totalPressureError > 1e-9) report.failures.push('Independent section total-pressure error.');
  }
} catch (error) { report.failures.push(error.message); report.error = error.stack; }
report.seconds = (performance.now() - started) / 1000;
report.changedSources = changedSources(report.sourceHashes);
if (report.changedSources.length) report.failures.push('Numerical source changed during calculation.');
report.status = report.failures.length ? 'unresolved section-speed coupled restart' : 'section-speed coupled root; physical accuracy pending';
save(); console.log(JSON.stringify({ output, status: report.status, seconds: report.seconds, families: report.result?.families, failures: report.failures }));
if (report.failures.length) process.exitCode = 1;
