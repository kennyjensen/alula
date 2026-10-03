// SPDX-License-Identifier: GPL-2.0-or-later
// Refine a retained coupled root, then solve the unchanged full equations.
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { refineCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-refinement.js';
import { directChannelConservation } from '../tests/oracles/streamtube.js';
import { numericalSourceHashes, changedSources, sha256 } from './validation/provenance.js';
const args = process.argv.slice(2);
const [path = 'docs/coupled-isentropic-radius-recovery.json', output = 'docs/coupled-nested-isentropic-root.json', iterations = '24'] = args.filter(s => !s.startsWith('--'));
const factor = key => Number(args.find(s => s.startsWith(`--${key}-factor=`))?.split('=')[1] ?? 2);
assert.notEqual(path, output);
const source = JSON.parse(fs.readFileSync(path)), r = source.result, started = performance.now();
const wallFactorArgument = args.find(s => s.startsWith('--wall-normal-factor='));
assert.ok(!wallFactorArgument || !args.some(s => s.startsWith('--normal-factor=')), 'Select uniform or wall-adjacent normal refinement.');
const refinementControls = { streamwiseFactor: factor('streamwise'), ...(wallFactorArgument
  ? { normalSubdivisions: source.input.weights.map((row, g) => row.map((_, j) =>
    (g > 0 && j === 0) || (g < source.input.bodies.length && j === row.length - 1) ? factor('wall-normal') : 1)) }
  : { normalFactor: factor('normal') }) };
const report = { date: new Date().toISOString(), source: { path, sha256: sha256(path) },
  sourceHashes: numericalSourceHashes(['scripts/validate-streamtube-nested-coupled.js', 'tests/oracles/streamtube.js']),
  scope: 'Nested refinement of surface/cut intervals and uniformly or selectively partitioned streamtube masses; full coupled smooth isentropic equations. Independent physical acceptance remains required.',
  physicalAcceptance: false, failures: [], refinementControls,
  controls: { maxIterations: Number(iterations), tolerance: 1e-10, stepMethod: 'dogleg', initialTrustRadius: 1 } };
assert.ok(Number.isInteger(report.controls.maxIterations) && report.controls.maxIterations >= 0 && report.controls.maxIterations <= 32);
const save = () => fs.writeFileSync(output, JSON.stringify(report, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v instanceof Map ? [...v] : v, 2) + '\n');
try {
  assert.equal(r.converged, true); assert.equal(source.input.streamwiseMode, 'isentropic');
  const parent = createCoupledStreamtubeBody(source.input, { ...source.options, initialEuler: r.flow, initialBL: r.x.slice(source.initial.ne) });
  assert.ok(Math.max(...Array.from(parent.residual(parent.initial), Math.abs)) < 1e-10);
  const refined = refineCoupledStreamtubeBody(source.input, parent, report.refinementControls), system = refined.system;
  report.input = refined.input; report.options = refined.options; report.refinement = refined.diagnostics;
  report.initial = { ne: system.ne, x: system.initial, families: system.evaluate(system.initial).families,
    initialEuler: refined.initialEuler, initialBL: refined.initialBL }; save();
  console.log(JSON.stringify({ initialization: report.refinement, families: report.initial.families }));
  report.result = solveCoupledStreamtubeBody(system, { ...report.controls, onIteration: h => console.log(JSON.stringify(h)) });
  report.conservation = report.result.flow.nodes.map((nodes, g) => directChannelConservation({ nodes,
    sections: report.result.flow.sections.map(row => row[g]), cells: report.result.flow.cells.map(row => row[g]) }, system.euler.conditions.gamma));
  const { gamma, mach, pInf } = system.euler.conditions, p0 = pInf * (1 + .5 * (gamma - 1) * mach ** 2) ** (gamma / (gamma - 1));
  report.totalPressureError = Math.max(...report.result.flow.sections.flat(2).map(s =>
    Math.abs(s.p * (1 + .5 * (gamma - 1) * s.q ** 2 / (gamma * s.p / s.rho)) ** (gamma / (gamma - 1)) / p0 - 1)));
  if (!report.result.converged) report.failures.push(`Coupled equations: ${report.result.reason}`);
  for (const c of report.conservation) for (const key of ['maxLocal', 'total', 'internalCancellation'])
    for (const i of key === 'internalCancellation' ? [0, 1, 2, 3] : [0, 3])
      if (Math.abs(c[key][i]) > 2e-9) report.failures.push(`Independent balance failed: ${key}[${i}]`);
  if (report.totalPressureError > 1e-9) report.failures.push('Independent total-pressure preservation failed.');
} catch (error) { report.failures.push(error.message); report.error = error.stack; report.cause = error.cause; }
report.seconds = (performance.now() - started) / 1000; report.changedSources = changedSources(report.sourceHashes);
if (report.changedSources.length) report.failures.push('Numerical source changed during calculation.');
report.status = report.failures.length ? 'unresolved nested coupled refinement' : 'nested coupled equation root; physical refinement gates pending'; save();
console.log(JSON.stringify({ output, status: report.status, seconds: report.seconds, families: report.result?.families, failures: report.failures }));
if (report.failures.length) process.exitCode = 1;
