// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { solveStreamtubeIses } from '../src/euler/streamtube-ises-update.js';
import { smoothStreamtubePassages, recoverStreamtubePassageGrid } from '../src/euler/streamtube-passage-startup.js';
import { assertConvexStreamtubeGrid } from '../src/geometry/streamtube-convex-step.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { buildReliabilityCase } from '../scripts/validation/solver-reliability-cases.js';
import { solveStreamtubeAssembly } from '../src/euler/streamtube-result.js';

test('passage coordinates smooth interiors while preserving all prescribed boundaries and inputs', () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 3 });
  const system = createStreamtubeBodySystem(input), initial = system.initial, nodes = system.decode(initial).nodes;
  const before = structuredClone({ input, nodes, initial });
  const r = smoothStreamtubePassages({ system, initial, nodes });
  assert.equal(r.converged, true, JSON.stringify(r.reports));
  assertConvexStreamtubeGrid(r.nodes);
  r.nodes.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => {
    if (!i || i === group.length - 1 || !j || j === row.length - 1) assert.deepEqual(p, nodes[g][i][j]);
  })));
  assert.deepEqual({ input, nodes, initial }, before);
  const stopped = { converged: false, history: [{ iteration: 1 }], lastRejectedStep: null };
  assert.equal(recoverStreamtubePassageGrid(stopped, { input }, { maxIterations: 40 }), stopped);
});

test('three-element inlet contact admits a compensated Newton step and an exact-replaying root', () => {
  const cp = JSON.parse(gunzipSync(readFileSync(new URL('./fixtures/three-element-inlet-contact.json.gz', import.meta.url))));
  const before = structuredClone(cp);
  const blocked = solveStreamtubeIses(undefined, { ...cp.continuation, resume: cp, maxIterations: 1 });
  assert.equal(blocked.converged, false);
  assert.match(blocked.reason, /inlet grid repair/);
  const resume = structuredClone(cp);
  resume.continuation.tangentialGridRecovery = true;
  const result = solveStreamtubeIses(undefined, { ...resume.continuation, resume, maxIterations: 8 });
  assert.equal(result.converged, true, result.reason);
  assert(result.diagnostics.residual < 1e-10);
  assert.equal(result.history[1].step, 1);
  assert(result.history[1].residual < 3e-4);
  assert.equal(result.history[1].coordinateRecovery.method, 'compensated-tangential-coordinates');
  assertConvexStreamtubeGrid(result.nodes);
  assert.equal(result.solverInput.mach, .2);
  assert.equal(result.solverInput.alpha, 2);
  assert.deepEqual(result.solverInput.bodies, cp.input.bodies);
  assert(result.linearDiagnostics.solves > result.history.length - 1);
  const replay = solveStreamtubeIses(undefined, { ...result.checkpoint.continuation,
    resume: result.checkpoint, maxIterations: 0 });
  assert.equal(replay.converged, true);
  assert.deepEqual(replay.checkpoint.residual, result.checkpoint.residual);
  assert.deepEqual(cp, before);
  assert.throws(() => solveStreamtubeIses(undefined, { ...result.checkpoint.continuation,
    tangentialGridRecovery: false, resume: result.checkpoint, maxIterations: 0 }), /controls do not match/);
});

test('default flap with automatic half-surface end counts escapes inlet contact and closes Euler', () => {
  const { caseData } = buildReliabilityCase({ preset: 'flap', mode: 'streamtube-bl' });
  const before = structuredClone(caseData);
  const result = solveStreamtubeAssembly(caseData, { maxIterations: 20 });
  assert.equal(result.status, 'research-converged', result.diagnostics.reason);
  assert.equal(result.mesh.cells.length, 3531);
  assert.equal(result.mesh.quality.valid, true);
  assert(result.diagnostics.equationResidual <= 1e-10);
  assert(result.flow.history.some(h => h.coordinateRecovery?.method === 'compensated-tangential-coordinates'
    && h.step > .1), 'The ordinary solve must use the available coordinate freedom at inlet contact.');
  for (const h of result.flow.history.slice(1)) {
    const merit = h.residualDecrease;
    assert(merit.converged || merit.afterSquaredNorm <= merit.allowedSquaredNorm);
  }
  // Rebuild the physical equations independently of the convergence flag.
  const system = createStreamtubeBodySystem(result.solverInput);
  const state = system.adoptGeometry(result.flow.x, result.flow.nodes);
  const replay = system.evaluate(state);
  assertConvexStreamtubeGrid(replay.nodes);
  assert(replay.diagnostics.residual <= 1e-10);
  assert.deepEqual(caseData, before);
});
