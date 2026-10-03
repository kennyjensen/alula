// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { buildReliabilityCase } from '../scripts/validation/solver-reliability-cases.js';
import { coupledConvergenceSatisfied } from '../src/euler/streamtube-coupled-convergence.js';
import { coupledRefinementPlan } from '../src/euler/streamtube-coupled-refinement-assembly.js';
import { assessSolverResult } from '../scripts/validation/solver-reliability-acceptance.js';

test('the app worker uses and exports change convergence on a cold default NACA solve', async () => {
  const { caseData } = buildReliabilityCase({ preset: 'single', mode: 'streamtube-bl' });
  const messages = [];
  globalThis.self = { postMessage: m => { if (['result', 'error'].includes(m.type)) messages.push(m); } };
  try {
    await import('../src/worker/solver.js');
    await self.onmessage({ data: { id: 1, task: 'solve', caseData } });
    const message = messages.at(-1);
    assert.equal(message?.type, 'result', message?.message);
    const r = message.result;
    assert.equal(r.converged, true, r.reason);
    assert.equal(r.reason, 'solution changes');
    assert.equal(r.residualConverged, false);
    assert.equal(r.solverSettings.convergence, 'mses');
    assert.equal(r.checkpoint.continuation.convergence, 'mses');
    assert.equal(coupledConvergenceSatisfied(r, 1e-10), true);
    assert.deepEqual(r.diagnostics.convergence, r.convergence);
    assert.equal(r.mesh.quality.valid, true);
    assert.ok([r.cl, r.cd, r.cm].every(Number.isFinite));
    assert.ok(coupledRefinementPlan(caseData, r).nodeCount > 0);
    const assessment = assessSolverResult(caseData, r);
    assert.equal(assessment.passed, true, JSON.stringify(assessment));
    assert.equal(assessment.residualConverged, false);
    // A success flag without the change evidence cannot qualify for refinement.
    assert.throws(() => coupledRefinementPlan(caseData, { ...r, convergence: undefined }), /convergence tolerance/);
  } finally { delete globalThis.self; }
});
