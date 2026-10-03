// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { pathToFileURL } from 'node:url';
import { resolve } from 'node:path';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { changeCoupledDissipationOrder, recoverCoupledDissipationOrder,
  coupledFirstOrderRecoveryEligible } from '../src/euler/tests/streamtube-coupled-dissipation-recovery.js';
import { createQuadSolveProgress } from '../src/ui/quad-solve-progress.js';

let cached;
function source() {
  cached ??= solveCoupledStreamtubeIses({ ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2, mach: .2 }),
    streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 },
    upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } }, {
    edgeMatching: 'section-velocity', blThermodynamics: 'historical-common-isentrope', transitionMode: 'automatic',
    maxIterations: 12, tolerance: 1e-10, stepAcceptance: 'admissible', iterationRecovery: false });
  assert.ok(cached.converged, cached.reason);
  return structuredClone(cached.checkpoint);
}

test('order switch preserves physical unknowns, BL, transition, geometry and input ownership', () => {
  const cp = source(), before = structuredClone(cp), first = changeCoupledDissipationOrder(cp, -1);
  assert.deepEqual(cp, before);
  assert.equal(first.restart.input.upwind.mucon, -1);
  assert.deepEqual(first.restart.initialEuler, cp.restart.initialEuler);
  assert.deepEqual(first.restart.initialBL, cp.restart.initialBL);
  assert.deepEqual(first.restart.options, cp.restart.options);
  const restored = changeCoupledDissipationOrder(first, 1);
  assert.deepEqual(restored, cp);
  assert.throws(() => changeCoupledDissipationOrder(cp, -.3), /stable/);
  assert.throws(() => changeCoupledDissipationOrder(cp, -1.2), /only the sign/);
  const corrupt = structuredClone(cp); corrupt.families.euler = 5;
  assert.throws(() => changeCoupledDissipationOrder(corrupt, -1), /replay/);
});

test('real same-state recovery returns only a rechecked second-order root and announces both phases', () => {
  const cp = source(), phases = [];
  const r = recoverCoupledDissipationOrder(cp, { maxIterations: 2, dissipationEnhancement: false,
    iterationRecovery: false, onPhase: info => phases.push(info.phase) });
  assert.deepEqual(phases, ['first-order-initialization', 'second-order-restoration']);
  assert.ok(r.accepted);
  assert.equal(r.result.checkpoint.restart.input.upwind.mucon, 1);
  assert.equal(r.result.checkpoint.restart.input.upwind.mcrit, .99);
  assert.ok(r.result.residual.every(v => Math.abs(v) <= 1e-10));
  const sentinel = Error('observer cancellation');
  assert.throws(() => recoverCoupledDissipationOrder(cp, { onPhase: () => { throw sentinel; } }), e => e === sentinel);
});

test('only failed admissible near-sonic second-order stages are eligible', () => {
  const t = { checkpoint: source(), converged: false, mesh: { quality: { valid: true } },
    flow: { diagnostics: { maxMach: 1.1 } }, residual: [1] };
  assert.ok(coupledFirstOrderRecoveryEligible(t));
  for (const patch of [{ converged: true }, { mesh: { quality: { valid: false } } },
    { flow: { diagnostics: { maxMach: .4 } } }, { residual: [NaN] }])
    assert.equal(coupledFirstOrderRecoveryEligible({ ...t, ...patch }), false);
  t.checkpoint.restart.input.upwind.mucon = -1;
  assert.equal(coupledFirstOrderRecoveryEligible(t), false);
});

// Exercise acceptance independently of whether a small numerical fixture
// happens to converge: a converged first-order root is never sufficient.
test('failed second-order restoration cannot escape as an accepted result', async () => {
  const file = resolve('src/euler/tests/streamtube-coupled-dissipation-recovery.js');
  let text = fs.readFileSync(file, 'utf8');
  const stub = `export function solveCoupledStreamtubeIses(_, o) {
    return { checkpoint:o.resume, converged:o.resume.restart.input.upwind.mucon < 0,
      conditions:{mach:o.resume.restart.input.mach}, mesh:{quality:{valid:true}}, residual:[1],
      history:[{},{}], linearDiagnostics:{solves:1}, reason:'restoration stalled' };
  }`;
  const url = 'data:text/javascript;base64,' + Buffer.from(stub).toString('base64');
  text = text.replace(/from '([^']+)'/g, (_, p) => `from '${p === '../streamtube-coupled-ises.js' ? url : new URL(p, pathToFileURL(file)).href}'`);
  const module = await import('data:text/javascript;base64,' + Buffer.from(text).toString('base64'));
  const r = module.recoverCoupledDissipationOrder(source(), { maxIterations: 2 });
  assert.equal(r.accepted, false); assert.equal(r.result, undefined);
  assert.equal(r.diagnostics.phases.length, 2);
});

test('progress distinguishes temporary first order, restoration, and ordinary Mach continuation', () => {
  const progress = createQuadSolveProgress({ mach: .74 });
  const a = progress.update({ stage: 'coupled-mach', mach: .72,
    dissipationRecovery: { phase: 'first-order-initialization' } }, { stageChange: true });
  assert.match(JSON.stringify(a), /Temporary first-order dissipation/);
  const b = progress.update({ stage: 'coupled-mach', mach: .72,
    dissipationRecovery: { phase: 'second-order-restoration' } }, { stageChange: true });
  assert.match(JSON.stringify(b), /Restore second-order dissipation/);
  assert.match(JSON.stringify(progress.update({ stage: 'coupled-mach', mach: .73 }, { stageChange: true })), /Mach continuation/);
});

test('automatic controller bounds unsuccessful recovery and retains the converged source', async () => {
  const file = resolve('src/euler/tests/streamtube-coupled-automatic.js'), base = pathToFileURL(file);
  const realSolver = new URL('../streamtube-coupled-ises.js', base).href;
  const solverStub = `import { solveCoupledStreamtubeIses as solve } from '${realSolver}';
    export function solveCoupledStreamtubeIses(input, options) {
      const r=solve(input, options);
      if(options.maxIterations > 0) { r.converged=false; r.flow.diagnostics.maxMach=1.1; }
      return r;
    }`;
  const recoveryStub = `export { coupledFirstOrderRecoveryEligible } from '${new URL('./streamtube-coupled-dissipation-recovery.js', base).href}';
    export function recoverCoupledDissipationOrder(cp, options) {
      options.onPhase({phase:'first-order-initialization'});
      return {accepted:false, diagnostics:{reason:'test restoration rejection'}};
    }`;
  const data = s => 'data:text/javascript;base64,' + Buffer.from(s).toString('base64');
  const text = fs.readFileSync(file, 'utf8').replace(/from '([^']+)'/g, (_, p) =>
    `from '${p === '../streamtube-coupled-ises.js' ? data(solverStub)
      : p === './streamtube-coupled-dissipation-recovery.js' ? data(recoveryStub) : new URL(p, base).href}'`);
  const module = await import(data(text)), cp = source(), before = structuredClone(cp);
  const options = { initialCheckpoint: cp, stageMaxIterations: 1, maxStages: 3, maxSubdivisions: 3,
    maxFirstOrderRecoveries: 1, iterationRecovery: false };
  const r = module.solveCoupledStreamtubeAutomatic(.21, options);
  assert.equal(r.converged, false);
  assert.equal(r.firstOrderRecovery.attempts.length, 1);
  assert.equal(r.continuation.attempts.length, 3);
  assert.deepEqual(r.checkpoint, before);
  assert.deepEqual(cp, before);
  const sentinel = Error('cancel recovery');
  assert.throws(() => module.solveCoupledStreamtubeAutomatic(.21, { ...options,
    onStage: e => { if (e.dissipationRecovery) throw sentinel; } }), e => e === sentinel);
});
