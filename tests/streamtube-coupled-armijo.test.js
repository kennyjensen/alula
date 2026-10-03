// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { pathToFileURL } from 'node:url';

// The opt-in draft can be tested before its runtime patch is applied.
const runtimeRoot = process.cwd(), draftRoot = process.env.MSES_ARMIJO_DRAFT_ROOT ?? runtimeRoot;
const runtime = p => pathToFileURL(path.resolve(runtimeRoot, p)).href;
const changed = p => pathToFileURL(path.resolve(draftRoot, p)).href;
const { requireCoupledResidualDecrease } = await import(changed('src/euler/streamtube-iteration-progress.js'));
const { createCoupledStreamtubeBody, coupledStreamtubeTripEvents, coupledStreamtubeResult } = await import(runtime('src/euler/streamtube-coupled.js'));
const { solveCoupledStreamtubeIses } = await import(runtime('src/euler/streamtube-coupled-ises.js'));
const { intrinsicBodyFixture } = await import(runtime('tests/fixtures/intrinsic-body.js'));
const { redistributeStreamtubeTangentially } = await import(runtime('src/geometry/streamtube-tangential-redistribution.js'));
const plain = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));

test('Armijo uses all existing scaled rows, exact allowed decrease and the actual limited step', () => {
  const current = [3, 4], candidate = [3, 3.999], before = structuredClone({ current, candidate });
  const accepted = requireCoupledResidualDecrease(current, candidate, { step: .5, tolerance: 1e-10 });
  assert.equal(accepted.beforeSquaredNorm, 25); assert.equal(accepted.afterSquaredNorm, 3 ** 2 + 3.999 ** 2);
  assert.equal(accepted.allowedSquaredNorm, (1 - 1e-4 * .5) * 25);
  assert.equal(accepted.coefficient, 1e-4); assert.equal(accepted.converged, false);
  assert.deepEqual({ current, candidate }, before);
  let failure;
  try { requireCoupledResidualDecrease(current, current, { step: .5, tolerance: 1e-10 }); } catch (e) { failure = e; }
  assert.equal(failure.code, 'COUPLED_RESIDUAL_DECREASE');
  assert.equal(failure.diagnostics.afterSquaredNorm, 25);
  assert.equal(failure.diagnostics.allowedSquaredNorm, accepted.allowedSquaredNorm);
  assert.match(failure.message, /Maintained coupled squared residual/);
});

test('maximum residual tolerance remains authoritative and malformed/nonfinite norms cannot pass', () => {
  const accepted = requireCoupledResidualDecrease([0, 0], [1e-12, -1e-12], { step: 1, tolerance: 1e-10 });
  assert.equal(accepted.converged, true);
  for (const [a, b, controls] of [[[], [], { step: 1, tolerance: 1e-10 }],
    [[1], [1, 2], { step: 1, tolerance: 1e-10 }], [[NaN], [0], { step: 1, tolerance: 1e-10 }],
    [[1], [Infinity], { step: 1, tolerance: 1e-10 }], [[1], [0], { step: 0, tolerance: 1e-10 }],
    [[1], [0], { step: 1.1, tolerance: 1e-10 }], [[1], [0], { step: .5, tolerance: 0 }]])
    assert.throws(() => requireCoupledResidualDecrease(a, b, controls), /Invalid coupled/);
  assert.throws(() => requireCoupledResidualDecrease([1e300], [1], { step: 1, tolerance: 1e-10 }), /Nonfinite coupled squared/);
});

let cached;
function fixture(policy = 'armijo') {
  if (!cached) {
    const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
      wakeGeometry: 'independent-banks', wakeOutlet: 'banks' };
    const result = solveCoupledStreamtubeIses(input, { edgeMatching: 'section-velocity', transitionMode: 'automatic',
      maxIterations: 0, stepAcceptance: 'admissible', linearOrdering: 'auto', iterationRecovery: false });
    assert(result.checkpoint, result.reason); cached = plain(result.checkpoint);
  }
  const checkpoint = structuredClone(cached), r = checkpoint.restart;
  const system = createCoupledStreamtubeBody(r.input, { ...r.options, initialEuler: r.initialEuler, initialBL: r.initialBL });
  assert(system.n < 500);
  checkpoint.continuation.stepAcceptance = policy;
  checkpoint.families = { euler: 1, boundaryLayer: 0, edgeMatching: 0 };
  const i = system.euler.layout.bodies[0].leadingIndex + 1;
  const lower = r.initialEuler.nodes[0][i].at(-1), upper = r.initialEuler.nodes[1][i][0];
  const spacing = .5 * Math.hypot(lower.x - upper.x, lower.y - upper.y);
  checkpoint.continuation.lastRedistributedStagnation[0] -= (.5 - 1e-8) * spacing;
  const direction = new Float64Array(system.n), column = system.euler.layout.globals.stagnation[0];
  direction[column] = .001 * spacing / system.euler.curves[0].length;
  return { checkpoint, direction, column };
}

// Keep real gas/BL/convex geometry checks and real SMOVE. Substitute only a
// prescribed direction and residual magnitudes to isolate acceptance order.
// These artificial merit values are not claimed to solve the flow equations.
async function driver(direction, magnitudes = [1, 2, 1.2, .9]) {
  const key = `armijo-test-${Math.random()}`, order = [], corrections = [], systems = [], solves = [], merits = [];
  globalThis[key] = {
    createCoupledStreamtubeBody(input, options) {
      const system = createCoupledStreamtubeBody(input, options), index = systems.length;
      const actualEvaluate = system.evaluate, actualAdmissibleValue = system.admissibleValue; systems.push(system);
      system.jacobian = () => { order.push('jacobian'); return { controlled: true }; };
      const controlled = value => {
        const residual = new Float64Array(system.n);
        residual[0] = magnitudes[Math.min(index, magnitudes.length - 1)];
        return { ...value, residual, families: { euler: residual[0], boundaryLayer: 0, edgeMatching: 0 } };
      };
      system.evaluate = x => controlled(actualEvaluate(x));
      system.admissibleValue = (x, controls) => {
        assert.equal(controls.requireConvex, true, 'Armijo retains the mandatory convex candidate domain.');
        const admitted = actualAdmissibleValue(x, controls); order.push(`admissible:${index}:${!!admitted}`);
        return admitted ? controlled(admitted) : false;
      };
      return system;
    },
    coupledStreamtubeTripEvents, coupledStreamtubeResult,
    solveSparseDirect(matrix) {
      assert.equal(matrix.controlled, true); solves.push(1);
      return { x: direction.slice(), relativeResidual: 0, refinements: 0, ordering: 'amd', pivotTolerance: .001, attempts: [] };
    },
    redistributeStreamtubeTangentially(nodes, options) {
      corrections.push(options.correctionScale); order.push('SMOVE'); return redistributeStreamtubeTangentially(nodes, options);
    },
    requireCoupledResidualDecrease(a, b, controls) {
      order.push('merit'); merits.push({ a: a[0], b: b[0], ...controls }); return requireCoupledResidualDecrease(a, b, controls);
    },
  };
  let text = fs.readFileSync(path.resolve(draftRoot, 'src/euler/streamtube-coupled-ises.js'), 'utf8');
  for (const statement of [
    "import { createCoupledStreamtubeBody, coupledStreamtubeTripEvents, coupledStreamtubeResult } from './streamtube-coupled.js';",
    "import { solveSparseDirect } from '../numerics/klu.js';",
    "import { redistributeStreamtubeTangentially } from '../geometry/streamtube-tangential-redistribution.js';",
    "import { requireCoupledResidualDecrease } from './streamtube-iteration-progress.js';",
  ]) {
    assert(text.includes(statement)); const names = statement.slice(statement.indexOf('{'), statement.indexOf('}') + 1);
    text = text.replace(statement, `const ${names} = globalThis[${JSON.stringify(key)}];`);
  }
  text = text.replace(/from '(\.[^']+)'/g, (_, relative) =>
    `from '${new URL(relative, runtime('src/euler/streamtube-coupled-ises.js')).href}'`);
  try {
    const imported = await import('data:text/javascript;base64,' + Buffer.from(text).toString('base64'));
    return { solve: imported.solveCoupledStreamtubeIses, order, corrections, systems, solves, merits };
  } finally { delete globalThis[key]; }
}

test('Armijo retries the final maintained state, retains full SMOVE and reuses one global direction', async () => {
  const { checkpoint, direction, column } = fixture(), before = structuredClone(checkpoint);
  const d = await driver(direction), published = [];
  const result = d.solve(undefined, { resume: checkpoint, stepAcceptance: 'armijo', maxIterations: 1, maxBacktracks: 2,
    onCheckpoint: value => published.push(plain(value)) });
  assert.equal(result.history.length, 2, result.reason); assert.equal(result.history[1].step, .25);
  assert.deepEqual(result.history[1].rejections.map(r => r.stage), ['residual decrease', 'residual decrease']);
  assert.deepEqual(d.merits.map(m => m.step), [1, .5, .25]);
  assert.deepEqual(d.corrections, Array(9).fill(1), 'A merit retry does not halve any SMOVE correction.');
  assert.equal(d.solves.length, 1); assert.equal(result.linearDiagnostics.solves, 1);
  assert.equal(result.x[column], checkpoint.restart.initialEuler.x[column] + .25 * direction[column]);
  for (let index = 1; index <= 3; index++) {
    const admission = d.order.indexOf(`admissible:${index}:true`);
    const merit = d.order.indexOf('merit', admission); assert(admission >= 0 && merit > admission);
  }
  assert.equal(result.history[1].residualDecrease.afterSquaredNorm, .9 ** 2);
  assert.equal(result.gridAcceptance, 'convex'); assert.equal(result.mesh.quality.valid, true);
  assert.equal(result.checkpoint.continuation.stepAcceptance, 'armijo');
  assert.equal(published.length, 2); assert.deepEqual(checkpoint, before);
  const resumedDriver = await driver(direction, [.9]);
  const resumed = resumedDriver.solve(undefined, { resume: result.checkpoint, stepAcceptance: 'armijo', maxIterations: 0 });
  assert.deepEqual(plain(resumed.checkpoint), plain(result.checkpoint));
  assert.equal(resumedDriver.solves.length, 0); assert.equal(resumedDriver.corrections.length, 0);
  assert.throws(() => resumedDriver.solve(undefined, { resume: result.checkpoint,
    stepAcceptance: 'admissible', maxIterations: 0 }), /update controls do not match/);
});

test('exhausted merit retries retain the complete accepted checkpoint and explicit failure evidence', async () => {
  const { checkpoint, direction } = fixture(), before = structuredClone(checkpoint);
  const d = await driver(direction, [1, 2]), published = [];
  const result = d.solve(undefined, { resume: checkpoint, stepAcceptance: 'armijo', maxIterations: 1, maxBacktracks: 1,
    onCheckpoint: cp => published.push(plain(cp)) });
  assert.equal(result.history.length, 1); assert.equal(result.lastRejectedStep.stage, 'residual decrease');
  assert.equal(result.lastRejectedStep.code, 'COUPLED_RESIDUAL_DECREASE');
  assert.equal(result.lastRejectedStep.diagnostics.afterSquaredNorm, 4);
  assert.equal(result.lastRejectedStep.diagnostics.allowedSquaredNorm, 1 - 1e-4 * .5);
  assert.deepEqual(plain(result.checkpoint), before); assert.deepEqual(checkpoint, before);
  assert.equal(published.length, 1); assert.equal(d.solves.length, 1);
});

test('legacy admissible keeps its first valid step and its checkpoint schema without a merit call', async () => {
  const { checkpoint, direction } = fixture('admissible'), d = await driver(direction, [1, 2]);
  const result = d.solve(undefined, { resume: checkpoint, stepAcceptance: 'admissible', maxIterations: 1 });
  assert.equal(result.history[1].step, 1); assert.equal(d.merits.length, 0);
  assert.equal(Object.hasOwn(result.history[1], 'residualDecrease'), false);
  assert.equal(result.checkpoint.continuation.stepAcceptance, 'admissible');
  assert.throws(() => d.solve(undefined, { resume: checkpoint, stepAcceptance: 'armijo', maxIterations: 0 }), /update controls do not match/);
});

test('plateau recovery backtracks the maintained candidate with one factorization per update', async () => {
  const { checkpoint, direction } = fixture('admissible');
  const d = await driver(direction, [1, 1, 1, 1, 1, 1, 1, 2, .8]);
  const result = d.solve(undefined, { resume: checkpoint, stepAcceptance: 'admissible',
    iterationRecovery: true, maxIterations: 7, maxBacktracks: 2 });
  assert.equal(result.history.length, 8, result.reason);
  assert.equal(result.history[6].progress.cause, 'residual-plateau');
  assert.equal(result.history[6].progress.action, 'damped-recovery');
  assert.equal(result.history[7].step, .25);
  assert.deepEqual(d.merits.map(m => m.step), [.5, .25]);
  assert.equal(d.solves.length, 7, 'line-search retry must not refactor');
  assert.equal(result.mesh.quality.valid, true);
});

test('motionless high-residual state stops after three updates without a false convergence', async () => {
  const { checkpoint, direction } = fixture('admissible'); direction.fill(0);
  // Keep the ordinary coordinate chart fixed for this no-motion fixture.
  checkpoint.continuation.lastRedistributedStagnation = createCoupledStreamtubeBody(checkpoint.restart.input,
    { ...checkpoint.restart.options, initialEuler: checkpoint.restart.initialEuler, initialBL: checkpoint.restart.initialBL })
    .evaluate(Float64Array.from([...checkpoint.restart.initialEuler.x, ...checkpoint.restart.initialBL])).outer.stagnation;
  const d = await driver(direction, [1]);
  const result = d.solve(undefined, { resume: checkpoint, stepAcceptance: 'admissible', iterationRecovery: true, maxIterations: 20 });
  assert.equal(result.converged, false);
  assert(result.linearDiagnostics.solves < 20, result.reason);
  assert.match(result.reason, /stalled|Recovery/);
});

test('rapid reduction grants only the explicitly bounded extra iterations', async () => {
  const { checkpoint, direction } = fixture('admissible');
  const d = await driver(direction, [1, .1, .01, .001, .0001, .00001, .000001]);
  const result = d.solve(undefined, { resume: checkpoint, stepAcceptance: 'admissible', iterationRecovery: true,
    maxIterations: 3, maxProgressExtraIterations: 2 });
  assert.equal(result.linearDiagnostics.solves, 5, result.reason);
  assert.equal(result.history[4].progress.extendedBudget, true);
  assert.equal(result.history[5].progress.extendedBudget, true);
  const limited = await driver(direction, [1, .1, .01, .001, .0001]);
  assert.equal(limited.solve(undefined, { resume: checkpoint, stepAcceptance: 'admissible', iterationRecovery: true,
    maxIterations: 3, maxProgressExtraIterations: 0 }).linearDiagnostics.solves, 3);
});

test('automatic exhausted merit search restores the ordinary policy without refactoring or admitting rejected geometry', async () => {
  const { checkpoint, direction } = fixture('admissible');
  const before = structuredClone(checkpoint);
  const d = await driver(direction, [1, 1, 1, 1, 1, 1, 1, 1.01]);
  const r = d.solve(undefined, { resume: checkpoint, stepAcceptance: 'admissible',
    iterationRecovery: true, maxIterations: 9, maxBacktracks: 2 });
  assert.equal(r.history.length, 10, r.reason);
  assert.equal(r.history[7].step, .5);
  assert.equal(r.history[7].progress.recoveryOutcome.kind, 'ordinary-policy-restored');
  assert.equal(r.history[7].progress.recoveryOutcome.rejectedTrials, 3);
  assert.equal(r.history[8].step, 1);
  assert.equal(r.lastRejectedStep, null);
  assert.equal(r.converged, false);
  assert.equal(r.reason, 'iteration limit');
  assert.equal(d.solves.length, 9);
  assert.equal(r.mesh.quality.valid, true);
  assert.deepEqual(checkpoint, before);
});

test('slow but moving damped recovery returns its remaining budget to ordinary Newton', async () => {
  const { checkpoint, direction } = fixture('admissible');
  const d = await driver(direction, [1, 1, 1, 1, 1, 1, 1, .999, .998, .997, .996]);
  const r = d.solve(undefined, { resume: checkpoint, stepAcceptance: 'admissible',
    iterationRecovery: true, maxIterations: 10 });
  assert.equal(r.history.length, 11, r.reason);
  assert.equal(r.history[9].progress.recoveryOutcome.reason, 'Damping did not accelerate residual reduction');
  assert.equal(r.history[10].step, 1);
  assert.equal(r.reason, 'iteration limit');
});
