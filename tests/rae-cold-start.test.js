// Cold app-route regressions. The remaining grid combinations are unresolved.
import test from 'node:test';
import assert from 'node:assert/strict';
import { buildReliabilityCase } from '../scripts/validation/solver-reliability-cases.js';
import { solveCoupledStreamtubeAssembly } from '../src/euler/streamtube-coupled-assembly.js';
import { solveStreamtubeAssembly } from '../src/euler/streamtube-result.js';

test('RAE 2822 8x7 standalone inviscid app route converges without driving the outlet sonic', () => {
  const { caseData: input } = buildReliabilityCase({
    preset: 'rae2822-mses', mode: 'streamtube-grid',
    changes: { mach: .74, alpha: 2.68, gridIntervals: 8, gridTubes: 7 },
  });
  const before = structuredClone(input);
  const result = solveStreamtubeAssembly(input, { maxIterations: 40 });
  assert.equal(result.status, 'research-converged', result.diagnostics.reason);
  assert.equal(result.solverSettings.stepAcceptance, 'armijo');
  assert.equal(result.mach, .74);
  assert.equal(result.alpha, 2.68);
  assert.equal(result.mesh.cells.length, 1260);
  assert.equal(result.mesh.quality.valid, true);
  assert.equal(result.diagnostics.equationResidual <= 1e-10, true);
  assert.deepEqual(result.solverSettings.upwind, { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } });
  assert.ok(result.diagnostics.iterations <= 25);
  for (const h of result.flow.history.slice(1)) {
    assert.ok(h.residualDecrease.afterSquaredNorm <= h.residualDecrease.allowedSquaredNorm);
  }
  assert.deepEqual(input, before);
});

test('RAE 2822 16x11 standalone Euler removes the premature double shock with adaptive MCRIT', () => {
  const { caseData: input } = buildReliabilityCase({ preset: 'rae2822-mses', mode: 'streamtube-grid',
    changes: { mach: .74, alpha: 2.68, gridIntervals: 16, gridTubes: 11 } });
  const before = structuredClone(input), progress = [];
  const result = solveStreamtubeAssembly(input, { maxIterations: 40, onIteration: h => progress.push(h) });
  assert.equal(result.status, 'research-converged', result.diagnostics.reason);
  assert.ok(result.diagnostics.equationResidual <= 1e-10);
  assert.equal(result.mesh.quality.valid, true);
  assert.equal(result.mesh.cells.length, 2772);
  assert.equal(result.mach, .74);
  assert.equal(result.alpha, 2.68);
  assert.equal(result.solverSettings.adaptiveMcrit, true);
  assert.deepEqual(result.flow.solverInput.upwind, { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } });
  assert.ok(progress.some(h => h.dissipation?.mcrit < .9));
  assert.equal(progress.at(-1).dissipation.mcrit, .99);
  assert.ok(result.diagnostics.iterations <= 20, 'The cold solve should recover full Newton convergence.');
  for (const h of result.flow.history.slice(1))
    assert.ok(h.residualDecrease.afterSquaredNorm <= h.residualDecrease.allowedSquaredNorm,
      'Residual decrease is checked with the dissipation law frozen.');
  // The failed iterate has a subsonic pocket at x/c≈.29 followed by another
  // supersonic acceleration. Check the resolved plateau, not a fitted Cp.
  const upper = result.mesh.flow.lines.find(line => line.group === 1 && line.tube === 0);
  for (let i = 0; i < upper.machNumbers.length; i++) {
    const x = .5 * (upper.points[i].x + upper.points[i + 1].x);
    if (x > .1 && x < .5) assert.ok(upper.machNumbers[i] > 1);
  }
  assert.deepEqual(input, before);
});

for (const preset of ['single', 'flap'])
test(`standalone ISMOM default preserves the subsonic ${preset} Euler solve`, () => {
  const { caseData: input } = buildReliabilityCase({ preset, mode: 'streamtube-grid',
    changes: { mach: .2, alpha: 4, gridIntervals: 16, gridTubes: 7 } });
  const before = structuredClone(input);
  const result = solveStreamtubeAssembly(input, { maxIterations: 40 });
  assert.equal(result.status, 'research-converged', result.diagnostics.reason);
  assert.equal(result.solverSettings.stepAcceptance, 'armijo');
  assert.equal(result.mach, .2);
  assert.equal(result.alpha, 4);
  assert.equal(result.mesh.quality.valid, true);
  assert.ok(result.diagnostics.equationResidual <= 1e-10);
  assert.deepEqual(input, before);
});

for (const gridTubes of [7, 9])
test(`RAE 2822 32x${gridTubes} standalone app worker converges after restoring second-order dissipation`, async () => {
  const { caseData: input } = buildReliabilityCase({ preset: 'rae2822-mses', mode: 'streamtube-grid',
    changes: { mach: .74, alpha: 2.68, gridIntervals: 32, gridTubes } });
  const before = structuredClone(input), progress = [], messages = [], previousSelf = globalThis.self;
  // Run the real message handler and its numerical imports. In particular,
  // the app must reach the startup policy without test-only solver options.
  globalThis.self = { postMessage(message) {
    if (message.type === 'iteration') progress.push(structuredClone(message.iteration));
    if (message.type === 'result' || message.type === 'error') messages.push(structuredClone(message));
  } };
  try {
    await import(`../src/worker/solver.js?cold-rae-tubes=${gridTubes}`);
    await self.onmessage({ data: { id: 15, task: 'solve', caseData: input } });
  } finally { globalThis.self = previousSelf; }
  assert.equal(messages.length, 1);
  assert.equal(messages[0].type, 'result', messages[0].message);
  assert.equal(messages[0].id, 15);
  const { result } = messages[0];
  assert.equal(result.status, 'research-converged', result.diagnostics.reason);
  assert.ok(result.diagnostics.equationResidual <= 1e-10);
  assert.equal(result.mesh.quality.valid, true);
  assert.equal(result.mesh.cells.length, 191 * 2 * (gridTubes + 3));
  assert.equal(result.mach, .74);
  assert.equal(result.alpha, 2.68);
  assert.equal(result.solverSettings.firstOrderStartup, true);
  assert.deepEqual(result.flow.solverInput.upwind, { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } });
  assert.ok(result.diagnostics.iterations <= 35, 'Both orders must converge within one 40-update app budget.');
  assert.ok(result.flow.linearDiagnostics.maxRelativeResidual <= 1e-10, 'Do not loosen the sparse solve accuracy gate.');
  const switchEvent = progress.find(h => h.dissipation?.restoredSecondOrder);
  assert.ok(switchEvent, 'The requested second-order operator must be restored.');
  assert.ok(switchEvent.dissipation.firstOrderResidual <= 1e-10);
  assert.ok(switchEvent.residual > 1e-10, 'The first-order root is not a root of the requested equations.');
  assert.ok(progress.at(-1).iteration > switchEvent.iteration);
  assert.equal(progress[0].dissipation.mucon, -1);
  assert.equal(progress.at(-1).dissipation.mucon, 1);
  assert.equal(progress.at(-1).dissipation.mcrit, .99);
  for (const h of result.flow.history.slice(1))
    assert.ok(h.residualDecrease.afterSquaredNorm <= h.residualDecrease.allowedSquaredNorm);
  assert.deepEqual(input, before);
});

for (const [gridIntervals, gridTubes] of [[8, 7], [8, 9], [8, 11], [16, 7], [16, 9], [16, 11], [32, 7], [32, 9], [32, 11]])
test(`RAE 2822 ${gridIntervals}x${gridTubes} reaches the requested viscous condition through the direct cold app route`, t => {
  const { caseData: input } = buildReliabilityCase({
    preset: 'rae2822-mses', mode: 'streamtube-bl',
    changes: { mach: .74, alpha: 2.68, reynolds: 2.7e6, ncrit: 4, gridIntervals, gridTubes },
  });
  const before = structuredClone(input);
  const result = solveCoupledStreamtubeAssembly(input, {
    direct: true, maxIterations: 40, tolerance: 1e-10,
  });
  t.diagnostic(JSON.stringify({ gridIntervals, gridTubes, converged: result.converged,
    reason: result.reason, iterations: result.history.length - 1, families: result.families,
    finalStep: result.history.at(-1)?.step, finalBacktracks: result.history.at(-1)?.backtracks }));
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.status, 'research-coupled-equations-converged');
  assert.equal(result.mach, .74);
  assert.equal(result.alpha, 2.68);
  assert.equal(result.referenceReynolds, 2.7e6);
  assert.equal(result.coupledOptions.ncrit, 4);
  assert.equal(result.coupledOptions.transitionMode, 'automatic');
  assert.ok(Object.values(result.families).every(value => value <= 1e-10));
  assert.equal(result.mesh.quality.valid, true);
  assert.equal(result.solverSettings.direct, true);
  assert.equal(result.solverSettings.maxStartupAttempts, 1);
  assert.equal(result.solverSettings.upwind.mcrit, .99);
  assert.equal(result.solverSettings.upwind.mucon, 1);
  assert.equal(result.solverSettings.eulerMaxIterations, 40, 'Match the direct app precursor budget.');
  assert.ok(result.initialization.euler.residual <= 1e-10,
    'The shock startup must restore and converge the requested Euler equations before this handoff.');
  assert.ok(result.initialization.euler.iterations <= 35, 'Both Euler orders share one budget.');
  assert.ok(result.history.length <= 35, 'The viscous solve should regain full Newton convergence.');
  assert.ok(result.history.slice(-3).every(h => h.step === 1 && h.backtracks === 0));
  assert.ok(result.flow.sections.flat(2).every(({ rho, p }) =>
    Number.isFinite(rho) && rho > 0 && Number.isFinite(p) && p > 0));
  assert.deepEqual(input, before);
});

for (const gridTubes of [7, 9, 11])
test(`RAE 2822 64x${gridTubes} standard startup converges through the app worker`, async t => {
  const { caseData: input } = buildReliabilityCase({ preset: 'rae2822-mses', mode: 'streamtube-grid',
    changes: { mach: .74, alpha: 2.68, gridIntervals: 64, gridTubes } });
  input.eulerStartup = 'standard';
  const before = structuredClone(input), progress = [], messages = [], previousSelf = globalThis.self;
  globalThis.self = { postMessage(message) {
    if (message.type === 'iteration') {
      progress.push(structuredClone(message.iteration));
      if (process.env.RAE_COLD_PROGRESS) {
        const h = message.iteration;
        console.log(JSON.stringify({ gridTubes, iteration: h.iteration, residual: h.residual,
          step: h.step, dissipation: h.dissipation }));
      }
    }
    if (message.type === 'result' || message.type === 'error') messages.push(message);
  } };
  try {
    await import(`../src/worker/solver.js?harmonic-cold-rae-tubes=${gridTubes}`);
    await self.onmessage({ data: { id: 64, task: 'solve', caseData: input } });
  } finally { globalThis.self = previousSelf; }
  assert.equal(messages.length, 1);
  assert.equal(messages[0].type, 'result', messages[0].message);
  const { result } = messages[0];
  assert.equal(result.status, 'research-converged', result.diagnostics.reason);
  assert.ok(result.diagnostics.equationResidual <= 1e-10);
  assert.equal(result.mach, .74);
  assert.equal(result.alpha, 2.68);
  assert.equal(result.mesh.quality.valid, true);
  assert.equal(result.mesh.cells.length, 380 * 2 * (gridTubes + 3));
  assert.equal(result.solverSettings.eulerStartup, 'harmonic-shock');
  assert.equal(result.solverSettings.maxIterations, 80);
  assert.deepEqual(result.flow.solverInput.upwind, { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } });
  assert.ok(result.flow.linearDiagnostics.maxRelativeResidual <= 1e-10);
  assert.equal(progress[0].dissipation.mucon, -2);
  assert.equal(progress[0].dissipation.mcrit, .75);
  const restored = progress.find(h => h.dissipation?.restoredSecondOrder);
  assert.ok(restored);
  assert.ok(restored.dissipation.firstOrderResidual <= 1e-10);
  assert.ok(restored.residual > 1e-10, 'The temporary root is not the requested-law root.');
  assert.ok(progress.at(-1).iteration > restored.iteration);
  assert.equal(progress.at(-1).dissipation.mucon, 1);
  assert.equal(progress.at(-1).dissipation.mcrit, .99);
  assert.deepEqual(input, before);
  t.diagnostic(JSON.stringify({ gridTubes, iterations: result.diagnostics.iterations,
    residual: result.diagnostics.equationResidual, seconds: messages[0].elapsed / 1000 }));
});
