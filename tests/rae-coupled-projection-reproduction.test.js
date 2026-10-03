// Exact transonic transition and grid-update replays, independent of cold
// grid construction. Keep the failing iterations and subsequent root checks
// reproducible without selecting an experimental app initialization policy.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { proposeCoupledXfoilBLUpdate } from '../src/euler/streamtube-coupled-xfoil-update.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { relaxCoupledWakeSeed } from '../src/euler/streamtube-coupled-initializer.js';

const fixture = name => JSON.parse(fs.readFileSync(new URL(
  `../docs/solver-reliability/rae-inviscid-audit/fixtures/${name}.json`, import.meta.url)));
const controls = { maxIterations: 1, tolerance: 1e-10 };

test('RAE32x7 first-wake contact preserves a usable raw Newton bound', () => {
  const cp = fixture('rae32x7-before-wake-backtracking-regression'), before = structuredClone(cp);
  const result = solveCoupledStreamtubeIses(undefined, { resume: cp, ...cp.continuation, ...controls });
  const h = result.history[1];
  assert.ok(h, result.reason);
  const rejection = h.rejections.find(r => r.code === 'streamtube-grid-step');
  assert.deepEqual(rejection.diagnostics.cell, { group: 1, i: 127, tube: 0, corner: 3 });
  assert.ok(rejection.diagnostics.stepFraction > .3 && rejection.diagnostics.stepFraction < .32);
  assert.equal(rejection.diagnostics.coordinateRepairBacktrack, undefined);
  assert.ok(Math.abs(h.step - rejection.step * rejection.diagnostics.stepFraction) < 1e-14,
    'A rejected pairing must not replace a usable convexity bound with a larger half step.');
  assert.ok(h.step > .16 && h.step < .17);
  assert.equal(h.wakeCoordinateRepair, undefined);
  assert.ok(h.residualDecrease.afterSquaredNorm < h.residualDecrease.beforeSquaredNorm);
  assert.equal(result.mesh.quality.valid, true);
  assert.deepEqual(result.checkpoint.restart.input, cp.restart.input);
  assert.deepEqual(result.checkpoint.continuation.fractions, cp.continuation.fractions);
  assert.equal(result.coupledOptions.ncrit, 4);
  assert.equal(result.coupledOptions.transitionMode, 'automatic');
  assert.deepEqual(cp, before);
});

test('RAE64x9 wake contact recovers through corrected trial steps and reaches the requested root', t => {
  const cp = fixture('rae64x9-before-wake-trial-backtracking'), before = structuredClone(cp);
  const policies = [];
  const result = solveCoupledStreamtubeIses(undefined, {
    resume: cp, ...cp.continuation, maxIterations: 30, tolerance: 1e-10,
    onCheckpoint: saved => policies.push(saved.continuation.wakeCoordinateRecovery === true),
  });
  t.diagnostic(JSON.stringify({ converged: result.converged, reason: result.reason,
    updates: result.history.length - 1, families: result.families,
    firstStep: result.history[1]?.step, lastStep: result.history.at(-1)?.step }));
  assert.equal(result.converged, true, result.reason);
  assert.ok(Object.values(result.families).every(r => r <= 1e-10));
  assert.equal(result.mesh.quality.valid, true);
  assert.ok(result.linearDiagnostics.maxRelativeResidual <= 1e-10);
  assert.ok(result.history[1].step > .01,
    'Reusing the uncorrected-grid bound skipped useful corrected trials and accepted only a 1e-6 step.');
  assert.ok(result.history.some(h => h.wakeCoordinateRepair?.method === 'paired-wake-trial'));
  assert.equal(policies[0], false);
  assert.ok(policies.slice(1, -1).some(Boolean), 'The active coordinate policy must survive a checkpoint.');
  assert.equal(policies.at(-1), false, 'A requested-law root releases the temporary coordinate recovery.');
  assert.equal(result.checkpoint.continuation.wakeCoordinateRecovery, undefined);
  assert.ok(result.history.slice(-3).every(h => h.step === 1 && h.backtracks === 0));
  assert.deepEqual(result.checkpoint.restart.input, cp.restart.input);
  assert.deepEqual(result.checkpoint.continuation.fractions, cp.continuation.fractions);
  const f = result.checkpoint.restart;
  assert.equal(f.options.reynolds, cp.restart.options.reynolds);
  assert.equal(f.options.ncrit, 4);
  assert.equal(f.options.transitionMode, 'automatic');
  const exact = createCoupledStreamtubeBody(f.input, {
    ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL,
  });
  assert.deepEqual(exact.evaluate(exact.initial).families, result.families);
  assert.deepEqual(cp, before);
  t.diagnostic(JSON.stringify({ families: result.families, updates: result.history.length - 1,
    firstStep: result.history[1].step }));
});

test('RAE32x7 cold wake correction owns its grid and preserves BL, gas and prescribed boundaries', () => {
  const cp = fixture('rae32x7-original-grid-cold-viscous-initial'), before = structuredClone(cp), f = cp.restart;
  const system = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  const state = system.initial.slice(), value = system.evaluate(state);
  assert.deepEqual(value.families, cp.families);
  const corrected = relaxCoupledWakeSeed(f.input, f.options, system, state, value);
  assert.equal(corrected?.report.accepted, true, corrected?.report.reason);
  assert.notEqual(corrected.system, system);
  assert.deepEqual(corrected.state, state, 'The correction changes the chart, not packed flow/BL unknowns.');
  assert.deepEqual(corrected.system.bl.snapshotActive(), system.bl.snapshotActive());
  assert.equal(corrected.value.families.boundaryLayer, value.families.boundaryLayer);
  assert.ok(corrected.report.afterSquaredNorm < corrected.report.beforeSquaredNorm);
  assert.ok(corrected.system.admissibleValue(corrected.state, { requireConvex: true }));
  let moved = false;
  const old = value.outer.nodes, next = corrected.value.outer.nodes, le = system.euler.layout.bodies[0].leadingIndex;
  for (let g = 0; g < old.length; g++) for (let i = 0; i < old[g].length; i++)
    for (let j = 0; j < old[g][i].length; j++) {
      if (i < le || i === old[g].length - 1 || j === 0 || j === old[g][i].length - 1)
        assert.deepEqual(next[g][i][j], old[g][i][j], 'Inlet, forward body and boundary nodes must be preserved.');
      else moved ||= next[g][i][j].x !== old[g][i][j].x || next[g][i][j].y !== old[g][i][j].y;
    }
  assert.ok(moved);
  assert.deepEqual(system.evaluate(state).residual, value.residual, 'Candidate evaluation cannot mutate its source.');
  assert.deepEqual(cp, before);
});

test('a BL-dominated RAE32x9 state does not trigger Euler wake smoothing', () => {
  const cp = fixture('rae32x9-before-cold-transition-step7'), f = cp.restart;
  const system = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  const state = system.initial.slice(), value = system.evaluate(state);
  assert.ok(value.families.boundaryLayer > value.families.euler);
  assert.equal(relaxCoupledWakeSeed(f.input, f.options, system, state, value), null);
  assert.deepEqual(system.evaluate(state).residual, value.residual);
});

test('checkpoint resumes cannot request another cold wake-grid initialization', () => {
  const cp = fixture('rae32x7-original-grid-cold-viscous-initial');
  assert.throws(() => solveCoupledStreamtubeIses(undefined, {
    resume: cp, ...cp.continuation, wakeGridInitialization: true, maxIterations: 0,
  }), /requires a fresh coupled solve/);
});

test('RAE32x9 cold upstream transition transfers shear without losing Newton convergence', t => {
  const cp = fixture('rae32x9-before-cold-transition-step7'), before = structuredClone(cp);
  const result = solveCoupledStreamtubeIses(undefined, {
    resume: cp, ...cp.continuation, maxIterations: 14, tolerance: 1e-10,
  });
  const crossing = result.history[1], change = crossing.changes.find(c => c.side === 'upper');
  assert.equal(crossing.step, 1);
  assert.equal(change.from, 41);
  assert.equal(change.to, 40);
  const transfer = change.converted.find(c => c.index === 40);
  assert.equal(transfer.shearSelection.method, 'mixed-and-following-shear-minimum');
  assert.ok(transfer.shearSelection.minimum.squaredResidual < .5 * transfer.shearSelection.initial.squaredResidual);
  assert.equal(transfer.shearInitialization, undefined, 'Do not misreport the selected shear as a single-row root.');
  assert.ok(crossing.boundaryLayer < .6, 'The old single-row transfer produced BL residual .695.');
  assert.equal(result.converged, true, result.reason);
  assert.ok(Object.values(result.families).every(v => v <= 1e-10));
  assert.ok(result.history.slice(-3).every(h => h.step === 1 && h.backtracks === 0));
  assert.equal(result.mesh.quality.valid, true);
  const f = result.checkpoint.restart;
  assert.equal(f.input.mach, .74);
  assert.equal(f.input.alpha, 2.68);
  assert.deepEqual(f.input.upwind, cp.restart.input.upwind);
  assert.equal(f.options.reynolds, cp.restart.options.reynolds);
  assert.equal(f.options.ncrit, 4);
  assert.equal(f.options.transitionMode, 'automatic');
  const exact = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  assert.ok(Math.abs(f.options.reynolds / exact.euler.conditions.lengthScale - 2.7e6) < 1e-8);
  assert.deepEqual(exact.evaluate(exact.initial).families, result.families);
  assert.deepEqual(cp, before);
  t.diagnostic(JSON.stringify({ families: result.families, updates: result.history.length - 1 }));
});

for (const [label, name, nx, tubes, root = false, probeSteps,
  rootSource = 'rae64x11-after-prescribed-euler-before-step3'] of [
  ['RAE32x24 late trust-region stall', 'rae32x24-late-trust-region-stall', 127, 27],
  ['RAE64x9 fresh cold linear stall', 'rae64x9-fresh-cold-before-linear-step23', 380, 12],
  ['RAE64x11 transition-feedback recovered root', 'rae64x11-transition-feedback-recovered-root', 380, 14, true, undefined,
    'rae64x11-after-wake-trial-linear-stall'],
  ['RAE64x9 late wake-corrected root', 'rae64x9-late-wake-corrected-root', 380, 12, true, undefined,
    'rae64x9-before-late-wake-merit-rejection'],
  ['RAE64x9 wake-trial backtracking root', 'rae64x9-wake-trial-backtracking-root', 380, 12, true, undefined,
    'rae64x9-before-wake-trial-backtracking'],
  ['RAE64x7 fresh localized cold repeat root', 'rae64x7-localized-repeat-cold-root', 380, 10, true, undefined,
    'rae64x7-localized-repeat-cold-initial'],
  ['RAE64x7 localized wake-corrected root', 'rae64x7-localized-wake-corrected-root', 380, 10, true, undefined,
    'rae64x7-localized-before-wake-contact'],
  ['RAE64x11 cold coordinate linear stall', 'rae64x11-cold-coordinate-linear-stall', 380, 14,
    false, [1e-4, 1e-5, 1e-6]],
  ['RAE64x11 repaired-wake cold-path root', 'rae64x11-repaired-wake-cold-path-root', 380, 14, true, undefined,
    'rae64x11-before-first-wake-collapse'],
  ['RAE64x11 first-wake collapse', 'rae64x11-before-first-wake-collapse', 380, 14, false, [1e-5, 3e-6, 1e-6]],
  // Larger coordinate perturbations cross this saved cell's narrow domain.
  ['RAE64x11 postshock grid contact', 'rae64x11-before-postshock-grid-contact', 380, 14, false, [1e-5, 3e-6, 1e-6]],
  ['RAE64x11 warm root on exact cold grid', 'rae64x11-warm-exact-cold-grid-root', 380, 14, true],
  ['RAE64x11 released BL before large update', 'rae64x11-after-prescribed-euler-before-step3', 380, 14],
  ['RAE64x11 tapered BL wake stall', 'rae64x11-tapered-bl-before-wake-step', 380, 14],
  ['RAE64x11 harmonic cold wake stall', 'rae64x11-harmonic-cold-coupled-wake-stall', 380, 14],
  ['RAE64x11 pre-pinch', 'rae64x11-before-coupled-wake-pinch', 380, 14],
  ['RAE64x11 pre-inlet-pinch', 'rae64x11-before-coupled-inlet-pinch', 380, 14],
  ['RAE64x11 pre-surface-projection', 'rae64x11-before-coupled-surface-projection', 380, 14],
  ['RAE32x7 pre-surface-projection', 'rae32x7-before-surface-projection-rejection', 191, 10],
  ['RAE32x7 cold wake-grid stall', 'rae32x7-before-cold-viscous-grid-stall', 191, 10],
]) test(`${label} coupled state replays and has a consistent full and block Jacobian`, t => {
  const cp = fixture(name), before = structuredClone(cp), f = cp.restart;
  const system = createCoupledStreamtubeBody(f.input, { ...f.options,
    initialEuler: f.initialEuler, initialBL: f.initialBL });
  const x = system.initial, value = system.admissibleValue(x, { requireConvex: true });
  assert.ok(value, 'The retained state before the grid-pinch trial must be admissible.');
  assert.deepEqual(value.families, cp.families);
  if (root) {
    assert.ok(Object.values(value.families).every(r => r <= 1e-10));
    assert.deepEqual(value.outer.nodes, f.initialEuler.nodes, 'The root physical grid must replay exactly.');
    const cold = fixture(rootSource);
    assert.deepEqual(f.input, cold.restart.input,
      'The diagnostic root preserves its source input, including surface fractions and mass weights.');
    assert.deepEqual(cp.continuation.fractions, cold.continuation.fractions,
      'The prescribed inlet spacing must also match the source solve.');
  }
  assert.equal(system.euler.layout.nx, nx);
  assert.deepEqual(system.euler.layout.tubes, [tubes, tubes]);
  assert.equal(f.input.mach, .74);
  assert.equal(f.input.alpha, 2.68);
  assert.equal(f.input.upwind.mucon, 1);
  assert.equal(f.input.upwind.mcrit, .99);
  assert.equal(f.options.ncrit, 4);
  assert.equal(f.options.transitionMode, 'automatic');
  assert.ok(Math.abs(f.options.reynolds / system.euler.conditions.lengthScale - 2.7e6) < 1e-8);
  const matrix = system.jacobian(x), peak = values => values.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
  const sweeps = ['Euler', 'BL', 'both'].map(block => {
    const direction = x.map((v, i) => Math.sin(1.7 * i + .4) * (i < system.ne
      ? block === 'BL' ? 0 : i < system.euler.layout.densityCount ? .02 : .0001
      : block === 'Euler' ? 0 : .001 * Math.max(1, Math.abs(v))));
    const product = sparseProduct(matrix, direction), scale = peak(product);
    const increments = probeSteps ?? (root ? [1e-3, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7] : [1e-3, 1e-4, 1e-5, 1e-6]);
    const errors = increments.map(h => {
      const plus = system.residual(x.map((v, i) => v + h * direction[i]));
      const minus = system.residual(x.map((v, i) => v - h * direction[i]));
      if (root) {
        // A resolved shock's narrow blending region gives large O(h²)
        // curvature error before small-step roundoff takes over. Fourth-order
        // centered differences separate that effect from a Jacobian defect;
        // the accuracy threshold and tested direction remain unchanged.
        const halfPlus = system.residual(x.map((v, i) => v + .5 * h * direction[i]));
        const halfMinus = system.residual(x.map((v, i) => v - .5 * h * direction[i]));
        return peak(product.map((v, i) => v - (4 * (halfPlus[i] - halfMinus[i]) / h
          - (plus[i] - minus[i]) / (2 * h)) / 3)) / scale;
      }
      return peak(product.map((v, i) => v - (plus[i] - minus[i]) / (2 * h))) / scale;
    });
    assert.ok(Math.min(...errors) < 1e-6, `${label} ${block} Jv: ${errors}`);
    return { block, errors };
  });
  assert.deepEqual(cp, before);
  t.diagnostic(JSON.stringify({ families: value.families, sweeps }));
});

for (const [label, name] of [
  ['RAE64x11 inlet-limited', 'rae64x11-before-coupled-inlet-pinch'],
  ['RAE32x7 cold wake-limited', 'rae32x7-before-cold-viscous-grid-stall'],
]) test(`${label} update preserves the physical domain and decreases coupled merit`, t => {
  const cp = fixture(name), before = structuredClone(cp);
  const result = solveCoupledStreamtubeIses(undefined, { resume: cp, ...cp.continuation, ...controls });
  const step = result.history[1];
  assert.ok(step?.step > 0, result.reason);
  assert.ok(step.residualDecrease.afterSquaredNorm < step.residualDecrease.beforeSquaredNorm);
  assert.equal(result.mesh.quality.valid, true);
  assert.ok(Object.values(result.families).every(Number.isFinite));
  assert.equal(result.solverInput.mach, .74);
  assert.equal(result.solverInput.alpha, 2.68);
  assert.deepEqual(result.solverInput.upwind, cp.restart.input.upwind);
  assert.equal(result.coupledOptions.reynolds, cp.restart.options.reynolds);
  assert.equal(result.coupledOptions.ncrit, 4);
  assert.equal(result.coupledOptions.transitionMode, 'automatic');
  assert.deepEqual(cp, before);
  // A larger safe step is welcome. Do not require the old tiny step or its
  // rejection sequence: this fixture should allow a genuine grid fix.
  t.diagnostic(JSON.stringify({ step: step.step, backtracks: step.backtracks,
    families: result.families, grid: result.mesh.quality,
    rejectedCells: step.rejections.map(r => r.diagnostics?.cell).filter(Boolean) }));
});

test('RAE64x11 Hk recovery preserves coordinates at displacement roundoff', () => {
  const source = fixture('rae64x11-before-hk-recovery-roundoff'), before = structuredClone(source), f = source.restart;
  const system = createCoupledStreamtubeBody(f.input, { ...f.options,
    initialEuler: f.initialEuler, initialBL: f.initialBL });
  const value = system.evaluate(system.initial);
  assert.deepEqual(value.families, source.families);
  // On this state, the former constructor rebased one node by 3.5e-18
  // and failed the strict recovery guard. Keep that guard exact.
  const result = solveCoupledStreamtubeIses(undefined, {
    resume: source, ...source.continuation, maxIterations: 1, tolerance: 1e-10,
  });
  const event = result.history.at(-1);
  assert.ok(event.hkProjectionRecovery, result.reason);
  assert.equal(event.step, 0);
  assert.equal(event.accepted, false);
  assert.equal(event.hkProjectionRecovery.equationsChanged, false);
  assert.ok(event.hkProjectionRecovery.rawMerit < event.hkProjectionRecovery.beforeMerit);
  assert.ok(event.hkProjectionRecovery.projectedMerit >= event.hkProjectionRecovery.beforeMerit);
  assert.equal(result.converged, false, 'A method change alone does not converge the flow.');
  const checkpoint = JSON.parse(JSON.stringify(result.checkpoint,
    (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
  assert.equal(checkpoint.restart.options.geometryReplay, 'preserve-undisplaced');
  assert.equal(checkpoint.restart.options.hkFloorLinearization, 'native');
  const next = checkpoint.restart;
  const replay = createCoupledStreamtubeBody(next.input, { ...next.options,
    initialEuler: next.initialEuler, initialBL: next.initialBL });
  const replayed = replay.evaluate(replay.initial);
  assert.deepEqual(replay.initial, system.initial);
  assert.deepEqual(replayed.outer.nodes, value.outer.nodes);
  assert.deepEqual(replayed.residual, value.residual);
  const resumed = solveCoupledStreamtubeIses(undefined, {
    resume: checkpoint, ...checkpoint.continuation, maxIterations: 0, tolerance: 1e-10,
  });
  assert.deepEqual(resumed.families, source.families);
  assert.deepEqual(source, before);
});

test('RAE32x7 targeted Hk recovery preserves the rejected state, reports its method and converges', t => {
  const source = fixture('rae32x7-before-surface-projection-rejection'), before = structuredClone(source);
  const ordinary = solveCoupledStreamtubeIses(undefined, {
    resume: source, ...source.continuation, maxIterations: 1, tolerance: 1e-10,
  });
  assert.equal(ordinary.converged, false);
  assert.equal(ordinary.lastRejectedStep.code, 'COUPLED_RESIDUAL_DECREASE');
  assert.equal(ordinary.checkpoint.continuation.hkProjectionRecovery, undefined,
    'An archived checkpoint must retain its original policy.');
  let methodCheckpoint;
  const result = solveCoupledStreamtubeIses(undefined, {
    resume: source, ...source.continuation, hkProjectionRecovery: true, maxIterations: 18, tolerance: 1e-10,
    onCheckpoint: (cp, details) => {
      if (details.history.at(-1).hkProjectionRecovery) methodCheckpoint = structuredClone(cp);
    },
  });
  const events = result.history.filter(h => h.hkProjectionRecovery);
  assert.equal(events.length, 1);
  const event = events[0];
  assert.equal(event.accepted, false);
  assert.equal(event.step, 0);
  assert.equal(event.hkFloorLinearization, 'native');
  assert.equal(event.hkProjectionRecovery.equationsChanged, false);
  assert.ok(event.hkProjectionRecovery.rawMerit < event.hkProjectionRecovery.beforeMerit);
  assert.ok(event.hkProjectionRecovery.projectedMerit >= event.hkProjectionRecovery.beforeMerit);
  assert.equal(methodCheckpoint.continuation.hkProjectionRecovery, true);
  assert.equal(methodCheckpoint.restart.options.hkFloorLinearization, 'native');
  assert.deepEqual(methodCheckpoint.families, source.families);
  assert.deepEqual(methodCheckpoint.restart.initialBL, Float64Array.from(source.restart.initialBL));
  assert.deepEqual(methodCheckpoint.restart.initialEuler.x, Float64Array.from(source.restart.initialEuler.x));
  assert.deepEqual(methodCheckpoint.restart.initialEuler.nodes, source.restart.initialEuler.nodes);
  // JSON is the browser/checkpoint persistence boundary. Replay must retain
  // the method without another event or any physical residual change.
  const serial = JSON.parse(JSON.stringify(methodCheckpoint, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
  const replay = solveCoupledStreamtubeIses(undefined, {
    resume: serial, ...serial.continuation, maxIterations: 0, tolerance: 1e-10,
  });
  assert.deepEqual(replay.families, source.families);
  assert.equal(replay.history[0].hkFloorLinearization, 'native');
  assert.equal(replay.checkpoint.continuation.hkProjectionRecovery, true);
  assert.equal(result.converged, true, result.reason);
  assert.ok(result.history.length <= 14);
  assert.ok(result.history.slice(-3).every(h => h.step === 1 && h.backtracks === 0));
  assert.equal(result.mesh.quality.valid, true);
  const f = result.checkpoint.restart;
  const exact = createCoupledStreamtubeBody(f.input, { ...f.options, hkFloorLinearization: 'exact',
    initialEuler: f.initialEuler, initialBL: f.initialBL });
  assert.deepEqual(exact.evaluate(exact.initial).families, result.families);
  assert.ok(Object.values(result.families).every(v => v <= 1e-10));
  assert.equal(f.input.mach, .74);
  assert.equal(f.input.alpha, 2.68);
  assert.deepEqual(f.input.upwind, source.restart.input.upwind);
  assert.equal(f.options.reynolds, source.restart.options.reynolds);
  assert.equal(f.options.ncrit, 4);
  assert.equal(f.options.transitionMode, 'automatic');
  assert.deepEqual(source, before);
  t.diagnostic(JSON.stringify({ events: events.map(h => h.hkProjectionRecovery), families: result.families,
    updates: result.history.length - 2 }));
});

test('RAE32x7 native Hk-floor linearization escapes the projected-step stall without changing residual equations', t => {
  const source = fixture('rae32x7-before-surface-projection-rejection'), before = structuredClone(source);
  const cp = structuredClone(source);
  cp.restart.options.hkFloorLinearization = 'native';
  const nativeBefore = structuredClone(cp);
  const result = solveCoupledStreamtubeIses(undefined, {
    resume: cp, ...cp.continuation, maxIterations: 18, tolerance: 1e-10,
  });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.mesh.quality.valid, true);
  assert.ok(Object.values(result.families).every(v => v <= 1e-10));
  assert.equal(result.solverInput.mach, .74);
  assert.equal(result.solverInput.alpha, 2.68);
  assert.deepEqual(result.solverInput.upwind, source.restart.input.upwind);
  assert.equal(result.coupledOptions.reynolds, source.restart.options.reynolds);
  assert.equal(result.coupledOptions.ncrit, 4);
  assert.equal(result.coupledOptions.transitionMode, 'automatic');
  assert.deepEqual(result.coupledOptions.tripFractions, [[1, 1]]);
  assert.ok(result.history.slice(-3).every(h => h.step === 1 && h.backtracks === 0));
  // Native floor sensitivities are a quasi-Newton direction, not a new
  // residual. Independently certify the root with the exact policy restored.
  const f = result.checkpoint.restart;
  const exact = createCoupledStreamtubeBody(f.input, { ...f.options,
    hkFloorLinearization: 'exact', initialEuler: f.initialEuler, initialBL: f.initialBL });
  const value = exact.admissibleValue(exact.initial, { requireConvex: true });
  assert.ok(value);
  assert.deepEqual(value.families, result.families);
  assert.ok(Math.abs(f.options.reynolds / exact.euler.conditions.lengthScale - 2.7e6) < 1e-8);
  const x = exact.initial, matrix = exact.jacobian(x);
  const peak = values => values.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
  const direction = x.map((v, i) => Math.sin(1.7 * i + .4) * (i < exact.ne
    ? i < exact.euler.layout.densityCount ? .02 : .0001
    : .001 * Math.max(1, Math.abs(v))));
  const product = sparseProduct(matrix, direction), scale = peak(product);
  const errors = [1e-3, 1e-4, 1e-5, 1e-6].map(h => {
    const plus = exact.residual(x.map((v, i) => v + h * direction[i]));
    const minus = exact.residual(x.map((v, i) => v - h * direction[i]));
    return peak(product.map((v, i) => v - (plus[i] - minus[i]) / (2 * h))) / scale;
  });
  assert.ok(Math.min(...errors) < 1e-6, `Exact root Jacobian: ${errors}`);
  assert.deepEqual(cp, nativeBefore);
  assert.deepEqual(source, before);
  t.diagnostic(JSON.stringify({ iterations: result.history.length - 1,
    families: value.families, exactRootJacobianErrors: errors }));
});

test('RAE8x11 terminal transition trial preserves the wake tail and recovers Newton convergence', t => {
  const cp = fixture('rae8x11-before-terminal-replay'), before = structuredClone(cp);
  const result = solveCoupledStreamtubeIses(undefined, {
    resume: cp, ...cp.continuation, maxIterations: 25, tolerance: 1e-10,
  });
  assert.equal(result.converged, true, result.reason);
  const first = result.history[1];
  assert.equal(first.activeChange ?? false, false);
  assert.ok(first.rejections.some(r => r.code === 'COUPLED_RESIDUAL_DECREASE'));
  assert.ok(first.step > .09 && first.step < .1,
    'Reject the old full terminal-event trial and accept the half step.');
  assert.ok(result.history.length <= 21);
  assert.ok(result.history.slice(-4).every(h => h.step === 1 && h.backtracks === 0));
  assert.ok(Object.values(result.families).every(v => v < 1e-10));
  assert.equal(result.mesh.quality.valid, true);
  assert.equal(result.solverInput.mach, .74);
  assert.equal(result.solverInput.alpha, 2.68);
  assert.deepEqual(result.solverInput.upwind, cp.restart.input.upwind);
  assert.equal(result.coupledOptions.reynolds, cp.restart.options.reynolds);
  assert.equal(result.coupledOptions.ncrit, 4);
  assert.equal(result.coupledOptions.transitionMode, 'automatic');
  assert.deepEqual(result.coupledOptions.tripFractions, [[1, 1]]);
  assert.deepEqual(cp, before);
  const root = result.checkpoint.restart;
  const system = createCoupledStreamtubeBody(root.input, { ...root.options,
    initialEuler: root.initialEuler, initialBL: root.initialBL });
  assert.deepEqual(system.euler.layout.tubes, [14, 14]);
  assert.equal(system.euler.layout.nx, 63);
  assert.ok(Math.abs(root.options.reynolds / system.euler.conditions.lengthScale - 2.7e6) < 1e-8);
  const x = system.initial, matrix = system.jacobian(x);
  const peak = values => values.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
  const sweeps = ['Euler', 'BL', 'both'].map(block => {
    const direction = x.map((v, i) => Math.sin(1.7 * i + .4) * (i < system.ne
      ? block === 'BL' ? 0 : i < system.euler.layout.densityCount ? .02 : .0001
      : block === 'Euler' ? 0 : .001 * Math.max(1, Math.abs(v))));
    const product = sparseProduct(matrix, direction), scale = peak(product);
    const errors = [1e-4, 1e-5, 1e-6].map(h => {
      const plus = system.residual(x.map((v, i) => v + h * direction[i]));
      const minus = system.residual(x.map((v, i) => v - h * direction[i]));
      return peak(product.map((v, i) => v - (plus[i] - minus[i]) / (2 * h))) / scale;
    });
    assert.ok(Math.min(...errors) < 1e-6, `${block} root Jv: ${errors}`);
    return { block, errors };
  });
  t.diagnostic(JSON.stringify({ firstStep: first.step, iterations: result.history.length - 1,
    families: result.families, sweeps }));
});

test('RAE32x11 crosses transition with the smaller admitted event and recovers full Newton convergence', t => {
  const cp = fixture('rae32x11-before-transition-event'), before = structuredClone(cp);
  const result = solveCoupledStreamtubeIses(undefined, {
    resume: cp, ...cp.continuation, maxIterations: 12, tolerance: 1e-10,
  });
  assert.equal(result.converged, true, result.reason);
  const event = result.history[1];
  assert.equal(event.transitionEventAcceptance.method, 'admissible-transition-event');
  assert.ok(event.step < event.rejections.find(r => r.stage === 'residual decrease').step / 4);
  assert.ok(event.residual < .05, 'The first large event trial raised the residual above one.');
  assert.deepEqual(event.changes.map(c => [c.side, c.from, c.to]), [['upper', 40, 39]]);
  assert.ok(event.changes[0].fraction > .99, 'Cross near the interval boundary.');
  // The six-row transition transfer can need one damped step after crossing.
  // Require decreasing maintained merit and terminal full Newton convergence,
  // not the earlier single-row initializer's exact damping sequence.
  for (const h of result.history.slice(2))
    assert.ok(h.residualDecrease.afterSquaredNorm <= h.residualDecrease.allowedSquaredNorm);
  assert.ok(result.history.slice(-3).every(h => h.step === 1 && h.backtracks === 0));
  assert.ok(result.history.length <= 7);
  assert.ok(Object.values(result.families).every(v => v < 1e-10));
  assert.equal(result.mesh.quality.valid, true);
  assert.equal(result.solverInput.mach, .74);
  assert.equal(result.solverInput.alpha, 2.68);
  assert.deepEqual(result.solverInput.upwind, cp.restart.input.upwind);
  assert.equal(result.coupledOptions.reynolds, cp.restart.options.reynolds);
  assert.equal(result.coupledOptions.ncrit, 4);
  assert.equal(result.coupledOptions.transitionMode, 'automatic');
  assert.deepEqual(result.coupledOptions.tripFractions, [[1, 1]]);
  assert.deepEqual(cp, before);
  t.diagnostic(JSON.stringify({ eventStep: event.step, eventResidual: event.residual,
    iterations: result.history.length - 1, families: result.families }));
});

test('RAE16x11 logarithmic shear reaches the target root from the saved ordinary BL seed', t => {
  const source = fixture('rae16x11-coupled-initial'), before = structuredClone(source), f = source.restart;
  const result = solveCoupledStreamtubeIses(f.input, { ...f.options,
    initialEuler: f.initialEuler, initialBL: f.initialBL,
    maxIterations: 40, tolerance: 1e-10, iterationGeometry: 'ises-sampled',
    stepAcceptance: 'event-armijo', blUpdate: 'xfoil', projectionGeometry: 'fixed',
    shearCoordinate: 'logarithmic' });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.solverInput.mach, .74);
  assert.equal(result.solverInput.alpha, 2.68);
  assert.equal(result.coupledOptions.ncrit, 4);
  assert.equal(result.coupledOptions.transitionMode, 'automatic');
  assert.deepEqual(result.coupledOptions.tripFractions, [[1, 1]]);
  assert.equal(result.coupledOptions.reynolds, f.options.reynolds);
  assert.deepEqual(result.solverInput.upwind, f.input.upwind);
  assert.equal(result.solverInput.upwind.mucon, 1);
  assert.equal(result.solverInput.upwind.mcrit, .99);
  assert.ok(Object.values(result.families).every(v => v < 1e-10));
  assert.equal(result.mesh.quality.valid, true);
  assert.ok(result.mesh.quality.minCornerSine > .1);
  assert.ok(result.history.slice(-3).every(h => h.step === 1 && h.backtracks === 0));
  assert.equal(result.checkpoint.continuation.shearCoordinate, 'logarithmic');
  assert.ok(result.history.slice(1).every(h => h.stepKind === 'coupled-density-newton-log-shear'
    && h.logarithmicShear.every(c => c.from > 0 && c.to > 0 && Number.isFinite(c.to))));
  assert.deepEqual(source, before);
  const cp = result.checkpoint.restart;
  const system = createCoupledStreamtubeBody(cp.input, { ...cp.options,
    initialEuler: cp.initialEuler, initialBL: cp.initialBL });
  assert.equal(system.euler.layout.nx, 99);
  assert.deepEqual(system.euler.layout.tubes, [14, 14]);
  assert.ok(Math.abs(cp.options.reynolds / system.euler.conditions.lengthScale - 2.7e6) < 1e-8);
  const x = system.initial, matrix = system.jacobian(x);
  const peak = values => values.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
  const sweeps = ['Euler', 'BL', 'both'].map(block => {
    const direction = x.map((v, i) => Math.sin(1.7 * i + .4) * (i < system.ne
      ? block === 'BL' ? 0 : i < system.euler.layout.densityCount ? .02 : .0001
      : block === 'Euler' ? 0 : .001 * Math.max(1, Math.abs(v))));
    const product = sparseProduct(matrix, direction), scale = peak(product);
    const errors = [1e-4, 1e-5, 1e-6].map(h => {
      const plus = system.residual(x.map((v, i) => v + h * direction[i]));
      const minus = system.residual(x.map((v, i) => v - h * direction[i]));
      return peak(product.map((v, i) => v - (plus[i] - minus[i]) / (2 * h))) / scale;
    });
    assert.ok(Math.min(...errors) < 1e-6, `${block} root Jv: ${errors}`);
    return { block, errors };
  });
  t.diagnostic(JSON.stringify({ iterations: result.history.length - 1, families: result.families, sweeps }));
});

test('RAE32x11 damps trial grid repair and converges beyond the first wake pinch', t => {
  const cp = fixture('rae32x11-coupled-grid-repair'), before = structuredClone(cp);
  const result = solveCoupledStreamtubeIses(undefined, { resume: cp, ...cp.continuation, ...controls, maxIterations: 25 });
  assert.equal(result.converged, true, result.reason);
  const step = result.history[1];
  assert.equal(step.rejections[0].stage, 'Newton grid step');
  assert.equal(step.rejections[1].code, 'COUPLED_RESIDUAL_DECREASE');
  assert.deepEqual(step.maintenance.triggeredBodies, []);
  assert.equal(step.maintenance.geometryRedistribution, true);
  assert.deepEqual(step.maintenance.passages.map(p => p.correctionScale), [.5, .5]);
  assert.equal(step.step, step.rejections[1].step, 'Damp the coordinate correction before shortening Newton.');
  assert.ok(step.step > 9e-3);
  assert.ok(step.residualDecrease.afterSquaredNorm < step.residualDecrease.beforeSquaredNorm);
  assert.ok(result.families.euler < cp.families.euler);
  assert.ok(result.families.boundaryLayer < cp.families.boundaryLayer);
  assert.equal(result.mesh.quality.valid, true);
  assert.equal(result.solverInput.mach, .74);
  assert.equal(result.solverInput.alpha, 2.68);
  assert.deepEqual(result.solverInput.upwind, cp.restart.input.upwind);
  assert.equal(result.coupledOptions.transitionMode, 'automatic');
  assert.equal(result.coupledOptions.ncrit, 4);
  assert.equal(result.coupledOptions.reynolds, cp.restart.options.reynolds);
  assert.ok(Object.values(result.families).every(v => v < 1e-10));
  assert.ok(result.history.slice(-3).every(h => h.step === 1 && h.backtracks === 0));
  assert.deepEqual(cp, before);
  t.diagnostic(JSON.stringify({ acceptedStep: step.step, iterations: result.history.length - 1, families: result.families }));
});

test('RAE32 stalled projection is continuous at zero step', () => {
  const cp = fixture('rae32-coupled-before-step7'), f = cp.restart;
  const system = createCoupledStreamtubeBody(f.input, { ...f.options,
    initialEuler: f.initialEuler, initialBL: f.initialBL });
  const state = system.initial, before = state.slice();
  assert.deepEqual(system.evaluate(state).families, cp.families);
  const proposal = proposeCoupledXfoilBLUpdate(system, state, new Float64Array(system.n));
  assert.ok(proposal.x.every((v, i) => Math.abs(v - state[i]) < 1e-13));
  assert.ok(proposal.projection.displacementChanges.every(c => Math.abs(c.correction) < 1e-13));
  assert.deepEqual(state, before);
});

test('RAE32 extra projection extension rejects a residual-decreasing viscous step', t => {
  const cp = fixture('rae32-coupled-before-step7'), before = structuredClone(cp);
  const extended = solveCoupledStreamtubeIses(undefined, { resume: cp, ...cp.continuation, ...controls });
  assert.equal(extended.history.length, 1);
  assert.equal(extended.lastRejectedStep.code, 'streamtube-interface-pressure');
  assert.deepEqual(extended.lastRejectedStep.diagnostics.cell, { i: 108, group: 1, tube: 10 });
  const fixed = structuredClone(cp);
  fixed.continuation.projectionGeometry = 'fixed';
  const ordinary = solveCoupledStreamtubeIses(undefined, { resume: fixed, ...fixed.continuation, ...controls });
  assert.equal(ordinary.history.length, 2, ordinary.reason);
  assert.equal(ordinary.mesh.quality.valid, true);
  assert.ok(ordinary.history[1].residualDecrease.afterSquaredNorm
    < ordinary.history[1].residualDecrease.beforeSquaredNorm);
  assert.deepEqual(cp, before);
  assert.equal(ordinary.solverInput.mach, .74);
  assert.equal(ordinary.solverInput.alpha, 2.68);
  t.diagnostic(JSON.stringify({ acceptedStep: ordinary.history[1].step, families: ordinary.families }));
});

test('RAE32 viscous startup from restored second-order Euler converges without extra grid extension', t => {
  const source = fixture('rae32-coupled-initial'), f = source.restart, before = structuredClone(source);
  const result = solveCoupledStreamtubeIses(f.input, { ...f.options, initialEuler: f.initialEuler,
    initialBL: f.initialBL, maxIterations: 40, tolerance: 1e-10, iterationGeometry: 'ises-sampled',
    stepAcceptance: 'event-armijo', blUpdate: 'xfoil', projectionGeometry: 'fixed', shearCoordinate: 'linear' });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.mesh.quality.valid, true);
  assert.ok(Object.values(result.families).every(v => v <= 1e-10));
  assert.equal(result.solverInput.mach, .74);
  assert.equal(result.solverInput.alpha, 2.68);
  assert.deepEqual(result.solverInput.upwind, f.input.upwind);
  assert.equal(result.solverInput.upwind.mucon, 1);
  assert.equal(result.solverInput.upwind.mcrit, .99);
  assert.equal(result.coupledOptions.ncrit, 4);
  assert.equal(result.coupledOptions.reynolds, f.options.reynolds);
  assert.equal(result.coupledOptions.transitionMode, 'automatic');
  assert.ok(result.flow.sections.flat(2).every(({ rho, p }) =>
    Number.isFinite(rho) && rho > 0 && Number.isFinite(p) && p > 0));
  assert.deepEqual(source, before);
  const cp = result.checkpoint.restart;
  const system = createCoupledStreamtubeBody(cp.input, { ...cp.options,
    initialEuler: cp.initialEuler, initialBL: cp.initialBL });
  assert.deepEqual(system.euler.layout.tubes, [12, 12]);
  assert.equal(system.euler.layout.nx, 191);
  assert.ok(Math.abs(result.coupledOptions.reynolds / system.euler.conditions.lengthScale - 2.7e6) < 1e-8,
    'Reference-chord Reynolds number must remain 2.7e6');
  const x = system.initial, jacobian = system.jacobian(x);
  const maximum = a => a.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
  const direction = x.map((v, i) => Math.sin(1.7 * i + .4) * (i < system.ne
    ? i < system.euler.layout.densityCount ? .02 : .0001 : .001 * Math.max(1, Math.abs(v))));
  const product = sparseProduct(jacobian, direction), scale = maximum(product);
  assert.ok(scale > 0);
  const errors = [1e-4, 1e-5, 1e-6].map(h => {
    const plus = system.residual(x.map((v, i) => v + h * direction[i]));
    const minus = system.residual(x.map((v, i) => v - h * direction[i]));
    return maximum(product.map((v, i) => v - (plus[i] - minus[i]) / (2 * h))) / scale;
  });
  assert.ok(Math.min(...errors) < 1e-6, `Converged transonic coupled Jv: ${errors}`);
  t.diagnostic(JSON.stringify({ iterations: result.history.length - 1, families: result.families,
    relativeDirectionalErrors: errors }));
});

test('RAE32x24 retained surface-floor state admits a decreasing constrained trust-region step', t => {
  // Exact user case: automatic spacing, 32 explicit inlet/outlet intervals.
  // This tests escape from the failed iteration, not convergence of the case.
  const cp = fixture('rae32x24-before-viscous-merit-rejection'), before = structuredClone(cp), f = cp.restart;
  const system = createCoupledStreamtubeBody(f.input, {
    ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL,
  });
  assert.deepEqual(system.evaluate(system.initial).families, cp.families);
  const phase = system.bl.snapshotActive();
  const result = solveCoupledStreamtubeBody(system, {
    stepMethod: 'dogleg', initialTrustRadius: 1, maxIterations: 1, tolerance: 1e-10,
  });
  const step = result.history[1];
  assert.ok(step, result.reason);
  assert.ok(step.actualReduction > 0);
  assert.ok(step.predictedReduction > 0);
  assert.ok(step.residual < result.history[0].residual);
  assert.equal(result.mesh.quality.valid, true);
  assert.equal(result.converged, false, 'One decreasing step is not a converged solution.');
  assert.deepEqual(system.bl.snapshotActive(), phase);
  assert.deepEqual(cp, before);
  t.diagnostic(JSON.stringify({ kind: step.stepKind, residual: step.residual,
    actualReduction: step.actualReduction, families: result.families }));
});
