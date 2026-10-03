// SPDX-License-Identifier: GPL-2.0-or-later
import { solveStreamtubeIses } from './streamtube-ises-update.js';
import { solveCoupledStreamtubeIses } from './streamtube-coupled-ises.js';
import { createCoupledStreamtubeBody } from './streamtube-coupled.js';
import { initializeStreamtubeWakeCorrespondence } from './streamtube-wake-correspondence.js';
import { requireCoupledResidualDecrease } from './streamtube-iteration-progress.js';

const same = (a, b) => {
  if (Object.is(a, b)) return true;
  if (!a || !b || typeof a !== 'object' || typeof b !== 'object') return false;
  const keys = Object.keys(a);
  return keys.length === Object.keys(b).length && keys.every(k => Object.hasOwn(b, k) && same(a[k], b[k]));
};

// An inviscid precursor is only a starting guess for the coupled equations.
// When its accepted terminal iterate is worse in both residual and geometry,
// use an earlier complete state. Do not reinterpret this as Euler convergence.
// The ordinary inviscid API and converged precursors retain their endpoints.
export function selectCoupledEulerStartup(precursor, { tolerance = 1e-10 } = {}) {
  const terminal = precursor?.flow, best = terminal?.bestCheckpoint;
  if (precursor?.status !== 'unconverged' || terminal?.converged !== false
    || terminal?.reason !== 'iteration limit' || terminal.lastRejectedStep != null
    || terminal.initialRedistribution?.accepted !== true || terminal.finalQuality?.valid !== true
    || !Number.isFinite(terminal.finalQuality.minCornerSine) || !best?.checkpoint
    || !Number.isInteger(best.iteration) || best.iteration < 0
    || !(best.iteration < terminal.history?.at(-1)?.iteration)
    || !Number.isFinite(best.residual) || !(best.residual < terminal.diagnostics?.residual)) return null;
  const { iterationGeometry, stepAcceptance, stagnationLimiter } = best.checkpoint.continuation;
  const adaptive = terminal.adaptiveMcrit === true, stagedOrder = terminal.firstOrderStartup === true;
  const comparableInput = adaptive || stagedOrder ? { ...terminal.solverInput,
    upwind: { ...terminal.solverInput.upwind,
      ...(adaptive ? { mcrit: best.checkpoint.input.upwind?.mcrit } : {}),
      ...(stagedOrder ? { mucon: best.checkpoint.input.upwind?.mucon } : {}) } } : terminal.solverInput;
  if (!same(best.checkpoint.input, comparableInput)
    || iterationGeometry !== terminal.iterationGeometry || stepAcceptance !== terminal.stepAcceptance
    || stagnationLimiter !== terminal.stagnationLimiter
    || adaptive !== (best.checkpoint.continuation.adaptiveMcrit ?? false)
    || adaptive && best.checkpoint.continuation.targetMcrit !== terminal.targetMcrit
    || stagedOrder !== (best.checkpoint.continuation.firstOrderStartup ?? false)
    || stagedOrder && best.checkpoint.continuation.targetMucon !== terminal.targetMucon)
    throw new Error('Retained Euler startup has different geometry, equations or update controls.');
  // Residuals from different temporary dissipation laws are not comparable.
  // Keep the terminal state instead of ranking or rejecting it on that basis.
  if (adaptive && best.checkpoint.input.upwind.mcrit !== terminal.solverInput.upwind.mcrit
    || stagedOrder && best.checkpoint.input.upwind.mucon !== terminal.solverInput.upwind.mucon) return null;
  // Exact complete-state replay also rechecks gas, primal and dual geometry.
  // maxIterations=0 performs no Jacobian, linear solve or grid maintenance.
  const flow = solveStreamtubeIses(undefined, { resume: best.checkpoint, maxIterations: 0,
    tolerance, ...best.checkpoint.continuation, retainCheckpoint: true });
  if (flow.diagnostics.residual !== best.residual)
    throw new Error('Retained Euler startup residual does not match its complete checkpoint.');
  if (flow.finalQuality.valid !== true || !(flow.finalQuality.minCornerSine >= terminal.finalQuality.minCornerSine)) return null;
  return { flow, solverInput: flow.solverInput, selection: {
    kind: 'earlier-admissible-inviscid-state', selectedIteration: best.iteration,
    selectedResidual: flow.diagnostics.residual, selectedMinCornerSine: flow.finalQuality.minCornerSine,
    terminalIteration: terminal.history.at(-1).iteration, terminalResidual: terminal.diagnostics.residual,
    terminalMinCornerSine: terminal.finalQuality.minCornerSine,
    stateConverged: flow.converged, reason: 'Earlier state improves residual without reducing the convex-corner margin.',
  } };
}

// Work selection only. No evaluation, Jacobian or solve.
// A decreasing tail is evidence to spend a bounded additional budget; it
// never substitutes for the unchanged complete coupled convergence criterion.
// Each extension is capped at 20 updates and the original stage size. The
// finite reserve also accounts for the resolved surface length. A moving
// transition can cross many intervals before full Newton steps become useful;
// refining the mesh must not leave it with a coarse mesh's work allowance.
const families = ['euler', 'boundaryLayer', 'edgeMatching'];
const vector = x => (Array.isArray(x) || ArrayBuffer.isView(x)) && x.length > 0 && x.every(Number.isFinite);
const positiveFamilies = f => f && families.every(k => Number.isFinite(f[k]) && f[k] > 0);
const sameVector = (a, b) => vector(a) && vector(b) && a.length === b.length && a.every((v, i) => Object.is(v, b[i]));
const sameFamilies = (a, b) => positiveFamilies(a) && positiveFamilies(b) && families.every(k => Object.is(a[k], b[k]));
const points = nodes => Array.isArray(nodes) && nodes.length > 0 && nodes.every(group => Array.isArray(group)
  && group.length > 1 && group.every(row => Array.isArray(row) && row.length > 1
    && row.every(p => p && Number.isFinite(p.x) && Number.isFinite(p.y))));
const sameNodes = (a, b) => points(a) && points(b) && a.length === b.length && a.every((g, i) => g.length === b[i].length
  && g.every((row, j) => row.length === b[i][j].length && row.every((p, k) => Object.is(p.x, b[i][j][k].x) && Object.is(p.y, b[i][j][k].y))));
const empty = value => value === undefined || Array.isArray(value) && value.length === 0;
const inactive = value => value === undefined || value === false;
const finiteTree = x => typeof x === 'number' ? Number.isFinite(x) : Array.isArray(x) && x.length > 0 && x.every(finiteTree);
// Trip location does not determine whether Newton progress earns more work.
// Preserve any valid prescribed trips, checking exact agreement with the
// retained checkpoint below, just as for a natural-transition-only run.
const trips = x => Array.isArray(x) && x.length > 0 && x.every(pair => Array.isArray(pair) && pair.length === 2
  && pair.every(v => Number.isFinite(v) && v > 0 && v <= 1));

export function coupledStartupExtensionPlan(result, { startupAttempt, maxIterations, remainingIterations = 20,
  extensionCount = 0, requestedConditions, coarseInitialization, tolerance = 1e-10 } = {}) {
  const cp = result?.checkpoint, f = cp?.restart, c = cp?.continuation, expected = requestedConditions, history = result?.history;
  const previous = result?.startupExtension, completed = history?.length - 1;
  const extended = completed - maxIterations;
  const surfaceIntervals = Array.isArray(f?.input?.bodies) ? f.input.bodies.reduce((longest, body) => {
    const { leadingIndex, trailingIndex } = body ?? {};
    return Number.isSafeInteger(leadingIndex) && leadingIndex >= 0 && Number.isSafeInteger(trailingIndex)
      && trailingIndex > leadingIndex ? Math.max(longest, trailingIndex - leadingIndex) : longest;
  }, 0) : 0;
  // Allow a phase move and a correction per interval on the longest surface.
  // Surfaces advance together, so neither element count nor farfield/wake
  // resolution increases this reserve. Every chunk still has to earn progress.
  const maximumExtensionIterations = 2 * Math.max(maxIterations, surfaceIntervals);
  if (startupAttempt !== 1 || !Number.isInteger(extensionCount) || extensionCount < 0
    || !Number.isInteger(maxIterations) || maxIterations < 3
    || (extensionCount === 0 ? previous !== undefined || completed !== maxIterations
      : previous?.extensionCount !== extensionCount || previous.initialIterations !== maxIterations
        || previous.iterations !== extended || extended <= 0)
    || extended >= maximumExtensionIterations
    || !Number.isInteger(remainingIterations) || remainingIterations < 1 || !Number.isFinite(tolerance) || tolerance <= 0
    || coarseInitialization || result?.automaticRefinement || result?.gridSequence || result?.refinement
    || result?.recoveryAttempt || result?.ncritContinuation || result?.machContinuation
    || result?.converged !== false || result?.reason !== 'iteration limit' || result.lastRejectedStep != null
    || result.initialRedistribution?.accepted !== true || result.mesh?.quality?.valid !== true
    || !Array.isArray(result.mesh.quality.invalidCells) || result.mesh.quality.invalidCells.length
    || cp?.version !== 1 || !f?.input || !f.options || !c || !positiveFamilies(result.families)
    || Math.max(...families.map(k => result.families[k])) <= tolerance || !sameFamilies(result.families, cp.families)
    || !expected || !Number.isFinite(expected.mach) || expected.mach <= 0
    || !Number.isFinite(expected.reynolds) || expected.reynolds <= 0 || !Number.isFinite(expected.ncrit) || expected.ncrit <= 0
    || expected.transitionMode !== 'automatic' || !trips(expected.tripFractions)
    || expected.ismom !== undefined && (!Number.isInteger(expected.ismom) || expected.ismom < 1 || expected.ismom > 4)
    || f.input.mach !== expected.mach || f.options.reynolds !== expected.reynolds || f.options.ncrit !== expected.ncrit
    || f.options.transitionMode !== expected.transitionMode || f.input.hybrid?.ismom !== expected.ismom
    || JSON.stringify(f.options.tripFractions) !== JSON.stringify(expected.tripFractions)
    || result.conditions?.mach !== expected.mach || result.conditions?.reynolds !== expected.reynolds
    || result.conditions?.ncrit !== expected.ncrit || result.conditions?.transitionMode !== expected.transitionMode
    || !Array.isArray(f.input.bodies) || f.input.bodies.length !== expected.tripFractions.length
    || !Array.isArray(f.options.transitionState) || f.options.transitionState.length !== 2 * f.input.bodies.length
    || !f.options.transitionState.every(x => Number.isInteger(x) && x >= 0)
    || !vector(f.initialEuler?.x) || !vector(f.initialBL) || !vector(result.x) || !vector(result.residual)
    || result.x.length !== f.initialEuler.x.length + f.initialBL.length || result.residual.length !== result.x.length
    || result.residual.reduce((m, v) => Math.max(m, Math.abs(v)), 0) !== Math.max(...families.map(k => result.families[k]))
    || !sameVector(result.x.slice(0, f.initialEuler.x.length), f.initialEuler.x)
    || !sameVector(result.x.slice(f.initialEuler.x.length), f.initialBL)
    || !sameNodes(result.flow?.nodes, f.initialEuler.nodes) || !sameNodes(result.flow?.undisplacedNodes, f.initialEuler.undisplacedNodes)
    || !finiteTree(c.fractions) || !vector(c.lastRedistributedStagnation) || c.lastRedistributedStagnation.length !== f.input.bodies.length
    || !['convex', 'ises-sampled'].includes(c.iterationGeometry) || !['admissible', 'armijo', 'event-armijo'].includes(c.stepAcceptance)
    || !['listing', 'prose'].includes(c.stagnationLimiter) || c.eventProfile !== undefined && c.eventProfile !== 'none'
    || !Array.isArray(history)
    || history.some((h, i) => h?.iteration !== i)) return null;
  const tail = history.slice(-3), last = tail.at(-1), maintenance = last.maintenance;
  const familyDecrease = families.every(k => tail[0][k] > tail[1][k] && tail[1][k] > tail[2][k]);
  // A shock can redistribute error between rows/families while Newton is
  // reducing the actual merit rapidly. Accept two consecutive full Armijo
  // steps as alternative evidence, with a useful net decrease, no events,
  // and the terminal merit checked against the complete residual vector.
  // Pressure coupling uses weighted rows, so retain its family-based rule.
  const recent = tail.slice(1), merits = recent.map(h => h.residualDecrease);
  const squaredNorm = result.residual.reduce((sum, v) => sum + v * v, 0);
  const close = (a, b) => Number.isFinite(a) && Number.isFinite(b)
    && Math.abs(a - b) <= 1e-12 * Math.max(Math.abs(a), Math.abs(b), Number.MIN_VALUE);
  const meritDecrease = f.options.edgeMatching === 'section-velocity'
    && recent.every((h, i) => h.step === 1 && h.backtracks === 0 && empty(h.rejections)
      && inactive(h.activeChange) && empty(h.changes) && h.eventProfile === undefined
      && merits[i]?.method === 'armijo-squared-residual' && merits[i].step === 1
      && Number.isFinite(merits[i].beforeSquaredNorm) && merits[i].beforeSquaredNorm > 0
      && Number.isFinite(merits[i].afterSquaredNorm) && merits[i].afterSquaredNorm > 0
      && merits[i].afterSquaredNorm <= merits[i].allowedSquaredNorm
      && merits[i].allowedSquaredNorm < merits[i].beforeSquaredNorm)
    && close(merits[0].afterSquaredNorm, merits[1].beforeSquaredNorm)
    && close(merits[1].afterSquaredNorm, squaredNorm)
    && merits[1].afterSquaredNorm <= .8 * merits[0].beforeSquaredNorm;
  const smoothProgress = (familyDecrease || meritDecrease)
    && last.step === 1 && last.backtracks === 0 && empty(last.rejections) && inactive(last.activeChange) && empty(last.changes)
    && last.eventProfile === undefined && maintenance && inactive(maintenance.geometryRedistribution)
    && empty(maintenance.triggeredBodies) && empty(maintenance.passages) && empty(maintenance.dekinkRepairs);
  // Compare only each accepted step's own unchanged equations. Never compare
  // N-amplification rows with shear rows across a transition event. This is
  // evidence to buy more work, not a claim of net decrease across events or
  // convergence. A finite reserve bounds even an event cycle with local descent.
  const window = history.slice(-Math.min(20, maxIterations));
  const descending = window.filter(h => {
    const d = h.residualDecrease;
    return h.accepted !== false && h.step > 0 && h.step <= 1 && inactive(h.activeChange) && empty(h.changes)
      && h.eventProfile === undefined && d?.method === 'armijo-squared-residual'
      && d.step === h.step && Number.isFinite(d.beforeSquaredNorm) && d.beforeSquaredNorm > 0
      && Number.isFinite(d.afterSquaredNorm) && d.afterSquaredNorm > 0
      && d.afterSquaredNorm <= d.allowedSquaredNorm && d.allowedSquaredNorm < d.beforeSquaredNorm;
  });
  const logReduction = descending.reduce((sum, h) => sum
    + Math.log(h.residualDecrease.afterSquaredNorm / h.residualDecrease.beforeSquaredNorm), 0);
  const measuredProgress = ['armijo', 'event-armijo'].includes(c.stepAcceptance)
    && !c.dissipationEnhancement && descending.length >= 2
    && descending.filter(h => h.step >= 1e-3).length >= 2 && logReduction <= Math.log(.8);
  if (tail.some(h => !positiveFamilies(h)) || !sameFamilies(last, result.families)
    || !(smoothProgress || measuredProgress)) return null;
  const additionalIterations = Math.min(20, maxIterations, remainingIterations, maximumExtensionIterations - extended);
  return { kind: 'coupled-startup-earned-extension', selection: smoothProgress
      ? familyDecrease ? 'one bounded work chunk after a full step and two strict all-family decreases'
        : 'one bounded work chunk after two full Armijo steps and a 20% squared-residual decrease'
      : 'one bounded work chunk after measured descent on unchanged equations',
    ...(measuredProgress ? { progress: { window: window.length, descendingSteps: descending.length,
      pairedMeritRatio: Math.exp(logReduction), comparesAcrossEvents: false } } : {}),
    startupAttempt: 1, extensionCount: extensionCount + 1, initialIterations: maxIterations,
    originalIterations: completed, additionalIterations, maximumExtensionIterations,
    maximumTotalIterations: completed + additionalIterations, tolerance, equationsChanged: false,
    actualNcrit: expected.ncrit, targetNcrit: expected.ncrit,
    resume: structuredClone(cp),
    // Kept separately even if a later trial fails: never overwrite this
    // accepted endpoint or discard its complete, unrenumbered first history.
    original: structuredClone({ reason: result.reason, families: result.families, checkpoint: cp,
      history, quality: result.mesh.quality, initialRedistribution: result.initialRedistribution,
      linearDiagnostics: result.linearDiagnostics, lastRejectedStep: result.lastRejectedStep }),
    accounting: { continuationInitialFrameIsAnUpdate: false, continuationIterationOffset: completed,
      originalHistoryLength: history.length } };
}

// Bounded correction for a retained, admissible state whose
// Newton trials crowd the first wake cell. The correspondence guess is never
// accepted just because its cells look better: coordinate damping or bounded
// coupled updates must lower the ORIGINAL squared residual before the caller
// receives a replacement checkpoint.
// This is a coordinate globalization experiment, not an MSES formula.
export function correctCoupledWakeCoordinates(checkpoint, { tolerance = 1e-10,
  maxCorrectorIterations = 6, maxCoordinateBacktracks = 0, onIteration } = {}) {
  if (!(Number.isFinite(tolerance) && tolerance > 0)
    || !Number.isInteger(maxCorrectorIterations) || maxCorrectorIterations < 0 || maxCorrectorIterations > 12
    || !Number.isInteger(maxCoordinateBacktracks) || maxCoordinateBacktracks < 0 || maxCoordinateBacktracks > 12)
    throw new Error('Invalid bounded wake-coordinate correction controls.');
  if (checkpoint?.continuation?.dissipationEnhancement)
    throw new Error('Wake-coordinate correction requires a fixed dissipation law for comparable residuals.');
  const create = cp => {
    const f = cp.restart;
    return createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  };
  const source = structuredClone(checkpoint), system = create(source), x = system.initial.slice();
  const value = system.admissibleValue(x, { requireConvex: true });
  if (!value || Object.keys(value.families).some(k => value.families[k] !== source.families[k]))
    throw new Error('Wake-coordinate correction requires an admissible, exactly replayed source.');
  const squared = residual => residual.reduce((sum, r) => sum + r * r, 0);
  const report = { method: 'bounded-wake-coordinate-corrector', accepted: false,
    beforeFamilies: value.families, beforeSquaredNorm: squared(value.residual), steps: [],
    linearDiagnostics: { solves: 0, refinements: 0, pivotRecoveries: 0, maxRelativeResidual: 0, iterations: [] },
    equationsChanged: false, physicalConditionsChanged: false };
  const massFractions = system.euler.decode(x.subarray(0, system.ne)).allocation.groups.map(group => {
    const total = group.reduce((sum, tube) => sum + tube.massFlow, 0), fractions = [0];
    for (const tube of group) fractions.push(fractions.at(-1) + tube.massFlow / total);
    fractions[fractions.length - 1] = 1;
    return fractions;
  });
  const pairing = initializeStreamtubeWakeCorrespondence({ nodes: value.outer.nodes,
    layout: system.euler.layout, massFractions });
  report.geometry = pairing.diagnostics;
  let current;
  report.coordinateTrials = [];
  for (let attempt = 0; attempt <= maxCoordinateBacktracks; attempt++) {
    const fraction = 2 ** -attempt, nodes = attempt === 0 ? pairing.nodes : value.outer.nodes.map((group, g) =>
      group.map((row, i) => row.map((p, j) => ({
        x: p.x + fraction * (pairing.nodes[g][i][j].x - p.x),
        y: p.y + fraction * (pairing.nodes[g][i][j].y - p.y),
      }))));
    const candidate = create(source), y = candidate.initial.slice(), trial = { fraction };
    report.coordinateTrials.push(trial);
    candidate.euler.setDisplacement(candidate.bl.thicknesses(y.subarray(candidate.ne)));
    y.set(candidate.euler.adoptGeometry(y.subarray(0, candidate.ne), nodes));
    const initial = candidate.admissibleValue(y, { requireConvex: true,
      onFailure: failure => { trial.admissibility = failure; } });
    if (!initial) continue;
    const next = { ...source, families: initial.families, restart: {
      input: source.restart.input,
      options: { ...source.restart.options, geometryReplay: 'preserve-undisplaced',
        transitionState: candidate.bl.snapshotActive() },
      initialEuler: { x: y.slice(0, candidate.ne), nodes: initial.outer.nodes,
        undisplacedNodes: initial.outer.undisplacedNodes }, initialBL: y.slice(candidate.ne),
    } };
    trial.families = initial.families;
    trial.squaredNorm = squared(initial.residual);
    if (attempt === 0) {
      current = next;
      report.uncommittedGuess = { families: initial.families, squaredNorm: trial.squaredNorm };
    }
    // First try reducing the coordinate movement itself. No BL, density,
    // transition or Newton state is reconstructed during this line search.
    try {
      report.decrease = requireCoupledResidualDecrease(value.residual, initial.residual, { step: fraction, tolerance });
      return { checkpoint: next, report: { ...report, accepted: true, coordinateFraction: fraction,
        reason: 'Wake-coordinate movement reduces the original coupled residual.', correctorIterations: 0,
        afterFamilies: initial.families, afterSquaredNorm: trial.squaredNorm } };
    } catch (error) {
      if (error.code !== 'COUPLED_RESIDUAL_DECREASE') throw error;
    }
  }
  if (!current) return { report: { ...report, reason: 'Wake-coordinate guess is not physically admissible.' } };
  for (let iteration = 1; iteration <= maxCorrectorIterations; iteration++) {
    const result = solveCoupledStreamtubeIses(undefined, { ...current.continuation, resume: current,
      maxIterations: 1, tolerance });
    report.linearDiagnostics = combinedCoupledLinearWork(report.linearDiagnostics, result.linearDiagnostics, iteration - 1);
    current = result.checkpoint;
    const step = result.history.at(-1);
    const replay = create(current), after = replay.admissibleValue(replay.initial, { requireConvex: true });
    if (!after || Object.keys(after.families).some(k => after.families[k] !== current.families[k]))
      throw new Error('Wake-coordinate corrector lost physical admissibility or exact replay.');
    const row = { iteration, families: result.families, squaredNorm: squared(after.residual),
      step: step.step, backtracks: step.backtracks, reason: result.reason,
      rejections: step.rejections, changes: step.changes };
    report.steps.push(row);
    onIteration?.(row);
    if (step.iteration > 0 && step.accepted !== false) {
      try {
        report.decrease = requireCoupledResidualDecrease(value.residual, after.residual, { step: 1, tolerance });
        return { checkpoint: current, report: { ...report, accepted: true,
          reason: 'Bounded correction reduces the original coupled residual.', correctorIterations: iteration,
          afterFamilies: after.families, afterSquaredNorm: squared(after.residual) } };
      } catch (error) {
        if (error.code !== 'COUPLED_RESIDUAL_DECREASE') throw error;
      }
    }
    if (result.converged || !(step.iteration > 0) || step.accepted === false)
      return { report: { ...report, reason: result.reason } };
  }
  return { report: { ...report, reason: 'Bounded correction did not reduce the original coupled residual.' } };
}

function combinedCoupledLinearWork(before, after, offset) {
  return { solves: before.solves + after.solves, refinements: before.refinements + after.refinements,
    pivotRecoveries: before.pivotRecoveries + after.pivotRecoveries,
    maxRelativeResidual: Math.max(before.maxRelativeResidual, after.maxRelativeResidual),
    iterations: [...before.iterations, ...after.iterations.map(h => ({ ...h, iteration: offset + h.iteration }))] };
}

// One correction inside the original work budget, selected by the observed
// paired-wake merit failure. All tentative coordinates stay private until
// the original residual decreases. A failed correction leaves the accepted
// source available for display and restart, without publishing its trial.
export function continueCoupledWakeCoordinates(before, { maxIterations, tolerance = 1e-10,
  onIteration, onCheckpoint, onMesh, onCorrection } = {}) {
  if (!Number.isInteger(maxIterations) || maxIterations < 0 || !Number.isFinite(tolerance) || tolerance <= 0)
    throw new Error('Invalid coupled wake continuation budget.');
  const cp = before?.checkpoint, rejection = before?.lastRejectedStep;
  const used = before?.history?.at(-1)?.iteration;
  if (before?.converged !== false || before.wakeCorrection || !Number.isInteger(used) || used < 0
    || used >= maxIterations || cp?.version !== 1 || cp.restart?.input?.wakeGeometry !== 'independent-banks'
    || cp.continuation?.dissipationEnhancement || !['armijo', 'event-armijo'].includes(cp.continuation?.stepAcceptance)
    || before.mesh?.quality?.valid !== true || before.mesh.quality.invalidCells?.length !== 0
    || rejection?.code !== 'COUPLED_RESIDUAL_DECREASE'
    || rejection.wakeCoordinateRepair?.method !== 'paired-wake-trial') return before;
  onCorrection?.({ phase: 'start', originalIterations: used, maxCorrectorIterations: Math.min(6, maxIterations - used) });
  const corrected = correctCoupledWakeCoordinates(cp, { tolerance,
    maxCorrectorIterations: Math.min(6, maxIterations - used),
    onIteration: h => onCorrection?.({ phase: 'trial', ...h, iteration: used + h.iteration }) });
  const work = corrected.report.steps.length;
  const evidence = { ...corrected.report, originalIterations: used, attemptedCorrectorIterations: work,
    maximumTotalIterations: maxIterations };
  onCorrection?.({ phase: corrected.report.accepted ? 'accepted' : 'rejected',
    originalIterations: used, correctorIterations: work, reason: evidence.reason });
  const correctedLinear = combinedCoupledLinearWork(before.linearDiagnostics, evidence.linearDiagnostics, used);
  if (!corrected.report.accepted) return { ...before, wakeCorrection: evidence, linearDiagnostics: correctedLinear };
  // Record work on uncommitted guesses against the retained physical state.
  // The new flow is observable only in the final, accepted correction frame.
  const correctionHistory = corrected.report.steps.map((h, i) => ({
    iteration: used + i + 1, step: 0, accepted: false, backtracks: h.backtracks,
    ...before.families, residual: Math.max(...Object.values(before.families)),
    wakeCorrection: { phase: 'trial', correctorIteration: h.iteration, squaredNorm: h.squaredNorm },
  }));
  const offset = used + work;
  let completeHistory;
  const result = solveCoupledStreamtubeIses(undefined, { ...corrected.checkpoint.continuation,
    resume: corrected.checkpoint, maxIterations: maxIterations - offset, tolerance,
    onIteration: h => {
      const entry = { ...h, iteration: offset + h.iteration,
        ...(h.iteration === 0 ? { accepted: true, wakeCorrection: { phase: 'accepted', correctorIterations: work } } : {}) };
      if (h.iteration === 0) {
        completeHistory = [...before.history, ...correctionHistory];
        // A pure coordinate correction consumes no Newton update.
        completeHistory[completeHistory.length - 1] = entry;
      } else completeHistory.push(entry);
      onIteration?.(entry);
    },
    onCheckpoint: (checkpoint, details) => onCheckpoint?.(checkpoint, { ...details,
      history: completeHistory, linearDiagnostics: combinedCoupledLinearWork(correctedLinear, details.linearDiagnostics, offset) }),
    onMesh: state => onMesh?.({ ...state, iteration: completeHistory.at(-1) }),
  });
  result.history = completeHistory;
  result.linearDiagnostics = combinedCoupledLinearWork(correctedLinear, result.linearDiagnostics, offset);
  result.wakeCorrection = { ...evidence, continuedIterations: result.history.at(-1).iteration - offset,
    converged: result.converged, families: { ...result.families } };
  return result;
}

// Section-speed matching supplies a forgiving viscous startup on an
// accelerated subsonic flow. Finish with pressure matching once the coupled
// grid/BL state is established: its stagnation-point limit remains regular.
// This is a fresh formulation at the SAME physical state, not a checkpoint
// resume with silently changed equations. Only the final formulation can
// certify convergence; both phases share the caller's Newton budget.
export function finishSubsonicPressureCoupling(before, { maxIterations, tolerance = 1e-10,
  linearOrdering, onIteration, onCheckpoint, onMesh, onCorrection } = {}) {
  if (!Number.isInteger(maxIterations) || maxIterations < 0 || !Number.isFinite(tolerance) || tolerance <= 0)
    throw new Error('Invalid subsonic coupling budget.');
  const cp = before?.checkpoint, used = before?.history?.at(-1)?.iteration;
  if (!cp || !Number.isInteger(used) || used >= maxIterations
    || cp.restart.options.edgeMatching !== 'section-velocity'
    || !(cp.restart.input.mach <= .3) || !(before.flow?.diagnostics?.maxMach < 1)
    || before.mesh?.quality?.valid !== true) return before;
  const source = cp.restart, options = { ...source.options, edgeMatching: 'pressure' };
  delete options.blThermodynamics;
  const evidence = { method: 'subsonic-pressure-handoff', startupIterations: used,
    maximumTotalIterations: maxIterations, startupFamilies: { ...before.families } };
  onCorrection?.(evidence);
  let completeHistory;
  const result = solveCoupledStreamtubeIses(source.input, { ...options,
    convergence: cp.continuation.convergence,
    initialEuler: source.initialEuler, initialBL: source.initialBL,
    maxIterations: maxIterations - used, tolerance,
    iterationGeometry: cp.continuation.iterationGeometry, stepAcceptance: cp.continuation.stepAcceptance,
    stagnationLimiter: cp.continuation.stagnationLimiter, blUpdate: cp.continuation.blUpdate,
    projectionGeometry: 'fixed', shearCoordinate: cp.continuation.shearCoordinate ?? 'linear',
    linearOrdering: linearOrdering ?? cp.continuation.linearOrdering ?? 'auto',
    onIteration: h => {
      const entry = { ...h, iteration: used + h.iteration, coupling: 'pressure' };
      if (h.iteration === 0) completeHistory = [...before.history.slice(0, -1), entry];
      else completeHistory.push(entry);
      onIteration?.(entry);
    },
    onCheckpoint: (checkpoint, details) => onCheckpoint?.(checkpoint, { ...details,
      history: completeHistory, linearDiagnostics: combinedCoupledLinearWork(before.linearDiagnostics, details.linearDiagnostics, used) }),
    onMesh: state => onMesh?.({ ...state, iteration: completeHistory.at(-1) }),
  });
  result.history = completeHistory;
  result.linearDiagnostics = combinedCoupledLinearWork(before.linearDiagnostics, result.linearDiagnostics, used);
  result.couplingStartup = { ...evidence, pressureIterations: result.history.at(-1).iteration - used,
    converged: result.converged, families: { ...result.families } };
  return result;
}
