// SPDX-License-Identifier: GPL-2.0-or-later
// Validation-only MSES §1.2.7 composition. Each flow call is an unchanged,
// complete-checkpoint ISES resume; no Newton/event/maintenance code is copied.
import { solveCoupledStreamtubeIses } from '../../src/euler/streamtube-coupled-ises.js';
import { temporaryShockMcrit, maximumAppliedDensityChange, rebaseCoupledMcrit } from './coupled-shock-broadening.js';
import { isDeepStrictEqual } from 'node:util';

const copy = value => value === undefined ? undefined : JSON.parse(JSON.stringify(value,
  (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const mcrit = cp => cp.restart.input.upwind.mcrit;
const maxResidual = residual => {
  if (!residual || !Number.isInteger(residual.length) || residual.length === 0)
    throw new Error('Missing shock-broadening residual.');
  let maximum = 0;
  for (const value of residual) {
    if (!Number.isFinite(value)) throw new Error('Nonfinite shock-broadening residual.');
    maximum = Math.max(maximum, Math.abs(value));
  }
  return maximum;
};
const errorData = error => ({ name: error?.name ?? 'Error', message: error?.message ?? String(error),
  ...copy(Object.fromEntries(Object.entries(error ?? {}))) });
const policy = { formula: 'mses-1.2.7-preferred', epsilon: .15, densityChange: 'accepted-physical-ratio' };

// The dependency factory permits state-machine tests without a flow solve.
// Ordinary callers use the production-backed export below.
export function createCoupledShockBroadeningDriver({ solve = solveCoupledStreamtubeIses,
  rebase = rebaseCoupledMcrit, threshold = temporaryShockMcrit,
  densityChange = maximumAppliedDensityChange } = {}) {
  return function solveCoupledWithShockBroadening(source, options = {}) {
    const { maxIterations = 4, tolerance = source?.kind === 'mses-shock-broadening'
      ? source.tolerance ?? 1e-10 : 1e-10, onRecord, onIteration } = options;
    if (!Number.isInteger(maxIterations) || maxIterations < 0 || !Number.isFinite(tolerance) || tolerance <= 0
      || onRecord !== undefined && typeof onRecord !== 'function'
      || onIteration !== undefined && typeof onIteration !== 'function')
      throw new Error('Invalid shock-broadening controller controls.');
    const saved = copy(source), wrapped = saved?.kind === 'mses-shock-broadening';
    let controller = wrapped ? saved : { version: 1, kind: 'mses-shock-broadening', checkpoint: saved,
      targetMcrit: saved?.restart?.input?.upwind?.mcrit, previousDensityChange: 0, acceptedUpdates: 0 };
    if (controller.version !== 1 || controller.kind !== 'mses-shock-broadening'
      || controller.checkpoint?.version !== 1 || !controller.checkpoint.restart || !controller.checkpoint.continuation
      || !Number.isInteger(controller.acceptedUpdates) || controller.acceptedUpdates < 0
      || !Number.isFinite(controller.previousDensityChange) || controller.previousDensityChange < 0
      || !Number.isInteger(controller.linearSolves ?? 0) || (controller.linearSolves ?? 0) < 0)
      throw new Error('Invalid shock-broadening controller checkpoint.');
    // Also validates the requested threshold even for a zero-update call.
    threshold({ targetMcrit: controller.targetMcrit, densityChange: controller.previousDensityChange });
    if (controller.policy !== undefined && !isDeepStrictEqual(controller.policy, policy)
      || controller.tolerance !== undefined && controller.tolerance !== tolerance
      || controller.effectiveMcrit !== undefined && controller.effectiveMcrit !== mcrit(controller.checkpoint)
      || controller.nextAction !== undefined && !['scheduled', 'target-check', 'complete'].includes(controller.nextAction))
      throw new Error('Shock-broadening controller history or controls do not match.');
    controller.linearSolves ??= 0;
    controller.policy = copy(policy); controller.tolerance = tolerance;
    controller.effectiveMcrit = mcrit(controller.checkpoint); controller.nextAction ??= 'scheduled';
    const history = [], rebases = [], calls = [], startUpdates = controller.acceptedUpdates,
      startSolves = controller.linearSolves;
    let result, failure, reason = 'iteration limit', stage = 'source', targetChecked = false,
      validated = false, requestedMcrit, callsWithUpdates = 0, zeroUpdateCalls = 0, observerFailed = false;
    const emit = (kind, value) => {
      const previousStage = stage; stage = `observer ${kind}`;
      try { onRecord?.(kind, copy(value)); } catch (error) { observerFailed = true; throw error; }
      stage = previousStage;
    };
    const finish = converged => ({ validationOnly: true, physicalAcceptance: false, converged, reason,
      checkpoint: copy(controller.checkpoint), controller: copy(controller), history: copy(history),
      result, ...(failure ? { error: failure } : {}),
      targetRestored: validated && mcrit(controller.checkpoint) === controller.targetMcrit, targetChecked,
      linearDiagnostics: { solves: controller.linearSolves - startSolves,
        cumulativeSolves: controller.linearSolves, calls: copy(calls) },
      operations: { acceptedUpdates: controller.acceptedUpdates - startUpdates, rebases: rebases.length,
        ordinaryCalls: callsWithUpdates, zeroUpdateCalls, linearSolves: controller.linearSolves - startSolves },
      rebases: copy(rebases), controls: { maxIterations, tolerance, maxBacktracks: 12,
        densityChange: 'max |expm1(z_after-z_before)| over accepted physical density columns; explicit port convention' } });
    const ordinary = (limit, densityCount) => {
      const before = copy(controller), c = before.checkpoint.continuation;
      let reportedSolves = 0, reportedUpdates = 0;
      const account = diagnostics => {
        const solves = diagnostics?.solves;
        if (!Number.isInteger(solves) || solves < reportedSolves || solves > limit)
          throw new Error('Unexpected ordinary ISES linear-solve count.');
        controller.linearSolves += solves - reportedSolves; reportedSolves = solves;
      };
      if (limit) callsWithUpdates++; else zeroUpdateCalls++;
      stage = limit ? 'ordinary Newton update' : 'target zero-update check';
      try {
        const next = solve(undefined, { resume: copy(before.checkpoint), maxIterations: limit,
          tolerance, maxBacktracks: 12, iterationGeometry: c.iterationGeometry,
          stepAcceptance: c.stepAcceptance, stagnationLimiter: c.stagnationLimiter, blUpdate: c.blUpdate ?? 'giles',
          onCheckpoint: (cp, details) => {
            account(details.linearDiagnostics);
            const entry = details.history.at(-1);
            if (entry.iteration === 0) return; // Complete resume publishes zero again.
            if (limit !== 1 || reportedUpdates !== 0 || entry.iteration !== 1)
              throw new Error('More than one ordinary accepted update was reported.');
            if (mcrit(cp) !== mcrit(before.checkpoint)) throw new Error('Ordinary update changed the frozen threshold.');
            const d = densityChange(before.checkpoint.restart.initialEuler.x, cp.restart.initialEuler.x, densityCount);
            controller = { ...controller, checkpoint: copy(cp), previousDensityChange: d,
              acceptedUpdates: before.acceptedUpdates + 1,
              nextAction: maxResidual(Object.values(cp.families)) <= tolerance
                || before.acceptedUpdates + 1 - startUpdates === maxIterations ? 'target-check' : 'scheduled' };
            reportedUpdates++;
            const record = { ...copy(entry), iteration: controller.acceptedUpdates,
              localIteration: entry.iteration, newIteration: controller.acceptedUpdates - startUpdates,
              effectiveMcrit: mcrit(cp), previousDensityChange: before.previousDensityChange, appliedDensityChange: d };
            history.push(record);
            emit('update', { controller, checkpoint: controller.checkpoint, iteration: record });
            stage = 'observer iteration';
            try { onIteration?.(copy(record)); } catch (error) { observerFailed = true; throw error; }
            stage = 'ordinary Newton update';
          } });
        account(next.linearDiagnostics);
        const accepted = next.history.filter(entry => entry.iteration > 0);
        if (accepted.length !== reportedUpdates) throw new Error('An accepted ISES update lacks its complete checkpoint.');
        return { result: next, updates: reportedUpdates };
      } finally {
        calls.push({ maxIterations: limit, acceptedUpdates: reportedUpdates, linearSolves: reportedSolves,
          effectiveMcrit: mcrit(before.checkpoint) });
      }
    };
    try {
      emit('source', { controller });
      // Check a possible existing root under its requested equations before
      // broadening it. A zero-budget call is also an exact target check.
      let checkTarget = controller.nextAction !== 'scheduled' || maxIterations === 0
        || maxResidual(Object.values(controller.checkpoint.families)) <= tolerance;
      for (;;) {
        const remaining = maxIterations - (controller.acceptedUpdates - startUpdates);
        if (remaining === 0) checkTarget = true;
        requestedMcrit = checkTarget ? controller.targetMcrit : threshold({ targetMcrit: controller.targetMcrit,
          densityChange: controller.previousDensityChange });
        controller.nextAction = checkTarget ? 'target-check' : 'scheduled';
        if (rebases.length >= 2 * maxIterations + 2) throw new Error('Shock-broadening threshold-event budget exhausted.');
        stage = 'threshold rebase';
        const prepared = rebase(copy(controller.checkpoint), requestedMcrit);
        if (mcrit(prepared.checkpoint) !== requestedMcrit) throw new Error('Rebase did not install the exact requested threshold.');
        controller = { ...controller, checkpoint: copy(prepared.checkpoint), effectiveMcrit: requestedMcrit }; validated = true;
        targetChecked = false;
        rebases.push(copy(prepared.diagnostics));
        emit('rebase', { controller, checkpoint: controller.checkpoint, diagnostics: prepared.diagnostics });
        const residual = maxResidual(prepared.residual);
        if (residual <= tolerance && requestedMcrit !== controller.targetMcrit) {
          checkTarget = true; continue; // Never factor a broadened root or reuse stale d indefinitely.
        }
        if (residual <= tolerance || remaining === 0) {
          // The unchanged result adapter enforces the final strict convexity
          // gate even when ordinary iterates permit positive-simple cells.
          const checked = ordinary(0, prepared.diagnostics.densityCount); result = checked.result;
          const cp = controller.checkpoint, r = cp.restart;
          if (!isDeepStrictEqual(copy(result.checkpoint), cp)
            || !isDeepStrictEqual(copy(result.x), [...r.initialEuler.x, ...r.initialBL])
            || !isDeepStrictEqual(copy(result.flow?.nodes), r.initialEuler.nodes)
            || !isDeepStrictEqual(copy(result.flow?.undisplacedNodes), r.initialEuler.undisplacedNodes)
            || !isDeepStrictEqual(copy(result.boundaryLayer?.transitionState), r.options.transitionState)
            || !isDeepStrictEqual(copy(result.residual), copy(prepared.residual)))
            throw new Error('Target finalization does not exactly match the retained controller state.');
          targetChecked = true;
          const converged = requestedMcrit === controller.targetMcrit && result.converged === true
            && result.mesh?.quality?.valid === true && maxResidual(result.residual) <= tolerance;
          reason = converged ? 'requested-threshold residual' : result.reason;
          controller.nextAction = converged ? 'complete' : 'scheduled';
          emit('target-check', { controller, checkpoint: controller.checkpoint, result });
          return finish(converged);
        }
        const next = ordinary(1, prepared.diagnostics.densityCount); result = next.result;
        emit('ordinary-result', { controller, checkpoint: controller.checkpoint, result });
        if (next.updates === 0) {
          reason = result.reason; failure = { stage: 'ordinary Newton update', message: result.reason,
            lastRejectedStep: copy(result.lastRejectedStep ?? null) };
          return finish(false);
        }
        // Intermediate nonconvex reports are retained under ises-sampled.
        // They do not become an extra rejection rule in this controller.
        checkTarget = maxResidual(result.residual) <= tolerance;
      }
    } catch (error) {
      if (observerFailed) throw error; // Preserve primitive/frozen cancellation values unchanged.
      failure = { ...errorData(error), stage, requestedMcrit };
      reason = `Shock-broadening ${stage} failed: ${error?.message ?? String(error)}`;
      return finish(false);
    }
  };
}

export const solveCoupledWithShockBroadening = createCoupledShockBroadeningDriver();
