import { shockRefinementRequest } from './streamtube-shock-refinement.js';
// SPDX-License-Identifier: GPL-2.0-or-later
import { observableFlow, observableBL } from '../streamtube-flow-preview.js';
import { transitionRecoveryPlan, recoverCoupledTransition } from '../streamtube-transition-recovery.js';
import { initializeCoupledStreamtubeFromFlow } from './streamtube-coupled-flow-restart.js';
import { streamtubeMeshSnapshot } from '../streamtube-mesh-preview.js';
import { solveCoupledStreamtubeIses } from '../streamtube-coupled-ises.js';

// Incidence continuation, optionally with adaptive Mach/alpha directions.
// Only converged roots advance the source; a bounded lower-Mach detour can
// cross an incidence barrier before restoring the requested Mach.
export function solveCoupledStreamtubeAlpha(targetAlpha, { initialCheckpoint, tolerance = 1e-10,
  maxAlphaStep = .1, minAlphaStep = Math.min(1e-4, maxAlphaStep), maxTransitionRecoveries = 2, maxStages = 80, maxSubdivisions = 6, stageMaxIterations = 20,
  onStage, onIteration, onCheckpoint, onMesh, onFlow, normalization,
  targetMach, maxMachStep = .05, stopAtShock = false, preferAlphaApproach = false } = {}) {
  if (typeof stopAtShock !== 'boolean' || typeof preferAlphaApproach !== 'boolean' || !Number.isFinite(targetAlpha) || !(Number.isFinite(maxAlphaStep) && maxAlphaStep > 0)
    || !(Number.isFinite(minAlphaStep) && minAlphaStep > 0 && minAlphaStep <= maxAlphaStep)
    || !Number.isInteger(maxTransitionRecoveries) || maxTransitionRecoveries < 0 || maxTransitionRecoveries > 4
    || !(Number.isFinite(tolerance) && tolerance > 0)
    || ![maxStages, maxSubdivisions, stageMaxIterations].every(Number.isInteger)
    || maxStages < 1 || maxSubdivisions < 0 || stageMaxIterations < 0)
    throw new Error('Invalid alpha continuation controls.');
  let mach = initialCheckpoint?.restart?.input?.mach;
  const combined = targetMach !== undefined;
  targetMach ??= mach;
  if (!(Number.isFinite(targetMach) && targetMach > 0 && targetMach < 1)
    || !(Number.isFinite(maxMachStep) && maxMachStep > 0 && maxMachStep <= 1))
    throw new Error('Invalid operating-point continuation controls.');
  const stage = combined ? 'coupled-operating-point' : 'coupled-alpha';
  const sourceMach = mach;
  let machStepCap = maxMachStep;
  let scale = 1, direction = 0, detourUsed = false, detourMach, detourPending = false;

  const source = initializeCoupledStreamtubeFromFlow(mach, initialCheckpoint, { tolerance });
  let cp = source.checkpoint;
  const controls = c => ({ iterationGeometry: c.continuation.iterationGeometry,
    stepAcceptance: c.continuation.stepAcceptance, stagnationLimiter: c.continuation.stagnationLimiter });
  let current = solveCoupledStreamtubeIses(undefined, { resume: cp, ...controls(cp), maxIterations: 0, tolerance });
  if (!current.converged || !current.mesh?.quality?.valid) throw new Error('Alpha source must replay as a converged valid state.');
  let alpha = cp.restart.input.alpha ?? 0, increment = Math.sign(targetAlpha - alpha) * maxAlphaStep, subdivisions = 0;
  const sourceAlpha = alpha, attempts = [], transitionRecoveries = [];
  let stopReason, refinementRequested;
  let acceptedFrame = onFlow ? { checkpoint: cp, flow: observableFlow(current.flow),
    bl: observableBL(source.system.bl), bodies: cp.restart.input.bodies, normalization,
    iteration: current.history.at(-1), mach, actualAlpha: alpha, targetAlpha, stage } : undefined;
  const retainFlow = info => { if (onFlow && acceptedFrame) emit(onFlow, { ...acceptedFrame, ...info,
    mach, actualMach: mach, targetMach, alpha, actualAlpha: alpha, retained: true, stage }); };
  const emit = (fn, ...args) => fn?.(...args.map(structuredCloneSafe));
  function structuredCloneSafe(v) { return structuredClone(v); }
  while ((alpha !== targetAlpha || mach !== targetMach) && attempts.length < maxStages) {
    if (stopAtShock && (refinementRequested = shockRefinementRequest(current))) {
      stopReason = refinementRequested.reason; break;
    }
    let next = Math.abs(targetAlpha - alpha) <= Math.abs(increment) ? targetAlpha : alpha + increment;
    let nextMach = mach;
    let stepMethod = 'alpha';
    if (combined) {
      const aimMach = detourPending ? detourMach : detourMach !== undefined && alpha !== targetAlpha ? mach : targetMach;
      const da = targetAlpha - alpha, dm = aimMach - mach;
      const directions = da !== 0 && dm !== 0 && !detourPending ? (preferAlphaApproach ? ['mach', 'combined', 'alpha'] : ['combined', 'alpha', 'mach'])
        : [dm !== 0 ? 'mach' : 'alpha'];
      stepMethod = directions[direction % directions.length];
      next = stepMethod === 'mach' ? alpha : alpha + Math.sign(da) * Math.min(Math.abs(da), maxAlphaStep * scale);
      nextMach = stepMethod === 'alpha' ? mach : mach + Math.sign(dm) * Math.min(Math.abs(dm), Math.min(machStepCap, maxMachStep * scale));
      // Assign exact endpoints, avoiding accumulation at the final step.
      if (stepMethod !== 'mach' && Math.abs(da) <= maxAlphaStep * scale) next = targetAlpha;
      if (stepMethod !== 'alpha' && Math.abs(dm) <= Math.min(machStepCap, maxMachStep * scale)) nextMach = aimMach;
    }
    if (next === alpha && nextMach === mach) { stopReason = 'floating-point step limit'; break; }
    if (!combined && Math.abs(increment) < minAlphaStep && Math.abs(targetAlpha - alpha) > minAlphaStep) {
      stopReason = 'minimum alpha step'; break;
    }
    const info = { stage, alpha: next, actualAlpha: next, targetAlpha, mach: nextMach, actualMach: nextMach,
      targetMach, stepMethod: detourPending ? 'lower-mach-detour' : stepMethod, lastAcceptedMach: mach, lastAcceptedAlpha: alpha, attempt: attempts.length + 1 };
    emit(onStage, info);
    let trial, error, trialFrame, observerFailed = false;
    try {
      const seed = initializeCoupledStreamtubeFromFlow(nextMach, cp, { tolerance, targetAlpha: next });
      let iterationCheckpoint;
      trial = solveCoupledStreamtubeIses(undefined, { resume: seed.checkpoint, ...controls(cp),
        maxIterations: stageMaxIterations, tolerance,
        onIteration: e => { try { emit(onIteration, { ...e, ...info }); } catch (e) { observerFailed = true; throw e; } },
        onCheckpoint: onFlow ? value => { iterationCheckpoint = value; } : undefined,
        onMesh: onMesh || onFlow ? snapshot => { try {
          if (onMesh) emit(onMesh, streamtubeMeshSnapshot(snapshot), 'solving', { ...info, iteration: snapshot.iteration });
          if (onFlow && iterationCheckpoint) {
            trialFrame = { checkpoint: iterationCheckpoint, flow: observableFlow(snapshot.flow),
              bl: observableBL(seed.system.bl), bodies: snapshot.system.layout.bodies, normalization,
              ...info, iteration: snapshot.iteration };
            emit(onFlow, trialFrame);
          }
        } catch (e) { observerFailed = true; throw e; } } : undefined });
    } catch (e) { if (observerFailed) throw e; error = e; }
    if (trial && !trial.converged && transitionRecoveries.length < maxTransitionRecoveries) {
      const plan = transitionRecoveryPlan(trial, { tolerance });
      if (plan) {
        const recovery = { alpha: next, reason: plan.reason, parentNx: plan.parentNx, refinedNx: plan.refinedNx, accepted: false };
        transitionRecoveries.push(recovery);
        const observe = fn => (...args) => { try { emit(fn, ...args); } catch (e) { observerFailed = true; throw e; } };
        try {
          const refined = recoverCoupledTransition(trial, { plan, tolerance, maxIterations: stageMaxIterations, normalization,
            onFlow: onFlow ? observe(frame => {
              trialFrame = { ...frame, ...info, transitionRecovery: plan }; emit(onFlow, trialFrame);
            }) : undefined,
            onStage: observe(e => emit(onStage, { ...info, ...e, stage, transitionRecovery: plan })),
            onIteration: observe(e => emit(onIteration, { ...e, ...info, transitionRecovery: plan })),
            onMesh: onMesh ? observe((mesh, phase) => emit(onMesh, mesh, phase, { ...info, transitionRecovery: plan })) : undefined });
          recovery.families = refined?.families;
          recovery.reason = refined?.reason ?? plan.reason;
          if (refined?.converged && refined.mesh?.quality?.valid
            && Object.values(refined.families).every(v => Number.isFinite(v) && v <= tolerance)) {
            trial = refined; recovery.accepted = true;
          } else observe(onMesh)(current.mesh, 'retained', { ...info, alpha, actualAlpha: alpha, retained: true });
        } catch (e) {
          if (observerFailed) throw e;
          recovery.reason = e.message;
          emit(onMesh, current.mesh, 'retained', { ...info, alpha, actualAlpha: alpha, retained: true });
        }
      }
    }
    const accepted = trial?.converged === true && trial.mesh?.quality?.valid === true
      && ['euler', 'boundaryLayer', 'edgeMatching'].every(k => Number.isFinite(trial.families?.[k])
        && trial.families[k] >= 0 && trial.families[k] <= tolerance)
      && (!combined || trial.checkpoint?.restart.input.mach === nextMach
        && (trial.checkpoint.restart.input.alpha ?? 0) === next);
    attempts.push({ mach: nextMach, stepMethod: info.stepMethod, alpha: next, increment: next - alpha, accepted, reason: error?.message ?? trial?.reason,
      families: trial?.families, progress: trial?.history?.at(-1)?.progress,
      transition: trial?.checkpoint?.restart?.options?.transitionState });
    if (accepted) {
      if (nextMach !== mach) machStepCap = Math.min(maxMachStep, machStepCap * 1.5);
      current = trial; cp = structuredClone(trial.checkpoint); alpha = next; mach = nextMach; subdivisions = 0;
      direction = 0; scale = Math.min(1, scale * 1.5);
      if (detourPending && mach === detourMach) detourPending = false;
      if (onFlow && trialFrame) acceptedFrame = trialFrame;
      emit(onCheckpoint, cp, { ...info, kind: 'accepted', reachedTarget: alpha === targetAlpha && mach === targetMach });
      increment = Math.sign(targetAlpha - alpha) * Math.min(maxAlphaStep, Math.abs(increment) * 1.5);
    } else {
      retainFlow(info);
      // An alpha-only success must not erase a rejected Mach step. Otherwise
      // an inadmissible Mach transfer repeats unchanged after every alpha step.
      if (combined && nextMach !== mach) machStepCap = Math.max(
        maxMachStep * minAlphaStep / maxAlphaStep,
        Math.min(machStepCap, Math.abs(nextMach - mach) * .5));
      if (combined && detourMach === undefined && !detourPending && alpha !== targetAlpha && mach !== targetMach && direction < 2) {
        direction++; continue;
      }
      direction = 0;
      if (++subdivisions > maxSubdivisions || combined && scale * .5 < minAlphaStep / maxAlphaStep) {
        if (combined && !detourUsed && alpha !== targetAlpha && mach > .3) {
          detourUsed = true; detourPending = true; detourMach = Math.max(.2, mach - Math.min(.02, maxMachStep));
          scale = 1; subdivisions = 0; continue;
        }
        stopReason = 'subdivision limit'; break;
      }
      scale *= .5;
      increment *= .5;
    }
  }
  const reachedTarget = alpha === targetAlpha && mach === targetMach;
  return { ...current, checkpoint: cp, mach, actualMach: mach, targetMach, alpha, actualAlpha: alpha, targetAlpha,
    ...(refinementRequested ? { refinementRequested } : {}),
    converged: reachedTarget && current.converged, stateConverged: current.converged,
    ...(combined ? { operatingPointContinuation: { sourceMach, sourceAlpha, targetMach, targetAlpha,
      actualMach: mach, actualAlpha: alpha, reachedTarget, attempts, detourUsed, stopReason } } : {}),
    alphaContinuation: { sourceAlpha, targetAlpha, actualAlpha: alpha, reachedTarget: alpha === targetAlpha, maxAlphaStep, minAlphaStep, attempts, transitionRecoveries,
      ...(!reachedTarget ? { stopReason: stopReason ?? 'stage limit', lastAttempt: attempts.at(-1) } : {}) },
    ...(!reachedTarget ? { reason: `Continuation retained Mach ${mach}, alpha ${alpha}°; requested Mach ${targetMach}, alpha ${targetAlpha}° was not reached (${stopReason ?? 'stage limit'}). Last attempt: ${attempts.at(-1)?.reason ?? 'unavailable'}.` } : {}) };
}
