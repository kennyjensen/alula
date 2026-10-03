import { shockRefinementRequest } from './streamtube-shock-refinement.js';
// SPDX-License-Identifier: GPL-2.0-or-later
// MPOLAR-style bounded Mach continuation from a complete accepted Euler/BL
// checkpoint. Each stage uses the unchanged simultaneous ISES update. Failed
// stages restart from the last converged state with a halved Mach increment.
import { initializeCoupledStreamtubeFromFlow } from './streamtube-coupled-flow-restart.js';
import { solveCoupledStreamtubeIses } from '../streamtube-coupled-ises.js';
import { streamtubeMeshSnapshot } from '../streamtube-mesh-preview.js';
import { coupledWakeRecoveryPlan, recoverCoupledWakeCorrespondence } from './streamtube-coupled-wake-recovery.js';

import { coupledFirstOrderRecoveryEligible, recoverCoupledDissipationOrder } from './streamtube-coupled-dissipation-recovery.js';

import { auditStreamtubeShocks } from '../streamtube-shock-audit.js';
import { coupledMachGrowth, compareShockAudits } from './streamtube-shock-continuation.js';

const maximum = values => values.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
// The evaluator's bodyPressure accessor closes over cells. Progress keeps
// those exact numeric cells, but never sends a live function to the caller.
const physicalFlowData = ({ bodyPressure, ...flow }) => flow;
const controlsFrom = checkpoint => {
  const { iterationGeometry, stepAcceptance, stagnationLimiter } = checkpoint.continuation;
  return { iterationGeometry, stepAcceptance, stagnationLimiter };
};

export function solveCoupledStreamtubeAutomatic(targetMach, { initialCheckpoint, tolerance = 1e-10,
  stageMaxIterations = 12, maxMachStep = 1, maxSubdivisions = 5, maxStages = 32, maxBacktracks = 12,
  maxFirstOrderRecoveries = 0, maxWakeRecoveries = 0, blPredictor = 'preserve', dissipationEnhancement, iterationRecovery, maxProgressExtraIterations = 0,
  stopAtShock = false, onStage, onIteration, onMesh, onCheckpoint, includeFlowState = false } = {}) {
  dissipationEnhancement ??= initialCheckpoint?.continuation?.dissipationEnhancement !== undefined;
  iterationRecovery ??= initialCheckpoint?.continuation?.iterationRecovery === true;
  if (typeof stopAtShock !== 'boolean' || !Number.isFinite(targetMach) || targetMach <= 0 || targetMach >= 1
    || !Number.isFinite(maxMachStep) || maxMachStep <= 0 || maxMachStep > 1
    || !Number.isFinite(tolerance) || tolerance <= 0
    || !Number.isInteger(stageMaxIterations) || stageMaxIterations < 0
    || !Number.isInteger(maxSubdivisions) || maxSubdivisions < 0 || maxSubdivisions > 20
    || !Number.isInteger(maxStages) || maxStages < 1 || maxStages > 256
    || !Number.isInteger(maxBacktracks) || maxBacktracks < 0 || maxBacktracks > 20
    || !Number.isInteger(maxFirstOrderRecoveries) || maxFirstOrderRecoveries < 0 || maxFirstOrderRecoveries > 4
    || !Number.isInteger(maxWakeRecoveries) || maxWakeRecoveries < 0 || maxWakeRecoveries > 4
    || !['preserve', 'xfoil-mrchdu'].includes(blPredictor)
    || typeof dissipationEnhancement !== 'boolean' || typeof includeFlowState !== 'boolean' || typeof iterationRecovery !== 'boolean'
    || !Number.isInteger(maxProgressExtraIterations) || maxProgressExtraIterations < 0 || maxProgressExtraIterations > 20
    || [onStage, onIteration, onMesh, onCheckpoint].some(fn => fn !== undefined && typeof fn !== 'function'))
    throw new Error('Invalid automatic coupled Mach-continuation controls.');
  const sourceMach = initialCheckpoint?.restart?.input?.mach;
  if (!Number.isFinite(sourceMach) || sourceMach <= 0 || sourceMach >= 1)
    throw new Error('Automatic coupled continuation requires a complete accepted source checkpoint.');
  // Validate source equations, packed state, BL phases, geometry and update
  // history before permitting any recovery. A bad source is an input error,
  // not a reason to try a sequence of unrelated operating points.
  let current, currentCheckpoint, currentMach = sourceMach;
  {
    const prepared = initializeCoupledStreamtubeFromFlow(sourceMach, initialCheckpoint, { tolerance });
    current = solveCoupledStreamtubeIses(undefined, { resume: prepared.checkpoint,
      ...controlsFrom(prepared.checkpoint), maxIterations: 0, maxBacktracks, tolerance });
    if (!current.converged || !current.mesh.quality.valid || maximum(current.residual) > tolerance)
      throw new Error('Automatic coupled source must replay as a converged state on a convex grid.');
    currentCheckpoint = structuredClone(prepared.checkpoint);
  }
  const attempts = [], wakeRecoveries = [], firstOrderRecoveries = [];
  const finish = (reachedTarget, reason) => ({ ...current,
    // Retain the exact checkpoint seed supplied by the complete transfer,
    // including its original undisplaced coordinates and ISES history.
    checkpoint: structuredClone(currentCheckpoint), converged: reachedTarget && current.converged,
    status: reachedTarget ? current.status : 'research-coupled-target-not-reached',
    stateConverged: true,
    ...(reason ? { reason } : {}),
    continuation: { method: 'freestream-mach', solver: 'coupled-ises', sourceMach, targetMach, currentMach,
      reachedTarget, attempts, ...(maxMachStep === 1 ? {} : { maxMachStep }), maxStages, maxSubdivisions, stageMaxIterations,
      coldInitializationSkipped: true, densityReinitializedOnWarmRestart: false,
      boundaryLayerReinitializedOnWarmRestart: attempts.some(a => a.transfer?.targetBLPrediction), initialRedistributionSkippedOnWarmRestart: true,
      ...(blPredictor === 'preserve' ? {} : { blPredictor }),
      analyticShockInformationUsed: false },
    ...(maxWakeRecoveries ? { wakeRecovery: { maximumAttempts:maxWakeRecoveries, attempts:structuredClone(wakeRecoveries) } } : {}),
    ...(maxFirstOrderRecoveries ? { firstOrderRecovery: { maximumAttempts: maxFirstOrderRecoveries, attempts: structuredClone(firstOrderRecoveries) } } : {}),
    physicalAcceptance: false, fullSolverComplete: false,
  });
  const refinementStop = () => {
    const request = stopAtShock && shockRefinementRequest(current);
    return request ? { ...finish(false, request.reason), refinementRequested: request } : null;
  };
  const initialRefinement = refinementStop();
  if (initialRefinement) return initialRefinement;
  if (sourceMach === targetMach) return finish(true);

  let increment = Math.sign(targetMach - sourceMach) * Math.min(maxMachStep, Math.abs(targetMach - sourceMach)), subdivisions = 0;
  while (currentMach !== targetMach && attempts.length < maxStages) {
    const distance = targetMach - currentMach;
    const mach = Math.abs(increment) >= Math.abs(distance) ? targetMach : currentMach + increment;
    if (mach === currentMach) return finish(false, 'Automatic coupled Mach increment cannot advance in floating-point arithmetic.');
    const info = { label: attempts.length === 0 ? 'target-warm' : 'mach-continuation',
      attempt: attempts.length + 1, mach, targetMach, sourceMach, lastAcceptedMach: currentMach,
      fraction: (mach - sourceMach) / (targetMach - sourceMach), subdivisions };
    // Callback arguments never alias a live solver or continuation state.
    let observerFailed = false, observerError;
    const observe = (fn, ...args) => {
      if (!fn) return;
      try { fn(...args.map(arg => structuredClone(arg))); }
      catch (error) { observerFailed = true; observerError = error; throw error; }
    };
    observe(onStage, info);
    const history = []; let stage = 'transfer', trial, transfer;
    try {
      {
        const prepared = initializeCoupledStreamtubeFromFlow(mach, currentCheckpoint, { tolerance, blPredictor });
        transfer = { checkpoint: prepared.checkpoint, diagnostics: prepared.diagnostics };
      }
      stage = 'newton';
      trial = solveCoupledStreamtubeIses(undefined, { resume: transfer.checkpoint,
        ...controlsFrom(transfer.checkpoint), maxIterations: stageMaxIterations, tolerance, maxBacktracks, dissipationEnhancement,
        iterationRecovery, maxProgressExtraIterations,
        onIteration: entry => {
          history.push(structuredClone(entry));
          observe(onIteration, { ...entry, ...info });
        },
        onMesh: onMesh ? snapshot => observe(onMesh, {
          mesh: streamtubeMeshSnapshot(snapshot), iteration: snapshot.iteration,
          families: snapshot.coupledFamilies, ...info }, ...(includeFlowState ? [{
            flow: physicalFlowData(snapshot.flow), iteration: snapshot.iteration, coupledFamilies: snapshot.coupledFamilies,
          }] : [])) : undefined,
        onCheckpoint: onCheckpoint ? (checkpoint, details) => observe(onCheckpoint, checkpoint,
          { ...details, ...info, kind: 'iterate' }) : undefined,
      });
      // Preserve cancellation even if an underlying numerical driver ever
      // catches an observer exception as a rejected step.
      if (observerFailed) throw observerError;
      if (stageMaxIterations > 0 && firstOrderRecoveries.length < maxFirstOrderRecoveries
        && coupledFirstOrderRecoveryEligible(trial)) {
        const record = { attempt: firstOrderRecoveries.length + 1, mach, afterMachAttempt: info.attempt,
          accepted: false, originalReason: trial.reason, originalResidual: maximum(trial.residual),
          originalLinearSolves: trial.linearDiagnostics.solves };
        firstOrderRecoveries.push(record);
        const recoveryInfo = phase => ({ ...info, stage: 'coupled-mach', dissipationRecovery: phase });
        try {
          const recovered = recoverCoupledDissipationOrder(trial.checkpoint, {
            tolerance, maxIterations: stageMaxIterations, maxBacktracks, dissipationEnhancement, iterationRecovery,
            onPhase: phase => observe(onStage, recoveryInfo(phase)),
            onIteration: (entry, phase) => observe(onIteration, { ...entry, ...recoveryInfo(phase) }),
            onMesh: onMesh ? (snapshot, phase) => observe(onMesh, {
              mesh: streamtubeMeshSnapshot(snapshot), iteration: snapshot.iteration,
              families: snapshot.coupledFamilies, ...recoveryInfo(phase) }, ...(includeFlowState ? [{
                flow: physicalFlowData(snapshot.flow), iteration: snapshot.iteration, coupledFamilies: snapshot.coupledFamilies,
              }] : [])) : undefined,
            onCheckpoint: onCheckpoint ? (cp, details) => observe(onCheckpoint, cp,
              { ...details, ...recoveryInfo(details.dissipationRecovery), kind: 'dissipation-initializer-iterate' }) : undefined,
          });
          if (observerFailed) throw observerError;
          Object.assign(record, recovered.diagnostics, { accepted: recovered.accepted });
          if (recovered.accepted) trial = recovered.result;
        } catch (error) {
          if (observerFailed) throw observerError;
          record.reason = error.message;
        }
        // Restore the accepted-source display checkpoint if the temporary
        // equations did not yield a requested-equation root.
        if (!record.accepted) observe(onCheckpoint, currentCheckpoint, { ...info, mach: currentMach, kind: 'retained' });
      }
      const accepted = trial.converged && trial.conditions.mach === mach
        && trial.mesh.quality.valid && maximum(trial.residual) <= tolerance;
      const dissipationAttempt = firstOrderRecoveries.find(r => r.afterMachAttempt === info.attempt);
      attempts.push({ ...info, converged: accepted, reason: trial.reason,
        iterations: Math.max(0, history.length - 1) + (dissipationAttempt?.phases ?? []).reduce((n, p) => n + p.iterations, 0),
        globalLinearSolves: (dissipationAttempt?.originalLinearSolves ?? trial.linearDiagnostics.solves)
          + (dissipationAttempt?.phases ?? []).reduce((n, p) => n + p.linearSolves, 0),
        ...(dissipationAttempt ? { firstOrderRecovery: structuredClone(dissipationAttempt) } : {}),
        maximumResidual: maximum(trial.residual), families: { ...trial.families },
        gridValid: trial.mesh.quality.valid, ...(trial.lastRejectedStep ? { lastRejectedStep: structuredClone(trial.lastRejectedStep) } : {}),
        ...(trial.progressControl ? { progressControl: trial.progressControl } : {}),
        ...(trial.dissipationEnhancement ? { dissipationEnhancement: trial.dissipationEnhancement } : {}), initialRedistributionSkipped: trial.initialRedistribution.resumed === true,
        transfer: structuredClone(transfer.diagnostics), history });
      if (accepted) {
        const usedIncrement = mach - currentMach;
        const audit = trial.flow.diagnostics.maxMach >= .9 ? auditStreamtubeShocks(trial.flow) : undefined;
        const previousAudit = audit && current.flow.diagnostics.maxMach >= .9 ? auditStreamtubeShocks(current.flow) : undefined;
        const shockMovement = compareShockAudits(previousAudit, audit);
        current = trial; currentMach = mach; currentCheckpoint = structuredClone(trial.checkpoint);
        observe(onCheckpoint, currentCheckpoint, { ...info, kind: 'accepted', families: trial.families,
          reachedTarget: currentMach === targetMach });
        const refine = refinementStop();
        if (refine) return refine;
        if (currentMach === targetMach) return finish(true);
        const growth = coupledMachGrowth({ iterations: attempts.at(-1).iterations, budget: stageMaxIterations,
          backtracks: history.reduce((sum, h) => sum + (h.backtracks ?? 0), 0),
          failedLargerStep: subdivisions > 0, shockAudit: audit, shockMovement });
        attempts.at(-1).nextIncrement = growth;
        increment = Math.sign(distance) * Math.min(maxMachStep, Math.abs(targetMach - currentMach), growth.factor * Math.abs(usedIncrement));
        subdivisions = 0;
        continue;
      }
      const recoveryPlan = wakeRecoveries.length < maxWakeRecoveries && coupledWakeRecoveryPlan(trial,currentCheckpoint,tolerance);
      if (recoveryPlan) {
        const recovery = { attempt:wakeRecoveries.length+1, afterMachAttempt:info.attempt,
          sourceMach:currentMach, failedMach:mach, plan:recoveryPlan, accepted:false };
        wakeRecoveries.push(recovery);
        const recoveryInfo = { ...info, mach:currentMach, fraction:(currentMach-sourceMach)/(targetMach-sourceMach),
          stage:'coupled-wake-grid',
          label:'wake-grid-reinitialization', wakeRecovery:recovery.attempt };
        observe(onStage,recoveryInfo);
        try {
          const recovered = recoverCoupledWakeCorrespondence(currentCheckpoint,{tolerance,
            maxIterations:Math.min(12,stageMaxIterations),maxBacktracks,
            onIteration:entry=>observe(onIteration,{...entry,...recoveryInfo}),
            onMesh:onMesh ? snapshot=>observe(onMesh,{
              mesh:streamtubeMeshSnapshot(snapshot),iteration:snapshot.iteration,
              families:snapshot.coupledFamilies,...recoveryInfo },...(includeFlowState ? [{
                flow:physicalFlowData(snapshot.flow),iteration:snapshot.iteration,coupledFamilies:snapshot.coupledFamilies }] : [])) : undefined,
            onCheckpoint:onCheckpoint ? (value,details)=>observe(onCheckpoint,value,
              {...details,...recoveryInfo,kind:'wake-initializer-iterate'}) : undefined });
          if (observerFailed) throw observerError;
          recovery.diagnostics = recovered.diagnostics;
          const root = recovered.result;
          if (root.converged && root.mesh.quality.valid && root.conditions.mach === currentMach
            && maximum(root.residual) <= tolerance) {
            current = root; currentCheckpoint = structuredClone(root.checkpoint); recovery.accepted = true;
            observe(onCheckpoint,currentCheckpoint,{...recoveryInfo,kind:'accepted',families:root.families,reachedTarget:false});
            // Retry the failed physical condition from the new, fully solved
            // same-Mach grid. The source Mach and continuation budget remain.
            increment = mach-currentMach;
            continue;
          }
          recovery.reason = root.reason;
        } catch(error) {
          if (observerFailed) throw observerError;
          recovery.reason = error.message;
        }
      }
    } catch (error) {
      if (observerFailed) throw observerError;
      attempts.push({ ...info, converged: false, reason: error.message, failureStage: stage,
        iterations: Math.max(0, history.length - 1), history,
        ...(error.code ? { code: error.code } : {}) });
    }
    if (subdivisions >= maxSubdivisions)
      return finish(false, `Automatic coupled Mach continuation exhausted ${maxSubdivisions} subdivisions at Mach ${currentMach}.`);
    increment = .5 * (mach - currentMach); subdivisions++;
  }
  return finish(false, `Automatic coupled Mach continuation reached its ${maxStages}-stage limit at Mach ${currentMach}.`);
}
