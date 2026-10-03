// SPDX-License-Identifier: GPL-2.0-or-later
// Presentation only: retain explicitly reported solve methods between updates.
const finite = Number.isFinite;
const labels = {
  euler: 'Euler precursor', initialization: 'Boundary-layer initialization',
  startup: 'Ordinary Euler/BL startup', lowerNcrit: 'Lower-Ncrit recovery',
  logarithmic: 'Logarithmic shear recovery', transition: 'Transition refinement',
  hkProjection: 'Hk thickness-projection recovery',
  grid: 'Grid refinement', ncrit: 'Ncrit continuation', mach: 'Mach continuation', alpha: 'Alpha continuation', operatingPoint: 'Mach–alpha continuation',
  wake: 'Wake-grid rebuild', cold: 'Cold Mach recovery', saved: 'Saved-flow initialization',
  firstOrder: 'Temporary first-order dissipation', secondOrder: 'Restore second-order dissipation',
  retained: 'Retained flow',
};

function method(event, state, previous) {
  const stage = state.stage;
  if (event.retained) return 'retained';
  if (event.hkProjectionRecovery) return 'hkProjection';
  if (event.dissipationRecovery?.phase === 'first-order-initialization') return 'firstOrder';
  if (event.dissipationRecovery?.phase === 'second-order-restoration') return 'secondOrder';
  if (stage === 'euler') return 'euler';
  if (stage === 'boundary-layer-initialization') return 'initialization';
  if (event.shearRecovery?.shearCoordinate === 'logarithmic') return 'logarithmic';
  if (stage === 'transition-refinement' || event.transitionRecovery) return 'transition';
  if (stage === 'coupled-refinement' || stage === 'coupled-grid-refinement' || event.refinement) return 'grid';
  if (stage?.includes('ncrit') || event.ncritStartup
    || previous && finite(event.actualNcrit) && event.actualNcrit !== previous.actualNcrit && stage === 'coupled') return 'ncrit';
  if (stage === 'coupled-operating-point') return 'operatingPoint';
  if (stage === 'coupled-alpha') return 'alpha';
  if (stage === 'coupled-mach') return 'mach';
  if (stage === 'coupled-wake-grid') return 'wake';
  if (stage === 'coupled-cold-recovery' || stage === 'euler-cold-startup') return 'cold';
  if (stage === 'hybrid-certification') return 'saved';
  if (stage === 'coupled-coarse-initialization') return 'initialization';
  if (stage === 'coupled') {
    if (previous && state.startupAttempt === previous.startupAttempt
      && ['logarithmic', 'hkProjection', 'transition', 'grid', 'ncrit', 'mach'].includes(previous.kind)) return previous.kind;
    return state.actualNcrit < state.targetNcrit ? 'lowerNcrit' : 'startup';
  }
  return previous?.kind;
}

function label(state) {
  const direct = state.startupStrategy === 'direct-requested-conditions';
  const parts = [direct ? (state.kind === 'euler' ? 'Step 1/2: inviscid at requested conditions'
    : state.kind === 'initialization' ? 'Initializing boundary layers at requested conditions'
    : 'Step 2/2: viscous at requested conditions') : labels[state.kind]];
  if (direct && state.kind === 'hkProjection') parts.push(labels.hkProjection);
  if (direct && state.kind === 'initialization' && state.eulerPreparation?.converged === false)
    parts.push('using partial inviscid state after iteration limit');
  if (state.alphaSeed) parts.push('lower-incidence seed for alpha continuation');
  if (state.startupStrategy === 'target-inviscid-then-viscous')
    parts.push(state.kind === 'euler' ? 'target-condition inviscid startup' : 'viscous solve at requested conditions');
  if (state.startupStrategy === 'low-mach-fallback') parts.push('low-Mach fallback');
  if (state.kind === 'initialization' && state.boundaryLayerInitialization) {
    const seed = state.boundaryLayerInitialization;
    if (seed.wakeInitialization === 'iset-linear-shape' || seed.method?.includes('ISET')) parts.push('ISET wake guess');
    else if (seed.method) parts.push(seed.method === 'uniform initial thickness backtracking' ? 'initial BL guess' : seed.method);
    if (finite(seed.thicknessFactor) && seed.thicknessFactor !== 1) parts.push(`thickness × ${seed.thicknessFactor}`);
  }
  if (state.startupAttempt > 0) parts.push(`attempt ${state.startupAttempt}`);
  if (state.kind !== 'euler') {
    if (state.shearCoordinate && state.kind !== 'logarithmic') parts.push(`${state.shearCoordinate} shear`);
    if (state.hkFloorLinearization) parts.push(`${state.hkFloorLinearization} Hk derivative`);
    if (finite(state.actualNcrit)) parts.push(`Ncrit ${state.actualNcrit}`
      + (state.actualNcrit !== state.targetNcrit && finite(state.targetNcrit) ? ` → ${state.targetNcrit}` : ''));
  }
  if (finite(state.actualAlpha) && finite(state.targetAlpha)) parts.push(`α ${state.actualAlpha.toFixed(2)}°`
    + (state.actualAlpha !== state.targetAlpha ? ` → ${state.targetAlpha.toFixed(2)}°` : ''));
  if (finite(state.mach)) parts.push(`Mach ${state.mach.toFixed(3)}`
    + (finite(state.targetMach) && state.mach !== state.targetMach ? ` → ${state.targetMach.toFixed(3)}` : ''));
  if (state.dissipation) parts.push(`MSES dissipation · MCRIT ${state.dissipation.mcrit.toFixed(3)} → ${state.dissipation.targetMcrit.toFixed(3)}`);
  if (state.dissipationDamping) parts.push('damping suppresses shock broadening');
  if (state.shockAudit?.weakMomentumCandidates) parts.push('shock momentum check flagged');
  if (state.shockAudit?.maximumSkewDegrees > 45) parts.push('high shock-region grid shear');
  if (state.stepMethod) parts.push(state.stepMethod.replaceAll('-', ' '));
  if (state.gridStrategy === 'shock-triggered-refinement') parts.push('shock-triggered refinement');
  if (state.gridStrategy === 'refine-retained-operating-point') parts.push('refine retained operating point');
  if (state.gridStrategy === 'target-mach-before-refinement') parts.push('coarse-to-fine continuation');
  if (state.gridLevel !== undefined) parts.push(`grid ${state.gridLevel}`);
  if (state.resolutionRecovery) {
    const r = state.resolutionRecovery;
    parts.push(r.finalNominalGridIntervals > r.requestedGridIntervals
      ? `adaptive refinement ${r.requestedGridIntervals} → ${r.finalNominalGridIntervals}`
      : 'coarse-to-fine grid recovery');
  }
  if (state.progress?.action === 'damped-recovery') parts.push('damped Newton recovery',
    ...(finite(state.progress.frozenMcrit) ? [`MCRIT held at ${state.progress.frozenMcrit.toFixed(3)}`] : []));
  if (state.progress?.cause) parts.push(state.progress.cause.replaceAll('-', ' '));
  if (state.progress?.recoveryOutcome?.kind === 'ordinary-policy-restored') parts.push('ordinary Newton resumed after damping trial');
  if (state.progress?.action === 'stop-stage') parts.push('stalled stage stopped');
  if (state.progress?.extendedBudget) parts.push('extra iteration: rapid residual decrease');
  if (state.residualContext?.equationsChanged) parts.push('dissipation equations changed');
  if (state.transitionChanged) parts.push('transition changed');
  if (state.wakeCorrection?.phase === 'start') parts.push('checking wake-grid correction');
  if (state.wakeCorrection?.phase === 'trial') parts.push(`checking wake-grid correction · trial ${state.wakeCorrection.iteration}`);
  if (state.wakeCorrection?.phase === 'accepted') parts.push('wake-grid correction accepted');
  if (state.wakeCorrection?.phase === 'rejected') parts.push('wake-grid correction rejected · retained previous flow');
  return parts.filter(Boolean).join(' · ');
}

export function createQuadSolveProgress(request = {}, parent) {
  const options = parent?.checkpoint?.restart?.options;
  let state = {
    targetNcrit: request.ncrit, actualNcrit: options?.ncrit ?? request.ncrit,
    targetMach: request.mach, mach: parent?.checkpoint?.restart?.input?.mach ?? request.mach,
    ...(parent?.checkpoint ? { shearCoordinate: parent.checkpoint.continuation?.shearCoordinate ?? 'linear',
      hkFloorLinearization: options?.hkFloorLinearization ?? 'exact' } : {}),
  };
  const history = [];
  return {
    update(event, { stageChange = false } = {}) {
      const previous = state;
      const fresh = stageChange && ['euler', 'boundary-layer-initialization', 'euler-cold-startup', 'coupled-cold-recovery'].includes(event.stage);
      if (stageChange || event.iteration !== undefined) state.wakeCorrection = event.wakeCorrection;
      state = { ...state, ...(fresh ? { shearCoordinate: undefined, hkFloorLinearization: undefined, kind: undefined,
        startupStrategy: undefined, boundaryLayerInitialization: undefined, dissipation: undefined, gridStrategy: undefined, resolutionRecovery: undefined } : {}) };
      if (stageChange) { delete state.alphaSeed; delete state.stepMethod; delete state.shockAudit; delete state.dissipationDamping; delete state.gridLevel; delete state.progress; delete state.residualContext; delete state.transitionChanged; }
      if (event.iteration !== undefined) {
        state.shockAudit = event.shockAudit;
        state.dissipationDamping = event.dissipation?.damping?.dampingSuppressesBroadening === true;
        state.progress = event.progress; state.residualContext = event.residualContext; state.transitionChanged = event.activeChange === true;
      }
      if (event.retained) delete state.resolutionRecovery;
      if (event.resolutionRecovery && typeof event.resolutionRecovery === 'object')
        state.resolutionRecovery = event.resolutionRecovery;
      for (const key of ['stage', 'startupAttempt', 'actualNcrit', 'targetNcrit', 'mach', 'targetMach', 'gridLevel', 'gridStrategy', 'startupStrategy', 'actualAlpha', 'targetAlpha', 'stepMethod', 'alphaSeed', 'eulerPreparation'])
        if (event[key] !== undefined) state[key] = event[key];
      if (event.dissipation && finite(event.dissipation.mcrit) && finite(event.dissipation.targetMcrit))
        state.dissipation = { mcrit: event.dissipation.mcrit, targetMcrit: event.dissipation.targetMcrit };
      const coordinate = event.shearCoordinate ?? event.shearRecovery?.shearCoordinate;
      if (['linear', 'logarithmic'].includes(coordinate)) state.shearCoordinate = coordinate;
      if (['exact', 'native'].includes(event.hkFloorLinearization)) state.hkFloorLinearization = event.hkFloorLinearization;
      if (event.boundaryLayerInitialization) {
        const { method, wakeInitialization, thicknessFactor } = event.boundaryLayerInitialization;
        state.boundaryLayerInitialization = { method, wakeInitialization, thicknessFactor };
      }
      state.kind = method(event, state, fresh ? null : previous);
      if (!state.kind) return { current: '', history: [] };
      const current = label(state);
      const key = JSON.stringify([state.kind, state.startupAttempt, state.actualNcrit, state.mach, state.gridLevel]);
      if (history.at(-1)?.key === key) history[history.length - 1] = { key, label: current };
      else history.push({ key, label: current });
      return { current, history: history.map(entry => entry.label) };
    },
  };
}

export function quadStartingGridLabel(mesh, phase) {
  if (mesh?.iteration?.stage?.startsWith('coupled')) return 'Starting Euler/BL grid · before the first Newton update';
  if (phase === 'smoothing') return 'SLOR initialization · Euler has not started';
  return mesh?.initialization?.flowSolved ? 'Euler converged on the starting grid · 0 updates'
    : 'Starting grid · no accepted Euler update yet';
}
