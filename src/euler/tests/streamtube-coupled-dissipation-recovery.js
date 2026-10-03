// SPDX-License-Identifier: GPL-2.0-or-later
// Bounded first-order initializer; only a re-solved target-order root may
// enter Mach continuation's accepted checkpoint chain.
import { createCoupledStreamtubeBody } from '../streamtube-coupled.js';
import { solveCoupledStreamtubeIses } from '../streamtube-coupled-ises.js';
import { gridReplayDeparture } from '../streamtube-coupled-grid-levels.js';

const maximum = values => values.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
const controls = cp => Object.fromEntries(['iterationGeometry', 'stepAcceptance', 'stagnationLimiter']
  .map(k => [k, cp.continuation[k]]));

export function coupledFirstOrderRecoveryEligible(trial) {
  const input = trial?.checkpoint?.restart?.input;
  return Boolean(!trial?.converged && trial?.mesh?.quality?.valid && trial?.checkpoint
    && trial.residual?.length && trial.residual.every(Number.isFinite)
    && trial.flow?.diagnostics?.maxMach >= .9
    && ['momentum', 'hybrid'].includes(input?.streamwiseMode) && input.hybrid?.ismom !== 2
    && input.upwind?.mucon >= .5 && input.upwind?.mcrit <= 1);
}

export function changeCoupledDissipationOrder(checkpoint, mucon) {
  if (!Number.isFinite(mucon) || Math.abs(mucon) < .5)
    throw new Error('Dissipation recovery requires a stable first-order MUCON magnitude.');
  const cp = structuredClone(checkpoint), f = cp.restart;
  if (!f?.input?.upwind || Math.abs(f.input.upwind.mucon) !== Math.abs(mucon)
    || !['momentum', 'hybrid'].includes(f.input.streamwiseMode))
    throw new Error('Dissipation recovery may change only the sign of MUCON.');
  const make = () => createCoupledStreamtubeBody(f.input, { ...f.options,
    initialEuler: f.initialEuler, initialBL: f.initialBL });
  const source = make(), before = source.admissibleValue(source.initial, { requireConvex: true });
  if (!before || Object.keys(before.families).some(k => before.families[k] !== cp.families[k]))
    throw new Error('Dissipation recovery source does not replay.');
  f.input.upwind.mucon = mucon;
  const target = make(), after = target.admissibleValue(target.initial, { requireConvex: true });
  if (!after || source.n !== target.n || source.initial.some((v, i) => v !== target.initial[i])
    || !gridReplayDeparture(after.outer.nodes, before.outer.nodes).equivalent
    || !gridReplayDeparture(after.outer.undisplacedNodes, before.outer.undisplacedNodes).equivalent)
    throw new Error('Dissipation order change did not preserve an admissible physical state.');
  cp.families = { ...after.families };
  return cp;
}

export function recoverCoupledDissipationOrder(checkpoint, { tolerance = 1e-10,
  maxIterations = 20, maxBacktracks = 12, dissipationEnhancement = true,
  iterationRecovery = true, onPhase, onIteration, onMesh, onCheckpoint } = {}) {
  if (!Number.isInteger(maxIterations) || maxIterations < 1)
    throw new Error('Dissipation recovery needs a positive iteration budget.');
  const target = structuredClone(checkpoint.restart.input.upwind), mach = checkpoint.restart.input.mach;
  target.mcrit = checkpoint.continuation.dissipationEnhancement?.targetMcrit ?? target.mcrit;
  if (!(target.mucon >= .5)) throw new Error('Dissipation recovery requires second-order target equations.');
  const phases = [];
  let bestFirstOrder, bestFirstOrderResidual = Infinity, initialFirstOrderResidual;
  const run = (resume, phase, budget) => {
    const info = { phase, mucon: resume.restart.input.upwind.mucon, targetMucon: target.mucon,
      provisional: true, targetEquationsRestored: phase === 'second-order-restoration' };
    onPhase?.(structuredClone(info));
    const result = solveCoupledStreamtubeIses(undefined, { resume, ...controls(resume),
      tolerance, maxIterations: budget, maxBacktracks, dissipationEnhancement, iterationRecovery,
      onIteration: h => onIteration?.(h, info), onMesh: s => onMesh?.(s, info),
      onCheckpoint: (cp, details) => {
        // Compare only identical target-MCRIT equations, not temporarily
        // broadened shocks with a different residual definition.
        if (phase === 'first-order-initialization' && cp.restart.input.upwind.mcrit === target.mcrit) {
          const merit = Math.max(...Object.values(cp.families));
          initialFirstOrderResidual ??= merit;
          if (merit < bestFirstOrderResidual) {
            bestFirstOrderResidual = merit; bestFirstOrder = structuredClone(cp);
          }
        }
        onCheckpoint?.(cp, { ...details, dissipationRecovery: info });
      } });
    phases.push({ ...info, converged: result.converged, residual: maximum(result.residual),
      iterations: result.history.length - 1, linearSolves: result.linearDiagnostics.solves, reason: result.reason });
    return result;
  };
  const first = run(changeCoupledDissipationOrder(checkpoint, -target.mucon), 'first-order-initialization', Math.min(8, maxIterations));
  if (!first.checkpoint || !first.mesh.quality.valid || !first.residual.every(Number.isFinite))
    return { accepted: false, diagnostics: { phases, reason: 'First-order initializer has no admissible checkpoint.' } };
  if (!first.converged && !(bestFirstOrderResidual < initialFirstOrderResidual))
    return { accepted: false, diagnostics: { phases, reason: 'First-order iterations did not improve the target-MCRIT residual.' } };
  let restoredCheckpoint;
  try { restoredCheckpoint = changeCoupledDissipationOrder(first.converged ? first.checkpoint : bestFirstOrder, target.mucon); }
  catch (error) {
    return { accepted: false, diagnostics: { phases, reason: `Second-order restoration rejected: ${error.message}` } };
  }
  const restored = run(restoredCheckpoint, 'second-order-restoration', maxIterations);
  const finalUpwind = restored.checkpoint?.restart.input.upwind;
  const accepted = restored.converged && restored.conditions.mach === mach && restored.mesh.quality.valid
    && maximum(restored.residual) <= tolerance && finalUpwind?.mucon === target.mucon
    && finalUpwind?.mcrit === target.mcrit && (!restored.dissipationEnhancement || restored.dissipationEnhancement.finalEquations);
  return { accepted, ...(accepted ? { result: restored } : {}), diagnostics: { phases, accepted,
    targetMucon: target.mucon, targetMcrit: target.mcrit, mach,
    reason: accepted ? 'Converged with requested second-order dissipation.' : restored.reason } };
}
