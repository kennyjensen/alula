// SPDX-License-Identifier: GPL-2.0-or-later
// Every visible GUI preset, including the three reported cold failures.
import { buildReliabilityCase } from './solver-reliability-cases.js';
import { assessSolverResult } from './solver-reliability-acceptance.js';

export const guiDefaultCoupledCases = Object.freeze([
  Object.freeze({ id: 'flap', preset: 'flap', label: 'Main element + flap' }),
  Object.freeze({ id: 'single', preset: 'single', label: 'Single NACA 0012' }),
  Object.freeze({ id: 'three', preset: 'three', label: 'Slat + main + flap' }),
  Object.freeze({ id: 'rae2822', preset: 'rae2822-mses', label: 'RAE 2822' }),
  Object.freeze({ id: 'nlr7301', preset: 'nlr7301', label: 'NLR 7301 + flap' }),
  Object.freeze({ id: '30p30n', preset: '30p30n', label: '30P30N' }),
]);
export const guiDefaultCoupledControls = Object.freeze({ maxIterations: 40, eulerMaxIterations: 20, tolerance: 1e-10, maxStartupAttempts: 2 });

export function buildGuiDefaultCoupledCase(id) {
  const spec = guiDefaultCoupledCases.find(c => c.id === id);
  if (!spec) throw new Error(`Unknown GUI default case: ${id}`);
  const { caseData } = buildReliabilityCase({ preset: spec.preset, mode: 'streamtube-bl', changes: { panels: 160, eulerIsmom: 4 } });
  return caseData;
}

export function guiDefaultGeometrySummary(input) {
  return input.elements.map(e => ({ name: e.name, points: e.points.length,
    distinctPoints: e.points.length - Number(e.points[0].x === e.points.at(-1).x && e.points[0].y === e.points.at(-1).y),
    trailingEdge: e.trailingEdge ?? { kind: 'sharp' },
    sourcePoints: e.sourcePoints?.length ?? null }));
}

export function assessGuiDefaultCoupledResult(input, raw, displayed) {
  const acceptance = assessSolverResult(input, displayed);
  const actualNcrit = raw.checkpoint?.restart?.options?.ncrit ?? raw.conditions?.ncrit ?? raw.actualNcrit;
  const expectedKernelReynolds = input.reynolds * raw.solverLength / input.referenceChord;
  const normalizedReynolds = [raw.solverLength, expectedKernelReynolds, raw.kernelReynolds].every(v => Number.isFinite(v) && v > 0)
    && Math.abs(raw.kernelReynolds / expectedKernelReynolds - 1) <= 8 * Number.EPSILON
    && raw.conditions?.reynolds === raw.kernelReynolds
    && raw.checkpoint?.restart?.options?.reynolds === raw.kernelReynolds;
  Object.assign(acceptance.checks, {
    requestedReynolds: raw.referenceReynolds === input.reynolds && raw.referenceChord === input.referenceChord,
    reynoldsNormalization: normalizedReynolds,
    requestedNcrit: actualNcrit === input.ncrit && raw.ncritContinuation?.reachedTarget !== false,
    selectedEquations: raw.solverInput?.hybrid?.ismom === input.eulerIsmom,
    requestedGrid: raw.gridSequence?.reachedTarget !== false,
    transitionMode: raw.conditions?.transitionMode === 'automatic',
    materialTrips: JSON.stringify(raw.materialTrips) === JSON.stringify(input.materialTrips),
  });
  acceptance.failures = Object.keys(acceptance.checks).filter(k => !acceptance.checks[k]);
  acceptance.passed = !acceptance.failures.length;
  return { ...acceptance, actualNcrit };
}

// Compact an existing observation only. Never evaluate the system or change a
// numerical step to collect diagnostics. Tiny accepted steps remain failures
// when the final requested equations have not converged.
export function compactGuiDefaultIteration(h) {
  const keys = ['stage', 'coarseStage', 'startupAttempt', 'attempt', 'gridLevel', 'mach', 'actualNcrit', 'targetNcrit', 'ncritContinuation',
    'iteration', 'residual', 'euler', 'boundaryLayer', 'edgeMatching', 'step', 'stepKind', 'backtracks',
    'limiter', 'viscousLimiter', 'viscousStep', 'stagnation', 'undampedUpdate', 'activeChange', 'changes', 'rejections', 'eventProfile'];
  return { ...Object.fromEntries(keys.filter(k => h[k] !== undefined).map(k => [k, h[k]])),
    ...(h.maintenance ? { maintenance: { triggeredBodies: h.maintenance.triggeredBodies,
      maxDisplacement: h.maintenance.maxDisplacement, accepted: h.maintenance.accepted,
      geometryRedistribution: h.maintenance.geometryRedistribution,
      passages: h.maintenance.passages?.map(p => ({ group: p.group, correctionScale: p.correctionScale,
        maxDisplacement: p.maxDisplacement, referenceBank: p.referenceBank, fixedBanks: p.fixedBanks })) } } : {}),
    ...(h.projection ? { projection: { method: h.projection.method, active: h.projection.active,
      displacementChanges: h.projection.displacementChanges, auxiliaryChanges: h.projection.auxiliaryChanges } } : {}),
    tinyAcceptedStep: Number.isFinite(h.step) && h.step > 0 && h.step < 1e-10 };
}
