// SPDX-License-Identifier: GPL-2.0-or-later
// Automatic startup and bounded pressure continuation for the verification
// channel. No analytic nozzle shape, sonic branch or shock position is used.
// Restarting from the last converged point and halving a failed parameter
// increment follows the MSES manual's MPOLAR procedure (§2.5, printed p.28).
import { createUpwindStreamtubeChannel } from './streamtube-upwind-channel.js';
import { initializeUpwindStreamtubeChannel } from './streamtube-channel-startup.js';
import { solveStreamtubeChannel } from './streamtube-channel.js';

export function solveUpwindStreamtubeChannelAutomatic(input, { tolerance = 1e-10,
  directMaxIterations = 12, stageMaxIterations = 20, maxSubdivisions = 5,
  initialFractionStep = .25, maxStages = 64, onIteration, onState, onStage, ...solverControls } = {}) {
  if (![directMaxIterations, stageMaxIterations, maxSubdivisions].every(v => Number.isInteger(v) && v >= 0)
    || !Number.isInteger(maxStages) || maxStages < 2 || maxSubdivisions > 20
    || !(initialFractionStep > 0 && initialFractionStep <= 1)) throw new Error('Invalid automatic channel controls.');
  // Case data must keep identical meanings across continuation stages, even
  // if a UI observer edits its own copy while receiving progress. Ignore
  // unrelated caller metadata/functions rather than cloning an entire oracle.
  const inputKeys = ['x', 'lower', 'upper', 'massFlows', 'stagnationEnthalpy', 'stagnationDensity',
    'gamma', 'referenceDensity', 'referencePressure', 'pressureCorrectionFactor', 'inletSlopes',
    'outletSlopes', 'outletPressure', 'streamwiseMode', 'upwind', 'hybrid'];
  input = Object.fromEntries(inputKeys.filter(key => Object.hasOwn(input, key)).map(key => [key, structuredClone(input[key])]));
  const targetSystem = createUpwindStreamtubeChannel(input), { conditions, nt } = targetSystem;
  const { gamma, stagnationEnthalpy: h0, stagnationDensity: rhoTotal } = conditions;
  const targetPressure = conditions.outletPressure.slice();
  const totalPressure = h0.map((h, j) => (gamma - 1) / gamma * rhoTotal[j] * h);
  if (targetPressure.some((p, j) => p >= totalPressure[j]))
    throw new Error('Positive reservoir-to-outlet flow requires outlet static pressure below inlet total pressure.');
  const startup = initializeUpwindStreamtubeChannel(targetSystem), attempts = [];
  const run = (system, initial, label, fraction, maxIterations) => {
    const pressures = system.conditions.outletPressure.slice();
    onStage?.({ label, fraction, outletPressure: pressures.slice() });
    const result = solveStreamtubeChannel(system, { ...solverControls, initial, tolerance, maxIterations,
      onIteration: h => onIteration?.({ ...h, stage: label, continuationFraction: fraction, outletPressure: pressures.slice() }),
      onState: snapshot => onState?.({ ...snapshot, stage: label, continuationFraction: fraction, outletPressure: pressures.slice() }) });
    attempts.push({ label, fraction, outletPressure: pressures, converged: result.converged,
      reason: result.reason, iterations: result.history.length - 1,
      maximumResidual: result.residual.reduce((peak, value) => Math.max(peak, Math.abs(value)), -Infinity), history: result.history.map(h => ({ ...h })) });
    return result;
  };
  const finish = (system, result, reachedTarget, reason) => ({ ...result, system,
    converged: reachedTarget && result.converged, ...(reason ? { reason } : {}),
    automaticInitialization: startup.diagnostics,
    continuation: { method: 'outlet-pressure', reachedTarget, targetPressure, attempts,
      preservesPhysicalDensityAndMass: true, analyticShockInformationUsed: false } });
  const direct = run(targetSystem, startup.initial, 'target-direct', 1, directMaxIterations);
  if (direct.converged) return finish(targetSystem, direct, true);

  // Build an easy subcritical pressure from one quarter of each tube's
  // geometry-derived sonic-capacity seed. This fraction only selects an
  // initial problem; the target pressure and all target equations stay fixed.
  // Constant reservoir density is retained in the starting iterate.
  const lowMassFraction = .25, lowInitial = startup.initial.slice();
  const startingPressure = Array.from({ length: nt }, (_, j) => {
    lowInitial[targetSystem.massIndex(j)] += Math.log(lowMassFraction);
    const end = startup.flow.sections.at(-1)[j];
    const area = startup.flow.massFlows[j] / (end.rho * end.q);
    const massFlux = lowMassFraction * startup.flow.massFlows[j] / area;
    const sonicSpeed = Math.sqrt(2 * (gamma - 1) / (gamma + 1) * h0[j]);
    let lo = 0, hi = sonicSpeed;
    for (let k = 0; k < 60; k++) {
      const q = .5 * (lo + hi), rho = rhoTotal[j] * (1 - .5 * q * q / h0[j]) ** (1 / (gamma - 1));
      if (rho * q < massFlux) lo = q; else hi = q;
    }
    const q = .5 * (lo + hi), h = h0[j] - .5 * q * q;
    const p = (gamma - 1) / gamma * rhoTotal[j] * (h / h0[j]) ** (1 / (gamma - 1)) * h;
    return Math.max(targetPressure[j], p);
  });
  let currentSystem = createUpwindStreamtubeChannel({ ...input, outletPressure: startingPressure });
  let current = run(currentSystem, lowInitial, 'subcritical-start', 0, stageMaxIterations);
  if (!current.converged) return finish(currentSystem, current, false, `Automatic subcritical startup failed: ${current.reason}`);
  let fraction = 0, step = initialFractionStep, subdivisions = 0;
  while (fraction < 1) {
    if (attempts.length >= maxStages) return finish(currentSystem, current, false, 'Automatic pressure-continuation stage limit.');
    const next = Math.min(1, fraction + step);
    const pressure = startingPressure.map((p, j) => next === 1 ? targetPressure[j] : p + next * (targetPressure[j] - p));
    const nextSystem = createUpwindStreamtubeChannel({ ...input, outletPressure: pressure });
    // Encoded densities, physical node offsets and reference tube masses
    // have identical meanings throughout the pressure continuation.
    const attempt = run(nextSystem, current.x.slice(), 'pressure-continuation', next, stageMaxIterations);
    if (attempt.converged) {
      current = attempt; currentSystem = nextSystem; fraction = next; subdivisions = 0;
      step = Math.min(initialFractionStep, step * 1.5);
    } else {
      if (subdivisions >= maxSubdivisions || step <= initialFractionStep / 2 ** maxSubdivisions) return finish(currentSystem, current, false,
        `Automatic pressure continuation stopped at fraction ${fraction}: ${attempt.reason}`);
      step *= .5; subdivisions++;
    }
  }
  return finish(currentSystem, current, true);
}
