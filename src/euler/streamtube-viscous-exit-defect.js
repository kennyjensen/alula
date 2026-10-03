// SPDX-License-Identifier: GPL-2.0-or-later
// MSES manual §2.7.2, printed p.31: extrapolate each viscous wake's
// measured momentum defect to freestream pressure, then sum the defects.
// Gas provenance is explicit; this helper does not infer or mix wake gas.

const positive = value => Number.isFinite(value) && value > 0;
const gasSources = ['historical-common-isentrope', 'physical-euler', 'measured'];
function logRatio(a, b) {
  const departure = (a - b) / b;
  return Math.abs(departure) < .5 ? Math.log1p(departure) : Math.log(a) - Math.log(b);
}
function sum(values) {
  let result = 0, correction = 0;
  for (const value of values) {
    const adjusted = value - correction, next = result + adjusted;
    correction = (next - result) - adjusted; result = next;
  }
  return result;
}
function gas(state, gamma, label) {
  const { density, speed, pressure } = state ?? {};
  if (![density, speed, pressure].every(positive))
    throw new Error(`${label} requires positive finite density, speed and static pressure.`);
  const dynamicPressureTwice = density * speed * speed;
  const machSquared = dynamicPressureTwice / (gamma * pressure);
  if (![dynamicPressureTwice, machSquared].every(positive))
    throw new Error(`${label} has nonfinite or nonpositive derived momentum flux or Mach squared.`);
  return { density, speed, pressure, machSquared, dynamicPressureTwice };
}

/**
 * Use consistent physical units. theta, deltaStar and referenceChord have
 * length units; defects and forces are per unit span. gasSource describes
 * every supplied wake edge gas, whose consistency with the thickness
 * convention is the caller's responsibility. No Euler wave defect is added.
 */
export function streamtubeViscousExitDefect({ wakes, freestream, gamma = 1.4,
  referenceChord = 1, gasSource } = {}) {
  if (!Number.isFinite(gamma) || gamma <= 1 || !positive(referenceChord))
    throw new Error('Viscous exit defect requires finite gamma > 1 and a positive finite reference chord.');
  if (!gasSources.includes(gasSource))
    throw new Error('Viscous exit defect requires an explicit gasSource: historical-common-isentrope, physical-euler or measured.');
  if (!Array.isArray(wakes) || !wakes.length)
    throw new Error('Viscous exit defect requires a nonempty array of wake exit states.');
  const free = gas(freestream, gamma, 'Freestream');
  const coefficientScale = .5 * free.dynamicPressureTwice * referenceChord;
  const farWakeShape = 1 + (gamma - 1) * free.machSquared;
  if (![coefficientScale, farWakeShape].every(positive))
    throw new Error('Viscous exit defect has nonpositive or nonfinite normalization or asymptotic wake shape.');
  const results = wakes.map((wake, index) => {
    const exit = gas(wake, gamma, `Wake ${index}`), { theta, deltaStar } = wake;
    if (![theta, deltaStar].every(positive))
      throw new Error(`Wake ${index} requires positive finite momentum and displacement thickness.`);
    const shape = deltaStar / theta, averageShape = .5 * (shape + farWakeShape);
    const pressureExponent = averageShape / (gamma * free.machSquared);
    const logPressureRecovery = logRatio(free.pressure, exit.pressure);
    const logCorrection = pressureExponent * logPressureRecovery;
    const correctionFactor = Math.exp(logCorrection);
    const exitMomentumDeficit = exit.dynamicPressureTwice * theta;
    const momentumDeficit = exitMomentumDeficit * correctionFactor;
    // Keep the small signed pressure correction visible even when adding
    // it to the finite exit defect would round back to the same number.
    const pressureCorrection = exitMomentumDeficit * Math.expm1(logCorrection);
    const pressureRatio = exit.pressure / free.pressure;
    if (![shape, averageShape, pressureExponent, correctionFactor, exitMomentumDeficit, momentumDeficit, pressureRatio].every(positive)
      || ![logPressureRecovery, logCorrection, pressureCorrection].every(Number.isFinite))
      throw new Error(`Wake ${index} has nonpositive or nonfinite extrapolated defect or pressure correction.`);
    return { index, gasSource, exit: { ...exit, theta, deltaStar, shape }, averageShape,
      pressureExponent, pressureRatio, logPressureRecovery, logCorrection, correctionFactor,
      exitMomentumDeficit, pressureCorrection, momentumDeficit,
      viscousDragCoefficient: momentumDeficit / coefficientScale };
  });
  const momentumDeficit = sum(results.map(w => w.momentumDeficit));
  const exitMomentumDeficit = sum(results.map(w => w.exitMomentumDeficit));
  const pressureCorrection = sum(results.map(w => w.pressureCorrection));
  const viscousDragCoefficient = momentumDeficit / coefficientScale;
  if (![momentumDeficit, exitMomentumDeficit, viscousDragCoefficient].every(positive)
    || !Number.isFinite(pressureCorrection) || results.some(w => !positive(w.viscousDragCoefficient)))
    throw new Error('Viscous exit defect has nonfinite or nonpositive integrated force or coefficient.');
  return { momentumDeficit, exitMomentumDeficit, pressureCorrection, viscousDragCoefficient,
    wakes: results, freestream: free, gamma, referenceChord, coefficientScale, farWakeShape, gasSource,
    maxRelativeExitPressureDeparture: Math.max(...results.map(w => Math.abs(w.pressureRatio - 1))),
    dragKind: 'viscous-wake-exit', includesEulerWaveDefect: false, includesSolidPressureIntegral: false,
    physicalAcceptance: false,
    method: 'MSES manual §2.7.2, p.31: extrapolate rho_e*u_e^2*theta by (p_infinity/p_exit)^[H_average/(gamma*Mach_infinity^2)], with H_average=(deltaStar/theta+1+(gamma-1)*Mach_infinity^2)/2. Sum all supplied element wakes. No gas replacement, clipping, wave-loss integration or total-drag claim.' };
}
