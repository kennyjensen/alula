// SPDX-License-Identifier: GPL-2.0-or-later
// MSES manual §2.7.2: isentropically extrapolate each physical Euler exit
// tube to p_infinity, then integrate (U_infinity - q_recovered) dm.
// This preserves each tube's own entropy. It is not total viscous drag
// or a finite-plane vector momentum / pressure-thrust integration.

const consistencyTolerance = 256 * Number.EPSILON;
const positive = value => Number.isFinite(value) && value > 0;
const consistent = (a, b) => Math.abs(a - b) <= consistencyTolerance * Math.max(Math.abs(a), Math.abs(b));

function logRatio(a, b) {
  const relative = (a - b) / b;
  return Math.abs(relative) < .5 ? Math.log1p(relative) : Math.log(a) - Math.log(b);
}

function direction(value, label) {
  if (value === undefined) return null;
  if (!Number.isFinite(value?.x) || !Number.isFinite(value?.y)) throw new Error(`${label} direction must be finite.`);
  const length = Math.hypot(value.x, value.y);
  if (!positive(length)) throw new Error(`${label} direction must have positive finite length.`);
  return { x: value.x / length, y: value.y / length };
}

function physicalGas(state, gamma, label) {
  const { pressure, density, speed, stagnationEnthalpy } = state ?? {};
  if (![pressure, density, speed, stagnationEnthalpy].every(positive))
    throw new Error(`${label} requires positive finite pressure, density, speed and stagnation enthalpy.`);
  const enthalpy = gamma / (gamma - 1) * (pressure / density), kineticEnergy = .5 * speed * speed;
  if (!positive(enthalpy) || !positive(stagnationEnthalpy - kineticEnergy)
    || !consistent(enthalpy + kineticEnergy, stagnationEnthalpy))
    throw new Error(`${label} must have positive static enthalpy and consistent ideal-gas total energy.`);
  const logTemperatureRatio = logRatio(stagnationEnthalpy, enthalpy);
  const stagnationDensity = Math.exp(Math.log(density) + logTemperatureRatio / (gamma - 1));
  const stagnationPressure = Math.exp(Math.log(pressure) + gamma / (gamma - 1) * logTemperatureRatio);
  const machSquared = speed * speed / ((gamma - 1) * enthalpy);
  if (![stagnationDensity, stagnationPressure, machSquared].every(positive))
    throw new Error(`${label} has nonfinite or nonpositive derived gas quantities.`);
  return { pressure, density, speed, enthalpy, stagnationEnthalpy, machSquared, stagnationDensity, stagnationPressure };
}

function sum(values) {
  // Keep small signed entropy defects visible when summing many tubes.
  let total = 0, correction = 0;
  for (const value of values) {
    const adjusted = value - correction, next = total + adjusted;
    correction = (next - total) - adjusted; total = next;
  }
  return total;
}

/**
 * All dimensional inputs must use consistent units. massFlow is physical
 * mass flow per unit span, including the tube's width (not mass flux).
 * Optional section/free-stream directions only measure the finite-exit
 * angle; the source formula assumes eventual alignment with the free stream.
 */
export function streamtubeExitDefect({ sections, freestream, gamma = 1.4, referenceChord = 1 } = {}) {
  if (!Number.isFinite(gamma) || gamma <= 1 || !positive(referenceChord))
    throw new Error('Exit defect requires finite gamma > 1 and a positive finite reference chord.');
  if (!Array.isArray(sections) || !sections.length)
    throw new Error('Exit defect requires a nonempty array of physical inviscid exit sections.');
  const free = physicalGas(freestream, gamma, 'Freestream');
  const freeDirection = direction(freestream.direction, 'Freestream');
  const coefficientScale = .5 * free.density * free.speed * free.speed * referenceChord;
  if (!positive(coefficientScale)) throw new Error('Exit defect has a nonpositive or nonfinite force normalization.');
  const exponent = (gamma - 1) / gamma;
  let measuredExitDirections = 0, maxAbsoluteExitAngleRadians = null, maxRelativeExitPressureDeparture = 0;
  const results = sections.map((section, index) => {
    const exit = physicalGas(section, gamma, `Exit section ${index}`), massFlow = section.massFlow;
    if (!positive(massFlow)) throw new Error(`Exit section ${index} requires positive finite physical mass flow.`);
    if (!consistent(exit.stagnationEnthalpy, free.stagnationEnthalpy))
      throw new Error('Exit defect requires a common adiabatic stagnation enthalpy; different-h0 flow is not modeled.');
    const logPressureRecovery = logRatio(free.pressure, exit.pressure);
    const enthalpyChange = exit.enthalpy * Math.expm1(exponent * logPressureRecovery);
    const recoveredEnthalpy = exit.enthalpy * Math.exp(exponent * logPressureRecovery);
    // Equivalent to 2*(h0 - h_recovered), with exact q recovery if p=p_inf
    // and improved accuracy when the pressure correction is small.
    const recoveredSpeedSquared = exit.speed * exit.speed - 2 * enthalpyChange;
    if (!positive(recoveredEnthalpy) || !(recoveredEnthalpy < exit.stagnationEnthalpy)
      || !positive(recoveredSpeedSquared))
      throw new Error(`Exit section ${index} cannot recover positive finite kinetic energy at freestream pressure.`);
    const recoveredSpeed = Math.sqrt(recoveredSpeedSquared);
    const recoveredDensity = exit.density * Math.exp(logPressureRecovery / gamma);
    const recoveredArea = massFlow / (recoveredDensity * recoveredSpeed);
    const exitArea = massFlow / (exit.density * exit.speed);
    const recoveredMachSquared = recoveredSpeedSquared / ((gamma - 1) * recoveredEnthalpy);
    if (![recoveredDensity, recoveredArea, exitArea, recoveredMachSquared].every(positive))
      throw new Error(`Exit section ${index} has nonpositive or nonfinite recovered mass capacity or gas.`);
    // Rationalize U_inf - q_recovered; never clip negative defects.
    const speedDeficit = ((free.speed - exit.speed) * (free.speed + exit.speed) + 2 * enthalpyChange)
      / (free.speed + recoveredSpeed);
    const momentumDeficit = massFlow * speedDeficit;
    const entropyOverR = (logRatio(exit.pressure, free.pressure) - gamma * logRatio(exit.density, free.density)) / (gamma - 1);
    const stagnationPressureRatio = Math.exp(-entropyOverR);
    const pressureRatio = exit.pressure / free.pressure;
    if (![speedDeficit, momentumDeficit, entropyOverR].every(Number.isFinite)
      || ![stagnationPressureRatio, pressureRatio].every(positive))
      throw new Error(`Exit section ${index} has a nonfinite defect or entropy diagnostic.`);
    const exitDirection = direction(section.direction, `Exit section ${index}`);
    let exitAngleRadians = null;
    if (exitDirection) {
      if (!freeDirection) throw new Error('Measured exit angles require an explicit freestream direction.');
      exitAngleRadians = Math.atan2(freeDirection.x * exitDirection.y - freeDirection.y * exitDirection.x,
        freeDirection.x * exitDirection.x + freeDirection.y * exitDirection.y);
      measuredExitDirections++;
      maxAbsoluteExitAngleRadians = Math.max(maxAbsoluteExitAngleRadians ?? 0, Math.abs(exitAngleRadians));
    }
    maxRelativeExitPressureDeparture = Math.max(maxRelativeExitPressureDeparture, Math.abs(pressureRatio - 1));
    return { index, massFlow, momentumDeficit, speedDeficit, entropyOverR, stagnationPressureRatio,
      stagnationDensity: exit.stagnationDensity, stagnationPressure: exit.stagnationPressure,
      pressureRatio, exitAngleRadians, exit: { ...exit, area: exitArea },
      recovered: { pressure: free.pressure, density: recoveredDensity, speed: recoveredSpeed,
        enthalpy: recoveredEnthalpy, machSquared: recoveredMachSquared, area: recoveredArea } };
  });
  const momentumDeficit = sum(results.map(s => s.momentumDeficit)), massFlow = sum(results.map(s => s.massFlow));
  const eulerWaveDragCoefficient = momentumDeficit / coefficientScale;
  if (!positive(massFlow) || ![momentumDeficit, eulerWaveDragCoefficient].every(Number.isFinite))
    throw new Error('Exit defect has a nonfinite integrated mass flow or force.');
  return { momentumDeficit, eulerWaveDragCoefficient, massFlow, sections: results,
    freestream: free, gamma, referenceChord, coefficientScale, measuredExitDirections,
    maxAbsoluteExitAngleRadians, maxRelativeExitPressureDeparture,
    dragKind: 'euler-exit-entropy', includesViscousWakeDefect: false, includesFinitePlanePressureThrust: false,
    physicalAcceptance: false,
    method: 'MSES manual §2.7.2: preserve each physical Euler exit tube entropy and common stagnation enthalpy, extrapolate to freestream static pressure and direction, then sum physical mass times speed deficit. All supplied inviscid tubes contribute; signed entropy errors are retained. No BL wake mixing, skin friction, or finite-plane vector force balance is included.' };
}
