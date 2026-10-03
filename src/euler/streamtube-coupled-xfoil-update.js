// SPDX-License-Identifier: GPL-2.0-or-later
// Variable-limiting/projection portion of XFOIL UPDATE, combined
// with the existing Euler density/stagnation update. This is not the full
// XFOIL algorithm: no negative-Ue island repair or MRCHDU is performed.
// Source xbl.f:1415–1485,1526–1542; final governing rows are unchanged.
import { proposeDensityNewton } from './streamtube-density-newton.js';

// The positive turbulent-shear trial map used by the optional log update.
export function logarithmicShearIncrement(current, direction, step) {
  if (!Number.isFinite(current) || current <= 0 || !Number.isFinite(direction)
    || !Number.isFinite(step) || step < 0 || step > 1)
    throw new Error('Invalid logarithmic shear increment.');
  if (step === 0 || direction === 0) return current;
  const ratio = direction / current;
  const change = Number.isFinite(ratio) ? step * ratio : step * direction / current;
  const value = Math.abs(change) <= 700
    ? current * Math.exp(change) : Math.exp(Math.log(current) + change);
  if (!Number.isFinite(change) || !Number.isFinite(value) || value <= 0)
    throw new Error('Logarithmic shear increment leaves the finite positive domain.');
  return value;
}

// UPDATE normalizes DUEDG by .25 and caps the normalized increase at 1.5.
// This port applies the rule to its physical packed BL speed; native UEI
// and physical U2 differ under XFOIL's compressibility transformation.
const edgeSpeedScale = .25, upperIncrementLimit = 1.5;
export const XFOIL_BL_EDGE_SPEED_LIMIT = edgeSpeedScale * upperIncrementLimit;

// Native DSLIM acts on fluid displacement and restores the physical gap.
// HKIN: original xblsys.f:2276–2286; DSLIM: xbl.f:1579–1589.
export function projectXfoilDisplacement({ theta, deltaStar, ue, wakeGap = 0, wake = false, mach = 0, gamma = 1.4 }) {
  if (![theta, deltaStar, ue, wakeGap, mach, gamma].every(Number.isFinite) || theta <= 0 || ue <= 0
    || wakeGap < 0 || mach < 0 || mach >= 1 || gamma <= 1 || typeof wake !== 'boolean' || !wake && wakeGap !== 0)
    throw new Error('Invalid XFOIL DSLIM physical proposal.');
  const hstinv = (gamma - 1) * mach ** 2 / (1 + .5 * (gamma - 1) * mach ** 2);
  const enthalpy = 1 - .5 * hstinv * ue ** 2;
  if (!(enthalpy > 0)) throw new Error('DSLIM cannot repair a nonpositive thermal state.');
  const machSquared = ue ** 2 * hstinv / ((gamma - 1) * enthalpy);
  const denominator = 1 + .113 * machSquared;
  const rawHk = ((deltaStar - wakeGap) / theta - .29 * machSquared) / denominator;
  const minimumHk = wake ? 1.00005 : 1.02;
  // The native formula is delta += theta*max(0,Hkmin-Hk)/(dHk/dH).
  // Keep inactive states bit-identical and preserve TOTAL wake displacement.
  const correction = theta * Math.max(0, minimumHk - rawHk) / (1 / denominator);
  const projected = correction === 0 ? deltaStar : (deltaStar - wakeGap + correction) + wakeGap;
  const actualCorrection = projected - deltaStar;
  return { deltaStar: projected, correction: actualCorrection, active: actualCorrection !== 0, rawHk,
    projectedRawHk: ((projected - wakeGap) / theta - .29 * machSquared) / denominator,
    minimumHk, machSquared, enthalpy, wakeGap, massDefectChange: ue * actualCorrection };
}

const laminar = station => station.regime === 'laminar' || station.regime === 'similarity';

// Only a candidate-domain failure can be recovered by reducing this trial's
// step. Bad accepted states, controls and programming errors must still stop.
export const XFOIL_CANDIDATE_DOMAIN = 'COUPLED_XFOIL_CANDIDATE_DOMAIN';
function candidateError(reason, details, cause) {
  return Object.assign(new Error(`XFOIL BL candidate is inadmissible: ${reason}`, cause ? { cause } : undefined), {
    code: XFOIL_CANDIDATE_DOMAIN, recoverable: true, step: details.step,
    diagnostics: { stage: 'xfoil-bl-proposal', reason, ...details,
      ...(cause ? { cause: { name: cause.name, message: cause.message, diagnostics: cause.diagnostics } } : {}) },
  });
}

const candidateGeometryErrors = new Set([
  'Collapsed BL wake interval.', 'Degenerate displacement-surface tangent.',
  'Contour parameter is outside the curve.', 'Contour arc parameter is outside the curve.',
  'Invalid stagnation/branch parameter.',
  'Invalid captured streamtube mass levels.', 'Nonfinite streamtube mass allocation.',
]);

function checkThermal(ue, mach, gamma) {
  const hstinv = (gamma - 1) * mach ** 2 / (1 + .5 * (gamma - 1) * mach ** 2);
  return 1 - .5 * hstinv * ue ** 2;
}

export function xfoilBLStepLimit(system, state, direction, { logarithmicShear = false } = {}) {
  if (state.length !== system.n || direction.length !== system.n || !state.every(Number.isFinite) || !direction.every(Number.isFinite))
    throw new Error('Invalid XFOIL coupled update state or direction.');
  let step = 1, limiter = { kind: 'full-step' };
  for (const station of system.bl.stations) {
    const k = system.ne + 4 * station.id, isLaminar = laminar(station);
    if (!(state[k + 1] > 0) || !(state[k + 2] > 0) || !(state[k + 3] > 0) || !isLaminar && !(state[k] > 0))
      throw new Error('Invalid accepted XFOIL BL update state.');
    // DDSTR is already the Newton displacement unknown in our state. It
    // equals XFOIL's linearized (DMASS - DSTR*DUEDG)/UEDG representation.
    const increments = [
      ...(logarithmicShear && !isLaminar ? [] : [['auxiliary', direction[k] / (isLaminar ? 10 : state[k])]]),
      ['theta', direction[k + 1] / state[k + 1]],
      ['delta-star', direction[k + 2] / state[k + 2]],
      ['edge-speed', Math.abs(direction[k + 3]) / edgeSpeedScale],
    ];
    for (const [variable, dn] of increments) {
      let next;
      if (step * dn > upperIncrementLimit) next = upperIncrementLimit / dn;
      if (step * dn < -.5) next = -.5 / dn;
      if (next !== undefined) {
        step = next; limiter = { kind: 'xfoil-bl-update', station: station.id, variable,
          normalizedIncrement: dn, bound: dn > 0 ? upperIncrementLimit : -.5 };
      }
    }
  }
  return { step, limiter };
}

function proposeXfoilUpdate(system, state, direction, controls, logarithmicShear) {
  if (system.bl.hasFiniteBase && !system.euler.layout.independentWakeBanks)
    throw new Error('Finite-gap XFOIL update requires independently solved wake banks.');
  const maximumStep = controls.maximumStep ?? 1;
  if (!Number.isFinite(maximumStep) || maximumStep <= 0 || maximumStep > 1) throw new Error('Invalid maximum XFOIL update step.');
  const { mach, gamma } = system.bl.kernel.parameters;
  if (![mach, gamma, system.bl.scale].every(Number.isFinite) || mach < 0 || mach >= 1 || gamma <= 1 || system.bl.scale <= 0)
    throw new Error('Invalid accepted XFOIL BL thermodynamic controls.');
  const viscous = xfoilBLStepLimit(system, state, direction, { logarithmicShear });
  for (const { id } of system.bl.stations) if (!(checkThermal(state[system.ne + 4 * id + 3], mach, gamma) > 0))
    throw new Error('Invalid accepted XFOIL BL thermal state.');
  system.euler.setDisplacement(system.bl.thicknesses(state.subarray(system.ne)));
  const density = proposeDensityNewton(system.euler, state.subarray(0, system.ne), direction.subarray(0, system.ne),
    { ...controls, maximumStep: Math.min(maximumStep, viscous.step) });
  const x = state.map((v, k) => v + density.step * direction[k]); x.set(density.x);
  const logarithmicShearChanges = [];
  if (logarithmicShear) for (const station of system.bl.stations) if (!laminar(station)) {
    const k = system.ne + 4 * station.id;
    try { x[k] = logarithmicShearIncrement(state[k], direction[k], density.step); }
    catch (cause) { throw candidateError('logarithmic shear coordinate', {
      station: station.id, current: state[k], direction: direction[k], step: density.step }, cause); }
    logarithmicShearChanges.push({ station: station.id, from: state[k], to: x[k],
      logIncrement: density.step * direction[k] / state[k] });
  }
  if (!x.every(Number.isFinite)) throw candidateError('nonfinite state', { step: density.step });
  for (const { id } of system.bl.stations) {
    const k = system.ne + 4 * id;
    if (!(x[k + 3] > 0)) throw candidateError('nonpositive edge speed', { station: id, ue: x[k + 3], step: density.step });
    if (!(x[k + 1] > 0 && x[k + 2] > 0))
      throw candidateError('nonpositive thickness', { station: id, theta: x[k + 1], deltaStar: x[k + 2], step: density.step });
    const enthalpy = checkThermal(x[k + 3], mach, gamma);
    if (!(enthalpy > 0)) throw candidateError('nonpositive thermal state', { station: id, ue: x[k + 3], enthalpy, step: density.step });
  }
  let geometry;
  if (system.bl.hasFiniteBase) {
    try {
      system.euler.setDisplacement(system.bl.thicknesses(x.subarray(system.ne)));
      geometry = system.bl.geometry(x.subarray(0, system.ne));
    } catch (cause) {
      if (candidateGeometryErrors.has(cause?.message))
        throw candidateError('candidate geometry', { step: density.step }, cause);
      throw cause;
    } finally { system.euler.setDisplacement(system.bl.thicknesses(state.subarray(system.ne))); }
  }
  const displacementChanges = [], auxiliaryChanges = [];
  for (const station of system.bl.stations) {
    const k = system.ne + 4 * station.id;
    if (!laminar(station) && x[k] > .25) {
      auxiliaryChanges.push({ ...station, before: x[k], after: .25, correction: .25 - x[k] });
      x[k] = .25;
    }
    const before = x[k + 2], wakeGap = (geometry?.coordinates[station.id].wakeGap ?? 0) / system.bl.scale;
    const p = projectXfoilDisplacement({ theta: x[k + 1], deltaStar: before, ue: x[k + 3], wakeGap,
      wake: station.kind === 'wake', mach, gamma });
    if (p.active) {
      x[k + 2] = p.deltaStar;
      displacementChanges.push({ ...station, ...p, beforeDeltaStar: before, thicknessUnits: 'BL packed theta/delta units',
        physicalNormalizedCorrection: p.correction * system.bl.scale,
        physicalNormalizedMassDefectChange: p.massDefectChange * system.bl.scale });
    }
  }
  const active = displacementChanges.length > 0 || auxiliaryChanges.length > 0;
  return { ...density, x, stepKind: logarithmicShear ? 'coupled-density-newton-log-shear' : 'coupled-density-newton-xfoil-bl-update',
    ...(logarithmicShear ? { logarithmicShear: logarithmicShearChanges } : {}), viscousStep: viscous.step, viscousLimiter: viscous.limiter,
    limiter: viscous.step <= maximumStep && density.step === viscous.step && viscous.step < 1 ? viscous.limiter : density.limiter,
    projection: { method: logarithmicShear ? 'Logarithmic Ctau Newton; original N/thickness/Ue limits and DSLIM'
      : 'XFOIL UPDATE variable limits, Ctau cap and DSLIM', active,
      displacementChanges, auxiliaryChanges, equationsChanged: false,
      preservedAfterGlobalStep: ['Euler unknowns', 'Ue', 'theta', 'material trips', 'physical wake gap'],
      massDefect: 'TOTAL displacement times Ue, including the native DSLIM change.',
      omitted: ['Giles H>1 proposal limiter', 'XFOIL negative-Ue island repair', 'MRCHDU'] },
    ...(active ? { meritComparable: false } : {}) };
}

export function proposeCoupledXfoilBLUpdate(system, state, direction, controls = {}) {
  return proposeXfoilUpdate(system, state, direction, controls, false);
}

export function proposeLogarithmicShearUpdate(system, state, direction, controls = {}) {
  return proposeXfoilUpdate(system, state, direction, controls, true);
}
