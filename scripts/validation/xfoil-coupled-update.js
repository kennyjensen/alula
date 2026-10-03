// Research variable-limiting/projection portion of XFOIL UPDATE, combined
// with the existing Euler density/stagnation update. This is not the full
// XFOIL algorithm: no negative-Ue island repair or MRCHDU is performed.
// Source xbl.f:1415–1485,1526–1542; final governing rows are unchanged.
import { proposeDensityNewton } from '../../src/euler/streamtube-density-newton.js';
import { projectXfoilDisplacement } from './xfoil-bl-update-projection.js';

const laminar = station => station.regime === 'laminar' || station.regime === 'similarity';

export function xfoilBLStepLimit(system, state, direction) {
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
      ['auxiliary', direction[k] / (isLaminar ? 10 : state[k])],
      ['theta', direction[k + 1] / state[k + 1]],
      ['delta-star', direction[k + 2] / state[k + 2]],
      ['edge-speed', Math.abs(direction[k + 3]) / .25],
    ];
    for (const [variable, dn] of increments) {
      let next;
      if (step * dn > 1.5) next = 1.5 / dn;
      if (step * dn < -.5) next = -.5 / dn;
      if (next !== undefined) {
        step = next; limiter = { kind: 'xfoil-bl-update', station: station.id, variable,
          normalizedIncrement: dn, bound: dn > 0 ? 1.5 : -.5 };
      }
    }
  }
  return { step, limiter };
}

export function proposeCoupledXfoilBLUpdate(system, state, direction, controls = {}) {
  if (system.bl.hasFiniteBase && !system.euler.layout.independentWakeBanks)
    throw new Error('Research finite-gap XFOIL update requires independently solved wake banks.');
  const maximumStep = controls.maximumStep ?? 1;
  if (!Number.isFinite(maximumStep) || maximumStep <= 0 || maximumStep > 1) throw new Error('Invalid maximum XFOIL update step.');
  const viscous = xfoilBLStepLimit(system, state, direction);
  system.euler.setDisplacement(system.bl.thicknesses(state.subarray(system.ne)));
  const density = proposeDensityNewton(system.euler, state.subarray(0, system.ne), direction.subarray(0, system.ne),
    { ...controls, maximumStep: Math.min(maximumStep, viscous.step) });
  const x = state.map((v, k) => v + density.step * direction[k]); x.set(density.x);
  let geometry;
  if (system.bl.hasFiniteBase) {
    try {
      system.euler.setDisplacement(system.bl.thicknesses(x.subarray(system.ne)));
      geometry = system.bl.geometry(x.subarray(0, system.ne));
    } finally { system.euler.setDisplacement(system.bl.thicknesses(state.subarray(system.ne))); }
  }
  const displacementChanges = [], auxiliaryChanges = [], { mach, gamma } = system.bl.kernel.parameters;
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
  return { ...density, x, stepKind: 'coupled-density-newton-xfoil-bl-update', viscousStep: viscous.step, viscousLimiter: viscous.limiter,
    limiter: viscous.step <= maximumStep && density.step === viscous.step && viscous.step < 1 ? viscous.limiter : density.limiter,
    projection: { method: 'XFOIL UPDATE variable limits, Ctau cap and DSLIM', active,
      displacementChanges, auxiliaryChanges, equationsChanged: false,
      preservedAfterGlobalStep: ['Euler unknowns', 'Ue', 'theta', 'material trips', 'physical wake gap'],
      massDefect: 'TOTAL displacement times Ue, including the native DSLIM change.',
      omitted: ['Giles H>1 proposal limiter', 'XFOIL negative-Ue island repair', 'MRCHDU'] },
    ...(active ? { meritComparable: false } : {}) };
}
