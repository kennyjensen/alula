// Research-only XFOIL UPDATE/DSLIM proposal control. No governing row or
// admissibility condition is changed. Original xbl.f:1530–1542,1579–1589;
// HKIN xblsys.f:2276–2286. Not imported by production.
import { proposeCoupledDensityNewton } from '../../src/euler/streamtube-density-newton.js';

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

export function proposeCoupledDensityNewtonWithXfoilDstrProjection(system, state, direction, controls = {}) {
  if (system.bl.hasFiniteBase && !system.euler.layout.independentWakeBanks)
    throw new Error('Research finite-gap DSLIM requires independently solved wake banks.');
  const proposal = proposeCoupledDensityNewton(system, state, direction, controls);
  let geometry;
  if (system.bl.hasFiniteBase) {
    // Independent banks and a fixed solid TE-center arc origin make WGAP
    // depend on the Euler coordinates, not on this delta-only correction.
    // Inspect the proposed coordinates, then restore the accepted context.
    try {
      system.euler.setDisplacement(system.bl.thicknesses(proposal.x.subarray(system.ne)));
      geometry = system.bl.geometry(proposal.x.subarray(0, system.ne));
    } finally { system.euler.setDisplacement(system.bl.thicknesses(state.subarray(system.ne))); }
  }
  const { mach, gamma } = system.bl.kernel.parameters, changes = [];
  for (const station of system.bl.stations) {
    const k = system.ne + 4 * station.id, before = proposal.x[k + 2], gap = (geometry?.coordinates[station.id].wakeGap ?? 0) / system.bl.scale;
    const result = projectXfoilDisplacement({ theta: proposal.x[k + 1], deltaStar: before, ue: proposal.x[k + 3],
      wakeGap: gap, wake: station.kind === 'wake', mach, gamma });
    if (result.active) {
      proposal.x[k + 2] = result.deltaStar;
      changes.push({ ...station, ...result, beforeDeltaStar: before,
        thicknessUnits: 'BL packed theta/delta units', physicalNormalizedCorrection: result.correction * system.bl.scale,
        physicalNormalizedMassDefectChange: result.massDefectChange * system.bl.scale });
    }
  }
  return { ...proposal, projection: { method: 'XFOIL UPDATE DSLIM', active: changes.length > 0, changes,
    equationsChanged: false, preservedFields: ['Euler unknowns', 'Ue', 'theta', 'auxiliary', 'material trips', 'physical wake gap'],
    massDefect: 'Recomputed implicitly as TOTAL deltaStar times Ue; it changes with the projected displacement.' },
    ...(changes.length ? { meritComparable: false } : {}) };
}
