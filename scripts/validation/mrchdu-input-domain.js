// Research input validation for one local MRCHDU source comparison only.
// Never use this as a coupled-state or native OUTPUT admissibility test.
export function validateMrchduMarchInput(kernel, station, regime) {
  try { return kernel.station(station, regime); }
  catch (error) {
    const d = error?.diagnostics;
    if (error?.code !== 'BL_EDGE_STATE_DOMAIN' || d?.condition !== 'raw-hk'
      || d.failedChecks.length !== 1 || d.failedChecks[0] !== 'raw-hk'
      || !Number.isFinite(d.rawHk) || d.rawHk > 1
      || ![d.enthalpyRatio, d.density, d.viscosity, d.reTheta, d.machSquared, d.physicalUe].every(Number.isFinite)
      || !(d.enthalpyRatio > 0 && d.density > 0 && d.viscosity > 0 && d.reTheta > 0
        && d.machSquared >= 0 && d.physicalUe > 0)) throw error;
    // MRCHDU receives the exact supplied primary state. Its own BLVAR and
    // DSLIM operations establish the local iterate; no profile remap here.
    return { marchInputOnly: true, rawShapeOutsideOutputDomain: true, ...d };
  }
}
