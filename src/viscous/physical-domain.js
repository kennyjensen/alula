// SPDX-License-Identifier: GPL-2.0-or-later
// Domain of XFOIL BLKIN/HKIN with physical edge velocity Ue/Uinf.
// This is a domain check, not a closure or a modification of BLVAR.
export function createPhysicalBLDomain({ mach = 0, gamma = 1.4 } = {}) {
  if (!Number.isFinite(mach) || mach < 0 || mach >= 1 || !Number.isFinite(gamma) || gamma <= 1)
    throw new Error('Invalid physical BL domain parameters.');
  const machFactor = mach ** 2 / (1 + .5 * (gamma - 1) * mach ** 2);
  const hstinv = (gamma - 1) * machFactor;
  return ({ theta, deltaStar, ue, wakeGap = 0 }) => {
    const enthalpy = 1 - .5 * hstinv * ue ** 2, enthalpyUe = -hstinv * ue;
    const gap = deltaStar - wakeGap - theta, compressibility = (.29 + .113) * machFactor * ue ** 2;
    // M_e^2 = machFactor*ue^2/enthalpy and
    // Hk = (H - .29*M_e^2)/(1 + .113*M_e^2).
    // Multiplication by positive theta*enthalpy removes the singular
    // denominator, so values remain evaluable outside the thermal domain.
    // Thicknesses may share any consistent units (including scaled BL
    // unknowns). Gradients below are ordered [theta, deltaStar, ue].
    return { enthalpy, shape: gap * enthalpy - theta * compressibility,
      enthalpyGradient: [0, 0, enthalpyUe],
      shapeGradient: [-enthalpy - compressibility, enthalpy,
        gap * enthalpyUe - 2 * theta * (.29 + .113) * machFactor * ue] };
  };
}
