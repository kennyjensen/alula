// SPDX-License-Identifier: GPL-2.0-or-later
// Bounded cold initialization policy, not a physical Mach limit or an MSES
// constant. Only typed gas-domain rejection before iteration can authorize
// another seed. Geometry failures and unconverged Newton solves cannot.
export function colderStreamtubeStartupMach(error, mach) {
  if (!Number.isFinite(mach) || mach <= .025) return null;
  const pressure = error?.code === 'streamtube-interface-pressure' && error.stage === 'gas-initialization';
  const capacity = error?.code === 'streamtube-sonic-capacity' && error.stage === 'gas-initialization'
    && error.diagnostics?.stage === 'gas-initialization';
  const gasDomainRejection = failure => failure == null
    || ['streamtube-interface-pressure', 'streamtube-static-enthalpy'].includes(failure.code);
  // Every attempted density guess must fail in a recognized gas domain.
  // A later geometry/unknown rejection remains terminal even when the
  // preceding stagnation-density guess failed only in pressure or enthalpy.
  // Reducing Mach cannot repair a folded cell or certify an unknown failure.
  const fallbacks = [error?.diagnostics?.stagnationDensityFallback,
    error?.diagnostics?.pressureDomainDensityFallback];
  if (!(pressure || capacity) || !fallbacks.every(gasDomainRejection)) return null;
  return Math.max(.025, .5 * mach);
}
