// SPDX-License-Identifier: GPL-2.0-or-later
// Explicit UI choices map to the existing Worker API. A hidden checkbox or
// a remembered viscous Mach number cannot override the selected solver.
export function flowModelRoute(mode, quadMach = 0.2) {
  if (mode === 'inviscid' || mode === 'coupled') return { flowModel: mode, mach: 0 };
  if (mode === 'streamtube-grid' || mode === 'streamtube-bl') {
    if (!Number.isFinite(quadMach) || quadMach <= 0) throw new Error('Enter a positive quad Euler Mach number.');
    return { flowModel: 'streamtube-grid', mach: quadMach, quadBoundaryLayers: mode === 'streamtube-bl' };
  }
  throw new Error('Select one of the four supported flow models.');
}
