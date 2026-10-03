// SPDX-License-Identifier: GPL-2.0-or-later
// Port of GRAPE RELAX's P+R stencil, NASA Ames/COSMIC ARC-11379,
// third_party/grape/original.src lines 3123--3131. The original uses one
// uniform RDX; here each one-sided difference uses its own interval.
// For +D*r_xi in the elliptic equation, D>=0 uses the forward difference.
export function grapePoissonDrift({ drift, leftDistance, rightDistance }) {
  if (![drift, leftDistance, rightDistance].every(Number.isFinite) || !(leftDistance > 0 && rightDistance > 0))
    throw new Error('GRAPE drift requires a finite coefficient and positive computational intervals.');
  const first = drift < 0 ? [-1 / leftDistance, 1 / leftDistance, 0] : [0, -1 / rightDistance, 1 / rightDistance];
  const coefficients = first.map(w => drift * w);
  if (!coefficients.every(Number.isFinite) || !first.every(Number.isFinite)) throw new Error('Unresolved GRAPE drift stencil.');
  // At zero drift this is the positive one-sided generalized derivative.
  // A central finite difference across this kink need not equal it.
  return { coefficients, first };
}
