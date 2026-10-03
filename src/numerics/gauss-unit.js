// SPDX-License-Identifier: GPL-2.0-or-later
// Gauss-Legendre rules on [0,1]. The order is fixed for each assembled system.
export function gaussUnitRule(order = 3) {
  if (order === 3) return { nodes: [.5 - Math.sqrt(3 / 5) / 2, .5, .5 + Math.sqrt(3 / 5) / 2],
    weights: [5 / 18, 4 / 9, 5 / 18] };
  if (order === 5) {
    const inner = Math.sqrt(5 - 2 * Math.sqrt(10 / 7)) / 6, outer = Math.sqrt(5 + 2 * Math.sqrt(10 / 7)) / 6;
    const wi = (322 + 13 * Math.sqrt(70)) / 1800, wo = (322 - 13 * Math.sqrt(70)) / 1800;
    return { nodes: [.5 - outer, .5 - inner, .5, .5 + inner, .5 + outer], weights: [wo, wi, 64 / 225, wi, wo] };
  }
  throw new Error('Gauss integration order must be 3 or 5.');
}
