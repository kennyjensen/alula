// SPDX-License-Identifier: GPL-2.0-or-later
// Richardson table for centered differences at h, h/2, h/4, ... .
// Undefined entries preserve invalid perturbations; they are never filled.
export function extrapolateCentral(differences, maximumOrder = 8) {
  if (!Number.isInteger(maximumOrder) || maximumOrder < 2 || maximumOrder % 2 || maximumOrder > 12)
    throw new Error('Expected an even extrapolation order from 2 through 12.');
  const levels = [{ order: 2, estimates: differences }];
  for (let order = 4; order <= maximumOrder; order += 2) {
    const previous = levels.at(-1).estimates, factor = 2 ** (order - 2);
    const estimates = previous.slice(0, -1).map((v, i) => {
      const finer = previous[i + 1]; if (!v || !finer) return undefined;
      if (v.length !== finer.length) throw new Error('Centered-difference vector sizes differ.');
      return v.map((q, row) => (factor * finer[row] - q) / (factor - 1));
    });
    levels.push({ order, estimates });
  }
  return levels;
}
