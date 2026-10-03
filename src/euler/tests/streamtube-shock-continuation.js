// SPDX-License-Identifier: GPL-2.0-or-later
// Test-only continuation policy for the historical automatic Mach experiment.
// Avoid regrowing immediately into a recently failed Mach increment. This
// controls the path only; all accepted stages use the original tolerance.
export function coupledMachGrowth({ iterations, budget, backtracks = 0, failedLargerStep = false, shockAudit, shockMovement }) {
  const difficult = failedLargerStep || iterations >= .75 * budget || backtracks >= 4
    || shockMovement?.maximumCellWidths > 2 || shockAudit?.maximumSkewDegrees > 45 || shockAudit?.weakMomentumCandidates > 0;
  return { factor: difficult ? 1 : 1.5,
    reason: difficult ? 'Hold Mach increment after difficult stage' : 'Grow Mach increment after inexpensive stage',
    iterations, budget, backtracks, failedLargerStep, shockMovement };
}

// Match the strongest retained compression in each tube. This measures
// candidate motion; topology changes/new candidates are reported separately.
export function compareShockAudits(before, after) {
  const strongest = audit => {
    const map = new Map();
    for (const c of audit?.candidates ?? []) {
      const key = `${c.group}/${c.tube}`;
      if (!map.has(key) || c.pressureRise > map.get(key).pressureRise) map.set(key, c);
    }
    return map;
  };
  const a = strongest(before), b = strongest(after);
  let matched = 0, maximumCellWidths = 0;
  for (const [key, c] of b) {
    const p = a.get(key); if (!p) continue;
    matched++;
    maximumCellWidths = Math.max(maximumCellWidths,
      Math.hypot(c.x - p.x, c.y - p.y) / (.5 * (c.streamwiseSpacing + p.streamwiseSpacing)));
  }
  return { matched, maximumCellWidths: matched ? maximumCellWidths : null,
    unmatchedCandidates: b.size - matched, interpretation: 'Motion of strongest retained compression per tube, not branch identification' };
}

