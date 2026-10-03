// SPDX-License-Identifier: GPL-2.0-or-later
// Unpromoted finite inverse of a sampled streamwise coordinate. Caller owns
// the scalar equation and every geometry/physics acceptance check.
export function invertSampledStreamtubeCoordinate(nodes, baseline, corrections, scale) {
  if (!Number.isFinite(scale) || scale < 0 || scale > 1) throw new Error('Invalid scale.');
  const n = nodes.length, width = nodes[0].length;
  if (baseline.length !== n || corrections.length !== n || corrections.some(row => row.length !== width))
    throw new Error('Coordinate dimensions differ.');
  const moved = structuredClone(nodes), targets = Array.from(baseline), labelsByRow = [];
  for (let i = 0; i < n; i++) {
    if (!Number.isFinite(targets[i]) || i && !(targets[i] > targets[i - 1])) throw new Error('Invalid baseline.');
    for (let j = 0; j < width; j++) if (!Number.isFinite(corrections[i][j])
      || !Number.isFinite(nodes[i][j]?.x) || !Number.isFinite(nodes[i][j]?.y)) throw new Error('Nonfinite input.');
  }
  for (let j = 0; j < width; j++) {
    if (corrections[0][j] !== 0 || corrections[n - 1][j] !== 0) throw new Error('Endpoints must be fixed.');
    const labels = targets.map((v, i) => v + scale * corrections[i][j]);
    for (let i = 1; i < n; i++) if (!(labels[i] > labels[i - 1]))
      throw Object.assign(new Error('Corrected coordinate is not strictly increasing.'),
        { stage: 'coordinate monotonicity', details: { i, j, increment: labels[i] - labels[i - 1] } });
    labelsByRow.push(labels);
    if (scale === 0 || corrections.every(row => row[j] === 0)) continue;
    let segment = 0;
    for (let i = 1; i < n - 1; i++) {
      while (segment + 1 < n - 1 && labels[segment + 1] < targets[i]) segment++;
      const fraction = (targets[i] - labels[segment]) / (labels[segment + 1] - labels[segment]);
      if (!(fraction >= 0 && fraction <= 1)) throw new Error('Inverse left its segment.');
      const a = nodes[segment][j], b = nodes[segment + 1][j];
      moved[i][j] = fraction === 0 ? { ...a } : fraction === 1 ? { ...b }
        : { ...nodes[i][j], x: a.x + fraction * (b.x - a.x), y: a.y + fraction * (b.y - a.y) };
    }
  }
  return { nodes: moved, labelsByRow };
}
