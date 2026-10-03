// SPDX-License-Identifier: GPL-2.0-or-later
// Prescribed station control from the two boundary polygons. Their edge
// lengths are measured exactly; they approximate arc on curved walls/cuts.
// If both normalized interval sequences obey h[i+1] <= R*h[i] and the
// reverse inequality, a convex combination obeys the same bounds. This is
// an arc-spacing property, not an orthogonality or positive-cell proof.
export function createBoundaryArcStationMap({ lower, upper }) {
  const normalized = points => {
    if (!Array.isArray(points) || points.length < 2 || points.some(p => !p || !Number.isFinite(p.x) || !Number.isFinite(p.y)))
      throw new Error('Boundary arc stations require finite polygon points.');
    const positions = [0], intervals = [];
    for (let i = 1; i < points.length; i++) {
      const h = Math.hypot(points[i].x - points[i - 1].x, points[i].y - points[i - 1].y);
      if (!(h > 0) || !Number.isFinite(h)) throw new Error('Boundary arc intervals must be positive and finite.');
      intervals.push(h); positions.push(positions.at(-1) + h);
    }
    const length = positions.at(-1);
    if (!Number.isFinite(length)) throw new Error('Boundary arc length must be finite.');
    const maximumAdjacentRatio = intervals.slice(1).reduce((ratio, h, i) => Math.max(ratio, h / intervals[i], intervals[i] / h), 1);
    return { length, positions: positions.map(s => s / length), maximumAdjacentRatio };
  };
  const a = normalized(lower), b = normalized(upper);
  if (lower.length !== upper.length) throw new Error('Boundary arc stations require matching station counts.');
  return {
    metric: 'normalized boundary polygon arc', lowerLength: a.length, upperLength: b.length,
    maximumAdjacentRatio: Math.max(a.maximumAdjacentRatio, b.maximumAdjacentRatio),
    at(fraction) {
      if (!Number.isFinite(fraction) || fraction < 0 || fraction > 1) throw new Error('Boundary arc blending fraction must lie in [0,1].');
      return a.positions.map((s, i) => (1 - fraction) * s + fraction * b.positions[i]);
    },
  };
}
