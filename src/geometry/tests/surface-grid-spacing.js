// SPDX-License-Identifier: GPL-2.0-or-later
// Geometry-only screen for a surface window intended to have nearly uniform
// spacing. This measures cell-edge chords, not a spline's exact arc length.
// Supply the geometric chord explicitly; the stagnation point is not the LE.
export function measureSurfaceGridSpacing({ points, leadingEdge, trailingEdge,
  range = [.2, .8], maximumLengthRatio = 1.5, maximumCoefficientVariation = .2,
  minimumIntervals = 4 }) {
  const finitePoint = p => p && Number.isFinite(p.x) && Number.isFinite(p.y);
  if (!Array.isArray(points) || points.length < 2 || !points.every(finitePoint)
    || !finitePoint(leadingEdge) || !finitePoint(trailingEdge)
    || !Array.isArray(range) || range.length !== 2 || !range.every(Number.isFinite)
    || !(0 <= range[0] && range[0] < range[1] && range[1] <= 1)
    || !Number.isFinite(maximumLengthRatio) || maximumLengthRatio < 1
    || !Number.isFinite(maximumCoefficientVariation) || maximumCoefficientVariation < 0
    || !Number.isInteger(minimumIntervals) || minimumIntervals < 2)
    throw new Error('Invalid surface-spacing screen inputs.');
  const dx = trailingEdge.x - leadingEdge.x, dy = trailingEdge.y - leadingEdge.y;
  const chord = Math.hypot(dx, dy);
  if (!(chord > 0 && Number.isFinite(chord))) throw new Error('A nonzero geometric chord is required.');
  const position = points.map(p => ((p.x - leadingEdge.x) / chord * dx / chord
    + (p.y - leadingEdge.y) / chord * dy / chord));
  if (!position.every(Number.isFinite)) throw new Error('Surface chord coordinates are not resolved.');
  const intervals = [];
  for (let i = 1; i < points.length; i++) {
    const a = points[i - 1], b = points[i], length = Math.hypot(b.x - a.x, b.y - a.y) / chord;
    if (!(length > 0 && Number.isFinite(length))) throw new Error('Surface nodes must be distinct and resolved.');
    // LE wraps can run briefly upstream. Only the selected middle window
    // requires increasing chordwise position; reject folds within it.
    if (Math.max(position[i - 1], position[i]) >= range[0]
      && Math.min(position[i - 1], position[i]) <= range[1] && !(position[i] > position[i - 1]))
      throw new Error('Surface nodes reverse chordwise direction in the selected window.');
    const tolerance = 32 * Number.EPSILON;
    if (position[i - 1] >= range[0] - tolerance && position[i] <= range[1] + tolerance)
      intervals.push({ from: i - 1, to: i, start: position[i - 1], end: position[i], length });
  }
  const limits = { maximumLengthRatio, maximumCoefficientVariation, minimumIntervals };
  const scope = 'Nearly uniform spacing screen in the specified geometric-chord window; engineering limits, not universal MSET or physical-flow requirements.';
  if (Math.min(...position) > range[0] || Math.max(...position) < range[1])
    return { status: 'insufficient-coverage', passed: false, range, limits, intervals, scope };
  if (intervals.length < minimumIntervals) return { status: 'insufficient-resolution', passed: false,
    range, limits, intervals, scope };
  const mean = intervals.reduce((sum, p) => sum + p.length, 0) / intervals.length;
  const smallest = intervals.reduce((a, b) => a.length < b.length ? a : b);
  const largest = intervals.reduce((a, b) => a.length > b.length ? a : b);
  const lengthRatio = largest.length / smallest.length;
  const coefficientVariation = Math.sqrt(intervals.reduce((sum, p) => sum + (p.length / mean - 1) ** 2, 0) / intervals.length);
  const maximumAdjacentRatio = Math.max(1, ...intervals.slice(1).map((p, i) =>
    Math.max(p.length / intervals[i].length, intervals[i].length / p.length)));
  const passed = lengthRatio <= maximumLengthRatio && coefficientVariation <= maximumCoefficientVariation;
  return { status: passed ? 'passed' : 'excessive-clustering', passed, range, limits, intervals,
    mean, smallest, largest, lengthRatio, coefficientVariation, maximumAdjacentRatio, scope };
}

// Extract wall nodes independently from the initializer's curve/maps. These
// are the vertices actually displayed by the quad mesh, including after SLOR.
export function streamtubeSurfaceNodes(mesh, { body, side }) {
  const { bodies } = mesh?.initialization?.potentialCrosslines ?? {};
  const { tubes, streamwiseSegments: nx } = mesh?.initialization ?? {};
  if (!Number.isInteger(body) || !bodies?.[body] || !['upper', 'lower'].includes(side)
    || !Number.isInteger(nx) || !Array.isArray(tubes) || tubes.length !== bodies.length + 1
    || mesh.cells.length !== nx * tubes.reduce((a, b) => a + b, 0))
    throw new Error('Recorded multielement wall connectivity is required.');
  const { leadingIndex, trailingIndex } = bodies[body], g = side === 'upper' ? body + 1 : body;
  if (!(Number.isInteger(leadingIndex) && Number.isInteger(trailingIndex)
    && leadingIndex > 0 && trailingIndex > leadingIndex && trailingIndex < nx))
    throw new Error('Invalid surface index range.');
  const nt = tubes[g], offset = nx * tubes.slice(0, g).reduce((a, b) => a + b, 0);
  return Array.from({ length: trailingIndex - leadingIndex + 1 }, (_, k) => {
    const cell = mesh.cells[offset + (leadingIndex + k) * nt + (side === 'upper' ? 0 : nt - 1)];
    if (cell.length !== 4) throw new Error('Surface spacing requires quadrilateral cells.');
    return mesh.vertices[cell[side === 'upper' ? 0 : 3]];
  });
}
