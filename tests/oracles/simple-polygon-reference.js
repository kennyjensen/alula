// SPDX-License-Identifier: GPL-2.0-or-later
// Signed-area and edge-intersection domain check. Concavity is permitted;
// crossed/touching nonadjacent edges, reversals and zero area are not.
export function requirePositiveSimplePolygon(vertices, label = 'finite-volume') {
  if (!Array.isArray(vertices) || vertices.length < 3 || vertices.some(p => !Number.isFinite(p?.x) || !Number.isFinite(p?.y)))
    throw new Error('A simple polygon requires at least three finite vertices.');
  const p = vertices.map(q => ({ x: q.x - vertices[0].x, y: q.y - vertices[0].y })), n = p.length;
  const orient = (a, b, c) => (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  const between = (a, b, c) => c.x >= Math.min(a.x, b.x) && c.x <= Math.max(a.x, b.x)
    && c.y >= Math.min(a.y, b.y) && c.y <= Math.max(a.y, b.y);
  const opposite = (a, b) => a < 0 && b > 0 || a > 0 && b < 0;
  const intersects = (a, b, c, d) => {
    const u = orient(a, b, c), v = orient(a, b, d), w = orient(c, d, a), z = orient(c, d, b);
    return opposite(u, v) && opposite(w, z) || u === 0 && between(a, b, c) || v === 0 && between(a, b, d)
      || w === 0 && between(c, d, a) || z === 0 && between(c, d, b);
  };
  let twiceArea = 0;
  for (let i = 0; i < n; i++) {
    const a = p[i], b = p[(i + 1) % n], c = p[(i + 2) % n];
    if (!(Math.hypot(b.x - a.x, b.y - a.y) > 0) || orient(a, b, c) === 0
      && (b.x - a.x) * (c.x - b.x) + (b.y - a.y) * (c.y - b.y) <= 0)
      throw new Error(`Degenerate or reversed ${label} edge.`);
    if (i > 0 && i < n - 1) twiceArea += orient(p[0], p[i], p[i + 1]);
    for (let j = i + 2; j < n; j++) if (!(i === 0 && j === n - 1)
      && intersects(a, b, p[j], p[(j + 1) % n]))
      throw new Error(`Crossed, reversed or degenerate ${label} polygon.`);
  }
  if (!(twiceArea > 0) || !Number.isFinite(twiceArea)) throw new Error(`Crossed, reversed or degenerate ${label} polygon.`);
  return .5 * twiceArea;
}
