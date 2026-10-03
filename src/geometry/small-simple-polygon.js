// SPDX-License-Identifier: GPL-2.0-or-later
// Fixed-size hot paths. The caller validates finite vertices. Translation,
// triangle summation and intersection order match the general polygon check.
const orient = (ax, ay, bx, by, cx, cy) => (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
const between = (ax, ay, bx, by, cx, cy) => cx >= Math.min(ax, bx) && cx <= Math.max(ax, bx)
  && cy >= Math.min(ay, by) && cy <= Math.max(ay, by);
const opposite = (a, b) => a < 0 && b > 0 || a > 0 && b < 0;
function intersects(ax, ay, bx, by, cx, cy, dx, dy) {
  const u = orient(ax, ay, bx, by, cx, cy), v = orient(ax, ay, bx, by, dx, dy);
  const w = orient(cx, cy, dx, dy, ax, ay), z = orient(cx, cy, dx, dy, bx, by);
  return opposite(u, v) && opposite(w, z) || u === 0 && between(ax, ay, bx, by, cx, cy)
    || v === 0 && between(ax, ay, bx, by, dx, dy) || w === 0 && between(cx, cy, dx, dy, ax, ay)
    || z === 0 && between(cx, cy, dx, dy, bx, by);
}
function edge(ax, ay, bx, by, cx, cy, label) {
  // Coordinate differences are finite on the fast path; no norm is needed
  // to distinguish a zero-length edge, including subnormal coordinates.
  if (bx === ax && by === ay || orient(ax, ay, bx, by, cx, cy) === 0
    && (bx - ax) * (cx - bx) + (by - ay) * (cy - by) <= 0)
    throw new Error(`Degenerate or reversed ${label} edge.`);
}
function crossed(label) { throw new Error(`Crossed, reversed or degenerate ${label} polygon.`); }

export function smallSimpleQuad(p, label) {
  const x0 = p[0].x - p[0].x, y0 = p[0].y - p[0].y;
  const x1 = p[1].x - p[0].x, y1 = p[1].y - p[0].y;
  const x2 = p[2].x - p[0].x, y2 = p[2].y - p[0].y;
  const x3 = p[3].x - p[0].x, y3 = p[3].y - p[0].y;
  if (!Number.isFinite(x1) || !Number.isFinite(y1) || !Number.isFinite(x2) || !Number.isFinite(y2) || !Number.isFinite(x3) || !Number.isFinite(y3)) return undefined;
  let twiceArea = 0;
  edge(x0, y0, x1, y1, x2, y2, label);
  if (intersects(x0, y0, x1, y1, x2, y2, x3, y3)) crossed(label);
  edge(x1, y1, x2, y2, x3, y3, label);
  twiceArea += orient(x0, y0, x1, y1, x2, y2);
  if (intersects(x1, y1, x2, y2, x3, y3, x0, y0)) crossed(label);
  edge(x2, y2, x3, y3, x0, y0, label);
  twiceArea += orient(x0, y0, x2, y2, x3, y3);
  edge(x3, y3, x0, y0, x1, y1, label);
  if (!(twiceArea > 0) || !Number.isFinite(twiceArea)) crossed(label);
  return .5 * twiceArea;
}

export function smallSimpleHexagon(p, label) {
  const x0 = p[0].x - p[0].x, y0 = p[0].y - p[0].y;
  const x1 = p[1].x - p[0].x, y1 = p[1].y - p[0].y;
  const x2 = p[2].x - p[0].x, y2 = p[2].y - p[0].y;
  const x3 = p[3].x - p[0].x, y3 = p[3].y - p[0].y;
  const x4 = p[4].x - p[0].x, y4 = p[4].y - p[0].y;
  const x5 = p[5].x - p[0].x, y5 = p[5].y - p[0].y;
  if (!Number.isFinite(x1) || !Number.isFinite(y1) || !Number.isFinite(x2) || !Number.isFinite(y2) || !Number.isFinite(x3) || !Number.isFinite(y3) || !Number.isFinite(x4) || !Number.isFinite(y4) || !Number.isFinite(x5) || !Number.isFinite(y5)) return undefined;
  let twiceArea = 0;
  edge(x0, y0, x1, y1, x2, y2, label);
  if (intersects(x0, y0, x1, y1, x2, y2, x3, y3)) crossed(label);
  if (intersects(x0, y0, x1, y1, x3, y3, x4, y4)) crossed(label);
  if (intersects(x0, y0, x1, y1, x4, y4, x5, y5)) crossed(label);
  edge(x1, y1, x2, y2, x3, y3, label);
  twiceArea += orient(x0, y0, x1, y1, x2, y2);
  if (intersects(x1, y1, x2, y2, x3, y3, x4, y4)) crossed(label);
  if (intersects(x1, y1, x2, y2, x4, y4, x5, y5)) crossed(label);
  if (intersects(x1, y1, x2, y2, x5, y5, x0, y0)) crossed(label);
  edge(x2, y2, x3, y3, x4, y4, label);
  twiceArea += orient(x0, y0, x2, y2, x3, y3);
  if (intersects(x2, y2, x3, y3, x4, y4, x5, y5)) crossed(label);
  if (intersects(x2, y2, x3, y3, x5, y5, x0, y0)) crossed(label);
  edge(x3, y3, x4, y4, x5, y5, label);
  twiceArea += orient(x0, y0, x3, y3, x4, y4);
  if (intersects(x3, y3, x4, y4, x5, y5, x0, y0)) crossed(label);
  edge(x4, y4, x5, y5, x0, y0, label);
  twiceArea += orient(x0, y0, x4, y4, x5, y5);
  edge(x5, y5, x0, y0, x1, y1, label);
  if (!(twiceArea > 0) || !Number.isFinite(twiceArea)) crossed(label);
  return .5 * twiceArea;
}
