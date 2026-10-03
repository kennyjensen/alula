// SPDX-License-Identifier: GPL-2.0-or-later
// Fixed curved-edge corrections to a bilinear reference grid. This constructs
// a geometry map; it does not certify its global Jacobian or transversality.
// Boundary callbacks receive (fractionalIndex, {interval}). The interval
// selects the one-sided derivative at a corner between two curve pieces.
const finite = p => Number.isFinite(p?.x) && Number.isFinite(p?.y);
const copy = p => ({ x: p.x, y: p.y });

export function createCoonsGridGeometry({ nodes, boundary = {} }) {
  const nx = nodes?.length - 1, nt = nodes?.[0]?.length - 1;
  if (!(nx >= 1 && nt >= 1) || !Array.isArray(nodes)
    || nodes.some(row => !Array.isArray(row) || row.length !== nt + 1 || row.some(p => !finite(p))))
    throw new Error('Invalid Coons reference grid.');
  if (!boundary || typeof boundary !== 'object' || Object.entries(boundary).some(([side, value]) =>
    !['bottom', 'top', 'left', 'right'].includes(side) || typeof value !== 'function')) throw new Error('Invalid Coons boundary callbacks.');
  const initial = nodes.map(row => row.map(copy));
  let scale = 0;
  for (const row of initial) for (const p of row) scale = Math.max(scale, Math.hypot(p.x - initial[0][0].x, p.y - initial[0][0].y));
  if (!(scale > 0) || !Number.isFinite(scale)) throw new Error('Degenerate Coons reference grid.');
  const count = side => ['bottom', 'top'].includes(side) ? nx : nt;
  const nodeAt = (side, k) => side === 'bottom' ? initial[k][0] : side === 'top' ? initial[k][nt]
    : side === 'left' ? initial[0][k] : initial[nx][k];
  const sample = (side, k, interval) => {
    const value = boundary[side](k, { interval });
    if (!finite(value?.point) || !finite(value?.derivative)) throw new Error(`Invalid ${side} Coons boundary sample.`);
    return { point: copy(value.point), derivative: copy(value.derivative) };
  };
  const ends = Object.fromEntries(Object.keys(boundary).map(side => [side,
    Array.from({ length: count(side) + 1 }, (_, k) => {
      const p = sample(side, k, Math.min(k, count(side) - 1)).point, q = nodeAt(side, k);
      if (Math.hypot(p.x - q.x, p.y - q.y) > 64 * Number.EPSILON * scale)
        throw new Error(`The ${side} curve does not match its prescribed grid nodes.`);
      return p;
    })]));
  // Each interval owns its endpoint tangents. A callback may change tangent
  // at a cut/wall junction, but both one-sided positions must still agree.
  for (const side of Object.keys(boundary)) for (let k = 1; k < count(side); k++) {
    const p = sample(side, k, k - 1).point, q = nodeAt(side, k);
    if (Math.hypot(p.x - q.x, p.y - q.y) > 64 * Number.EPSILON * scale)
      throw new Error(`The ${side} one-sided curve does not match its prescribed grid nodes.`);
  }
  const edge = (side, k, u) => {
    const a = ends[side][k], b = ends[side][k + 1], value = sample(side, k + u, k);
    const dx = b.x - a.x, dy = b.y - a.y;
    return { point: { x: (value.point.x - a.x) - u * dx, y: (value.point.y - a.y) - u * dy },
      derivative: { x: value.derivative.x - dx, y: value.derivative.y - dy } };
  };
  const validIndex = (i, j, s, t) => {
    if (!Number.isInteger(i) || i < 0 || i >= nx || !Number.isInteger(j) || j < 0 || j >= nt
      || !Number.isFinite(s) || !Number.isFinite(t) || s < 0 || s > 1 || t < 0 || t > 1)
      throw new Error('Coons cell coordinates are outside the reference grid.');
  };
  const correction = (i, j, s, t) => {
    validIndex(i, j, s, t);
    const point = { x: 0, y: 0 }, ds = { x: 0, y: 0 }, dt = { x: 0, y: 0 };
    for (const [side, active, k, u, weight, sign] of [
      ['bottom', j === 0, i, s, 1 - t, -1], ['top', j === nt - 1, i, s, t, 1],
      ['left', i === 0, j, t, 1 - s, -1], ['right', i === nx - 1, j, t, s, 1],
    ]) if (active && boundary[side]) {
      const e = edge(side, k, u), horizontal = side === 'bottom' || side === 'top';
      for (const key of ['x', 'y']) {
        point[key] += weight * e.point[key];
        ds[key] += horizontal ? weight * e.derivative[key] : sign * e.point[key];
        dt[key] += horizontal ? sign * e.point[key] : weight * e.derivative[key];
      }
    }
    return { point, ds, dt };
  };
  const at = (i, j, s, t) => {
    const r = correction(i, j, s, t), p = [initial[i][j], initial[i + 1][j], initial[i + 1][j + 1], initial[i][j + 1]];
    const shape = [(1 - s) * (1 - t), s * (1 - t), s * t, (1 - s) * t];
    const ds = [t - 1, 1 - t, t, -t], dt = [s - 1, -s, s, 1 - s];
    r.point.x += p[0].x; r.point.y += p[0].y;
    for (let k = 1; k < 4; k++) for (const key of ['x', 'y']) {
      const d = p[k][key] - p[0][key];
      r.point[key] += shape[k] * d; r.ds[key] += ds[k] * d; r.dt[key] += dt[k] * d;
    }
    return r;
  };
  return { nx, nt, initial, at, correction,
    scope: 'Fixed physical boundary curves, bilinear interior edges, Coons edge corrections. Geometry construction only; no global Jacobian certificate.' };
}
