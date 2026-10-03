// SPDX-License-Identifier: GPL-2.0-or-later
import { velocityAt } from './linear-vortex.js';
import { pointInside, segmentsTouch } from '../geometry/airfoil.js';

// Display traces only: adaptive midpoint integration of the panel velocity.
// These are not an MSES computational grid or a boundary-layer wake.
export function traceStreamlines(result, bounds, count = 25) {
  const lines = [];
  const step = (bounds.xMax - bounds.xMin) / 180;
  if (!(step > 0) || !Object.values(bounds).every(Number.isFinite) || !(bounds.yMax > bounds.yMin)
    || !Number.isInteger(count) || count < 1 || count > 100) throw new Error('Invalid streamline domain.');
  const inside = p => result.elements.some(e => pointInside(p, e.points));
  const direction = p => {
    const v = velocityAt(p, result.field); const length = Math.hypot(v.u, v.v);
    return length < 1e-7 || !Number.isFinite(length) ? null : { x: v.u / length, y: v.v / length };
  };
  for (let i = 0; i < count; i++) {
    let p = { x: bounds.xMin, y: bounds.yMin + (i + 0.5) / count * (bounds.yMax - bounds.yMin) };
    if (inside(p)) { lines.push([]); continue; }
    const points = [p];
    for (let j = 0; j < 400; j++) {
      const d = direction(p); if (!d) break;
      let next = null;
      for (let h = step; h >= step / 256; h /= 2) {
        const mid = { x: p.x + d.x * h / 2, y: p.y + d.y * h / 2 };
        if (inside(mid)) continue;
        const dm = direction(mid); if (!dm) continue;
        const candidate = { x: p.x + h * dm.x, y: p.y + h * dm.y };
        if (inside(candidate) || result.field.panels.some(panel => segmentsTouch(p, candidate, panel.a, panel.b, step * 1e-9))) continue;
        // A curved segment must still align with the field at its actual
        // midpoint. Reduce its length near stagnation points and tight gaps.
        const actual = direction({ x: (p.x + candidate.x) / 2, y: (p.y + candidate.y) / 2 });
        if (!actual || Math.abs(dm.x * actual.y - dm.y * actual.x) > 0.001 || dm.x * actual.x + dm.y * actual.y <= 0) continue;
        next = candidate; break;
      }
      if (!next) break;
      points.push(next); p = next;
      if (p.x > bounds.xMax || p.x < bounds.xMin - step || p.y > bounds.yMax || p.y < bounds.yMin) break;
    }
    lines.push(points);
  }
  return lines;
}
