// SPDX-License-Identifier: GPL-2.0-or-later
// Analytic continuation from the upper TE to a finite-base wake-center seed.
// The integration path may change; the field, solid and endpoints never do.
import { potentialDifference } from './streamfunction.js';
import { pointInside } from '../geometry/airfoil.js';
const same = (a, b) => a.x === b.x && a.y === b.y;
const finite = p => Number.isFinite(p?.x) && Number.isFinite(p?.y);
const cross = (a, b, p) => (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
const side = (a, b, p) => {
  const dx = b.x - a.x, dy = b.y - a.y, px = p.x - a.x, py = p.y - a.y;
  const error = 16 * Number.EPSILON * (Math.abs(dx * py) + Math.abs(dy * px)
    + Math.max(Math.abs(a.x), Math.abs(a.y), Math.abs(b.x), Math.abs(b.y), Math.abs(p.x), Math.abs(p.y)) * Math.hypot(dx, dy));
  const z = cross(a, b, p); return z > error ? 1 : z < -error ? -1 : 0;
};
const on = (a, b, p) => side(a, b, p) === 0 && p.x >= Math.min(a.x, b.x) && p.x <= Math.max(a.x, b.x)
  && p.y >= Math.min(a.y, b.y) && p.y <= Math.max(a.y, b.y);

function exteriorPath(points, field, anchor) {
  const sheets = [...field.panels, ...(field.basePanels ?? [])], bodies = new Map();
  for (const p of sheets) { if (!bodies.has(p.element)) bodies.set(p.element, []); bodies.get(p.element).push(p); }
  const contours = Array.from(bodies.values(), panels => {
    if (panels.some((p, i) => i && !same(panels[i - 1].b, p.a)) || !same(panels.at(-1).b, panels[0].a))
      throw new Error('Exterior finite-base potential continuation requires complete ordered solid contours.');
    return [...panels.map(p => p.a), panels.at(-1).b];
  });
  for (let i = 1; i < points.length; i++) {
    const a = points[i - 1], b = points[i];
    if (same(a, b)) continue;
    for (const p of sheets) {
      const proper = side(a, b, p.a) * side(a, b, p.b) === -1 && side(p.a, p.b, a) * side(p.a, p.b, b) === -1;
      const forbidden = (q, c, d) => on(c, d, q) && !(same(q, anchor) && (same(c, anchor) || same(d, anchor)));
      if (proper || forbidden(p.a, a, b) || forbidden(p.b, a, b) || forbidden(a, p.a, p.b) || forbidden(b, p.a, p.b))
        throw new Error('Exterior finite-base potential continuation touches or crosses a physical sheet.');
    }
    for (const q of [a, { x: .5 * (a.x + b.x), y: .5 * (a.y + b.y) }, b])
      if (!same(q, anchor) && contours.some(c => pointInside(q, c)))
        throw new Error('Exterior finite-base potential continuation enters a solid body.');
  }
}

export function finiteBaseWakePotentialDifference({ trailingEdge: a, wakeSeed: b, direction, element, field }) {
  if (!finite(a) || !finite(b) || !finite(direction) || !(Math.hypot(direction.x, direction.y) > 0)
    || !Number.isInteger(element) || !Array.isArray(field?.panels) || !Array.isArray(field?.basePanels))
    throw new Error('Invalid finite-base wake potential continuation.');
  const panels = field.basePanels.filter(p => p.element === element);
  if (!panels.length || !same(a, panels.at(-1).b))
    throw new Error('Finite-base wake potential must start at the retained upper TE corner.');
  try {
    const value = potentialDifference(a, b, field);
    exteriorPath([a, b], field, a);
    return { value, path: [a, b].map(p => ({ ...p })), detour: false };
  }
  catch (error) { if (error.message !== 'Potential-difference path crosses a finite-base sheet.') throw error; }
  // A short TE-to-center chord can enter the measured cap. Translate both
  // endpoints along the prescribed downstream direction to a common exterior
  // region, connect there, then return along the seed ray. This is a sufficient
  // local chart for caps whose segment normals face downstream, not a general
  // obstacle router. Every leg is checked against every actual physical sheet.
  const length = Math.hypot(direction.x, direction.y), t = { x: direction.x / length, y: direction.y / length };
  const scale = Math.max(...panels.flatMap(p => [p.length, Math.abs(p.a.x), Math.abs(p.a.y), Math.abs(p.b.x), Math.abs(p.b.y)]),
    Math.abs(a.x), Math.abs(a.y), Math.abs(b.x), Math.abs(b.y));
  const roundoffClearance = 128 * Number.EPSILON * scale;
  let translation = 0;
  for (const panel of panels) {
    const nx = panel.ty, ny = -panel.tx, projection = nx * t.x + ny * t.y;
    if (!(projection > 0) || !Number.isFinite(projection))
      throw new Error('Finite-base wake potential needs an exterior chart for a cap that does not face downstream.');
    for (const p of [a, b]) translation = Math.max(translation,
      (roundoffClearance - (p.x - panel.a.x) * nx - (p.y - panel.a.y) * ny) / projection);
  }
  if (!(translation > 0) || !Number.isFinite(translation)) throw new Error('No representable finite-base exterior potential chart.');
  const path = [a, { x: a.x + translation * t.x, y: a.y + translation * t.y },
    { x: b.x + translation * t.x, y: b.y + translation * t.y }, b].map(p => ({ ...p }));
  exteriorPath(path, field, a);
  let value = 0;
  for (let i = 1; i < path.length; i++) value += potentialDifference(path[i - 1], path[i], field);
  return { value, path, detour: true, translation, roundoffClearance,
    method: 'Analytic potential continuation on three exterior segments; unchanged physical sheets and endpoints.' };
}
