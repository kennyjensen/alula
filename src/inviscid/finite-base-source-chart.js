// SPDX-License-Identifier: GPL-2.0-or-later
import { finiteBaseSourceChart } from './finite-base-influence.js';

const cross = (a, b) => a.x * b.y - a.y * b.x;
const subtract = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });

// A source's logarithm needs one common branch ray. A ray through another
// solid makes that body's constant-streamfunction boundary condition jump.
// Select a visible downstream sector from the geometry, identically in the
// inviscid and displacement operators. This changes the chart, not TECALC's
// physical tangent or source/vortex coefficients.
export function assemblyFiniteBaseSourceChart(topology, preferred, otherSolids) {
  const initial = finiteBaseSourceChart(topology, preferred);
  const blocked = chart => otherSolids.some(points => points.some((a, i) => {
    const b = points[(i + 1) % points.length], edge = subtract(b, a), offset = subtract(a, chart.origin);
    const denominator = cross(chart.direction, edge);
    const tolerance = 64 * Number.EPSILON * Math.max(Math.hypot(offset.x, offset.y), Math.hypot(edge.x, edge.y));
    if (Math.abs(denominator) <= tolerance) return Math.abs(cross(offset, chart.direction)) <= tolerance
      && Math.max(offset.x * chart.direction.x + offset.y * chart.direction.y,
        (b.x - chart.origin.x) * chart.direction.x + (b.y - chart.origin.y) * chart.direction.y) >= 0;
    const ray = cross(offset, edge) / denominator, segment = cross(offset, chart.direction) / denominator;
    return ray >= -tolerance && segment >= -64 * Number.EPSILON && segment <= 1 + 64 * Number.EPSILON;
  }));
  if (!blocked(initial)) return initial;
  const angle = Math.atan2(preferred.y, preferred.x), center = topology.trailingEdge.center;
  const angles = [-Math.PI / 2, Math.PI / 2];
  for (const points of otherSolids) for (const p of points) {
    const a = Math.atan2(cross(preferred, subtract(p, center)),
      preferred.x * (p.x - center.x) + preferred.y * (p.y - center.y));
    if (a > -Math.PI / 2 && a < Math.PI / 2) angles.push(a);
  }
  angles.sort((a, b) => a - b);
  const chartAt = a => {
    try { return finiteBaseSourceChart(topology, { x: Math.cos(angle + a), y: Math.sin(angle + a) }); }
    catch { return null; } // An oblique ray may not admit a common origin inside the base.
  };
  const sectors = [];
  for (let i = 1; i < angles.length; i++) {
    if (!(angles[i] > angles[i - 1])) continue;
    const chart = chartAt(.5 * (angles[i] + angles[i - 1]));
    if (!chart || blocked(chart)) continue;
    if (sectors.at(-1)?.[1] === angles[i - 1]) sectors.at(-1)[1] = angles[i];
    else sectors.push([angles[i - 1], angles[i]]);
  }
  // Prefer the nearest free sector and stay in its interior, away from
  // tangency. Recheck the actual shifted chart origin before using it.
  sectors.sort((a, b) => Math.abs(a[0] + a[1]) - Math.abs(b[0] + b[1]));
  for (const [lo, hi] of sectors) {
    const chart = chartAt(.5 * (lo + hi));
    if (chart && !blocked(chart)) return chart;
  }
  throw new Error('No solid-free downstream branch ray for the finite-base source chart.');
}
