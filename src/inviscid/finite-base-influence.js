// SPDX-License-Identifier: GPL-2.0-or-later
// Constant source/vortex influences on every retained finite-base segment.
// Source logs share one downstream cut from the element's TE center. The
// paths joining segment midpoints to that center stay inside the base cap.
import { pointInside } from '../geometry/airfoil.js';
const moments = Array.from({ length: 49 }, (_, k) => k % 2 ? 0 : 1 / ((k + 1) * 2 ** k));
const same = (a, b) => a.x === b.x && a.y === b.y;
const cross = (a, b) => a.x * b.y - a.y * b.x;
const dot = (a, b) => a.x * b.x + a.y * b.y;
const sub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });
const side = (a, b, p) => {
  const d = sub(b, a), r = sub(p, a), value = cross(d, r);
  const error = 32 * Number.EPSILON * (Math.abs(d.x * r.y) + Math.abs(d.y * r.x)
    + Math.max(Math.abs(a.x), Math.abs(a.y), Math.abs(b.x), Math.abs(b.y), Math.abs(p.x), Math.abs(p.y)) * Math.hypot(d.x, d.y));
  return value > error ? 1 : value < -error ? -1 : 0;
};

// A source chart changes no coordinate or physical influence. If roundoff in
// a retained nearly straight base prevents visibility from its chord center,
// choose the first interior point on the upstream continuation of the SAME
// downstream ray. The segment half-planes give this distance analytically.
export function finiteBaseSourceChart(topology, direction) {
  const center = topology.trailingEdge.center;
  const roundoff = 32 * Number.EPSILON * Math.max(topology.trailingEdge.gapLength, Math.abs(center.x), Math.abs(center.y));
  let distance = 0;
  for (const panel of topology.base.panels) {
    const midpoint = { x: .5 * (panel.start.x + panel.end.x), y: .5 * (panel.start.y + panel.end.y) };
    const departure = dot(sub(center, midpoint), panel.outwardNormal);
    if (departure > roundoff) {
      const projection = dot(direction, panel.outwardNormal);
      if (!(projection > 0)) throw new Error('No upstream finite-base source chart exists for this segment orientation.');
      distance = Math.max(distance, (departure + roundoff) / projection);
    }
  }
  const origin = { x: center.x - distance * direction.x, y: center.y - distance * direction.y };
  // When a shift is needed, prove all connectors stay inside the actual
  // solid, not merely the base half-planes. Unshifted straight bases can have
  // connectors on the solid boundary and retain their original exact gauge.
  if (distance > 0 && !pointInside(origin, topology.points)) throw new Error('Finite-base source-chart origin lies outside the retained solid.');
  for (const panel of topology.base.panels) {
    const midpoint = { x: .5 * (panel.start.x + panel.end.x), y: .5 * (panel.start.y + panel.end.y) };
    if (dot(sub(origin, midpoint), panel.outwardNormal) > roundoff)
      throw new Error('Finite-base source-chart origin cannot satisfy all retained base half-planes.');
    for (let i = 1; i < topology.points.length; i++) {
      const a = topology.points[i - 1], b = topology.points[i];
      // Strict crossing only: the connector terminates ON its base panel,
      // and native straight-base connectors can lie on that boundary.
      if (side(origin, midpoint, a) * side(origin, midpoint, b) === -1 && side(a, b, origin) * side(a, b, midpoint) === -1)
        throw new Error('Finite-base source-chart connector crosses the retained solid boundary.');
    }
    if (distance > 0 && !pointInside({ x: .5 * (origin.x + midpoint.x), y: .5 * (origin.y + midpoint.y) }, topology.points))
      throw new Error('Finite-base source-chart connector lies outside the retained solid.');
  }
  return { origin, direction: { ...direction }, upstreamShift: distance, originalCenter: { ...center },
    geometryUnchanged: true, definition: 'Minimum upstream half-plane shift of logarithm gauge; common exterior ray unchanged; connectors checked against actual solid.' };
}

// Integral log(z-s) ds in panel-local coordinates, with exact endpoint
// limits and a convergent far-field series to avoid short-panel cancellation.
function logarithmIntegral(point, panel) {
  const d = sub(point, panel.a), endpoint = same(point, panel.a) || same(point, panel.b);
  const x = same(point, panel.a) ? 0 : same(point, panel.b) ? 1 : (d.x * panel.tx + d.y * panel.ty) / panel.length;
  const y = endpoint ? 0 : (-d.x * panel.ty + d.y * panel.tx) / panel.length;
  const zx = x - .5, r = Math.hypot(zx, y);
  let real, imaginary;
  if (r >= 1) {
    const r2 = zx * zx + y * y, ir = zx / r2, ii = -y / r2;
    real = Math.log(panel.length) + Math.log(r); imaginary = Math.atan2(y, zx);
    let re = ir, im = ii;
    for (let k = 1; k <= 48; k++) {
      real -= re * moments[k] / k; imaginary -= im * moments[k] / k;
      const next = re * ir - im * ii; im = re * ii + im * ir; re = next;
    }
  } else {
    const ay = Math.abs(y), primitive = z => {
      const r2 = z * z + y * y;
      return (r2 === 0 ? 0 : z * Math.log(r2)) - 2 * z + (ay === 0 ? 0 : 2 * ay * Math.atan(z / ay));
    };
    real = Math.log(panel.length) + .5 * (primitive(x) - primitive(x - 1));
    imaginary = x * Math.atan2(y, x) - (x - 1) * Math.atan2(y, x - 1)
      + (y === 0 ? 0 : .5 * y * Math.log((x * x + y * y) / ((x - 1) ** 2 + y * y)));
  }
  return { real: panel.length * real, imaginary: panel.length * imaginary, x, y,
    midpointAngle: Math.atan2(y, zx) };
}

export function basePanelOnSourceCut(point, panel) {
  if (same(point, panel.a) || same(point, panel.b)) return false;
  const origin = panel.cutOrigin ?? { x: .5 * (panel.a.x + panel.b.x), y: .5 * (panel.a.y + panel.b.y) };
  const p = sub(point, origin), t = panel.cutDirection;
  const roundoff = 16 * Number.EPSILON * Math.max(panel.length, Math.abs(point.x), Math.abs(point.y), Math.abs(origin.x), Math.abs(origin.y));
  return dot(t, p) >= 0 && Math.abs(cross(t, p)) <= roundoff;
}

export function basePanelStreamfunctionBasis(point, panel) {
  const z = logarithmIntegral(point, panel), t = panel.cutDirection;
  const relative = { x: point.x - .5 * (panel.a.x + panel.b.x), y: point.y - .5 * (panel.a.y + panel.b.y) };
  const origin = panel.cutOrigin ?? { x: .5 * (panel.a.x + panel.b.x), y: .5 * (panel.a.y + panel.b.y) };
  const common = sub(point, origin);
  let angle = Math.atan2(cross(t, common), dot(t, common));
  if (angle < 0) angle += 2 * Math.PI;
  // This integral branch is continuous outside the sheet and the connecting
  // segment from its midpoint to the common origin. That connecting segment
  // is inside the solid base cap; only the common downstream ray is exterior.
  angle += Math.atan2(cross(common, relative), dot(common, relative));
  const source = (z.imaginary + panel.length * (angle - z.midpointAngle)) / (2 * Math.PI);
  return { source, vortex: -z.real / (2 * Math.PI), onSourceCut: basePanelOnSourceCut(point, panel) };
}

export function basePanelStreamfunction(point, panel) {
  const b = basePanelStreamfunctionBasis(point, panel);
  if (panel.sourceStrength !== 0 && b.onSourceCut) throw new Error('Streamfunction point lies on the finite-base downstream source-cut ray.');
  return panel.sourceStrength * b.source + panel.vortexStrength * b.vortex;
}

export function basePanelPotentialDifference(a, b, panel) {
  const za = logarithmIntegral(a, panel), zb = logarithmIntegral(b, panel);
  if (za.y * zb.y < 0) {
    const t = -za.y / (zb.y - za.y), x = za.x + t * (zb.x - za.x);
    if (x >= 0 && x <= 1) throw new Error('Potential-difference path crosses a finite-base sheet.');
  }
  const ax = za.x - .5, bx = zb.x - .5;
  const raw = zb.midpointAngle - za.midpointAngle;
  const unwrapped = Math.atan2(ax * zb.y - za.y * bx, ax * bx + za.y * zb.y);
  const turns = Math.round((unwrapped - raw) / (2 * Math.PI));
  return (panel.sourceStrength * (zb.real - za.real)
    + panel.vortexStrength * (zb.imaginary - za.imaginary + turns * 2 * Math.PI * panel.length)) / (2 * Math.PI);
}

// Literal SPLIND(...,-999,-999): zero third derivative at both endpoints.
// This TE derivative definition does not change or fit the solid geometry.
export function xfoilSurfaceDerivatives(points) {
  const s = [0];
  for (let i = 1; i < points.length; i++) s.push(s.at(-1) + Math.hypot(points[i].x - points[i - 1].x, points[i].y - points[i - 1].y));
  const n = points.length;
  if (n < 3 || s.some((v, i) => i && !(v > s[i - 1]))) throw new Error('Finite-base surface spline requires distinct ordered points.');
  const solve = key => {
    const a = new Float64Array(n), b = new Float64Array(n), c = new Float64Array(n), r = new Float64Array(n);
    for (let i = 1; i < n - 1; i++) {
      const before = s[i] - s[i - 1], after = s[i + 1] - s[i];
      b[i] = after; a[i] = 2 * (before + after); c[i] = before;
      r[i] = 3 * ((points[i + 1][key] - points[i][key]) * before / after
        + (points[i][key] - points[i - 1][key]) * after / before);
    }
    a[0] = c[0] = 1; r[0] = 2 * (points[1][key] - points[0][key]) / (s[1] - s[0]);
    b[n - 1] = a[n - 1] = 1; r[n - 1] = 2 * (points[n - 1][key] - points[n - 2][key]) / (s[n - 1] - s[n - 2]);
    for (let i = 1; i < n; i++) { const factor = b[i] / a[i - 1]; a[i] -= factor * c[i - 1]; r[i] -= factor * r[i - 1]; }
    r[n - 1] /= a[n - 1];
    for (let i = n - 2; i >= 0; i--) r[i] = (r[i] - c[i] * r[i + 1]) / a[i];
    return r;
  };
  const x = solve('x'), y = solve('y');
  return { s, derivatives: points.map((_, i) => ({ x: x[i], y: y[i] })) };
}
