// SPDX-License-Identifier: GPL-2.0-or-later
// Single-valued streamfunction of the existing linear vortex panel field.
// u = d(psi)/dy, v = -d(psi)/dx; a positive CCW point vortex contributes
// -Gamma log(r)/(2*pi). This is an incompressible grid-initialization tool.
import { basePanelStreamfunction, basePanelPotentialDifference } from './finite-base-influence.js';
const moments = Array.from({ length: 50 }, (_, k) => k % 2 ? 0 : 1 / ((k + 1) * 2 ** k));
const xlogx = x => x === 0 ? 0 : x * Math.log(x);

export function vortexStreamfunctionBasis(point, panel) {
  const { length, tx, ty } = panel;
  if (![point.x, point.y, panel.a.x, panel.a.y, length, tx, ty].every(Number.isFinite) || !(length > 0))
    throw new Error('Invalid vortex streamfunction geometry.');
  const dx = point.x - panel.a.x, dy = point.y - panel.a.y;
  const x = (dx * tx + dy * ty) / length, y = (-dx * ty + dy * tx) / length;
  const zx = x - .5, radius = Math.hypot(zx, y);
  if (radius >= 1) {
    // log|z-eta| = log|z| - Re sum (eta/z)^k/k, eta in [-1/2,1/2].
    // Far-field primitive subtraction loses the tiny first moment of short
    // panels. Here |eta/z| <= 1/2; 48 terms bound the omitted integral below
    // roundoff at the branch boundary, and converge faster farther away.
    const r2 = zx * zx + y * y, inverse = { re: zx / r2, im: -y / r2 };
    let re = inverse.re, im = inverse.im;
    const integral = [.5 * (Math.log(length) + Math.log(radius)), .5 * (Math.log(length) + Math.log(radius))];
    for (let k = 1; k <= 48; k++) {
      integral[0] -= re / k * (.5 * moments[k] - moments[k + 1]);
      integral[1] -= re / k * (.5 * moments[k] + moments[k + 1]);
      const nextRe = re * inverse.re - im * inverse.im;
      im = re * inverse.im + im * inverse.re; re = nextRe;
    }
    return integral.map(v => -length * v / (2 * Math.PI));
  }
  // Exact log-distance primitives on and near the panel. Their finite
  // endpoint limits include the logarithmic singularity without epsilon
  // offsets, clipping the distance, or a special arbitrary streamfunction.
  const ay = Math.abs(y), primitive = z => {
    const r2 = z * z + y * y;
    return (r2 === 0 ? 0 : z * Math.log(r2)) - 2 * z + (ay === 0 ? 0 : 2 * ay * Math.atan(z / ay));
  };
  const zeroth = primitive(x) - primitive(x - 1), a2 = x * x + y * y, b2 = (x - 1) ** 2 + y * y;
  const first = x * zeroth - .5 * (xlogx(a2) - a2 - xlogx(b2) + b2);
  return [zeroth - first, first].map(v => -length * (Math.log(length) + v) / (4 * Math.PI));
}

export function streamfunctionAt(point, field) {
  let value = field.u * point.y - field.v * point.x, correction = 0;
  for (const panel of field.panels) {
    const basis = vortexStreamfunctionBasis(point, panel);
    for (let k = 0; k < 2; k++) {
      const term = field.gamma[panel.node + k] * basis[k] - correction, next = value + term;
      correction = (next - value) - term; value = next;
    }
  }
  for (const panel of field.basePanels ?? []) value += basePanelStreamfunction(point, panel);
  if (!Number.isFinite(value)) throw new Error('Nonfinite panel streamfunction.');
  return value;
}

// Velocity potential is multi-valued around a lifting body. Expose shape
// functions and path-local differences, not an arbitrary global scalar field.
export function vortexPotentialBasis(point, panel) {
  const { length, tx, ty } = panel;
  if (![point.x, point.y, panel.a.x, panel.a.y, length, tx, ty].every(Number.isFinite) || !(length > 0))
    throw new Error('Invalid vortex potential geometry.');
  const dx = point.x - panel.a.x, dy = point.y - panel.a.y;
  const atStart = point.x === panel.a.x && point.y === panel.a.y;
  const atEnd = point.x === panel.b.x && point.y === panel.b.y;
  // Use the same canonical branch at exact endpoints as potentialDifference
  // uses for angle unwrapping. For an upstream-pointing oblique panel, the
  // projection at its start can otherwise produce -0: atan2(-0, negative)
  // differs by 2*pi from atan2(+0, negative). That introduces a whole-panel
  // circulation jump at a contour knot, even though subsequent increments
  // cancel the error over a complete path. These are exact endpoint limits,
  // not displaced evaluation points or changes to off-sheet values.
  const x = atStart ? 0 : atEnd ? 1 : (dx * tx + dy * ty) / length;
  const y = atStart || atEnd ? 0 : (-dx * ty + dy * tx) / length, zx = x - .5;
  if (Math.hypot(zx, y) >= 1) {
    const r2 = zx * zx + y * y, inverse = { re: zx / r2, im: -y / r2 };
    let re = inverse.re, im = inverse.im;
    const integral = Array(2).fill(.5 * Math.atan2(y, zx));
    for (let k = 1; k <= 48; k++) {
      integral[0] -= im / k * (.5 * moments[k] - moments[k + 1]);
      integral[1] -= im / k * (.5 * moments[k] + moments[k + 1]);
      const nextRe = re * inverse.re - im * inverse.im;
      im = re * inverse.im + im * inverse.re; re = nextRe;
    }
    return integral.map(v => length * v / (2 * Math.PI));
  }
  const a2 = x * x + y * y, b2 = (x - 1) ** 2 + y * y, theta0 = Math.atan2(y, x), theta1 = Math.atan2(y, x - 1);
  const zeroth = x * theta0 - (x - 1) * theta1 + (y === 0 ? 0 : .5 * y * Math.log(a2 / b2));
  const first = x * zeroth - .5 * (a2 * theta0 - b2 * theta1 + y);
  return [zeroth - first, first].map(v => length * v / (2 * Math.PI));
}

export function potentialDifference(a, b, field) {
  let value = field.u * (b.x - a.x) + field.v * (b.y - a.y), correction = 0;
  for (const panel of field.panels) {
    const ax = (a.x - panel.a.x) * panel.tx + (a.y - panel.a.y) * panel.ty;
    const endpoint = p => (p.x === panel.a.x && p.y === panel.a.y) || (p.x === panel.b.x && p.y === panel.b.y);
    // Exact geometric endpoints have zero normal coordinate. A floating
    // tangent projection of an oblique endpoint can otherwise have a tiny
    // wrong sign and falsely classify an outgoing TE path as a sheet crossing.
    const ay = endpoint(a) ? 0 : -(a.x - panel.a.x) * panel.ty + (a.y - panel.a.y) * panel.tx;
    const bx = (b.x - panel.a.x) * panel.tx + (b.y - panel.a.y) * panel.ty;
    const by = endpoint(b) ? 0 : -(b.x - panel.a.x) * panel.ty + (b.y - panel.a.y) * panel.tx;
    if (ay * by < 0) {
      const t = -ay / (by - ay), x = ax + t * (bx - ax);
      if (x >= 0 && x <= panel.length) throw new Error('Potential-difference path crosses a vortex sheet.');
    }
    const amx = ax - .5 * panel.length, bmx = bx - .5 * panel.length;
    const raw = Math.atan2(by, bmx) - Math.atan2(ay, amx);
    const unwrapped = Math.atan2(amx * by - ay * bmx, amx * bmx + ay * by);
    const turns = Math.round((unwrapped - raw) / (2 * Math.PI));
    const pa = vortexPotentialBasis(a, panel), pb = vortexPotentialBasis(b, panel);
    for (let k = 0; k < 2; k++) {
      const term = field.gamma[panel.node + k] * (pb[k] - pa[k] + turns * panel.length / 2) - correction;
      const next = value + term; correction = (next - value) - term; value = next;
    }
  }
  for (const panel of field.basePanels ?? []) value += basePanelPotentialDifference(a, b, panel);
  if (!Number.isFinite(value)) throw new Error('Nonfinite panel potential difference.');
  return value;
}

// Change of source-log sheet along an oriented exterior path. Subtract this
// jump to unwrap sampled psi, or add it to a body's constant to express the
// same dividing streamline in the inlet chart. Velocity is single-valued;
// the jump is exactly the base source flux, not an integration-drift estimate.
export function streamfunctionBranchIncrement(path, field) {
  let increment = 0;
  for (const panel of field.basePanels ?? []) {
    if (!panel.sourceStrength) continue;
    const d = panel.cutDirection, o = panel.cutOrigin ?? { x: .5 * (panel.a.x + panel.b.x), y: .5 * (panel.a.y + panel.b.y) };
    for (let i = 1; i < path.length; i++) {
      const a = path[i - 1], b = path[i];
      const sa = d.x * (a.y - o.y) - d.y * (a.x - o.x);
      const sb = d.x * (b.y - o.y) - d.y * (b.x - o.x);
      if (!((sa < 0 && sb >= 0) || (sb < 0 && sa >= 0))) continue;
      const f = sa / (sa - sb), x = a.x + f * (b.x - a.x), y = a.y + f * (b.y - a.y);
      if (d.x * (x - o.x) + d.y * (y - o.y) <= 0) continue;
      increment += (sa < 0 ? -1 : 1) * panel.sourceStrength * panel.length;
    }
  }
  return increment;
}
