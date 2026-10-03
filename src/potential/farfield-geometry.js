// SPDX-License-Identifier: GPL-2.0-or-later
// Spatial derivatives of the manual's velocity expansion, by an exact
// two-coordinate forward chain rule. The value implementation remains an
// independent reference in farfield.js; no finite differences are used here.
import { multipoleBasis } from './farfield.js';
const add = (a, b) => a.map((v, i) => v + b[i]);
const sub = (a, b) => a.map((v, i) => v - b[i]);
const scale = (a, s) => a.map(v => v * s);
const mul = (a, b) => [a[0] * b[0], a[1] * b[0] + a[0] * b[1], a[2] * b[0] + a[0] * b[2]];
const div = (a, b) => { const q = a[0] / b[0]; return [q, (a[1] - q * b[1]) / b[0], (a[2] - q * b[2]) / b[0]]; };
const log = a => [Math.log(a[0]), a[1] / a[0], a[2] / a[0]];

export function multipoleGeometryDerivatives(point, settings = {}) {
  const value = multipoleBasis(point, settings), { center = { x: 0, y: 0 }, alpha = 0, mach = 0, gamma = 1.4 } = settings;
  const angle = alpha * Math.PI / 180, c = Math.cos(angle), s = Math.sin(angle), beta = Math.sqrt(1 - mach * mach), twoPi = 2 * Math.PI;
  const x = [((point.x - center.x) * c + (point.y - center.y) * s) / beta, c / beta, s / beta];
  const y = [-(point.x - center.x) * s + (point.y - center.y) * c, -s, c];
  const xx = mul(x, x), yy = mul(y, y), xy = mul(x, y), r2 = add(xx, yy), r4 = mul(r2, r2), r6 = mul(r4, r2);
  const logr = scale(log(r2), .5), cubic = mul(x, sub(xx, scale(yy, 3)));
  const a = .25 * ((3 - gamma) / beta + (gamma + 1) / beta ** 3), b = (gamma + 1) / 16 * (1 / beta - 1 / beta ** 3);
  const gradients = [
    [scale(div(y, r2), -1 / twoPi), scale(div(x, r2), 1 / twoPi)],
    [scale(div(x, r2), 1 / twoPi), scale(div(y, r2), 1 / twoPi)],
    [scale(div(sub(yy, xx), r4), 1 / twoPi), scale(div(xy, r4), -2 / twoPi)],
    [scale(div(xy, r4), -2 / twoPi), scale(div(sub(xx, yy), r4), 1 / twoPi)],
    [scale(add(scale(div(add(xx, mul(logr, sub(yy, xx))), r4), a),
      scale(sub(div(scale(sub(xx, yy), 3), r4), div(scale(mul(x, cubic), 4), r6)), b)), (mach / twoPi) ** 2),
    scale(add(scale(div(mul(xy, sub([1, 0, 0], scale(logr, 2))), r4), a),
      scale(sub(div(scale(xy, -6), r4), div(scale(mul(y, cubic), 4), r6)), b)), (mach / twoPi) ** 2)] ];
  const derivatives = gradients.map(([gx, gy]) => {
    const vx = sub(scale(gx, c / beta), scale(gy, s)), vy = add(scale(gx, s / beta), scale(gy, c));
    return [[vx[1], vx[2]], [vy[1], vy[2]]];
  });
  if (!derivatives.flat(2).every(Number.isFinite)) throw new Error('Nonfinite farfield geometry derivatives.');
  return { ...value, velocityDerivatives: derivatives };
}
