// SPDX-License-Identifier: GPL-2.0-or-later
// Analytic geometry/fields only. No production solver or gas-kernel imports.
import { ringlebPoint, ringlebAtPoint } from './ringleb.js';

const average = points => ({ x: points.reduce((s, p) => s + p.x, 0) / points.length,
  y: points.reduce((s, p) => s + p.y, 0) / points.length });

// Reversing all velocities is an exact steady Euler symmetry. Swapping banks
// at the same time preserves positive cell orientation and positive tube mass.
export function gradedRinglebTube({ mach, h, widthRatio, skew = false, direction = 1, grading = 8, spacing = 'smooth' }) {
  if (!(mach > 0 && h > 0 && widthRatio > 0) || ![1, -1].includes(direction)
    || !Number.isFinite(grading) || !['smooth', 'fixed-ratio'].includes(spacing)) throw new Error('Invalid graded Ringleb controls.');
  const qCenter = mach / Math.sqrt(1 + .2 * mach * mach), psiCenter = skew ? .985 / qCenter : .7;
  const qs = [3, 2, 1, 0, -1].map(k => qCenter + direction * (spacing === 'smooth'
    ? k * h + grading * (k * h) ** 2 / qCenter : k * h * (k > 0 ? 1.5 : 1)));
  const psiWidth = widthRatio * h, lowerPsi = psiCenter - direction * psiWidth / 2, upperPsi = psiCenter + direction * psiWidth / 2;
  if (qs.some((q, i) => i && direction * (qs[i - 1] - q) <= 0)) throw new Error('Nonmonotone Ringleb stations.');
  const lower = qs.map(q => ringlebPoint(q, lowerPsi)), upper = qs.map(q => ringlebPoint(q, upperPsi));
  const exactAt = (p, guess = { q: qCenter, psi: psiCenter }) => {
    const s = ringlebAtPoint(p.x, p.y, { guess });
    return { ...s, u: direction * s.u, v: direction * s.v };
  };
  const exactSections = qs.slice(1).map((q, i) => exactAt(average([lower[i], lower[i + 1], upper[i], upper[i + 1]]),
    { q: (q + qs[i]) / 2, psi: psiCenter }));
  return { lower, upper, densities: exactSections.map(s => s.rho), exactSections, exactAt,
    massFlow: psiWidth, gamma: 1.4, stagnationEnthalpy: 2.5,
    center: ringlebPoint(qCenter, psiCenter), controls: { mach, h, widthRatio, skew, direction, grading, spacing, qCenter, psiCenter, psiWidth } };
}

export function vortexField({ x, y }, { gamma = 1.4, mach = .3 } = {}) {
  const r = Math.hypot(x, y), q = 1 / r, t = 1 + .5 * (gamma - 1) * mach ** 2 * (1 - q * q);
  if (!(r > 0 && t > 0)) throw new Error('Vortex field outside positive-enthalpy region.');
  const rho = t ** (1 / (gamma - 1)), p = t ** (gamma / (gamma - 1)) / (gamma * mach * mach);
  const enthalpy = t / ((gamma - 1) * mach * mach);
  return { x, y, r, q, u: -q * y / r, v: q * x / r, rho, p, enthalpy,
    mach: q / Math.sqrt((gamma - 1) * enthalpy), gamma, h0: .5 + 1 / ((gamma - 1) * mach * mach) };
}

export function simpsonIntegral(f, a, b, intervals = 64) {
  if (!Number.isInteger(intervals) || intervals <= 0 || intervals % 2) throw new Error('Even positive integration count required.');
  let sum = 0;
  for (let i = 0; i <= intervals; i++) sum += (i === 0 || i === intervals ? 1 : i % 2 ? 4 : 2) * f(a + (b - a) * i / intervals);
  return sum * (b - a) / (3 * intervals);
}

export function gradedVortexTube({ h, widthRatio = 1, radius = 1, center = .37, grading = .4 }) {
  const theta = t => .7 * t + grading * t * t, ts = [-h, 0, h].map(t => theta(center + t));
  const inner = radius - widthRatio * h / 2, outer = radius + widthRatio * h / 2;
  const point = (r, t) => ({ x: r * Math.cos(t), y: r * Math.sin(t) });
  const lower = ts.map(t => point(outer, t)), upper = ts.map(t => point(inner, t));
  const exactSections = [0, 1].map(i => vortexField(average([lower[i], lower[i + 1], upper[i], upper[i + 1]])));
  const mass = intervals => simpsonIntegral(r => { const s = vortexField({ x: r, y: 0 }); return s.rho * s.q; }, inner, outer, intervals);
  return { lower, upper, densities: exactSections.map(s => s.rho), exactSections, exactAt: vortexField,
    massFlow: mass(128), massQuadratureDifference: Math.abs(mass(128) - mass(64)), gamma: 1.4,
    stagnationEnthalpy: exactSections[0].h0, angles: ts, controls: { h, widthRatio, radius, center, grading } };
}

// Independent physical-face integration. Positive polygon orientation gives
// outward normal measure (dy,-dx); exact field velocities may cross chord walls.
export function exactEulerPolygonFlux(points, exactAt, intervals = 32) {
  if (!Number.isInteger(intervals) || intervals <= 0 || intervals % 2) throw new Error('Even positive integration count required.');
  const faces = points.map((a, i) => {
    const b = points[(i + 1) % points.length], nx = b.y - a.y, ny = a.x - b.x;
    const value = t => {
      const s = exactAt({ x: a.x + t * (b.x - a.x), y: a.y + t * (b.y - a.y) });
      const mass = s.rho * (s.u * nx + s.v * ny);
      return [mass, mass * s.u + s.p * nx, mass * s.v + s.p * ny, mass * (s.enthalpy + .5 * (s.u ** 2 + s.v ** 2))];
    };
    const flux = [0, 0, 0, 0];
    for (let k = 0; k <= intervals; k++) {
      const weight = (k === 0 || k === intervals ? 1 : k % 2 ? 4 : 2) / (3 * intervals);
      value(k / intervals).forEach((v, j) => { flux[j] += weight * v; });
    }
    return flux;
  });
  return { faces, total: [0, 1, 2, 3].map(k => faces.reduce((s, f) => s + f[k], 0)) };
}

// Radial equilibrium gives dp/dr=rho/r^3. For any oriented face a -> b,
// the integrated outward momentum is G(a)-G(b), G=p(r)*(-y,x).
// Unlike the quadrature oracle above, this identity uses endpoints only.
export function exactVortexMomentumFace(a, b) {
  const pa = vortexField(a).p, pb = vortexField(b).p;
  return { x: -pa * a.y + pb * b.y, y: pa * a.x - pb * b.x };
}
