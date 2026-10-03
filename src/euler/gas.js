// SPDX-License-Identifier: GPL-2.0-or-later
// Calorically perfect gas. Coordinates, density and velocity use c_ref, rho_inf, U_inf.
export function freestream({ mach = 0.3, alpha = 0, gamma = 1.4 } = {}) {
  if (!Number.isFinite(mach) || mach <= 0 || !Number.isFinite(alpha) || !Number.isFinite(gamma) || gamma <= 1) throw new Error('Euler requires positive Mach and gamma > 1.');
  const angle = alpha * Math.PI / 180;
  return { rho: 1, u: Math.cos(angle), v: Math.sin(angle), p: 1 / (gamma * mach ** 2), gamma };
}
export function conserved(s, gamma = s.gamma ?? 1.4) {
  if (![s.rho, s.u, s.v, s.p, gamma].every(Number.isFinite) || s.rho <= 0 || s.p <= 0 || gamma <= 1) throw new Error('Nonphysical Euler state.');
  return [s.rho, s.rho * s.u, s.rho * s.v, s.p / (gamma - 1) + 0.5 * s.rho * (s.u ** 2 + s.v ** 2)];
}
export function totalConditions(s, gamma = s.gamma ?? 1.4) {
  const speed2 = s.u ** 2 + s.v ** 2;
  const a2 = gamma * s.p / s.rho;
  const temperatureRatio = 1 + 0.5 * (gamma - 1) * speed2 / a2;
  return { h0: a2 / (gamma - 1) + speed2 / 2,
    p0: s.p * temperatureRatio ** (gamma / (gamma - 1)), entropy: Math.log(s.p) - gamma * Math.log(s.rho) };
}
export function physicalFlux(s, nx, ny, gamma = s.gamma ?? 1.4) {
  const q = conserved(s, gamma); const un = s.u * nx + s.v * ny;
  return [q[0] * un, q[1] * un + s.p * nx, q[2] * un + s.p * ny, (q[3] + s.p) * un];
}

// HLLC: Toro, Spruce & Speares (1994), with Davis min/max acoustic bounds.
// A single oriented flux is shared by both cells; it is never evaluated twice.
export function hllc(left, right, nx, ny, gamma = 1.4) {
  if (!Number.isFinite(nx) || !Number.isFinite(ny) || Math.abs(Math.hypot(nx, ny) - 1) > 1e-10) throw new Error('Flux needs a unit face normal.');
  const l = conserved(left, gamma); const r = conserved(right, gamma);
  const ul = left.u * nx + left.v * ny; const ur = right.u * nx + right.v * ny;
  const al = Math.sqrt(gamma * left.p / left.rho); const ar = Math.sqrt(gamma * right.p / right.rho);
  const sl = Math.min(ul - al, ur - ar); const sr = Math.max(ul + al, ur + ar);
  const fl = physicalFlux(left, nx, ny, gamma); const fr = physicalFlux(right, nx, ny, gamma);
  if (sl >= 0) return fl;
  if (sr <= 0) return fr;
  const dl = left.rho * (sl - ul); const dr = right.rho * (sr - ur);
  const contact = (right.p - left.p + dl * ul - dr * ur) / (dl - dr);
  const s = contact >= 0 ? left : right; const q = contact >= 0 ? l : r;
  const f = contact >= 0 ? fl : fr; const speed = contact >= 0 ? sl : sr; const un = contact >= 0 ? ul : ur;
  const rhoStar = s.rho * (speed - un) / (speed - contact);
  const star = [rhoStar, rhoStar * (s.u + (contact - un) * nx), rhoStar * (s.v + (contact - un) * ny),
    rhoStar * (q[3] / s.rho + (contact - un) * (contact + s.p / (s.rho * (speed - un))))];
  const flux = f.map((v, k) => v + speed * (star[k] - q[k]));
  const internalEnergy = star[3] - (star[1] ** 2 + star[2] ** 2) / (2 * rhoStar);
  if (!flux.every(Number.isFinite) || rhoStar <= 0 || !(internalEnergy > 0)) throw new Error('Inadmissible HLLC intermediate state.');
  return flux;
}
