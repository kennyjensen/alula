// SPDX-License-Identifier: GPL-2.0-or-later
// Shared momentum speeds along one physical tube, with local directional
// derivatives. Section k occupies the bank segment k -> k+1. Body cuts,
// wall intervals and wakes do not interrupt this downstream sequence.
import { streamtubeCellGeometry, streamtubeSection } from './streamtube-cell.js';
import { evaluateStreamtubeSpeedUpwind, linearizeStreamtubeSpeedUpwind } from './streamtube-speed-upwind.js';

const zero = Object.freeze({ x: 0, y: 0 });
const zeroPair = Object.freeze([zero, zero]);
const sub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });
const mean = (a, b) => ({ x: .5 * (a.x + b.x), y: .5 * (a.y + b.y) });
const dot = (a, b) => a.x * b.x + a.y * b.y;
const cross = (a, b) => a.x * b.y - a.y * b.x;
const point = p => p && Number.isFinite(p.x) && Number.isFinite(p.y);

export function prepareStreamtubeTransportChain({ lower, upper, densities, massFlow,
  stagnationEnthalpy, gamma = 1.4, geometryDomain = 'convex', upwind }, { linearize = false } = {}) {
  if (!Array.isArray(lower) || lower.length < 3 || !Array.isArray(upper) || upper.length !== lower.length
    || !lower.every(point) || !upper.every(point) || !densities || densities.length !== lower.length - 1)
    throw new Error('Supply matching physical banks and section densities for the transport chain.');
  if (!upwind || upwind.boundary?.kind !== 'unfiltered-first-two'
    || Object.keys(upwind.boundary).some(key => key !== 'kind'))
    throw new Error('Transport chain requires an explicit unfiltered-first-two inlet closure.');
  const { mucon = 1, mcrit = .99 } = upwind;
  if (![mucon, mcrit, gamma].every(Number.isFinite) || mcrit < 0 || mcrit > 1 || gamma <= 1)
    throw new Error('Invalid transport-chain upwind controls.');
  lower = lower.map(p => ({ ...p })); upper = upper.map(p => ({ ...p }));
  const n = densities.length, sections = Array(n), sectionGeometry = Array(n);
  // Do not evaluate undissipated cell pressures here: their admissibility
  // is not a prerequisite for the actual biased momentum equations.
  for (let i = 1; i < n; i++) {
    let geometry;
    try { geometry = streamtubeCellGeometry(lower.slice(i - 1, i + 2), upper.slice(i - 1, i + 2), { geometryDomain }); }
    catch (error) {
      // A stencil covers two primal intervals and their staggered volume.
      // Preserve the actual rejected six points; a later accepted grid cannot
      // identify which of those trial polygons failed. No geometry is repaired.
      throw Object.assign(error, { code: error.code ?? 'streamtube-transport-geometry', diagnostics: {
        ...error.diagnostics, stencil: { i, bankStations: [i - 1, i, i + 1], geometryDomain,
          lower: lower.slice(i - 1, i + 2).map(p => ({ ...p })),
          upper: upper.slice(i - 1, i + 2).map(p => ({ ...p })) },
      } });
    }
    for (const slot of i === 1 ? [0, 1] : [1]) {
      const k = i - 1 + slot;
      try { sections[k] = streamtubeSection({ density: densities[k], massFlow,
        normalArea: geometry.normalAreas[slot], stagnationEnthalpy, gamma }); }
      catch (error) {
        if (error.code === 'streamtube-static-enthalpy') error.diagnostics = { ...error.diagnostics,
          section: { i: k, bankStations: [k, k + 1] } };
        throw error;
      }
      sectionGeometry[k] = { normalArea: geometry.normalAreas[slot], length: geometry.streamwiseLengths[slot],
        direction: geometry.directions[slot], gap: geometry.sections[slot] };
    }
  }
  let arc = 0;
  const sectionArc = sectionGeometry.map(g => { const s = arc + .5 * g.length; arc += g.length; return s; });
  const transportSpeeds = sections.map(s => s.q), coefficients = Array(n).fill(0), corrections = Array(n).fill(0);
  const linears = linearize ? Array(n).fill(null) : null;
  for (let k = 2; k < n; k++) {
    const parameters = { speeds: sections.slice(k - 2, k + 1).map(s => s.q),
      machSquared: sections.slice(k - 1, k + 1).map(s => s.machSquared),
      // These are differences of midpoint arc coordinates. Form them
      // locally so long upstream arcs cannot erase a short interval.
      spacing: [.5 * (sectionGeometry[k - 2].length + sectionGeometry[k - 1].length),
        .5 * (sectionGeometry[k - 1].length + sectionGeometry[k].length)], mucon, mcrit, gamma };
    const local = linearize ? linearizeStreamtubeSpeedUpwind(parameters) : null;
    const value = local ? local.value : evaluateStreamtubeSpeedUpwind(parameters);
    transportSpeeds[k] = value.speed; coefficients[k] = value.coefficient; corrections[k] = value.correction;
    if (linears) linears[k] = local;
  }
  const result = { sections, sectionGeometry, sectionArc, transportSpeeds,
    upwind: { coefficients, corrections, filtered: sections.map((_, k) => k >= 2),
      secondOrder: mucon >= 0, boundary: 'unfiltered-first-two' } };
  if (!linearize) return result;
  const index = k => { if (!Number.isInteger(k) || k < 0 || k >= n) throw new Error('Invalid transport-chain section index.'); };
  const sectionTangent = (k, { lower: dl = zeroPair, upper: du = zeroPair, density: rho = 0,
    massFlow: dm = 0, stagnationEnthalpy: dh0 = 0 } = {}) => {
    index(k);
    if (![dl, du].every(row => Array.isArray(row) && row.length === 2 && row.every(point))
      || ![rho, dm, dh0].every(Number.isFinite)) throw new Error('Invalid transport-chain section tangent.');
    const g = sectionGeometry[k], s = sections[k];
    const dv = mean(sub(dl[1], dl[0]), sub(du[1], du[0]));
    const dg = mean(sub(du[0], dl[0]), sub(du[1], dl[1]));
    const length = dot(g.direction, dv);
    const dt = { x: (dv.x - g.direction.x * length) / g.length, y: (dv.y - g.direction.y * length) / g.length };
    const normalArea = cross(dt, g.gap) + cross(g.direction, dg);
    const chain = d => d.density * rho + d.massFlow * dm + d.normalArea * normalArea + d.stagnationEnthalpy * dh0;
    const q = chain(s.derivatives.q), p = chain(s.derivatives.p), enthalpy = dh0 - s.q * q;
    const machSquared = (2 * s.q * q - s.machSquared * (gamma - 1) * enthalpy) / ((gamma - 1) * s.enthalpy);
    const tangent = { rho, q, p, enthalpy, machSquared, length, normalArea };
    if (!Object.values(tangent).every(Number.isFinite)) throw new Error('Nonfinite transport-chain section tangent.');
    return tangent;
  };
  const transportTangent = (k, sectionTangentAt) => {
    index(k);
    if (typeof sectionTangentAt !== 'function') throw new Error('Supply a local section-tangent callback.');
    if (k < 2) {
      const q = sectionTangentAt(k).q;
      if (!Number.isFinite(q)) throw new Error('Nonfinite transport-chain speed tangent.');
      return q;
    }
    const a = sectionTangentAt(k - 2), b = sectionTangentAt(k - 1), c = sectionTangentAt(k);
    return linears[k].apply({ speeds: [a.q, b.q, c.q], machSquared: [b.machSquared, c.machSquared],
      spacing: [.5 * (a.length + b.length), .5 * (b.length + c.length)] }).speed;
  };
  return { ...result, sectionTangent, transportTangent };
}
