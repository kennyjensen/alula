// SPDX-License-Identifier: GPL-2.0-or-later
// Drela (1986), Eqs. 6.61–6.63: surface-streamtube speed interpolation
// with the signed geometry sawtooth correction. The scalar cross-product
// convention is defined on thesis p.24; the bars do NOT denote absolute value.
// The existing geometry uses successive node-to-node streamline vectors.
// Optional distance weighting is our nonuniform-grid reconstruction, not
// the literal thesis mean. It leaves the printed sawtooth correction intact.
const weights = (cell, interpolation) => {
  if (interpolation === 'arithmetic') return [.5, .5];
  if (interpolation !== 'distance-weighted') throw new Error('Unknown edge-speed interpolation.');
  const lengths = cell.geometry.streamwiseLengths;
  if (!Array.isArray(lengths) || lengths.length !== 2 || !lengths.every(v => Number.isFinite(v) && v > 0))
    throw new Error('Distance-weighted edge speed requires positive streamwise lengths.');
  const [left, right] = lengths, total = left + right;
  return [right / total, left / total];
};
export function streamtubeEdgeVelocity(cell, correctionFactor = .2, interpolation = 'arithmetic') {
  const [a, b] = cell.states, q = .5 * (a.q + b.q);
  const w = weights(cell, interpolation), interpolated = interpolation === 'arithmetic' ? q : w[0] * a.q + w[1] * b.q;
  const machSquared = .5 * ((a.machSquared ?? 0) + (b.machSquared ?? 0));
  const curvature = cell.geometry.pressureCurvature;
  if (![q, machSquared, curvature, correctionFactor].every(Number.isFinite)
    || q <= 0 || machSquared < 0 || correctionFactor < 0)
    throw new Error('Invalid section-speed matching state.');
  const correction = machSquared < 1 ? correctionFactor * q * machSquared * (machSquared - 1) * curvature : 0;
  const ue = interpolated + correction;
  if (!(ue > 0) || !Number.isFinite(ue)) throw new Error('Nonpositive corrected section edge speed.');
  return { ue, meanSpeed: q, meanMachSquared: machSquared, correction, curvature,
    ...(interpolation === 'arithmetic' ? {} : { interpolation, interpolatedSpeed: interpolated, weights: w }) };
}

// Exact chain through both section speeds/Mach numbers and every geometry
// coordinate. Sonic switching is tested one-sided, not differentiated across.
export function streamtubeEdgeVelocityTangent(cell, tangent, correctionFactor = .2, interpolation = 'arithmetic') {
  const w = weights(cell, interpolation);
  const q = .5 * (cell.states[0].q + cell.states[1].q);
  const m = .5 * ((cell.states[0].machSquared ?? 0) + (cell.states[1].machSquared ?? 0));
  const dq = .5 * (tangent.states[0].q + tangent.states[1].q);
  const dm = .5 * ((tangent.states[0].machSquared ?? 0) + (tangent.states[1].machSquared ?? 0));
  let interpolatedDerivative = dq;
  if (interpolation === 'distance-weighted') {
    const [dl, dr] = tangent.geometry.streamwiseLengths, [left, right] = cell.geometry.streamwiseLengths;
    const dw = (w[1] * dr - w[0] * dl) / (left + right);
    interpolatedDerivative = w[0] * tangent.states[0].q + w[1] * tangent.states[1].q
      + dw * (cell.states[0].q - cell.states[1].q);
  }
  const c = cell.geometry.pressureCurvature, dc = tangent.geometry.pressureCurvature;
  return interpolatedDerivative + (m < 1 ? correctionFactor * (dq * m * (m - 1) * c
    + q * (2 * m - 1) * dm * c + q * m * (m - 1) * dc) : 0);
}

// Pressure matching to a smooth, uniform-entropy BL edge. This explicitly
// excludes entropy jumps/shocks; it is not an entropy-aware post-shock map.
export function streamtubeEdgeState(pressure, { flowModel, mach, gamma, pInf }) {
  if (!Number.isFinite(pressure)) throw new Error('Nonfinite BL edge pressure.');
  let q2, rho, temperature, machSquared;
  if (flowModel === 'incompressible') {
    q2 = -2 * pressure; rho = 1; machSquared = 0;
  } else {
    if (!(pressure > 0)) throw new Error('Nonpositive BL edge pressure.');
    const logarithm = Math.log1p((pressure - pInf) / pInf);
    temperature = Math.exp((gamma - 1) / gamma * logarithm);
    q2 = 1 - 2 * Math.expm1((gamma - 1) / gamma * logarithm) / ((gamma - 1) * mach * mach);
    rho = Math.exp(logarithm / gamma); machSquared = mach * mach * q2 / temperature;
  }
  if (!(q2 > 0) || !Number.isFinite(q2) || !(machSquared < 1))
    throw new Error('BL edge pressure implies stagnating, reversed or sonic flow at a resolved station.');
  const ue = Math.sqrt(q2);
  return { ue, rho, machSquared, pressure, velocityPressureDerivative: -1 / (rho * ue) };
}

// Direct isentropic pressure matching avoids taking sqrt(q^2) of an
// off-solution Euler pressure above stagnation pressure. The BL speed is
// already a positive unknown; its thermal/subsonic domain is still checked.
export function streamtubeEdgePressure(ue, { flowModel, mach, gamma, pInf }) {
  if (!(ue > 0) || !Number.isFinite(ue)) throw new Error('Invalid BL edge speed.');
  if (flowModel === 'incompressible') return { pressure: -.5 * ue * ue, derivative: -ue };
  const temperature = 1 + .5 * (gamma - 1) * mach * mach * (1 - ue * ue);
  if (!(temperature > 0) || !(mach * mach * ue * ue < temperature))
    throw new Error('BL edge speed is sonic or outside its thermal domain.');
  const rho = temperature ** (1 / (gamma - 1));
  return { pressure: pInf * temperature ** (gamma / (gamma - 1)), derivative: -rho * ue };
}
