// SPDX-License-Identifier: GPL-2.0-or-later
// MSES manual §§1.2.5–1.2.6: speed bias on three successive sections.
// Coordinates increase downstream; spacing contains actual arc increments.
// This local filter does not choose boundary stencils or compute a shock.

function inputs({ speeds, machSquared, spacing, mucon = 1, mcrit = .99, gamma = 1.4 }) {
  if (!Array.isArray(speeds) || speeds.length !== 3 || !speeds.every(q => Number.isFinite(q) && q >= 0)
    || !Array.isArray(machSquared) || machSquared.length !== 2 || !machSquared.every(m => Number.isFinite(m) && m >= 0)
    || !Array.isArray(spacing) || spacing.length !== 2 || !spacing.every(s => Number.isFinite(s) && s > 0)
    || ![mucon, mcrit, gamma].every(Number.isFinite) || mcrit < 0 || mcrit > 1 || gamma <= 1)
    throw new Error('Invalid streamtube speed-upwind inputs.');
  return { speeds: [...speeds], machSquared: [...machSquared], spacing: [...spacing], mucon, mcrit, gamma };
}

function activation(x, e) {
  // The x=0 limit is flat. Avoid both division by zero and 0*Infinity.
  if (x === 0) return { base: 0, dx: 0, de: 0 };
  const a = 1 - 1 / x;
  if (e === 0) return { base: Math.max(0, a), dx: x > 1 ? 1 / x / x : 0, de: 0 };
  const z = a / e, az = Math.abs(z), tail = Math.exp(-az);
  const sigmoid = z >= 0 ? 1 / (1 + tail) : tail / (1 + tail);
  // e*softplus(a/e), without overflow as e approaches zero.
  const base = Math.max(a, 0) + e * Math.log1p(tail);
  // softplus(z)-z*sigmoid(z) is even. This form avoids subtracting
  // nearly equal large numbers in the derivative with respect to e.
  const de = Math.log1p(tail) + (Number.isFinite(az) ? az * tail / (1 + tail) : 0);
  return { base, dx: sigmoid === 0 ? 0 : sigmoid / x / x, de };
}

function evaluate(p) {
  const [q0, q1, q2] = p.speeds, [d0, d1] = p.spacing;
  const secondOrder = p.mucon >= 0, ratio = d1 / d0;
  const x = .5 * p.machSquared[0] + .5 * p.machSquared[1];
  const weight = activation(x, 1 - p.mcrit);
  const coefficient = Math.abs(p.mucon) / p.gamma * weight.base;
  const difference = -(q2 - q1) + (secondOrder ? ratio * (q1 - q0) : 0);
  const correction = coefficient * difference, speed = q2 + correction;
  if (![speed, coefficient, correction, ratio, difference].every(Number.isFinite))
    throw new Error('Nonfinite streamtube speed-upwind result.');
  return { value: { speed, coefficient, correction, secondOrder }, x, weight, ratio, difference };
}

export function evaluateStreamtubeSpeedUpwind(parameters) {
  return evaluate(inputs(parameters)).value;
}

export function linearizeStreamtubeSpeedUpwind(parameters) {
  const p = inputs(parameters), { value, x, weight, ratio, difference } = evaluate(p);
  if (p.mcrit === 1 && x === 1)
    throw new Error('Speed-upwind activation is not differentiable at sonic Mach with MCRIT=1.');
  const apply = ({ speeds = [0, 0, 0], machSquared = [0, 0], spacing = [0, 0],
    mucon = 0, mcrit = 0, gamma = 0 } = {}) => {
    if (![[speeds, 3], [machSquared, 2], [spacing, 2]].every(([a, n]) => Array.isArray(a) && a.length === n && a.every(Number.isFinite))
      || ![mucon, mcrit, gamma].every(Number.isFinite)) throw new Error('Invalid streamtube speed-upwind tangent.');
    if (p.mcrit === 1 && mcrit !== 0)
      throw new Error('MCRIT tangent must vanish on the MCRIT=1 boundary.');
    if (p.mucon === 0 && mucon !== 0)
      throw new Error('MUCON tangent is not differentiable at the order-switch boundary.');
    const dx = .5 * machSquared[0] + .5 * machSquared[1];
    const dc = Math.sign(p.mucon) * mucon / p.gamma - Math.abs(p.mucon) * gamma / (p.gamma * p.gamma);
    const coefficient = dc * weight.base + Math.abs(p.mucon) / p.gamma * (weight.dx * dx - weight.de * mcrit);
    const dr = (spacing[1] - ratio * spacing[0]) / p.spacing[0];
    const dd = -(speeds[2] - speeds[1]) + (value.secondOrder
      ? ratio * (speeds[1] - speeds[0]) + (p.speeds[1] - p.speeds[0]) * dr : 0);
    const correction = coefficient * difference + value.coefficient * dd;
    const speed = speeds[2] + correction;
    if (![speed, coefficient, correction].every(Number.isFinite)) throw new Error('Nonfinite streamtube speed-upwind derivative.');
    return { speed, coefficient, correction };
  };
  return { value, apply };
}
