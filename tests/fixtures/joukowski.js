// SPDX-License-Identifier: GPL-2.0-or-later
// Independent conformal-map oracle: z = zeta + 1/zeta, circle centered at
// (-epsilon,0), radius 1+epsilon, and a cusp at zeta=1. U_inf=1.
export function joukowski(panels = 160, epsilon = 0.08, alpha = 4) {
  const radius = 1 + epsilon;
  const leading = -1 - 2 * epsilon;
  const xLE = leading + 1 / leading;
  const chord = 2 - xLE;
  const points = [];
  for (let i = 0; i <= panels; i++) {
    const t = 2 * Math.PI * i / panels;
    const x = -epsilon + radius * Math.cos(t); const y = radius * Math.sin(t);
    const d = x * x + y * y;
    points.push({ x: (x + x / d - xLE) / chord, y: (y - y / d) / chord });
  }
  points[panels] = { ...points[0] };
  const a = alpha * Math.PI / 180;
  // Circle surface speed divided by |dz/dzeta|. Map scaling changes neither
  // nondimensional velocity nor Cp when the whole physical flow is scaled.
  const cpAtAngle = t => {
    const x = -epsilon + radius * Math.cos(t); const y = radius * Math.sin(t);
    const r2 = x * x + y * y;
    const derivative = Math.hypot(1 - (x * x - y * y) / r2 ** 2, 2 * x * y / r2 ** 2);
    return 1 - ((-2 * Math.sin(t - a) - 2 * Math.sin(a)) / derivative) ** 2;
  };
  return { points, cl: 8 * Math.PI * radius * Math.sin(a) / chord, cpAtAngle };
}
