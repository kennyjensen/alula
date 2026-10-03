// SPDX-License-Identifier: GPL-2.0-or-later
// Steinbrenner & Chawner, Gridgen's Implementation of PDE-Based Structured
// Grid Generation Methods, Eqs. (1), (4): Δxi = |grad xi|² F, with boundary
// F = -(r_xixi . r_xi)/(r_xi . r_xi), then interpolate into the interior.
// We retain harmonic mass eta, so only the streamwise control is used.
// This restriction and linear interpolation in mass fraction are explicit
// reconstruction choices, not recovered modern MSET X controls.
export function createBoundaryStretchControl({ nodes, xi, eta, discretization = 'giles-1985', metric = 'vector' }) {
  const ordered = a => Array.isArray(a) && a.length >= 3 && a[0] === 0 && a.at(-1) === 1
    && a.every((v, i) => Number.isFinite(v) && (!i || v > a[i - 1]));
  if (!ordered(xi) || !ordered(eta) || !['giles-1985', 'quadratic'].includes(discretization) || !['vector', 'polygon-arc'].includes(metric)
    || !Array.isArray(nodes) || nodes.length !== xi.length
    || nodes.some(row => !Array.isArray(row) || row.length !== eta.length || row.some(p => !Number.isFinite(p?.x) || !Number.isFinite(p?.y))))
    throw new Error('Boundary stretch requires a finite grid and ordered normalized coordinates.');
  const boundary = j => nodes.map((row, i) => {
    // No PDE is solved at inlet/outlet nodes; these unused values are zero.
    if (i === 0 || i === nodes.length - 1) return 0;
    const a = xi[i] - xi[i - 1], b = xi[i + 1] - xi[i];
    if (metric === 'polygon-arc') {
      // On a smooth curve r_xixi.r_xi / |r_xi|² = s_xixi/s_xi.
      // Measuring s directly also retains constant speed at a boundary
      // corner, where the centered vector derivative is not a tangent.
      // These interval lengths are polygon arc, not exact spline arc.
      const left = Math.hypot(row[j].x - nodes[i - 1][j].x, row[j].y - nodes[i - 1][j].y);
      const right = Math.hypot(nodes[i + 1][j].x - row[j].x, nodes[i + 1][j].y - row[j].y);
      const first = discretization === 'giles-1985' ? (left + right) / (a + b) : (b * left / a + a * right / b) / (a + b);
      const second = 2 * (right / b - left / a) / (a + b), f = -second / first;
      if (!(left > 0 && right > 0 && first > 0) || !Number.isFinite(f)) throw new Error(`Singular boundary arc stretch at station ${i}, streamline ${j}.`);
      return f;
    }
    const derivatives = ['x', 'y'].map(k => {
      const left = row[j][k] - nodes[i - 1][j][k], right = nodes[i + 1][j][k] - row[j][k];
      return { first: discretization === 'giles-1985' ? (left + right) / (a + b) : (b * left / a + a * right / b) / (a + b),
        second: 2 * (right / b - left / a) / (a + b) };
    });
    const gamma = derivatives.reduce((s, d) => s + d.first * d.first, 0);
    const f = -derivatives.reduce((s, d) => s + d.first * d.second, 0) / gamma;
    if (!(gamma > 0) || !Number.isFinite(f)) throw new Error(`Singular boundary stretch at station ${i}, streamline ${j}.`);
    return f;
  });
  const lower = boundary(0), upper = boundary(eta.length - 1);
  return { lower, upper, values: lower.map((f, i) => eta.map(e => (1 - e) * f + e * upper[i])), metric,
    source: 'Thomas-Middlecoff streamwise boundary formula, interpolated in mass fraction; transverse control zero',
    exactMsetSpacingLaw: false };
}
