// SPDX-License-Identifier: GPL-2.0-or-later
// Recover a candidate optimum from a nearly converged chain barrier iterate.
// Active homogeneous adjacent rows define disjoint one-parameter chains.
// Their original least-squares objective and the one total-length equality
// can be solved without forming I + A' diag(z/s) A near zero slack.
// The candidate is accepted ONLY by the original, unscaled KKT tolerances.
export function polishChainQuadratic({ constraints, slack, dual, equalityMultiplier, tolerance }) {
  const equality = constraints[0], rows = constraints.slice(1), n = equality.terms.length;
  const activeTolerance = Math.sqrt(tolerance);
  const active = rows.map((_, k) => slack[k] <= activeTolerance);
  const edges = Array(n - 1).fill(null), anchors = Array.from({ length: n }, () => []);
  for (let k = 0; k < rows.length; k++) if (active[k]) {
    const c = rows[k];
    if (c.terms.length === 1) {
      const [i, a] = c.terms[0];
      if (!a) return null;
      anchors[i].push({ k, a, rhs: c.rhs });
    } else {
      const [[i, a], [j, b]] = c.terms;
      // The spacing problem has homogeneous growth rows with opposite signs.
      // Other valid chain QPs continue through their normal barrier path.
      if (c.rhs !== 0 || !(a * b < 0)) return null;
      const lo = Math.min(i, j), edge = { k, a: i === lo ? a : b, b: i === lo ? b : a };
      if (edges[lo]) return null;
      edges[lo] = edge;
    }
  }
  const e = Array(n).fill(0);
  equality.terms.forEach(([i, a]) => { e[i] = a; });
  const chains = [];
  for (let first = 0; first < n;) {
    let last = first;
    const logs = [0];
    while (last < n - 1 && edges[last]) {
      logs.push(logs.at(-1) + Math.log(Math.abs(edges[last].a)) - Math.log(Math.abs(edges[last].b)));
      last++;
    }
    const peak = Math.max(...logs), v = logs.map(log => Math.exp(log - peak));
    if (v.some(w => !(w > 0))) return null;
    const fixed = [];
    for (let i = first; i <= last; i++) for (const a of anchors[i]) fixed.push({ ...a, i });
    const vv = v.reduce((s, w) => s + w * w, 0), one = v.reduce((s, w) => s + w, 0);
    const ev = v.reduce((s, w, j) => s + e[first + j] * w, 0);
    const value = fixed.length ? fixed[0].rhs / (fixed[0].a * v[fixed[0].i - first]) : undefined;
    if (fixed.length && !Number.isFinite(value)) return null;
    chains.push({ first, last, v, fixed, vv, one, ev, value });
    first = last + 1;
  }
  let target = -equality.rhs, schur = 0;
  for (const c of chains) {
    if (c.fixed.length) target += c.ev * c.value;
    else { target += c.ev * c.one / c.vv; schur += c.ev * c.ev / c.vv; }
  }
  const lambda = schur > 0 ? target / schur : equalityMultiplier;
  if (!Number.isFinite(lambda)) return null;
  const x = Array(n), z = Array(rows.length).fill(0);
  for (const c of chains) {
    const value = c.fixed.length ? c.value : (c.one - lambda * c.ev) / c.vv;
    c.v.forEach((v, j) => { x[c.first + j] = v * value; });
  }
  const residual = x.map((v, i) => v - 1 + lambda * e[i]);
  for (const c of chains) {
    const anchor = c.fixed[0], root = anchor?.i ?? c.last;
    // Redundant endpoint anchors retain their barrier multipliers; the
    // remaining tree equations determine corrections to an independent set.
    for (const a of c.fixed.slice(1)) { z[a.k] = dual[a.k]; residual[a.i] += a.a * z[a.k]; }
    for (let i = c.first; i < root; i++) {
      const edge = edges[i]; z[edge.k] = -residual[i] / edge.a;
      residual[i + 1] += edge.b * z[edge.k];
    }
    for (let i = c.last; i > root; i--) {
      const edge = edges[i - 1]; z[edge.k] = -residual[i] / edge.b;
      residual[i - 1] += edge.a * z[edge.k];
    }
    if (anchor) z[anchor.k] = -residual[root] / anchor.a;
  }
  if (![...x, ...z].every(Number.isFinite) || z.some(v => v < 0)) return null;
  const dot = (c, values) => c.terms.reduce((sum, [i, a]) => sum + a * values[i], 0);
  const stationarity = x.map((v, i) => v - 1 + lambda * e[i]);
  let primalResidual = Math.abs(dot(equality, x) - equality.rhs), complementarity = 0;
  rows.forEach((c, k) => {
    const r = dot(c, x) - c.rhs;
    primalResidual = Math.max(primalResidual, r);
    complementarity = Math.max(complementarity, Math.abs(r * z[k]));
    c.terms.forEach(([i, a]) => { stationarity[i] += a * z[k]; });
  });
  const dualResidual = Math.max(...stationarity.map(Math.abs));
  if (!(Math.max(primalResidual, dualResidual, complementarity) <= tolerance)) return null;
  return { x, equalityMultiplier: lambda, inequalityMultipliers: z, primalResidual, dualResidual,
    complementarity, polished: true, activeRows: active.flatMap((yes, k) => yes ? [k] : []) };
}
