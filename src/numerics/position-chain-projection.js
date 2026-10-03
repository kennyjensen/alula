// SPDX-License-Identifier: GPL-2.0-or-later
// Internal exact equality projection for active adjacent growth rows.
export function projectPositionChains({ active, rows, reference, scales, natural }) {
  const count = reference.length - 1, m = count - 2;
  const edges = Array(m - 1).fill(null), anchors = [];
  for (let k = 0; k < active.length; k++) {
    const row = rows[active[k]], entry = { ...row, activeIndex: k };
    if (row.intervals.length === 1) anchors.push(entry);
    else {
      const i = row.intervals[0][0];
      if (row.intervals[1][0] !== i + 1 || edges[i]) throw new Error('Dependent chain growth row.');
      edges[i] = entry;
    }
  }
  const chains = [];
  for (let a = 0; a < m;) {
    let b = a; while (b < m - 1 && edges[b]) b++;
    const logs = [0];
    for (let i = a; i < b; i++) logs.push(logs.at(-1) + Math.log(-edges[i].intervals[0][1] / edges[i].intervals[1][1]));
    const max = Math.max(...logs), v = logs.map(x => Math.exp(x - max));
    if (!v.every(x => x > 0 && Number.isFinite(x))) throw new Error('Unresolved chain growth representation.');
    const sum = v.reduce((s, x) => s + x, 0), f = [0];
    for (const x of v) f.push(f.at(-1) + x / sum);
    f[f.length - 1] = 1;
    if (f.some((x, i) => i && !(x > f[i - 1]))) throw new Error('Unresolved chain station fractions.');
    chains.push({ a, b, v, sum, f, anchors: anchors.filter(r => r.intervals[0][0] >= a && r.intervals[0][0] <= b) });
    a = b + 1;
  }
  const nc = chains.length, boundary = [1, ...chains.map(c => c.b + 2)];
  const sigma = boundary.map(i => .5 * (natural[i - 1] + natural[i]));
  const fixed = new Map([[0, 0], [nc, 0]]);
  const fix = (k, p) => {
    const e = (p - reference[boundary[k]]) / sigma[k];
    if (fixed.has(k) && Math.abs(fixed.get(k) - e) > 1e-10) throw new Error('Incompatible anchored chain.');
    fixed.set(k, e);
  };
  for (let c = 0; c < nc; c++) for (const row of chains[c].anchors) {
    const [i, a] = row.intervals[0], chain = chains[c], length = row.intervalRhs * chain.sum / (a * chain.v[i - chain.a]);
    if (i === 0) fix(c + 1, reference[1] + length);
    else if (i === m - 1) fix(c, reference[count - 1] - length);
    else throw new Error('Unexpected interior chain anchor.');
  }
  if (chains.every(c => c.anchors.length)) throw new Error('Dependent fully anchored chain system.');
  const diagonal = new Float64Array(nc + 1), off = new Float64Array(nc), rhs = new Float64Array(nc + 1);
  const terms = [];
  for (let c = 0; c < nc; c++) {
    const chain = chains[c];
    for (let j = 0; j < chain.v.length; j++) {
      const i = chain.a + 1 + j; if (i < 2 || i > count - 2) continue;
      const f = chain.f[j], s = scales[i - 2], a = (1 - f) * sigma[c] / s, b = f * sigma[c + 1] / s;
      const target = (reference[i] - ((1 - f) * reference[boundary[c]] + f * reference[boundary[c + 1]])) / s;
      diagonal[c] += a * a; diagonal[c + 1] += b * b; off[c] += a * b;
      rhs[c] += a * target; rhs[c + 1] += b * target;
      terms.push({ i, c, a, b, target });
    }
  }
  const e = new Float64Array(nc + 1);
  for (const [k, v] of fixed) {
    e[k] = v;
    if (k) rhs[k - 1] -= off[k - 1] * v;
    if (k < nc) rhs[k + 1] -= off[k] * v;
  }
  for (let a = 0; a <= nc;) {
    if (fixed.has(a)) { a++; continue; }
    let b = a; while (b < nc && !fixed.has(b + 1)) b++;
    for (let k = a; k <= b; k++) {
      if (k > a) { const q = off[k - 1] / diagonal[k - 1]; diagonal[k] -= q * off[k - 1]; rhs[k] -= q * rhs[k - 1]; }
      if (!(diagonal[k] > 0 && Number.isFinite(diagonal[k]))) throw new Error('Unresolved chain endpoint Hessian.');
    }
    for (let k = b; k >= a; k--) e[k] = (rhs[k] - (k < b ? off[k] * e[k + 1] : 0)) / diagonal[k];
    a = b + 1;
  }
  const target = scales.map(() => 0), gradient = new Float64Array(count + 1);
  for (const t of terms) target[t.i - 2] = t.a * e[t.c] + t.b * e[t.c + 1] - t.target;
  // Respect the actual normalized floating-point rows after the physical
  // endpoint reduction. This solves their original tridiagonal equalities;
  // no constraint or certificate tolerance is enlarged.
  for (const row of anchors) {
    if (row.terms.length !== 1) throw new Error('Unexpected anchor row stencil.');
    const [i, a] = row.terms[0]; target[i] = row.rhs / a;
  }
  for (const c of chains) {
    const n = c.b - c.a; if (!n) continue;
    const d = new Float64Array(n), l = new Float64Array(n), u = new Float64Array(n), b = new Float64Array(n);
    for (let k = 0; k < n; k++) {
      const center = c.a + k, row = edges[c.a + k]; b[k] = row.rhs;
      for (const [i, a] of row.terms) {
        if (i === center) d[k] = a;
        else if (i === center - 1 && k) l[k] = a;
        else if (i === center + 1 && k < n - 1) u[k] = a;
        else b[k] -= a * target[i];
      }
      if (k) { const q = l[k] / d[k - 1]; d[k] -= q * u[k - 1]; b[k] -= q * b[k - 1]; }
      if (!(Math.abs(d[k]) > 0 && Number.isFinite(d[k]))) throw new Error('Unresolved normalized chain equalities.');
    }
    for (let k = n - 1; k >= 0; k--) target[c.a + k] = (b[k] - (k < n - 1 ? u[k] * target[c.a + k + 1] : 0)) / d[k];
  }
  for (let i = 2; i <= count - 2; i++) gradient[i] = target[i - 2] / scales[i - 2];
  const g = new Float64Array(m); let suffix = 0;
  for (let j = count - 2; j >= 1; j--) { suffix += gradient[j + 1]; g[j - 1] = suffix; }
  let lambdaNumerator = 0, lambdaDenominator = 0;
  for (const c of chains) if (!c.anchors.length) {
    for (let k = c.a; k <= c.b; k++) { lambdaNumerator -= c.v[k - c.a] * g[k]; lambdaDenominator += c.v[k - c.a]; }
  }
  const lambda = lambdaNumerator / lambdaDenominator, multipliers = new Float64Array(active.length);
  const residual = g.map(x => x + lambda);
  for (const c of chains) {
    const anchor = c.anchors[0];
    if (c.anchors.length > 1) throw new Error('Dependent doubly anchored chain.');
    if (anchor?.intervals[0][0] === c.a) {
      for (let i = c.b; i > c.a; i--) {
        const row = edges[i - 1], mu = -residual[i] / row.intervals[1][1];
        multipliers[row.activeIndex] = mu * row.norm; residual[i - 1] += mu * row.intervals[0][1];
      }
    } else {
      for (let i = c.a; i < c.b; i++) {
        const row = edges[i], mu = -residual[i] / row.intervals[0][1];
        multipliers[row.activeIndex] = mu * row.norm; residual[i + 1] += mu * row.intervals[1][1];
      }
    }
    if (anchor) {
      const [i, a] = anchor.intervals[0]; multipliers[anchor.activeIndex] = -residual[i] / a * anchor.norm;
    }
  }
  return { target, multipliers: Array.from(multipliers) };
}
