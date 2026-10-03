// SPDX-License-Identifier: GPL-2.0-or-later
// One filtered momentum speed per physical streamtube section. Neighboring
// cells must share this output; filtering is never recursively applied to
// already filtered speeds. Physical mass/enthalpy states remain caller data.
import { evaluateStreamtubeSpeedUpwind, linearizeStreamtubeSpeedUpwind } from '../streamtube-speed-upwind.js';

const sequence = value => Array.isArray(value) || (ArrayBuffer.isView(value) && !(value instanceof DataView));
function numbers(value, n, label, nonnegative = false) {
  if (!sequence(value) || value.length !== n
    || !Array.from(value).every(v => Number.isFinite(v) && (!nonnegative || v >= 0)))
    throw new Error(`Invalid streamtube speed-stencil ${label}.`);
  return Array.from(value);
}
function increasing(arc) {
  for (let i = 1; i < arc.length; i++) {
    const ds = arc[i] - arc[i - 1];
    if (!(ds > 0) || !Number.isFinite(ds)) throw new Error('Streamtube section arcs must increase with finite positive increments.');
  }
}
function inputs({ speeds, machSquared, sectionArc, boundary, mucon = 1, mcrit = .99, gamma = 1.4 }) {
  if (!sequence(speeds) || speeds.length < 1) throw new Error('Supply at least one physical streamtube section.');
  const n = speeds.length;
  const q = numbers(speeds, n, 'speeds', true), m2 = numbers(machSquared, n, 'Mach squared', true);
  const arc = numbers(sectionArc, n, 'section arcs'); increasing(arc);
  if (![mucon, mcrit, gamma].every(Number.isFinite) || mcrit < 0 || mcrit > 1 || gamma <= 1)
    throw new Error('Invalid streamtube speed-stencil controls.');
  if (!boundary || !['upstream-history', 'unfiltered-first-two'].includes(boundary.kind))
    throw new Error('Choose an explicit streamtube speed-stencil boundary closure.');
  const history = boundary.kind === 'upstream-history';
  if (history) {
    q.unshift(...numbers(boundary.speeds, 2, 'upstream speeds', true));
    m2.unshift(...numbers(boundary.machSquared, 2, 'upstream Mach squared', true));
    arc.unshift(...numbers(boundary.sectionArc, 2, 'upstream section arcs')); increasing(arc);
  } else if (['speeds', 'machSquared', 'sectionArc'].some(key => Object.hasOwn(boundary, key))) {
    throw new Error('The unfiltered-first-two closure does not accept upstream history.');
  }
  return { n, q, m2, arc, offset: history ? 2 : 0, first: history ? 0 : Math.min(2, n),
    kind: boundary.kind, mucon, mcrit, gamma };
}
function build(parameters, linearize) {
  const p = inputs(parameters), { n, q, m2, arc, offset, first, mucon, mcrit, gamma } = p;
  const value = { speeds: q.slice(offset), coefficients: Array(n).fill(0), corrections: Array(n).fill(0),
    filtered: Array(n).fill(false), secondOrder: mucon >= 0, boundary: p.kind };
  const locals = Array(n).fill(null);
  for (let i = first; i < n; i++) {
    const j = i + offset;
    const localInput = { speeds: q.slice(j - 2, j + 1), machSquared: m2.slice(j - 1, j + 1),
      spacing: [arc[j - 1] - arc[j - 2], arc[j] - arc[j - 1]], mucon, mcrit, gamma };
    const local = linearize ? linearizeStreamtubeSpeedUpwind(localInput) : null;
    const filtered = local ? local.value : evaluateStreamtubeSpeedUpwind(localInput);
    value.speeds[i] = filtered.speed; value.coefficients[i] = filtered.coefficient;
    value.corrections[i] = filtered.correction; value.filtered[i] = true; locals[i] = local;
  }
  return { p, value, locals };
}

// The manual gives the interior stencil, not its inlet history. The caller
// must either supply two physical upstream states at strictly earlier arcs,
// or explicitly select qtilde_0=q_0, qtilde_1=q_1. The latter is an experimental
// boundary closure, not a claim about the inlet implementation in MSES.
export function evaluateStreamtubeSpeedStencil(parameters) {
  return build(parameters, false).value;
}

// O(N) construction and each application; no global basis vectors or dense
// Jacobian are formed. Row i depends only on physical q[i-2..i], M²[i-1..i],
// arc[i-2..i] and controls, with negative indices supplied by explicit history.
export function linearizeStreamtubeSpeedStencil(parameters) {
  const { p, value, locals } = build(parameters, true);
  const apply = ({ speeds, machSquared, sectionArc, boundary, mucon = 0, mcrit = 0, gamma = 0 } = {}) => {
    const tangent = (row, n, label) => row === undefined ? Array(n).fill(0) : numbers(row, n, label);
    const dq = tangent(speeds, p.n, 'speed tangent'), dm = tangent(machSquared, p.n, 'Mach-squared tangent');
    const ds = tangent(sectionArc, p.n, 'section-arc tangent');
    if (![mucon, mcrit, gamma].every(Number.isFinite)) throw new Error('Invalid streamtube speed-stencil control tangent.');
    if (boundary !== undefined && (!boundary || typeof boundary !== 'object' || Array.isArray(boundary)))
      throw new Error('Invalid streamtube speed-stencil boundary tangent.');
    if (p.offset) {
      dq.unshift(...tangent(boundary?.speeds, 2, 'upstream speed tangent'));
      dm.unshift(...tangent(boundary?.machSquared, 2, 'upstream Mach-squared tangent'));
      ds.unshift(...tangent(boundary?.sectionArc, 2, 'upstream section-arc tangent'));
    } else if (boundary && Object.keys(boundary).length) {
      throw new Error('The unfiltered-first-two closure has no upstream-history tangent.');
    }
    const result = { speeds: dq.slice(p.offset), coefficients: Array(p.n).fill(0), corrections: Array(p.n).fill(0) };
    for (let i = p.first; i < p.n; i++) {
      const j = i + p.offset;
      const d = locals[i].apply({ speeds: dq.slice(j - 2, j + 1), machSquared: dm.slice(j - 1, j + 1),
        spacing: [ds[j - 1] - ds[j - 2], ds[j] - ds[j - 1]], mucon, mcrit, gamma });
      result.speeds[i] = d.speed; result.coefficients[i] = d.coefficient; result.corrections[i] = d.correction;
    }
    return result;
  };
  return { value, apply };
}
