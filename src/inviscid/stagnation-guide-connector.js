// SPDX-License-Identifier: GPL-2.0-or-later
// An explicit geometric INITIALIZER connector, not a panel streamline.
// Giles OUTLIN uses geometric incoming cuts. This local C1 construction is
// our extension where the panel and material stagnation references differ.
import { velocityAt } from './linear-vortex.js';
import { potentialDifference, streamfunctionAt } from './streamfunction.js';
import { intervalPoint as ip, intervalAdd as ia, intervalSub as is,
  intervalMul as im } from '../geometry/bernstein-cell-certificate.js';

const finite = p => p && Number.isFinite(p.x) && Number.isFinite(p.y);
const lerp = (a, b, t) => ({ x: a.x + t * (b.x - a.x), y: a.y + t * (b.y - a.y) });
const evaluate = (controls, t) => {
  if (t === 0) return { ...controls[0] };
  if (t === 1) return { ...controls[3] };
  const a = controls.slice(1).map((p, k) => lerp(controls[k], p, t));
  return lerp(lerp(a[0], a[1], t), lerp(a[1], a[2], t), t);
};
const mixInterval = (a, b) => ({ x: im(ia(a.x, b.x), ip(.5)), y: im(ia(a.y, b.y), ip(.5)) });
const split = c => {
  const a = [mixInterval(c[0], c[1]), mixInterval(c[1], c[2]), mixInterval(c[2], c[3])];
  const b = [mixInterval(a[0], a[1]), mixInterval(a[1], a[2])], m = mixInterval(b[0], b[1]);
  return [[c[0], a[0], b[0], m], [m, b[1], a[2], c[3]]];
};
// Every point of a Bezier segment lies in its control hull. Outward-rounded
// half-plane/projection separation proves that hull misses each panel.
function separated(c, panel) {
  const ax = ip(panel.a.x), ay = ip(panel.a.y), dx = is(ip(panel.b.x), ax), dy = is(ip(panel.b.y), ay);
  const cross = c.map(p => is(im(dx, is(p.y, ay)), im(dy, is(p.x, ax))));
  if (cross.every(v => v[0] > 0) || cross.every(v => v[1] < 0)) return true;
  const along = c.map(p => ia(im(dx, is(p.x, ax)), im(dy, is(p.y, ay))));
  const lengthSquared = ia(im(dx, dx), im(dy, dy));
  return along.every(v => v[1] < 0) || along.every(v => v[0] > lengthSquared[1]);
}

export function createStagnationGuideConnector({ join, anchor, anchorDirection, stagnationPotential,
  field, tolerance = 2e-9, streamfunctionLevel, admissible = () => true,
  admissibleSegment = () => true } = {}) {
  if (![join, anchor, anchorDirection].every(finite) || !Number.isFinite(join.potential)
    || !Number.isFinite(stagnationPotential) || !Number.isFinite(tolerance) || tolerance <= 0
    || !(anchor.x > join.x) || !(anchorDirection.x > 0) || !Array.isArray(field?.panels)
    || typeof admissible !== 'function' || typeof admissibleSegment !== 'function'
    || !admissible(join) || !admissible(anchor)) throw new Error('Invalid geometric stagnation-connector endpoints or field.');
  const q0 = velocityAt(join, field);
  if (!(q0.u > 0) || !Number.isFinite(q0.v)) throw new Error('The incoming connector needs a forward panel-flow tangent.');
  const dx = anchor.x - join.x, m0 = q0.v / q0.u, m1 = anchorDirection.y / anchorDirection.x;
  const controls = [{ x: join.x, y: join.y }, { x: join.x + dx / 3, y: join.y + dx * m0 / 3 },
    { x: anchor.x - dx / 3, y: anchor.y - dx * m1 / 3 }, { ...anchor }];
  if (!controls.every(finite)) throw new Error('Nonfinite geometric stagnation connector.');
  const sheets = [...field.panels, ...(field.basePanels ?? [])], breaks = [0];
  let certificateLeaves = 0, maximumDepth = 0;
  const certify = (c, a, b, candidates, depth) => {
    const pending = candidates.filter(panel => !separated(c, panel));
    if (!pending.length) { breaks.push(b); certificateLeaves++; maximumDepth = Math.max(maximumDepth, depth); return; }
    if (depth >= 20 || certificateLeaves >= 4096)
      throw new Error('Geometric stagnation connector cannot be certified clear of every physical sheet.');
    const [left, right] = split(c), mid = .5 * (a + b);
    certify(left, a, mid, pending, depth + 1); certify(right, mid, b, pending, depth + 1);
  };
  certify(controls.map(p => ({ x: ip(p.x), y: ip(p.y) })), 0, 1, sheets, 0);
  const vertices = breaks.map(t => evaluate(controls, t)), potentials = Array(breaks.length);
  potentials[potentials.length - 1] = stagnationPotential;
  for (let k = vertices.length - 2; k >= 0; k--) {
    if (!admissible(vertices[k]) || !admissibleSegment(vertices[k], vertices[k + 1]))
      throw new Error('Geometric stagnation connector failed its physical segment guard.');
    potentials[k] = potentials[k + 1] - potentialDifference(vertices[k], vertices[k + 1], field);
  }
  const joinPotentialError = potentials[0] - join.potential;
  if (Math.abs(joinPotentialError) > 8 * tolerance)
    throw new Error('Geometric stagnation connector does not match the existing common potential gauge.');
  const atParameter = t => {
    if (!Number.isFinite(t) || t < 0 || t > 1) throw new Error('Geometric stagnation query is outside its source-panel interval.');
    const p = evaluate(controls, t);
    let k = breaks.findIndex(v => v >= t); if (k < 0) k = breaks.length - 1;
    const potential = t === breaks[k] ? potentials[k] : potentials[k] - potentialDifference(p, vertices[k], field);
    const dy = 3 * ((1 - t) ** 2 * (controls[1].y - controls[0].y)
      + 2 * t * (1 - t) * (controls[2].y - controls[1].y) + t * t * (controls[3].y - controls[2].y));
    const q = velocityAt(p, field), derivative = q.u + q.v * dy / dx;
    if (!admissible(p) || !Number.isFinite(potential) || !Number.isFinite(derivative) || !(derivative > 0))
      throw new Error('Geometric stagnation connector is outside the fluid or has nonincreasing panel potential.');
    return { ...p, potential, derivative };
  };
  const samples = [...new Set([...breaks, ...Array.from({ length: 65 }, (_, k) => k / 64)])].sort((a, b) => a - b);
  const checked = samples.map(atParameter);
  if (checked.some((p, k) => k && !(p.potential > checked[k - 1].potential)))
    throw new Error('Geometric stagnation connector potential is not resolved and increasing.');
  const atX = x => {
    if (x === anchor.x) return atParameter(1);
    if (x === join.x) return atParameter(0);
    return atParameter((x - join.x) / dx);
  };
  const atPotential = potential => {
    if (!Number.isFinite(potential) || potential < potentials[0] || potential > stagnationPotential)
      throw new Error('Potential query is outside the geometric stagnation connector.');
    if (potential === potentials[0]) return atParameter(0);
    if (potential === stagnationPotential) return atParameter(1);
    let lo = 0, hi = 1, t = (potential - potentials[0]) / (stagnationPotential - potentials[0]);
    for (let iteration = 0; iteration < 56; iteration++) {
      const p = atParameter(t), error = p.potential - potential;
      if (Math.abs(error) <= tolerance && (Math.abs(error) / p.derivative <= tolerance || (hi - lo) * dx <= tolerance)) return p;
      if (error < 0) lo = t; else hi = t;
      const next = t - error / (dx * p.derivative);
      t = next > lo && next < hi ? next : .5 * (lo + hi);
    }
    throw new Error('Geometric stagnation potential inverse did not converge.');
  };
  const points = checked.map(({ derivative, ...p }) => p);
  return { points, atX, atPotential, join: points[0], anchor: points.at(-1), anchorDirection: { ...anchorDirection },
    diagnostics: { method: 'Source-panel-sized C1 geometric incoming connector',
      interpretation: 'Initialization geometry, not an exact constant-streamfunction panel path; governing Euler/BL equations unchanged.',
      sourceRelationship: 'Geometric incoming cuts follow the role of Giles OUTLIN; this local C1 join is a declared extension, not recovered MSET code.',
      controls, joinPotentialError, potentialTolerance: tolerance,
      geometryCertificate: 'Outward-rounded Bezier control-hull separation from every retained physical sheet, with exterior endpoint and chord checks.',
      certificateLeaves, maximumDepth, monotonicity: 'Strict physical x; sampled positive analytic dphi/dx and ordered potential, with every query checked.',
      samples: checked.length, minimumPotentialDerivative: Math.min(...checked.map(p => p.derivative)),
      ...(Number.isFinite(streamfunctionLevel) ? { maximumStreamfunctionDefect: Math.max(...points.map(p => Math.abs(streamfunctionAt(p, field) - streamfunctionLevel))) } : {}) } };
}
