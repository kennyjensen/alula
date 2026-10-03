// SPDX-License-Identifier: GPL-2.0-or-later
import { prepareContour, validateAssembly, pointInside } from '../geometry/airfoil.js';
import { createContourCurve } from '../geometry/contour-curve.js';
import { normInf, solveLinear, linearResidual } from '../numerics/linear.js';
import { makePanel, sourceVelocity } from './panel.js';
import { vortexStreamfunctionBasis } from './streamfunction.js';
import { solveFiniteBaseInviscid } from './finite-base-panel.js';
import { prepareAirfoilElement } from '../geometry/airfoil-element.js';

// Analytic velocity from linear vortex basis functions (1-s/L) and s/L.
// Circulation density is positive counterclockwise. Panel orientation is CCW;
// the exterior velocity limit is taken on its right side.
export function vortexBasis(point, panel, self = false) {
  const source = sourceVelocity(point, panel, self);
  const x = ((point.x - panel.a.x) * panel.tx + (point.y - panel.a.y) * panel.ty) / panel.length;
  const y = (-(point.x - panel.a.x) * panel.ty + (point.y - panel.a.y) * panel.tx) / panel.length;
  const sx = source.u * panel.tx + source.v * panel.ty;
  const sy = -source.u * panel.ty + source.v * panel.tx;
  const endX = x * sx + y * sy - 1 / (2 * Math.PI);
  const endY = x * sy - y * sx;
  return [[sx - endX, sy - endY], [endX, endY]].map(([u, v]) => ({
    u: -v * panel.tx - u * panel.ty, v: -v * panel.ty + u * panel.tx,
  }));
}

export function solveInviscid({ elements, alpha = 0, mach = 0, viscous = false,
  referenceChord = 1, momentReference = { x: 0.25, y: 0 }, boundaryCondition = 'normal-velocity' }) {
  if (Array.isArray(elements)) elements = elements.map(prepareAirfoilElement);
  if (elements?.some(e => e.trailingEdge?.kind === 'finite-base')) return solveFiniteBaseInviscid({
    elements, alpha, mach, viscous, referenceChord, momentReference, boundaryCondition }, vortexBasis);
  if (!['normal-velocity', 'streamfunction'].includes(boundaryCondition)) throw new Error('Unknown panel boundary condition.');
  if (mach !== 0 || viscous) throw new Error('Only incompressible inviscid flow (Mach 0) is implemented.');
  if (!Number.isFinite(alpha) || Math.abs(alpha) > 90 || !(referenceChord > 0) || !Number.isFinite(referenceChord)
    || ![momentReference.x, momentReference.y].every(Number.isFinite)) throw new Error('Invalid flow or reference parameters.');
  if (!Array.isArray(elements) || !elements.length || elements.length > 6) throw new Error('Supply between one and six elements.');
  const contours = elements.map(e => prepareContour(e.points));
  validateAssembly(contours);
  const panels = []; const ranges = [];
  let size = 0;
  contours.forEach((points, element) => {
    const first = panels.length;
    const start = size;
    for (let i = 0; i < points.length - 1; i++) {
      panels.push({ ...makePanel(points[i], points[i + 1], element), node: start + i });
    }
    size += points.length;
    ranges.push({ first, last: panels.length - 1, start, end: size - 1 });
  });
  const n = panels.length;
  if (n > 700) throw new Error('The reference solver supports at most 700 total panels.');
  const strengthCount = size, streamfunction = boundaryCondition === 'streamfunction', teProbes = [];
  if (streamfunction) size += elements.length;
  const a = new Float64Array(size * size); const tangential = new Float64Array(n * strengthCount);
  const rhs = new Float64Array(size);
  const radians = alpha * Math.PI / 180;
  const u = Math.cos(radians); const v = Math.sin(radians);
  for (let i = 0; i < n; i++) {
    const target = panels[i];
    rhs[i] = streamfunction ? (-u * target.a.y + v * target.a.x) / referenceChord : -u * target.nx - v * target.ny;
    if (streamfunction) a[i * size + strengthCount + target.element] = -1;
    for (let j = 0; j < n; j++) {
      const source = panels[j];
      const basis = vortexBasis(target, source, i === j);
      const psi = streamfunction ? vortexStreamfunctionBasis(target.a, source) : null;
      for (let k = 0; k < 2; k++) {
        const col = source.node + k;
        a[i * size + col] += streamfunction ? psi[k] / referenceChord : basis[k].u * target.nx + basis[k].v * target.ny;
        tangential[i * strengthCount + col] += basis[k].u * target.tx + basis[k].v * target.ty;
      }
    }
  }
  // Independent top/bottom TE values are constrained by one Kutta condition
  // per element. Contiguous surface nodes elsewhere share their strengths.
  ranges.forEach(({ start, end }, e) => { a[(n + e) * size + start] = 1; a[(n + e) * size + end] = 1; });
  if (streamfunction) ranges.forEach(({ first, last }, e) => {
    // A closed sharp TE repeats a geometric node. Its second psi equation
    // would be identical, leaving a null mode. As in XFOIL's GGCALC, replace
    // it with zero interior velocity along the TE bisector, 0.1 of the
    // shorter adjacent panel ahead of the TE. Kutta remains a separate row.
    const curve = createContourCurve(contours[e]), top = curve.evaluate(0).derivative, bottom = curve.evaluate(curve.length).derivative;
    const tx = -top.x / Math.hypot(top.x, top.y) + bottom.x / Math.hypot(bottom.x, bottom.y);
    const ty = -top.y / Math.hypot(top.x, top.y) + bottom.y / Math.hypot(bottom.x, bottom.y), length = Math.hypot(tx, ty);
    if (!(length > 0)) throw new Error('Degenerate sharp-TE bisector.');
    const tangent = { x: tx / length, y: ty / length }, distance = .1 * Math.min(panels[first].length, panels[last].length);
    const point = { x: contours[e][0].x - distance * tangent.x, y: contours[e][0].y - distance * tangent.y };
    if (!pointInside(point, contours[e])) throw new Error('Sharp-TE velocity control point must lie inside the body.');
    const row = n + elements.length + e; rhs[row] = -u * tangent.x - v * tangent.y;
    for (const source of panels) vortexBasis(point, source).forEach((q, k) => { a[row * size + source.node + k] += q.u * tangent.x + q.v * tangent.y; });
    teProbes.push({ point, tangent });
  });
  const solution = solveLinear(a, rhs), gamma = solution.slice(0, strengthCount);
  const residual = linearResidual(a, solution, rhs);
  const perElement = elements.map((e, i) => ({ name: e.name ?? `Element ${i + 1}`, cx: 0, cy: 0, cm: 0,
    circulation: 0, cp: [], points: contours[i] }));
  for (let i = 0; i < n; i++) {
    const p = panels[i];
    let qt = u * p.tx + v * p.ty;
    for (let col = 0; col < strengthCount; col++) qt += tangential[i * strengthCount + col] * gamma[col];
    const cp = 1 - qt * qt;
    const fx = -cp * p.nx * p.length / referenceChord;
    const fy = -cp * p.ny * p.length / referenceChord;
    const e = perElement[p.element];
    e.cx += fx; e.cy += fy;
    e.cm -= ((p.x - momentReference.x) * fy - (p.y - momentReference.y) * fx) / referenceChord;
    e.circulation += (gamma[p.node] + gamma[p.node + 1]) * p.length / 2;
    e.cp.push({ x: p.x, y: p.y, cp, qt, length: p.length });
  }
  for (const e of perElement) {
    e.cl = -e.cx * v + e.cy * u;
    e.pressureDrag = e.cx * u + e.cy * v;
  }
  const sum = key => perElement.reduce((total, e) => total + e[key], 0);
  const circulationLift = -2 * sum('circulation') / referenceChord;
  const diagnostics = { linearResidual: normInf(residual),
    ...(streamfunction ? { nodalStreamfunctionResidual: normInf(residual.subarray(0, n)), interiorTEResidual: normInf(residual.subarray(n + elements.length)),
      surfaceStreamfunctions: Array.from(solution.slice(strengthCount), psi => referenceChord * psi), teProbes }
      : { normalVelocityResidual: normInf(residual.subarray(0, n)) }),
    kuttaResidual: normInf(residual.subarray(n, n + elements.length)), circulationLift, liftMismatch: Math.abs(sum('cl') - circulationLift),
    pressureDrag: sum('pressureDrag') };
  const warnings = [];
  if (Math.abs(diagnostics.pressureDrag) > 0.005) warnings.push('Pressure-force drag error is large; refine the panels.');
  if (diagnostics.liftMismatch > 0.01) warnings.push('Pressure and circulation lift disagree; refine the panels.');
  if (Math.abs(alpha) > 12) warnings.push('High incidence: potential flow cannot predict separation or stall.');
  return { model: 'incompressible-linear-vortex', status: diagnostics.linearResidual < 1e-9 ? 'solved' : 'failed',
    alpha, mach: 0, boundaryCondition, referenceChord, momentReference, panelCount: n, cl: sum('cl'), cm: sum('cm'), cd: null,
    elements: perElement, diagnostics, warnings, field: { panels, gamma, u, v } };
}

export function velocityAt(point, field) {
  let { u, v } = field;
  for (const p of field.panels) {
    // Same arithmetic and accumulation order as vortexBasis, without its
    // per-panel arrays, callback and two basis objects. Streamline tracing
    // evaluates this field millions of panel/point pairs.
    const source = sourceVelocity(point, p);
    const x = ((point.x - p.a.x) * p.tx + (point.y - p.a.y) * p.ty) / p.length;
    const y = (-(point.x - p.a.x) * p.ty + (point.y - p.a.y) * p.tx) / p.length;
    const sx = source.u * p.tx + source.v * p.ty;
    const sy = -source.u * p.ty + source.v * p.tx;
    const endX = x * sx + y * sy - 1 / (2 * Math.PI), endY = x * sy - y * sx;
    const startX = sx - endX, startY = sy - endY;
    u += field.gamma[p.node] * (-startY * p.tx - startX * p.ty);
    v += field.gamma[p.node] * (-startY * p.ty + startX * p.tx);
    u += field.gamma[p.node + 1] * (-endY * p.tx - endX * p.ty);
    v += field.gamma[p.node + 1] * (-endY * p.ty + endX * p.tx);
  }
  for (const p of field.basePanels ?? []) {
    const source = sourceVelocity(point, p);
    u += p.sourceStrength * source.u - p.vortexStrength * source.v;
    v += p.sourceStrength * source.v + p.vortexStrength * source.u;
  }
  return { u, v };
}
