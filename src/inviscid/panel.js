// SPDX-License-Identifier: GPL-2.0-or-later
import { prepareContour, validateAssembly } from '../geometry/airfoil.js';
import { solveLinear, linearResidual, normInf } from '../numerics/linear.js';

export function makePanel(a, b, element = 0) {
  const length = Math.hypot(b.x - a.x, b.y - a.y);
  if (!(length > 0)) throw new Error('Zero-length panel.');
  const tx = (b.x - a.x) / length; const ty = (b.y - a.y) / length;
  return { a, b, x: (a.x + b.x) / 2, y: (a.y + b.y) / 2, length, tx, ty, nx: ty, ny: -tx, element };
}

// Integral of a unit source sheet; +90° rotation gives a unit CCW vortex
// sheet. The analytic exterior self limit must be supplied at collocation.
export function sourceVelocity(point, panel, self = false) {
  if (self) return { u: panel.nx / 2, v: panel.ny / 2 };
  const dx = point.x - panel.a.x; const dy = point.y - panel.a.y;
  const x = dx * panel.tx + dy * panel.ty;
  const y = -dx * panel.ty + dy * panel.tx;
  const r1 = x * x + y * y;
  const r2 = (x - panel.length) ** 2 + y * y;
  if (r1 === 0 || r2 === 0) throw new Error('Velocity is singular at a panel endpoint.');
  const tangent = Math.log(r1 / r2) / (4 * Math.PI);
  const left = Math.atan2(y * panel.length, x * (x - panel.length) + y * y) / (2 * Math.PI);
  return { u: tangent * panel.tx - left * panel.ty, v: tangent * panel.ty + left * panel.tx };
}

// Hess–Smith baseline: N constant source strengths plus one circulation
// density per lifting element. N no-penetration rows + one Kutta row/element.
// This is a potential-flow initializer/reference, not the MSES Euler solver.
export function solveInviscid({ elements, alpha = 0, mach = 0, viscous = false,
  referenceChord = 1, momentReference = { x: 0.25, y: 0 } }) {
  if (mach !== 0 || viscous) throw new Error('Only incompressible inviscid flow (Mach 0) is implemented.');
  if (!Number.isFinite(alpha) || Math.abs(alpha) > 90 || !(referenceChord > 0) || !Number.isFinite(referenceChord)
    || ![momentReference.x, momentReference.y].every(Number.isFinite)) throw new Error('Invalid flow or reference parameters.');
  if (!Array.isArray(elements) || !elements.length || elements.length > 6) throw new Error('Supply between one and six elements.');
  const contours = elements.map(e => prepareContour(e.points, { lifting: e.lifting !== false }));
  validateAssembly(contours);
  const panels = []; const ranges = []; const lifting = [];
  contours.forEach((points, e) => {
    const first = panels.length;
    for (let i = 0; i < points.length - 1; i++) panels.push(makePanel(points[i], points[i + 1], e));
    ranges.push({ first, last: panels.length - 1 });
    if (elements[e].lifting !== false) lifting.push(e);
  });
  const n = panels.length; const size = n + lifting.length;
  if (n > 700) throw new Error('The reference solver supports at most 700 total panels.');
  const a = new Float64Array(size * size); const rhs = new Float64Array(size);
  const tangential = new Float64Array(n * size);
  const radians = alpha * Math.PI / 180;
  const u = Math.cos(radians); const v = Math.sin(radians);
  const vortexColumn = new Map(lifting.map((e, k) => [e, n + k]));
  for (let i = 0; i < n; i++) {
    const target = panels[i];
    rhs[i] = -(u * target.nx + v * target.ny);
    for (let j = 0; j < n; j++) {
      const source = sourceVelocity(target, panels[j], i === j);
      a[i * size + j] = source.u * target.nx + source.v * target.ny;
      tangential[i * size + j] = source.u * target.tx + source.v * target.ty;
      const col = vortexColumn.get(panels[j].element);
      if (col !== undefined) {
        a[i * size + col] += -source.v * target.nx + source.u * target.ny;
        tangential[i * size + col] += -source.v * target.tx + source.u * target.ty;
      }
    }
  }
  lifting.forEach((element, k) => {
    const { first, last } = ranges[element];
    for (let col = 0; col < size; col++) a[(n + k) * size + col] = tangential[first * size + col] + tangential[last * size + col];
    rhs[n + k] = -u * (panels[first].tx + panels[last].tx) - v * (panels[first].ty + panels[last].ty);
  });
  const strengths = solveLinear(a, rhs);
  const residual = linearResidual(a, strengths, rhs);
  const perElement = elements.map((e, i) => ({ name: e.name ?? `Element ${i + 1}`, cx: 0, cy: 0, cm: 0,
    sourceFlux: 0, circulation: 0, cp: [], points: contours[i] }));
  for (let i = 0; i < n; i++) {
    const p = panels[i];
    let qt = u * p.tx + v * p.ty;
    for (let col = 0; col < size; col++) qt += tangential[i * size + col] * strengths[col];
    const cp = 1 - qt * qt;
    const fx = -cp * p.nx * p.length / referenceChord;
    const fy = -cp * p.ny * p.length / referenceChord;
    const result = perElement[p.element];
    result.cx += fx; result.cy += fy;
    // x points downstream; positive nose-up pitching moment is clockwise.
    result.cm -= ((p.x - momentReference.x) * fy - (p.y - momentReference.y) * fx) / referenceChord;
    result.sourceFlux += strengths[i] * p.length;
    result.circulation += (strengths[vortexColumn.get(p.element)] ?? 0) * p.length;
    result.cp.push({ x: p.x, y: p.y, cp, qt, length: p.length });
  }
  for (const e of perElement) {
    e.cl = -e.cx * v + e.cy * u;
    e.pressureDrag = e.cx * u + e.cy * v;
  }
  const total = key => perElement.reduce((sum, e) => sum + e[key], 0);
  const circulationLift = -2 * total('circulation') / referenceChord;
  const diagnostics = {
    linearResidual: normInf(residual),
    normalVelocityResidual: normInf(residual.subarray(0, n)),
    kuttaResidual: normInf(residual.subarray(n)),
    sourceFlux: total('sourceFlux') / referenceChord,
    maxElementSourceFlux: Math.max(...perElement.map(e => Math.abs(e.sourceFlux))) / referenceChord,
    circulationLift,
    liftMismatch: Math.abs(total('cl') - circulationLift),
    pressureDrag: total('pressureDrag'),
  };
  const warnings = [];
  if (Math.abs(diagnostics.pressureDrag) > 0.005) warnings.push('Pressure-force drag error is large; refine the panels.');
  if (diagnostics.liftMismatch > 0.01) warnings.push('Pressure and circulation lift disagree; refine the panels.');
  if (diagnostics.maxElementSourceFlux > 0.005) warnings.push('Source-flux imbalance is large; refine the panels.');
  if (Math.abs(alpha) > 12) warnings.push('High incidence: this potential-flow solution cannot predict separation or stall.');
  return { model: 'incompressible-hess-smith', status: diagnostics.linearResidual < 1e-9 ? 'solved' : 'failed',
    alpha, mach: 0, referenceChord, momentReference, panelCount: n,
    cl: total('cl'), cm: total('cm'), cd: null, elements: perElement, diagnostics, warnings,
    // Needed for off-body field queries; all serializable and browser compatible.
    field: { panels, sources: strengths.slice(0, n), vortices: elements.map((_, e) => strengths[vortexColumn.get(e)] ?? 0), u, v },
  };
}

export function velocityAt(point, field) {
  let { u, v } = field;
  for (let i = 0; i < field.panels.length; i++) {
    const panel = field.panels[i];
    const source = sourceVelocity(point, panel);
    const gamma = field.vortices[panel.element];
    u += field.sources[i] * source.u - gamma * source.v;
    v += field.sources[i] * source.v + gamma * source.u;
  }
  return { u, v };
}
