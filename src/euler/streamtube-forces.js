// SPDX-License-Identifier: GPL-2.0-or-later

// Actual Euler interface pressure, integrated on the solid contour.
// Body normalization is rho_inf = q_inf = 1, hence Cp = 2*(p - p_inf).
// This is a pressure contribution only, with no shear or exit-defect model.
import { pressureForces } from '../inviscid/pressure-forces.js';
import { createContourTopology } from '../geometry/contour-topology.js';

export function streamtubeSolidPressureForces({ flow, layout, conditions,
  referenceChord = 1, momentReference = { x: referenceChord / 4, y: 0 } }) {
  const { pInf, alpha, flowModel } = conditions ?? {}, { nx, tubes, bodies } = layout ?? {};
  if (flowModel !== 'compressible' || ![pInf, alpha, referenceChord, momentReference?.x, momentReference?.y].every(Number.isFinite)
    || !(pInf > 0) || !(referenceChord > 0))
    throw new Error('Solid pressure loads require positive compressible reference pressure and chord, and finite reference conditions.');
  if (!Number.isInteger(nx) || nx < 3 || !Array.isArray(bodies) || !bodies.length
    || !Array.isArray(tubes) || tubes.length !== bodies.length + 1
    || !tubes.every(n => Number.isInteger(n) && n > 0)
    || bodies.some(b => !Number.isInteger(b.leadingIndex) || !Number.isInteger(b.trailingIndex)
      || b.leadingIndex < 1 || b.trailingIndex >= nx || b.trailingIndex - b.leadingIndex < 2))
    throw new Error('Solid pressure loads require complete body stations and passage dimensions.');
  const finitePoint = p => p && Number.isFinite(p.x) && Number.isFinite(p.y);
  const validNodes = nodes => Array.isArray(nodes) && nodes.length === tubes.length && nodes.every((group, g) =>
    Array.isArray(group) && group.length === nx + 1 && group.every(row =>
      Array.isArray(row) && row.length === tubes[g] + 1 && row.every(finitePoint)));
  const displaced = Boolean(layout.displacedBoundaries || flow?.displacement || flow?.undisplacedNodes !== undefined);
  if (displaced && !validNodes(flow?.undisplacedNodes))
    throw new Error('Displaced Euler flow requires the complete undisplaced solid grid for pressure loads.');
  if (!validNodes(flow?.nodes) || !Array.isArray(flow.cells) || flow.cells.length !== nx - 1
    || flow.cells.some(row => !Array.isArray(row) || row.length !== tubes.length
      || row.some((group, g) => !Array.isArray(group) || group.length !== tubes[g])))
    throw new Error('Solid pressure loads require complete finite flow nodes and cells.');
  const solid = displaced ? flow.undisplacedNodes : flow.nodes;
  const point = (body, i, side) => side === 'upper' ? solid[body + 1][i][0] : solid[body][i].at(-1);
  const pressure = (body, i, side) => {
    const p = side === 'upper' ? flow.cells[i - 1][body + 1][0]?.interfacePressure?.lower
      : flow.cells[i - 1][body].at(-1)?.interfacePressure?.upper;
    if (!Number.isFinite(p) || !(p > 0))
      throw new Error(`Nonpositive or nonfinite solid interface pressure: body ${body}, station ${i}, ${side}.`);
    return p;
  };
  const coefficients = { cl: 0, cm: 0, cx: 0, cy: 0, pressureIntegralDrag: 0 };
  const perBody = bodies.map((body, b) => {
    const base = body.trailingEdge?.kind === 'finite-base' ? createContourTopology(body.points, body).base : null;
    for (const i of [body.leadingIndex, body.trailingIndex]) {
      if (base && i === body.trailingIndex) {
        if (point(b, i, 'upper').x !== base.points.at(-1).x || point(b, i, 'upper').y !== base.points.at(-1).y
          || point(b, i, 'lower').x !== base.points[0].x || point(b, i, 'lower').y !== base.points[0].y)
          throw new Error('Finite-base pressure loads must retain both original solid TE corners.');
        continue;
      }
      const a = point(b, i, 'upper'), z = point(b, i, 'lower');
      const tolerance = 64 * Number.EPSILON * Math.max(referenceChord, Math.abs(a.x), Math.abs(a.y), Math.abs(z.x), Math.abs(z.y));
      if (Math.hypot(a.x - z.x, a.y - z.y) > tolerance)
        throw new Error('Solid pressure integration requires joined leading and trailing edges; no base-pressure model is supplied.');
    }
    const leadingPressures = { upper: pressure(b, body.leadingIndex, 'upper'), lower: pressure(b, body.leadingIndex, 'lower') };
    // The LE Kutta equation makes these equal at a root. During iteration
    // their mean defines one pressure at the shared stagnation vertex.
    const stagnationPressure = .5 * (leadingPressures.upper + leadingPressures.lower);
    const vertices = [], pressures = [];
    for (let i = body.trailingIndex; i > body.leadingIndex; i--) {
      vertices.push(point(b, i, 'upper')); pressures.push(pressure(b, i, 'upper'));
    }
    vertices.push(point(b, body.leadingIndex, 'upper')); pressures.push(stagnationPressure);
    for (let i = body.leadingIndex + 1; i <= body.trailingIndex; i++) {
      vertices.push(point(b, i, 'lower')); pressures.push(pressure(b, i, 'lower'));
    }
    if (base) {
      // As in XFOIL's closing pressure panel, interpolate between the two
      // TE pressures. They agree when Kutta converges. Keep every actual
      // base segment; the base pressure is a model, not a resolved cavity.
      const lowerPressure = pressures.at(-1), upperPressure = pressures[0];
      let arc = 0;
      for (let k = 1; k < base.points.length; k++) {
        arc += base.panels[k - 1].length;
        const fraction = k === base.points.length - 1 ? 1 : arc / base.length;
        vertices.push(base.points[k]); pressures.push((1 - fraction) * lowerPressure + fraction * upperPressure);
      }
    }
    // The upper-TE -> LE -> lower-TE order must enclose positive area.
    // Use an origin shift so a distant coordinate origin does not corrupt
    // this orientation/degeneracy check through cancellation.
    let twiceArea = 0;
    const origin = vertices[0];
    for (let k = 1; k < vertices.length; k++) {
      const a = vertices[k - 1], z = vertices[k];
      if (!(Math.hypot(z.x - a.x, z.y - a.y) > 0)) throw new Error('Degenerate solid pressure contour segment.');
      twiceArea += (a.x - origin.x) * (z.y - origin.y) - (a.y - origin.y) * (z.x - origin.x);
    }
    if (!(twiceArea > 0) || !Number.isFinite(twiceArea)) throw new Error('Solid pressure contour must enclose positive counterclockwise area.');
    const cp = pressures.map(p => 2 * (p - pInf));
    const forces = pressureForces(vertices.map(p => ({ x: p.x / referenceChord, y: p.y / referenceChord })), cp,
      { alpha, momentOrigin: { x: momentReference.x / referenceChord, y: momentReference.y / referenceChord } });
    const result = { body: b, ...(body.element === undefined ? {} : { element: body.element }),
      cl: forces.cl, cm: forces.cm, cx: forces.cx, cy: forces.cy, pressureIntegralDrag: forces.cd,
      stagnationPressure, leadingPressures, leadingPressureMismatch: leadingPressures.upper - leadingPressures.lower,
      surfaceStations: body.trailingIndex - body.leadingIndex + 1, enclosedSolidArea: .5 * twiceArea };
    if (base) result.basePressureModel = { kind: 'TE-pressure interpolation on retained base',
      panelCount: base.panels.length, upper: pressures[0], lower: pressure(b, body.trailingIndex, 'lower') };
    for (const key of Object.keys(coefficients)) {
      if (!Number.isFinite(result[key])) throw new Error('Nonfinite integrated solid pressure load.');
      coefficients[key] += result[key];
    }
    return result;
  });
  if (!Object.values(coefficients).every(Number.isFinite)) throw new Error('Nonfinite total solid pressure load.');
  return { ...coefficients, perBody, referenceChord, momentReference: { ...momentReference }, alpha,
    geometry: displaced ? 'undisplaced-solid-grid' : 'solid-grid',
    method: 'Linearly interpolated physical Euler interface pressure on the closed solid polygon; LE pressure is the two-bank mean during iteration; positive nose-up moment.',
    dragKind: 'pressure-contribution', includesSkinFriction: false, includesExitDefects: false, physicalAcceptance: false };
}

// Convert the solver's actual dual-volume wall tractions to the user's
// reference chord and nose-up moment origin. This does not certify a root.
export function streamtubeForceCoefficients(flow, conditions, {
  referenceChord = 1, momentReference = { x: referenceChord / 4, y: 0 },
} = {}) {
  const { lengthScale: length, center, alpha } = conditions;
  if (!(referenceChord > 0) || ![referenceChord, momentReference?.x, momentReference?.y,
    length, center?.x, center?.y, alpha].every(Number.isFinite) || !(length > 0))
    throw new Error('Invalid force reference chord, moment origin or solver normalization.');
  if (flow.displacement && !(flow.inviscidBaseWake && flow.displacement.surfaces.every(s =>
    [...s.upper, ...s.lower].every(d => d === 0)))) throw new Error('Displaced-wall tractions require a solid-surface load calculation.');
  if (!flow.diagnosticForces?.length || !flow.diagnosticForces.every(f => [f.cx, f.cy, f.cm].every(Number.isFinite)))
    throw new Error('No finite wall pressure loads are available.');
  const angle = alpha * Math.PI / 180;
  const perBody = flow.diagnosticForces.map(f => {
    const cx = f.cx * length / referenceChord, cy = f.cy * length / referenceChord;
    const cm = (f.cm * length ** 2 + (momentReference.x - center.x) * f.cy * length
      - (momentReference.y - center.y) * f.cx * length) / referenceChord ** 2;
    return { cx, cy, cm, cl: cy * Math.cos(angle) - cx * Math.sin(angle),
      cd: cx * Math.cos(angle) + cy * Math.sin(angle) };
  });
  const total = Object.fromEntries(['cx', 'cy', 'cl', 'cd', 'cm'].map(k => [k, perBody.reduce((s, f) => s + f[k], 0)]));
  return { ...total, perBody, referenceChord, momentReference: { ...momentReference },
    dragKind: 'pressure', method: 'Dual-volume wall pressure integration; no skin friction',
    physicalValidation: false };
}
