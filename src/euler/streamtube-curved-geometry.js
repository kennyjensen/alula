// SPDX-License-Identifier: GPL-2.0-or-later
// Cubic polynomial boundary approximations for the already prepared grid.
// No tracing, redistribution, clipping, repairs or BIE-condition changes.
import { createPolynomialGridGeometry } from '../geometry/polynomial-grid-geometry.js';
import { createTransverseHarmonicGrid, smoothTransverseHarmonicGrid } from '../geometry/transverse-harmonic-grid.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';
import { streamtubeMeshConnectivity } from '../geometry/streamtube-mesh-connectivity.js';

const copy = p => ({ x: p.x, y: p.y });
const finite = p => Number.isFinite(p?.x) && Number.isFinite(p?.y);
const linear = (a, b, t) => ({ x: a.x + t * (b.x - a.x), y: a.y + t * (b.y - a.y) });
const chord = (a, b) => [copy(a), linear(a, b, 1 / 3), linear(a, b, 2 / 3), copy(b)];
const hermite = (a, b, da, db) => [copy(a), { x: a.x + da.x / 3, y: a.y + da.y / 3 },
  { x: b.x - db.x / 3, y: b.y - db.y / 3 }, copy(b)];
const times = (p, s) => ({ x: p.x * s, y: p.y * s });

function polynomialControls(nodes, horizontal, vertical) {
  const nx = nodes.length - 1, nt = nodes[0].length - 1;
  return Array.from({ length: nx }, (_, i) => Array.from({ length: nt }, (_, j) => {
    const p = [nodes[i][j], nodes[i + 1][j], nodes[i + 1][j + 1], nodes[i][j + 1]];
    const bottom = horizontal[i][j], top = horizontal[i][j + 1], left = vertical[i][j], right = vertical[i + 1][j];
    return Array.from({ length: 4 }, (_, a) => Array.from({ length: 4 }, (_, b) => {
      // Copy the canonical complete faces, including exact original corners.
      if (!a) return copy(left[b]); if (a === 3) return copy(right[b]);
      if (!b) return copy(bottom[a]); if (b === 3) return copy(top[a]);
      const s = a / 3, t = b / 3, q = linear(linear(p[0], p[1], s), linear(p[3], p[2], s), t);
      for (const [edge, base, weight] of [[bottom[a], linear(p[0], p[1], s), 1 - t],
        [top[a], linear(p[3], p[2], s), t], [left[b], linear(p[0], p[3], t), 1 - s],
        [right[b], linear(p[1], p[2], t), s]]) {
        q.x += weight * (edge.x - base.x); q.y += weight * (edge.y - base.y);
      }
      return q;
    }));
  }));
}

export function createCurvedStreamtubeRegions({ input, system, initial = system?.initial, nodes, guideField, diagnostics }) {
  if (!system?.layout || !Array.isArray(input?.bodies) || !Array.isArray(nodes)
    || nodes.length !== input.bodies.length + 1 || typeof guideField?.velocityAt !== 'function')
    throw new Error('Curved streamtube geometry needs a prepared body system, nodes and its original panel velocity field.');
  const { nx, tubes } = system.layout, decoded = system.decode(initial), velocities = new Map();
  if (!nodes.every((group, g) => Array.isArray(group) && group.length === nx + 1
    && group.every(row => Array.isArray(row) && row.length === tubes[g] + 1 && row.every(finite))))
    throw new Error('Invalid prepared curved-streamtube node dimensions.');
  const velocity = p => {
    const key = `${p.x},${p.y}`;
    if (!velocities.has(key)) {
      const value = guideField.velocityAt(p);
      if (!Number.isFinite(value?.u) || !Number.isFinite(value?.v)) throw new Error('Nonfinite panel velocity at a curved boundary node.');
      velocities.set(key, { x: value.u, y: value.v });
    }
    return velocities.get(key);
  };
  const xDerivative = (direction, dx, location) => {
    const speed = Math.hypot(direction.x, direction.y);
    if (!(speed > 0) || !(direction.x > 1e-10 * speed) || !(dx > 0))
      throw new Error(`Unresolved streamwise-x boundary parameterization at ${location}.`);
    return { x: dx, y: dx * direction.y / direction.x };
  };
  const normalAtStagnation = (body, a, b) => {
    const d = system.curves[body].evaluate(decoded.stagnation[body]).derivative;
    const normal = { x: d.y, y: -d.x };
    if (normal.x * (b.x - a.x) + normal.y * (b.y - a.y) < 0) { normal.x = -normal.x; normal.y = -normal.y; }
    return normal;
  };
  const wakeBisector = body => {
    const curve = system.curves[body], upper = curve.evaluate(0).derivative, lower = curve.evaluate(curve.length).derivative;
    const u = Math.hypot(upper.x, upper.y), l = Math.hypot(lower.x, lower.y);
    if (!(u > 0 && l > 0)) throw new Error('Unresolved trailing-edge wall tangents.');
    return { x: -upper.x / u + lower.x / l, y: -upper.y / u + lower.y / l };
  };
  return nodes.map((group, g) => {
    const nt = tubes[g], massFlows = decoded.allocation.groups[g].map(tube => tube.massFlow);
    const metadata = { method: 'Q3 Hermite boundaries with fixed Coons bows and Q1 interior edges.',
      wallGeometry: 'Endpoint derivatives from the prescribed C2 contour. An interval spanning multiple spline knots is a cubic approximation.',
      nonwallGeometry: 'Same prepared panel field: streamwise-x cut/farfield derivatives; streamfunction-parametrized straight end-edge chords.',
      stagnationCondition: 'The incoming dividing cut meets a regular slip-wall stagnation point along the prescribed wall normal.',
      wakeOriginCondition: 'At the sharp trailing edge, the outgoing Kutta direction bisects the two outgoing wall tangents; panel endpoint velocity is singular.',
      panelBoundaryCondition: diagnostics?.panelBoundaryCondition ?? 'unchanged supplied field',
      wallEdges: 0, wallEdgesSpanningSplineKnots: 0, fieldEdges: 0, verticalEndEdges: 0,
      incomingStagnationNormals: 0, outgoingKuttaBisectors: 0, maximumWallEndpointMismatch: 0,
      physicsValidated: false, globalOverlapCertified: false };
    const streamwiseEdge = (i, j) => {
      const a = group[i][j], b = group[i + 1][j];
      const body = j === 0 && g > 0 ? g - 1 : j === nt && g < input.bodies.length ? g : null;
      const side = j === 0 ? 'upper' : 'lower', range = body === null ? null : input.bodies[body];
      if (range && i >= range.leadingIndex && i < range.trailingIndex) {
        const curve = system.curves[body], stag = decoded.stagnation[body], fractions = system.fractions[body][side];
        const k = i - range.leadingIndex, f0 = fractions[k], f1 = fractions[k + 1];
        const p = curve.branch(side, f0, stag), q = curve.branch(side, f1, stag);
        const factor = (f1 - f0) * (side === 'upper' ? -stag : curve.length - stag);
        metadata.wallEdges++;
        if (curve.knots.some(s => s > Math.min(p.parameter, q.parameter) && s < Math.max(p.parameter, q.parameter)))
          metadata.wallEdgesSpanningSplineKnots++;
        metadata.maximumWallEndpointMismatch = Math.max(metadata.maximumWallEndpointMismatch,
          Math.hypot(a.x - p.point.x, a.y - p.point.y), Math.hypot(b.x - q.point.x, b.y - q.point.y));
        return hermite(a, b, times(p.derivative, factor), times(q.derivative, factor));
      }
      metadata.fieldEdges++;
      const dx = b.x - a.x;
      const firstDirection = range && i === range.trailingIndex ? wakeBisector(body) : velocity(a);
      if (range && i === range.trailingIndex) metadata.outgoingKuttaBisectors++;
      const da = xDerivative(firstDirection, dx, `region ${g}, edge ${i}, boundary ${j}`);
      let lastDirection;
      if (range && i + 1 === range.leadingIndex) {
        lastDirection = normalAtStagnation(body, a, b); metadata.incomingStagnationNormals++;
      } else lastDirection = velocity(b);
      const db = xDerivative(lastDirection, dx, `region ${g}, edge ${i + 1}, boundary ${j}`);
      return hermite(a, b, da, db);
    };
    const endEdge = (i, j) => {
      const a = group[i][j], b = group[i][j + 1], mass = massFlows[j], dx = b.x - a.x, dy = b.y - a.y;
      const derivative = p => {
        const v = velocity(p), flux = v.x * dy - v.y * dx;
        if (!(flux > 1e-10 * Math.hypot(v.x, v.y) * Math.hypot(dx, dy)))
          throw new Error('Unresolved streamfunction parameterization on an end-edge chord.');
        return { x: mass * dx / flux, y: mass * dy / flux };
      };
      metadata.verticalEndEdges++;
      return hermite(a, b, derivative(a), derivative(b));
    };
    const horizontal = Array.from({ length: nx }, (_, i) => Array.from({ length: nt + 1 }, (_, j) =>
      !j || j === nt ? streamwiseEdge(i, j) : chord(group[i][j], group[i + 1][j])));
    const vertical = Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt }, (_, j) =>
      !i || i === nx ? endEdge(i, j) : chord(group[i][j], group[i][j + 1])));
    const tangent = (a, b) => { const x = b.x - a.x, y = b.y - a.y, d = Math.hypot(x, y);
      if (!(d > 0)) throw new Error('Degenerate prepared streamwise edge.'); return { x: x / d, y: y / d }; };
    const directions = group.map((row, i) => row.map((p, j) => {
      const a = i ? tangent(group[i - 1][j], p) : { x: 0, y: 0 };
      const b = i < nx ? tangent(p, group[i + 1][j]) : { x: 0, y: 0 };
      const x = a.x + b.x, y = a.y + b.y, length = Math.hypot(x, y);
      if (!(length > 1e-12)) throw new Error(`Unresolved transverse tangent bisector at region ${g}, station ${i}, streamline ${j}.`);
      return { x: -y / length, y: x / length };
    }));
    return { geometry: createPolynomialGridGeometry({ nodes: group, controlPoints: polynomialControls(group, horizontal, vertical) }),
      directions, massFlows, boundaryApproximation: metadata };
  });
}

// Inspection-only adapter for the certified curved harmonic research kernel.
// The candidate is published independently of physical accuracy acceptance.
function edgeCurves(layout, nodes, regions) {
  const { indices } = streamtubeMeshConnectivity(layout, nodes), curves = [], seen = new Set();
  for (let g = 0; g < nodes.length; g++) {
    const map = regions[g].geometry.onGrid(nodes[g]), nx = nodes[g].length - 1, nt = nodes[g][0].length - 1;
    const add = (a, b, first, last, derivative) => {
      const key = a < b ? `${a}:${b}` : `${b}:${a}`;
      if (seen.has(key)) return; seen.add(key);
      const p = first.point, q = last.point, d = first[derivative], e = last[derivative];
      curves.push({ a, b, c1: { x: p.x + d.x / 3, y: p.y + d.y / 3 }, c2: { x: q.x - e.x / 3, y: q.y - e.y / 3 } });
    };
    for (let i = 0; i < nx; i++) for (let j = 0; j <= nt; j++) {
      const cell = Math.min(j, nt - 1), t = j === nt ? 1 : 0;
      add(indices[g][i][j], indices[g][i + 1][j], map.at(i, cell, 0, t), map.at(i, cell, 1, t), 'ds');
    }
    for (let i = 0; i <= nx; i++) for (let j = 0; j < nt; j++) {
      const cell = Math.min(i, nx - 1), s = i === nx ? 1 : 0;
      add(indices[g][i][j], indices[g][i][j + 1], map.at(cell, j, s, 0), map.at(cell, j, s, 1), 'dt');
    }
  }
  return curves;
}

export function prepareCurvedStreamtubePreview(prepared, { onMesh, maxIterations = 20 } = {}) {
  const { system, guideField } = prepared;
  const nodes = structuredClone(prepared.nodes);
  const report = { attempted: true, method: 'curved harmonic', experimental: true, converged: false,
    physicalAcceptance: false, previewOnly: true, regions: [],
    scope: 'Curved harmonic geometry preview; polynomial boundary approximation and physical refinement remain unvalidated. Euler/BL is not run.' };
  const snapshot = (stage, curves) => {
    const mesh = streamtubeMeshSnapshot({ ...prepared, nodes, diagnostics: { ...prepared.diagnostics, gridSmoothing: report } });
    if (curves) mesh.edgeCurves = curves;
    onMesh?.(mesh, stage); return mesh;
  };
  let regions;
  try {
    regions = createCurvedStreamtubeRegions(prepared);
    for (let g = 0; g < regions.length; g++) {
      const region = regions[g], entry = { region: g, boundaryApproximation: region.boundaryApproximation };
      report.regions.push(entry);
      try {
        const chart = createTransverseHarmonicGrid({ nodes: nodes[g], massFlows: region.massFlows,
          directions: region.directions, curvedGeometry: region.geometry });
        const result = smoothTransverseHarmonicGrid(chart, { maxIterations,
          onIteration: (history, current) => {
            nodes[g] = current; report.current = { region: g, ...history };
            snapshot('smoothing');
          } });
        nodes[g] = result.nodes;
        const { nodes: omitted, ...details } = result; Object.assign(entry, details);
      } catch (error) {
        entry.converged = false; entry.reason = error.message;
        // Preserve the rejected geometry for inspection, with its actual
        // whole-cell evidence. It never becomes an Euler initial state.
        entry.quality = region.geometry.quality(nodes[g], region.directions);
      }
    }
    report.converged = report.regions.every(r => r.converged);
    report.geometryCertified = report.regions.every(r => r.quality?.valid);
    if (!report.converged) report.reason = report.regions.filter(r => !r.converged).map(r => `Region ${r.region + 1}: ${r.reason}`).join('; ');
    let maximumLocalIntervals = 0, maximumInteriorLocalIntervals = 0, bodyIntrusions = 0, maximumMove = 0;
    for (let g = 0; g < nodes.length; g++) {
      const weights = prepared.input.weights[g], total = weights.reduce((s, v) => s + v, 0), labels = [prepared.input.captureLevels[g]];
      const span = prepared.input.captureLevels[g + 1] - labels[0], widths = weights.map(w => w / total * span);
      for (const width of widths) labels.push(labels.at(-1) + width);
      nodes[g].forEach((row, i) => row.forEach((p, j) => {
        const width = Math.min(j ? widths[j - 1] : Infinity, j < widths.length ? widths[j] : Infinity);
        const error = Math.abs(guideField.streamfunctionAt(p) - labels[j]) / width;
        if (!Number.isFinite(error)) throw new Error('Panel comparison is nonfinite on the curved candidate.');
        maximumLocalIntervals = Math.max(maximumLocalIntervals, error);
        const interior = i > 0 && i < nodes[g].length - 1 && j > 0 && j < row.length - 1;
        if (interior) {
          maximumInteriorLocalIntervals = Math.max(maximumInteriorLocalIntervals, error);
          if (!guideField.admissibleNode(p)) bodyIntrusions++;
        }
        const q = prepared.nodes[g][i][j]; maximumMove = Math.max(maximumMove, Math.hypot(p.x - q.x, p.y - q.y));
      }));
    }
    report.accuracy = { maximumLocalIntervals, maximumInteriorLocalIntervals, limit: .05, passed: maximumLocalIntervals <= .05,
      scope: 'Nodal streamfunction comparison against the same panel field; physical refinement is still required.' };
    report.bodyIntrusions = bodyIntrusions; report.maximumMove = maximumMove;
  } catch (error) { report.converged = false; report.reason = error.message; }
  delete report.current;
  let curves;
  if (regions) {
    try { curves = edgeCurves(system.layout, nodes, regions); }
    catch (error) { report.displayError = error.message; }
  }
  const mesh = snapshot('experimental', curves);
  return { ...prepared, nodes, mesh, diagnostics: { ...prepared.diagnostics, gridSmoothing: report },
    status: 'experimental curved mesh preview; physical acceptance incomplete' };
}
