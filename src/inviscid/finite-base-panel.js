// SPDX-License-Identifier: GPL-2.0-or-later
// XFOIL-style nodal streamfunction/Kutta system with explicit finite bases.
// Segment-wise TECALC projections retain the actual supplied base polyline.
import { createContourTopology } from '../geometry/contour-topology.js';
import { createContourCurve } from '../geometry/contour-curve.js';
import { validateAssembly, pointInside } from '../geometry/airfoil.js';
import { makePanel } from './panel.js';
import { vortexStreamfunctionBasis } from './streamfunction.js';
import { basePanelStreamfunctionBasis, xfoilSurfaceDerivatives } from './finite-base-influence.js';
import { assemblyFiniteBaseSourceChart } from './finite-base-source-chart.js';
import { solveLinear, linearResidual, normInf } from '../numerics/linear.js';

export function solveFiniteBaseInviscid({ elements, alpha = 0, mach = 0, viscous = false,
  referenceChord = 1, momentReference = { x: .25, y: 0 }, boundaryCondition = 'streamfunction' }, vortexBasis) {
  if (mach !== 0 || viscous) throw new Error('Finite-base panel initialization is incompressible and inviscid.');
  if (!Number.isFinite(alpha) || Math.abs(alpha) > 90 || !Number.isFinite(referenceChord) || referenceChord <= 0
    || ![momentReference.x, momentReference.y].every(Number.isFinite) || !Array.isArray(elements)
    || !elements.length || elements.length > 6) throw new Error('Invalid finite-base panel input.');
  if (!['normal-velocity', 'streamfunction'].includes(boundaryCondition)) throw new Error('Unknown panel boundary condition.');
  const topologies = elements.map(e => createContourTopology(e.points, { trailingEdge: e.trailingEdge }));
  validateAssembly(topologies.map(t => t.points));
  const panels = [], basePanels = [], ranges = [];
  let strengthCount = 0;
  for (let element = 0; element < elements.length; element++) {
    const topology = topologies[element], points = topology.surface.points, start = strengthCount, first = panels.length;
    for (let i = 1; i < points.length; i++) panels.push({ ...makePanel(points[i - 1], points[i], element),
      node: start + i - 1, sourcePanelIndex: topology.surface.panels[i - 1].sourcePanelIndex });
    strengthCount += points.length;
    const range = { start, end: strengthCount - 1, first, last: panels.length - 1, topology,
      baseFirst: basePanels.length, teDerivative: null, sourceChart: null };
    if (topology.kind === 'finite-base') {
      const d = xfoilSurfaceDerivatives(points).derivatives;
      const t = { x: .5 * (-d[0].x + d.at(-1).x), y: .5 * (-d[0].y + d.at(-1).y) }, length = Math.hypot(t.x, t.y);
      if (!(length > 0)) throw new Error('Degenerate finite-base XFOIL TE derivative.');
      range.sourceChart = assemblyFiniteBaseSourceChart(topology, { x: t.x / length, y: t.y / length },
        topologies.filter((_, k) => k !== element).map(t => t.points));
      const cutDirection = range.sourceChart.direction;
      range.teDerivative = { vector: t, cutDirection,
        definition: 'SPLIND zero-third-derivative surface endpoints; half sum of downstream derivatives, without unit renormalization.' };
      for (const base of topology.base.panels) {
        const panel = makePanel(base.start, base.end, element);
        // Our CCW circulation sign is opposite XFOIL's GAM=QINV sign.
        const sourceCoefficient = -.5 * (t.x * panel.ty - t.y * panel.tx);
        const vortexCoefficient = -.5 * (t.x * panel.tx + t.y * panel.ty);
        basePanels.push({ ...panel, sourcePanelIndex: base.sourcePanelIndex, upperNode: start, lowerNode: range.end,
          sourceCoefficient, vortexCoefficient, cutDirection, cutOrigin: { ...range.sourceChart.origin } });
      }
    }
    ranges.push(range);
  }
  // Known base influences add no matrix unknowns. Preserve the existing
  // 700 wetted-panel budget while reporting all retained geometry panels.
  if (panels.length > 700) throw new Error('The finite-base reference solver supports at most 700 wetted-surface panels.');
  const size = strengthCount + elements.length, a = new Float64Array(size * size), rhs = new Float64Array(size);
  const radians = alpha * Math.PI / 180, u = Math.cos(radians), v = Math.sin(radians);
  const psiRows = [], kuttaRows = [], teProbes = [];
  const assemblePsi = (row, point, element) => {
    rhs[row] = (-u * point.y + v * point.x) / referenceChord;
    a[row * size + strengthCount + element] = -1;
    for (const panel of panels) {
      const b = vortexStreamfunctionBasis(point, panel);
      for (let k = 0; k < 2; k++) a[row * size + panel.node + k] += b[k] / referenceChord;
    }
    for (const panel of basePanels) {
      const b = basePanelStreamfunctionBasis(point, panel);
      if (b.onSourceCut && panel.sourceCoefficient !== 0) throw new Error('A surface node lies on another finite-base downstream source-cut ray.');
      const coefficient = (panel.sourceCoefficient * b.source + panel.vortexCoefficient * b.vortex) / referenceChord;
      a[row * size + panel.upperNode] += coefficient;
      a[row * size + panel.lowerNode] -= coefficient;
    }
    psiRows.push(row);
  };
  for (let element = 0; element < ranges.length; element++) {
    const range = ranges[element], topology = range.topology;
    for (let i = 0; i < topology.surface.points.length; i++) {
      const row = range.start + i;
      if (topology.kind === 'sharp' && i === topology.surface.points.length - 1) {
        const curve = createContourCurve(topology.surface.points), upper = curve.evaluate(0).derivative, lower = curve.evaluate(curve.length).derivative;
        const tx = -upper.x / Math.hypot(upper.x, upper.y) + lower.x / Math.hypot(lower.x, lower.y);
        const ty = -upper.y / Math.hypot(upper.x, upper.y) + lower.y / Math.hypot(lower.x, lower.y), length = Math.hypot(tx, ty);
        const tangent = { x: tx / length, y: ty / length }, distance = .1 * Math.min(panels[range.first].length, panels[range.last].length);
        const te = topology.surface.points[0], point = { x: te.x - distance * tangent.x, y: te.y - distance * tangent.y };
        if (!(length > 0) || !pointInside(point, topology.points)) throw new Error('Invalid sharp-TE interior control in finite-base assembly.');
        rhs[row] = -u * tangent.x - v * tangent.y;
        for (const panel of panels) vortexBasis(point, panel).forEach((q, k) => {
          a[row * size + panel.node + k] += q.u * tangent.x + q.v * tangent.y;
        });
        // The other elements' finite-base fields also reach this control.
        for (const panel of basePanels) {
          const q = vortexBasis(point, panel), vortex = { u: q[0].u + q[1].u, v: q[0].v + q[1].v };
          const source = { u: vortex.v, v: -vortex.u };
          const coefficient = panel.sourceCoefficient * (source.u * tangent.x + source.v * tangent.y)
            + panel.vortexCoefficient * (vortex.u * tangent.x + vortex.v * tangent.y);
          a[row * size + panel.upperNode] += coefficient; a[row * size + panel.lowerNode] -= coefficient;
        }
        teProbes.push({ element, row, point, tangent });
      } else assemblePsi(row, topology.surface.points[i], element);
    }
    const row = strengthCount + element;
    a[row * size + range.start] = 1; a[row * size + range.end] = 1; kuttaRows.push(row);
  }
  const solution = solveLinear(a, rhs), residual = linearResidual(a, solution, rhs), gamma = solution.slice(0, strengthCount);
  for (const panel of basePanels) {
    const difference = gamma[panel.upperNode] - gamma[panel.lowerNode];
    panel.sourceStrength = panel.sourceCoefficient * difference;
    panel.vortexStrength = panel.vortexCoefficient * difference;
  }
  const outputElements = ranges.map((range, element) => {
    const topology = range.topology, surfaceGamma = Array.from(gamma.slice(range.start, range.end + 1));
    const nodalCp = surfaceGamma.map(g => 1 - g * g), cp = Array(topology.sourcePanelCount);
    let cx = 0, cy = 0, cm = 0, circulation = 0, baseSourceFlux = 0, baseCirculation = 0, baseCx = 0, baseCy = 0;
    const integrate = (panel, cp0, cp1, details) => {
      const pressure = .5 * (cp0 + cp1), fx = -pressure * panel.nx * panel.length / referenceChord,
        fy = -pressure * panel.ny * panel.length / referenceChord;
      cx += fx; cy += fy;
      cm -= ((panel.x - momentReference.x) * fy - (panel.y - momentReference.y) * fx) / referenceChord
        + (cp1 - cp0) * panel.length ** 2 / (12 * referenceChord ** 2);
      cp[panel.sourcePanelIndex] = { x: panel.x, y: panel.y, cp: pressure, length: panel.length, ...details };
      return { fx, fy };
    };
    for (let i = range.first; i <= range.last; i++) {
      const panel = panels[i], k = panel.node - range.start;
      integrate(panel, nodalCp[k], nodalCp[k + 1], { qt: .5 * (surfaceGamma[k] + surfaceGamma[k + 1]), base: false });
      circulation += .5 * (surfaceGamma[k] + surfaceGamma[k + 1]) * panel.length;
    }
    const baseCp = .5 * (nodalCp[0] + nodalCp.at(-1));
    for (const panel of basePanels.filter(p => p.element === element)) {
      const force = integrate(panel, baseCp, baseCp, { qt: null, base: true, pressureModel: 'equal-endpoint-TE-pressure' });
      baseCx += force.fx; baseCy += force.fy;
      baseSourceFlux += panel.sourceStrength * panel.length;
      baseCirculation += panel.vortexStrength * panel.length;
    }
    circulation += baseCirculation;
    return { name: elements[element].name ?? `Element ${element + 1}`, points: topology.points, cp, cx, cy, cm,
      cl: -cx * v + cy * u, pressureDrag: cx * u + cy * v, circulation,
      surfaceGamma, nodalCp, surfacePoints: topology.surface.points,
      finiteBase: { kind: topology.kind, retainedSegments: topology.base.panels.length, sourceFlux: baseSourceFlux,
        circulation: baseCirculation, pressure: baseCp, pressureModel: 'XFOIL equal endpoint TE Cp; modeled closure, not measured recirculation pressure',
        force: { cx: baseCx, cy: baseCy }, teDerivative: range.teDerivative,
        surfaceStreamfunction: referenceChord * solution[strengthCount + element],
        sourceCutJump: baseSourceFlux, cutDirection: range.teDerivative?.cutDirection,
        cutOrigin: range.sourceChart ? { ...range.sourceChart.origin } : undefined,
        sourceChart: range.sourceChart } };
  });
  const total = key => outputElements.reduce((sum, e) => sum + e[key], 0);
  const diagnostics = { linearResidual: normInf(residual), nodalStreamfunctionResidual: Math.max(...psiRows.map(i => Math.abs(residual[i]))),
    kuttaResidual: Math.max(...kuttaRows.map(i => Math.abs(residual[i]))),
    interiorTEResidual: teProbes.length ? Math.max(...teProbes.map(p => Math.abs(residual[p.row]))) : 0,
    surfaceStreamfunctions: Array.from(solution.slice(strengthCount), p => referenceChord * p), teProbes,
    circulationLift: -2 * total('circulation') / referenceChord, liftMismatch: Math.abs(total('cl') + 2 * total('circulation') / referenceChord),
    pressureDrag: total('pressureDrag'), baseSourceFlux: basePanels.reduce((s, p) => s + p.sourceStrength * p.length, 0),
    baseCirculation: basePanels.reduce((s, p) => s + p.vortexStrength * p.length, 0),
    unknowns: size, strengthCount, streamfunctionConstants: elements.length, wettedPanelCount: panels.length,
    basePanelCount: basePanels.length, surfaceRows: psiRows.length, kuttaRows: kuttaRows.length,
    finiteBaseConvention: 'Segment-wise TECALC projections on retained geometry; SPLIND endpoint derivatives. Source streamfunction has explicit downstream cuts; nonzero base source flux is intentional.' };
  return { model: 'incompressible-linear-vortex', finiteBaseModel: 'xfoil-segmented-base', status: diagnostics.linearResidual < 1e-9 ? 'solved' : 'failed',
    alpha, mach: 0, boundaryCondition: 'streamfunction', requestedBoundaryCondition: boundaryCondition,
    referenceChord, momentReference, panelCount: panels.length + basePanels.length,
    cl: total('cl'), cm: total('cm'), cd: null, elements: outputElements, diagnostics,
    warnings: ['Finite-base inviscid panel closure; base pressure is modeled from TE pressure. This is not a viscous/base-drag or NLR experimental validation.'],
    field: { panels, basePanels, gamma, u, v } };
}
