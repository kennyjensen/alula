// SPDX-License-Identifier: GPL-2.0-or-later
// Research initializer for a SUPPLIED body/cut topology. Optional potential
// spacing matches its cross-line counts; automatic block topology is pending.
// It traces a simultaneous incompressible panel field and distributes each
// cross-line by ordered, locally unwrapped velocity potential, or seeds a
// linear interior for the companion fixed-boundary elliptic SLOR stage.
import { createStreamtubeBodySystem } from './streamtube-body.js';
import { initializeStreamtubeDensities } from './streamtube-initial-state.js';
import { solveInviscid, velocityAt } from '../inviscid/linear-vortex.js';
import { streamfunctionAt, potentialDifference, streamfunctionBranchIncrement } from '../inviscid/streamfunction.js';
import { finiteBaseWakePotentialDifference } from '../inviscid/finite-base-wake-potential.js';
import { createContourPotential } from '../inviscid/contour-potential.js';
import { createSheetSurfaceCoordinate } from '../inviscid/sheet-surface-coordinate.js';
import { createStagnationGuideConnector } from '../inviscid/stagnation-guide-connector.js';
import { refineSurfaceStagnation } from '../inviscid/surface-stagnation.js';
import { tracePotentialCurve } from '../inviscid/potential-curve.js';
import { traceArcCurve } from '../inviscid/arc-curve.js';
import { createContourCurve } from '../geometry/contour-curve.js';
import { createContourTopology, createSurfaceContourCurve } from '../geometry/contour-topology.js';
import { createStreamtubeDisplacement } from './streamtube-displacement.js';
import { recoverFiniteBaseWakeTangency, recoverFiniteBaseWakeInterior } from './streamtube-finite-base-wake-initializer.js';
import { prepareContour, pointInside, segmentsTouch } from '../geometry/airfoil.js';
import { matchPanelPotentialGuides } from './streamtube-potential-guides.js';
import { stagnationStreamtubeMass, resolveNormalMassWeights } from '../geometry/streamtube-normal-spacing.js';
import { fitNormalStreamtubeFraction } from '../geometry/streamtube-spacing-fit.js';
import { distributeCurvatureSurface, distributeSourceSurface } from '../geometry/curvature-surface-spacing.js';
import { createContourArc } from '../geometry/contour-arc.js';
import { createBoundaryArcStationMap } from '../geometry/boundary-arc-stations.js';
import { sampleUniqueGuideX } from '../geometry/guide-x-bracket.js';
import { validatePanelTopologySeed } from './streamtube-panel-topology.js';

// Stagnation-normal allocation belongs to the panel grid initializer: it
// consumes its panel profiles and guide curves and has no independent caller.
function matchStagnationNormalSpacing({ input, profiles, guides, outer, aspectRatio, pointAt, expandOnly = false, maximumAdjacentRatio = Infinity, maximumTubes = 63 }) {
  if (!(maximumAdjacentRatio > 1) || !Number.isInteger(maximumTubes) || maximumTubes < 3)
    throw new Error('Invalid normal-spacing refinement controls.');
  const groups = input.weights.map((original, g) => {
    const lower = g === 0 ? outer[0] : guides[g - 1].upper;
    const upper = g === profiles.length ? outer[1] : guides[g].lower;
    const mass = input.captureLevels[g + 1] - input.captureLevels[g], targets = {};
    const total = original.reduce((a, b) => a + b, 0);
    for (const [end, b, side, wall] of [['first', g - 1, 'upper', lower], ['last', g, 'lower', upper]]) {
      if (b < 0 || b >= profiles.length) continue;
      const i = input.bodies[b].leadingIndex, along = Math.hypot(wall[i + 1].x - wall[i].x, wall[i + 1].y - wall[i].y);
      const potentialSlope = (upper[i].potential - lower[i].potential) / mass;
      if (expandOnly) {
        const fraction = (end === 'first' ? original[0] : original.at(-1)) / total;
        const eta = end === 'first' ? fraction : 1 - fraction;
        const p = pointAt(g, i, eta, (1 - eta) * lower[i].potential + eta * upper[i].potential);
        const distance = Math.hypot(p.x - wall[i].x, p.y - wall[i].y);
        if (distance >= aspectRatio * along) {
          targets[end] = { body: b, side, along, fraction, distance, retained: true };
          continue;
        }
      }
      const target = stagnationStreamtubeMass({ strainRate: profiles[b].strainRate, distance: aspectRatio * along, potentialSlope });
      let fit;
      try {
        fit = fitNormalStreamtubeFraction({ initialFraction: target / mass, targetDistance: aspectRatio * along,
          distanceAtFraction: fraction => {
            const eta = end === 'first' ? fraction : 1 - fraction;
            const p = pointAt(g, i, eta, (1 - eta) * lower[i].potential + eta * upper[i].potential);
            return Math.hypot(p.x - wall[i].x, p.y - wall[i].y);
          } });
      } catch (error) { throw new Error(`Body ${b + 1} ${side} normal-spacing fit: ${error.message}`); }
      targets[end] = { body: b, side, along, strainRate: profiles[b].strainRate, potentialSlope,
        linearMass: target, mass: mass * fit.fraction, fraction: fit.fraction, fit };
    }
    const distribution = expandOnly && Object.values(targets).every(t => t.retained)
      ? { weights: original.slice(), requestedTubes: original.length, addedTubes: 0, retained: true }
      : resolveNormalMassWeights({ count: original.length, first: targets.first?.fraction, last: targets.last?.fraction,
      maximumAdjacentRatio, maximumTubes });
    return { mass, targets, ...distribution };
  });
  input.weights = groups.map(g => g.weights.slice());
  return { model: 'panel distance fit initialized by linear stagnation flow', aspectRatio, expandOnly, groups };
}

export function createPanelStreamtubeGrid(supplied, { relativeTolerance = 2e-9,
  seedDistance = .002, maxSpatialStep = .1, crosslinePlacement = 'supplied-x', crosslineGrowth = .25,
  maxCrosslineIntervals = 2000, outerCrosslineSpread = 0, normalSpacing = 'supplied', stagnationAspectRatio = 2.5,
  panelBoundaryCondition = 'normal-velocity', panelSolution, surfaceSpacing = 'supplied', curvatureSpacing = {}, interiorInitialization = 'traced', resolveSurfaceTurning = false,
  surfaceCrosslineMetric = 'potential', cutStationSpacing = 'potential', surfaceStationPlacement = 'common-rank', passageIntervalCounts,
  passageCountPlanning = false, interiorStationPlacement = 'potential', recordPotentialCoordinates = false, onProgress } = {}) {
  if (![relativeTolerance, seedDistance, maxSpatialStep].every(Number.isFinite)
    || Math.min(relativeTolerance, seedDistance, maxSpatialStep) <= 0)
    throw new Error('Invalid panel streamtube initialization controls.');
  if (!['supplied-x', 'potential'].includes(crosslinePlacement)) throw new Error('Unknown panel cross-line placement.');
  if (!['potential', 'arc'].includes(surfaceCrosslineMetric)) throw new Error('Unknown surface cross-line metric.');
  if (!['potential', 'physical-x'].includes(cutStationSpacing) || cutStationSpacing === 'physical-x' && crosslinePlacement !== 'potential')
    throw new Error('Physical cut station spacing requires potential-block crosslines.');
  if (cutStationSpacing === 'physical-x' && surfaceCrosslineMetric !== 'arc')
    throw new Error('Physical cut station spacing requires the contour-arc surface metric.');
  if (!['common-rank', 'local-density', 'local-density-joined', 'passage-density'].includes(surfaceStationPlacement)
    || surfaceStationPlacement !== 'common-rank' && (crosslinePlacement !== 'potential' || surfaceCrosslineMetric !== 'arc'))
    throw new Error('Local surface density placement requires contour-arc potential-block maps.');
  if (!Number.isFinite(crosslineGrowth) || crosslineGrowth <= 0 || !Number.isInteger(maxCrosslineIntervals) || maxCrosslineIntervals < 4)
    throw new Error('Invalid potential cross-line spacing controls.');
  if (!Number.isFinite(outerCrosslineSpread) || outerCrosslineSpread < 0 || outerCrosslineSpread > 1)
    throw new Error('Outer cross-line spreading must lie between zero and one.');
  if (!['supplied', 'stagnation', 'automatic'].includes(normalSpacing) || !Number.isFinite(stagnationAspectRatio) || stagnationAspectRatio <= 0)
    throw new Error('Invalid normal-spacing controls.');
  if (!['supplied', 'curvature', 'source'].includes(surfaceSpacing) || !curvatureSpacing || typeof curvatureSpacing !== 'object' || Array.isArray(curvatureSpacing))
    throw new Error('Invalid surface-spacing controls.');
  if (!['traced', 'linear'].includes(interiorInitialization)) throw new Error('Unknown interior-grid initialization.');
  if (!['potential', 'boundary-arc'].includes(interiorStationPlacement)
    || interiorStationPlacement === 'boundary-arc' && interiorInitialization !== 'traced')
    throw new Error('Boundary-arc station placement requires traced interior initialization.');
  if (supplied.potentialTopologyRequired && (crosslinePlacement !== 'potential' || panelSolution === undefined))
    throw new Error('Reentrant provisional topology requires its shared panel potential chart.');
  // Validate the supplied topology/geometry before running the panel solve.
  createStreamtubeBodySystem(supplied);
  const input = structuredClone(supplied);
  const bodyTopologies = input.bodies.map(b => b.trailingEdge?.kind === 'finite-base' ? createContourTopology(b.points, b) : null);
  input.bodies.forEach((b, k) => { b.points = bodyTopologies[k]?.points ?? prepareContour(b.points); });
  if (bodyTopologies.some(Boolean)) panelBoundaryCondition = 'streamfunction';
  if (input.bodies.some(b => b.roundedLeadingEdge === false))
    throw new Error('Panel streamtube initialization currently requires movable rounded leading edges.');
  const primary = input.bodies[input.primaryBody ?? 0].points;
  const scale = Math.max(...primary.map(p => Math.hypot(p.x - primary[0].x, p.y - primary[0].y)));
  const tolerance = relativeTolerance * scale, psiTolerance = Math.min(1e-11 * scale, tolerance);
  const xmin = input.outerLower[0].x, xmax = input.outerLower.at(-1).x, reference = input.outerLower[0];
  const paths = [input.outerLower, input.outerUpper, ...input.cutPaths];
  if (!(xmax > xmin) || paths.some(path => Math.abs(path[0].x - xmin) > tolerance || Math.abs(path.at(-1).x - xmax) > tolerance
    || path.some((p, i) => i && p.x <= path[i - 1].x)))
    throw new Error('Panel streamtube initialization requires ordered x paths and common vertical inlet/outlet planes.');
  const panel = panelSolution === undefined
    ? solveInviscid({ elements: input.bodies.map(b => ({ points: b.points, ...(b.trailingEdge ? { trailingEdge: b.trailingEdge } : {}) })), alpha: input.alpha, boundaryCondition: panelBoundaryCondition })
    : validatePanelTopologySeed(panelSolution, input);
  if (panelSolution !== undefined && panel.boundaryCondition !== panelBoundaryCondition)
    throw new Error('Prepared panel topology has a different boundary condition.');
  if (panel.status !== 'solved') throw new Error('Panel initializer did not solve its simultaneous boundary conditions.');
  const { field } = panel, psi = p => streamfunctionAt(p, field), velocity = p => velocityAt(p, field);
  const fluid = p => !input.bodies.some(b => pointInside(p, b.points));
  const boundaryPanels = [...field.panels, ...(field.basePanels ?? [])];
  const segment = (a, b) => !boundaryPanels.some(p => segmentsTouch(a, b, p.a, p.b, 1e-12 * scale));
  const diagnostics = { panelResidual: panel.diagnostics.linearResidual, panelBoundaryCondition, relativeTolerance, seedDistance, crosslinePlacement, cutStationSpacing, normalSpacing,
    ...(input.gridSpacing ? { gridSpacing: input.gridSpacing } : {}),
    traces: 0, acceptedSteps: 0, rejectedSteps: 0, cutTraceExtensions: 0,
    maxStreamfunctionDrift: 0, maxSurfaceStreamfunctionDefect: 0, profiles: [] };
  const geometricCuts = new WeakMap();
  const trace = (seed, initialPotential, endPotential, plane) => {
    const r = tracePotentialCurve({ seed, initialPotential, endPotential, plane, velocity, tolerance,
      maxStep: .05 * scale, maxSpatialStep: maxSpatialStep * scale, minStep: 1e-12 * scale,
      admissible: fluid, admissibleSegment: segment });
    diagnostics.traces++; diagnostics.acceptedSteps += r.accepted; diagnostics.rejectedSteps += r.rejected;
    if (!r.converged) throw Object.assign(new Error(`Panel streamline trace: ${r.reason}.`),
      { code: 'streamtube-guide-trace', trace: { seed, initialPotential, endPotential, plane, result: r } });
    return r.points;
  };
  const traceArc = (seed, initialPotential, initialArc, endArc, plane) => {
    const r = traceArcCurve({ seed, initialPotential, initialArc, endArc, plane, velocity, tolerance, potentialTolerance: tolerance,
      maxStep: maxSpatialStep * scale, minStep: 1e-12 * scale, admissible: fluid, admissibleSegment: segment });
    diagnostics.traces++; diagnostics.acceptedSteps += r.accepted; diagnostics.rejectedSteps += r.rejected;
    diagnostics.arcTraces = (diagnostics.arcTraces ?? 0) + 1;
    if (!r.converged) throw new Error(`Panel streamline arc trace: ${r.reason}.`);
    return r.points;
  };
  const sampleArc = (path, arc) => {
    if (!Number.isFinite(arc) || arc < path[0].arc || arc > path.at(-1).arc)
      throw new Error('Cross-line arc lies outside its traced streamline.');
    if (arc === path[0].arc) return { ...path[0] };
    if (arc === path.at(-1).arc) return { ...path.at(-1) };
    let lo = 0, hi = path.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (path[mid].arc < arc) lo = mid; else hi = mid; }
    const a = path[lo], b = path[hi];
    if (arc === a.arc) return { ...a }; if (arc === b.arc) return { ...b };
    return traceArc(a, a.potential, a.arc, arc).at(-1);
  };
  const sample = (path, potential) => {
    if (potential < path[0].potential - 8 * tolerance || potential > path.at(-1).potential + 8 * tolerance)
      throw new Error('Cross-line potential lies outside its traced streamline.');
    // Plane events have finite integration tolerance; use the actual traced
    // endpoint when a requested endpoint differs only within that tolerance.
    if (potential <= path[0].potential) return { ...path[0] };
    if (potential >= path.at(-1).potential) return { ...path.at(-1) };
    const geometric = geometricCuts.get(path);
    if (geometric && potential >= geometric.join.potential) return geometric.atPotential(potential);
    let lo = 0, hi = path.length - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (path[m].potential < potential) lo = m; else hi = m; }
    const a = path[lo], b = path[hi];
    if (potential === a.potential) return { ...a }; if (potential === b.potential) return { ...b };
    return trace(a, a.potential, potential).at(-1);
  };
  const sampleCut = (path, potential) => {
    // New potential stations can be closer to an edge than the original
    // x-based seed. Extend the actual streamline toward that station using
    // the same error and body-intersection checks; never extrapolate it.
    if (potential < path[0].potential - 8 * tolerance) {
      const r = trace(path[0], path[0].potential, potential).reverse();
      path.unshift(...r.slice(0, -1)); diagnostics.cutTraceExtensions++;
    } else if (potential > path.at(-1).potential + 8 * tolerance) {
      const r = trace(path.at(-1), path.at(-1).potential, potential);
      path.push(...r.slice(1)); diagnostics.cutTraceExtensions++;
    }
    return sample(path, potential);
  };
  const atX = (path, x) => {
    const geometric = geometricCuts.get(path);
    if (geometric && x >= geometric.join.x) return geometric.atX(x);
    if (path.some((p, i) => i && p.x <= path[i - 1].x)) {
      const point = sampleUniqueGuideX(path, x, { sample, velocity, tolerance });
      diagnostics.uniqueGuideIntersections ??= { count: 0, method: 'Unique resolved x-plane bracket on a potential-ordered guide',
        ambiguityGuard: 'Reject multiple intersections and vertical intervals at the requested x; retain every original guide point.' };
      diagnostics.uniqueGuideIntersections.count++;
      return point;
    }
    if (Math.abs(x - path[0].x) <= tolerance) return { ...path[0] };
    if (Math.abs(x - path.at(-1).x) <= tolerance) return { ...path.at(-1) };
    if (x < path[0].x || x > path.at(-1).x) throw new Error('Requested x lies outside its traced guide.');
    let lo = 0, hi = path.length - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (path[m].x < x) lo = m; else hi = m; }
    let lower = path[lo].potential, upper = path[hi].potential;
    let potential = lower + (upper - lower) * (x - path[lo].x) / (path[hi].x - path[lo].x);
    for (let k = 0; k < 40; k++) {
      const p = sample(path, potential), error = p.x - x; if (Math.abs(error) <= tolerance) return p;
      if (error < 0) lower = potential; else upper = potential;
      const q = velocity(p), next = potential - error * (q.u * q.u + q.v * q.v) / q.u;
      potential = next > lower && next < upper ? next : .5 * (lower + upper);
    }
    throw new Error('Streamline/x-plane intersection did not converge.');
  };
  const sampleCutX = (path, x) => {
    const geometric = geometricCuts.get(path);
    if (geometric && x >= geometric.join.x) return geometric.atX(x);
    // Extend an actual traced path to a requested physical x-plane. This
    // also serves continuous-boundary reference sampling below; no chord
    // interpolation or positive potential slope is imposed near stagnation.
    try {
    // A potential-ordered incoming guide can turn back in x near the
    // panel/material stagnation mismatch. Its last point is then not its
    // maximum x. Resolve existing brackets before attempting an extension;
    // multiple intersections remain an explicit error.
    if (path.some((p, i) => i && p.x <= path[i - 1].x)
      && x >= Math.min(...path.map(p => p.x)) && x <= Math.max(...path.map(p => p.x)))
      return atX(path, x);
    if (x > path.at(-1).x + tolerance) {
      const p = path.at(-1), r = trace(p, p.potential, p.potential + potentialBudget, { normal: { x: 1, y: 0 }, offset: x });
      path.push(...r.slice(1)); diagnostics.cutTraceExtensions++;
    } else if (x < path[0].x - tolerance) {
      const p = path[0], r = trace(p, p.potential, p.potential - potentialBudget, { normal: { x: 1, y: 0 }, offset: x }).reverse();
      path.unshift(...r.slice(0, -1)); diagnostics.cutTraceExtensions++;
    }
    return atX(path, x);
    } catch (error) {
      // A different C2 material anchor and panel-Dirichlet contour need not
      // share an exact incoming streamline. Only this proven, unresolved
      // near-anchor case can request a geometric initializer connector.
      const body = profiles.findIndex(p => p.upstream === path);
      const profile = profiles[body], anchor = profile?.curve.evaluate(profile.stag).point;
      const sourcePanel = profile ? field.panels.filter(p => p.element === body)[profile.stagnationSegment] : null;
      // Fine stations can enter the guide's double-valued part before
      // reaching its maximum x. Use the same source-panel-sized connector
      // for that mismatch; never recover unrelated or remote ambiguity.
      if (sourcePanel && !geometric
        && ['streamtube-guide-trace', 'GUIDE_X_AMBIGUOUS', 'GUIDE_X_UNRESOLVED'].includes(error.code)
        && Math.max(...path.map(p => p.x)) < anchor.x - tolerance
        && x > Math.min(sourcePanel.a.x, sourcePanel.b.x) + tolerance && x < anchor.x - tolerance)
        throw Object.assign(new Error('The panel incoming guide cannot reach its separate material stagnation reference.'),
          { code: 'STAGNATION_GUIDE_REFERENCE_MISMATCH', body, requestedX: x, cause: error });
      throw error;
    }
  };
  const onPsi = (point, target, normal) => {
    let p = { ...point };
    if (!fluid(p)) throw new Error('Streamfunction seed lies inside a body.');
    for (let k = 0; k < 40; k++) {
      const f = psi(p) - target; if (Math.abs(f) <= psiTolerance) return p;
      const q = velocity(p), d = normal ?? { x: -q.v, y: q.u }, df = -q.v * d.x + q.u * d.y;
      if (!Number.isFinite(df) || df === 0) throw new Error('Singular streamfunction seed correction.');
      const step = -f / df; let accepted = false;
      for (let factor = Math.min(1, .1 * scale / (Math.abs(step) * Math.hypot(d.x, d.y))); factor >= 2 ** -24; factor /= 2) {
        const next = { x: p.x + factor * step * d.x, y: p.y + factor * step * d.y };
        if (fluid(next) && segment(p, next) && Math.abs(psi(next) - target) < Math.abs(f)) { p = next; accepted = true; break; }
      }
      if (!accepted) throw new Error('Streamfunction seed correction could not remain in the fluid.');
    }
    throw new Error('Streamfunction seed correction did not converge.');
  };
  const inlet = level => onPsi({ x: xmin, y: reference.y + level - psi(reference) }, level, { x: 0, y: 1 });
  const potentialBudget = 4 * (xmax - xmin);
  const profiles = input.bodies.map((body, b) => {
    const curve = bodyTopologies[b] ? createSurfaceContourCurve(body.points, body) : createContourCurve(body.points);
    const cp = bodyTopologies[b] ? panel.elements[b].cp.filter(p => !p.base) : panel.elements[b].cp;
    const chord = Math.max(...body.points.map(p => Math.hypot(p.x - body.points[0].x, p.y - body.points[0].y)));
    const stations = cp.map((p, i) => ({ s: .5 * (curve.knots[i] + curve.knots[i + 1]), qt: p.qt }));
    stations.unshift({ s: 0, qt: stations[0].qt }); stations.push({ s: curve.length, qt: stations.at(-1).qt });
    const roots = [];
    for (let i = 1; i < stations.length; i++) if (stations[i - 1].qt < 0 && stations[i].qt >= 0) {
      const a = stations[i - 1], c = stations[i];
      roots.push({ s: a.s - (c.s - a.s) * a.qt / (c.qt - a.qt), derivative: (c.qt - a.qt) / (c.s - a.s), left: a.s, right: c.s });
    }
    if (roots.length !== 1) throw new Error('Panel initializer requires one resolved incoming stagnation point per body.');
    const buildProfile = (rejectedContourPath = null) => {
    // A concave spline wall need not lie outside its straight-panel polygon.
    // If that physical path is incompatible, use XFOIL's nodal sheet-speed
    // coordinate for initialization, with its own consistent stagnation zero.
    // Keep the measured wall and the strict exterior-potential guard intact.
    let sheetCoordinate = null;
    if (rejectedContourPath) {
      const source = { points: body.points, knots: curve.knots, field, element: b, trailingEdge: body.trailingEdge };
      try { sheetCoordinate = createSheetSurfaceCoordinate(source); }
      catch (error) {
        if (!input.potentialTopologyRequired || error.code !== 'SHEET_SURFACE_NONMONOTONE') throw error;
        sheetCoordinate = createSheetSurfaceCoordinate({ ...source, reconstruction: 'conservative-panel' });
        diagnostics.conservativeSurfaceCoordinates ??= [];
        diagnostics.conservativeSurfaceCoordinates.push({ body: b, ...sheetCoordinate.diagnostics });
      }
    }
    let stag = roots[0].s, stagnationRefinement = null;
    // Potential-block intersections require a monotone analytic surface
    // potential. Away from panel vertices, refine the actual tangential-
    // velocity zero on the C2 wall. At a nonsmooth panel junction only a
    // tightly resolved potential-minimum bracket may exist; record that
    // distinction instead of inventing a finite endpoint velocity.
    if (sheetCoordinate) {
      stag = sheetCoordinate.stagnationParameter;
      stagnationRefinement = { parameter: stag, kind: 'sheet-strength-zero',
        derivative: sheetCoordinate.derivative(stag), derivativeStatus: 'nodal surface-speed initializer' };
    } else if (crosslinePlacement === 'potential') {
      stagnationRefinement = refineSurfaceStagnation({ curve, field, left: roots[0].left, right: roots[0].right,
        chord, relativeTolerance });
      stag = stagnationRefinement.parameter;
    }
    if (!sheetCoordinate && stations.some(p => p.s < stag ? p.qt >= 0 : p.s > stag && p.qt <= 0))
      throw new Error('Panel surface potential is not monotone on both stagnation-to-TE branches.');
    // Surface potential is an initializer approximation: integrate the
    // piecewise-linear collocation velocity, with independent TE branches.
    // Neither the spline wall nor this profile is an exact panel streamline.
    const integral = s => {
      let value = 0;
      for (let i = 1; i < stations.length && s > stations[i - 1].s; i++) {
        const a = stations[i - 1], c = stations[i], d = Math.min(s, c.s) - a.s;
        value += d * a.qt + .5 * d * d * (c.qt - a.qt) / (c.s - a.s);
      }
      return value;
    };
    const value = curve.evaluate(stag), length = Math.hypot(value.derivative.x, value.derivative.y);
    if (sheetCoordinate) {
      const p = value.point;
      const vertex = field.panels.some(panel => [panel.a, panel.b].some(q => p.x === q.x && p.y === q.y));
      if (!vertex && !fluid(p)) throw new Error('The sheet-speed stagnation anchor lies inside the panel body.');
      if (!vertex && field.panels.some(panel => {
        const dx = panel.b.x - panel.a.x, dy = panel.b.y - panel.a.y;
        const px = p.x - panel.a.x, py = p.y - panel.a.y, dot = px * dx + py * dy;
        return px * dy - py * dx === 0 && dot > 0 && dot < dx * dx + dy * dy;
      })) throw new Error('The sheet-speed stagnation anchor needs an explicit one-sided panel-interior limit.');
    }
    const strainRate = (sheetCoordinate?.diagnostics.stagnationSlope ?? roots[0].derivative) / length;
    const normal = { x: value.derivative.y / length, y: -value.derivative.x / length };
    // Traced dividing/wake paths must cover the nearest requested cut
    // cross-lines. A fixed fraction of chord can start beyond them when
    // LE/TE intervals are refined. Keep the seed inside both cut intervals.
    const upstreamGap = value.point.x - input.cutPaths[b][body.leadingIndex - 1].x;
    const downstreamGap = input.cutPaths[b][body.trailingIndex + 1].x - body.points[0].x;
    const epsilon = Math.min(seedDistance * chord, .2 * upstreamGap, .2 * downstreamGap);
    const sampledStagnationStreamfunction = psi(value.point);
    // A streamfunction BIE solves one Dirichlet constant per body. Preserve
    // that solved value: the interpolated stagnation point on the C2 curve
    // is not an exact collocation point of the straight-panel field.
    const solvedSurfaceStreamfunction = panel.diagnostics.surfaceStreamfunctions?.[b] ?? null;
    const level = panelBoundaryCondition === 'streamfunction' ? solvedSurfaceStreamfunction : sampledStagnationStreamfunction;
    if (!Number.isFinite(level)) throw new Error('The panel field did not supply a finite body streamfunction level.');
    if (!(epsilon > 0)) throw new Error('Dividing/wake seed requires positive adjacent cut intervals.');
    const seed = onPsi({ x: value.point.x + epsilon * normal.x, y: value.point.y + epsilon * normal.y }, level, { x: -normal.y, y: normal.x });
    const reverse = trace(seed, 0, -potentialBudget, { normal: { x: 1, y: 0 }, offset: xmin }), start = reverse.at(-1);
    const offset = potentialDifference(reference, start, field) - start.potential;
    let phiStag;
    try { phiStag = offset - potentialDifference(value.point, seed, field); }
    catch (error) { throw new Error(`Body ${b + 1} stagnation-to-seed potential: ${error.message}`); }
    const upstream = reverse.map(p => ({ ...p, potential: p.potential + offset })).reverse();
    body.stagnationParameter = stag;
    const contourPotential = crosslinePlacement === 'potential'
      ? sheetCoordinate ? createSheetSurfaceCoordinate({ points: body.points, knots: curve.knots, field,
        element: b, stagnationPotential: phiStag, trailingEdge: body.trailingEdge,
        reconstruction: sheetCoordinate.diagnostics.reconstruction ?? 'nodal-linear' })
        : createContourPotential({ curve, stagnation: stag, stagnationPotential: phiStag, field }) : null;
    if (contourPotential && !contourPotential.diagnostics.sampledMonotone)
      throw new Error('Analytic panel surface potential is not monotone on the prescribed contour.');
    diagnostics.profiles.push({ stagnationParameter: stag, collocationStagnationParameter: roots[0].s, streamfunction: level,
      streamfunctionSource: panelBoundaryCondition === 'streamfunction' ? 'solved body Dirichlet constant' : 'sampled stagnation point',
      solvedSurfaceStreamfunction, sampledStagnationStreamfunction,
      stagnationStreamfunctionDefect: sampledStagnationStreamfunction - level,
      stagnationPotential: phiStag, seedDistance: epsilon / chord, strainRate,
      strainRateSource: sheetCoordinate ? 'nodal sheet-strength slope' : 'interpolated panel collocation slope',
      ...(rejectedContourPath ? { rejectedContourPath } : {}),
      ...(stagnationRefinement ? { stagnationRefinement } : {}),
      ...(contourPotential ? { surfacePotential: contourPotential.diagnostics } : {}) });
    onProgress?.({ phase: 'dividing-streamline', body: b });
    return { curve, stag, level, inletLevel: level + streamfunctionBranchIncrement(reverse, field), phiStag, upstream, epsilon, strainRate,
      // Both analytic-contour and sheet-speed anchors belong to a source
      // panel interval. The local connector must use that same interval.
      stagnationSegment: sheetCoordinate?.diagnostics.stagnation.segment ?? curve.knots.findIndex((s, i) => i < curve.knots.length - 1 && s <= stag && stag < curve.knots[i + 1]),
      ...(sheetCoordinate ? { sheetStagnation: sheetCoordinate.diagnostics.stagnation } : {}),
      phase: contourPotential?.phase ?? (s => phiStag + integral(s) - integral(stag)) };
    };
    try { return buildProfile(); }
    catch (error) {
      if (crosslinePlacement !== 'potential' || panelBoundaryCondition !== 'streamfunction'
        || error.code !== 'CONTOUR_POTENTIAL_PATH_INCOMPATIBLE') throw error;
      return buildProfile({ code: error.code, reason: error.message, interval: error.interval });
    }
  });
  // The provisional x-union supplies connectivity only. Reusing its foreign
  // element samples as surface-spacing demands feeds near coincidences back
  // into the matcher and creates artificial refinement patches.
  const independentCounts = surfaceCrosslineMetric === 'arc' && crosslinePlacement === 'potential'
    ? input.bodies.map(b => input.gridSpacing?.requestedSurfaceIntervalsByElement?.find(q => q.element === b.element)?.intervals
      ?? b.trailingIndex - b.leadingIndex) : null;
  let surfaceDemands = independentCounts?.map(count => Object.fromEntries(['upper', 'lower'].map(side => [side,
    Array.from({ length: count + 1 }, (_, i) => .5 * (1 - Math.cos(Math.PI * i / count)))])));
  if (surfaceSpacing === 'source') {
    diagnostics.surfaceSpacing = surfaceSpacing;
    input.bodies.forEach((body, b) => {
      const demand = Object.fromEntries(['upper', 'lower'].map(side => [side,
        distributeSourceSurface({ curve: profiles[b].curve, side, stagnation: profiles[b].stag,
          count: (independentCounts?.[b] ?? body.trailingIndex - body.leadingIndex) + 1 })]));
      if (surfaceDemands) surfaceDemands[b] = demand; else body.surfaceFractions = demand;
    });
  }
  if (surfaceSpacing === 'curvature') {
    diagnostics.surfaceSpacing = surfaceSpacing; diagnostics.surfaceDistributions = [];
    input.bodies.forEach((body, b) => {
      const demand = Object.fromEntries(['upper', 'lower'].map(side => {
        const distribution = distributeCurvatureSurface({ exponent: .5,
          ...(curvatureSpacing.curvatureWeight === undefined ? { leadingSpacingRatio: .2, trailingSpacingRatio: .4 } : {}),
          ...curvatureSpacing, curve: profiles[b].curve, side, stagnation: profiles[b].stag,
          count: (independentCounts?.[b] ?? body.trailingIndex - body.leadingIndex) + 1 });
        diagnostics.surfaceDistributions.push({ body: b, ...distribution.diagnostics });
        return [side, distribution.fractions];
      }));
      if (surfaceDemands) surfaceDemands[b] = demand; else body.surfaceFractions = demand;
    });
  }
  input.captureLevels = [psi(reference), ...profiles.map(p => p.inletLevel), psi(input.outerUpper[0])];
  let system = createStreamtubeBodySystem(input), { nx, elements } = system.layout;
  profiles.forEach((profile, b) => {
    const te = input.bodies[b].points[0];
    const base = system.baseGeometry?.[b];
    let wakeSeed;
    if (base) {
      // The finite-base panel field has nonzero source flux and a downstream
      // streamfunction cut. Trace a velocity-based center guess, then apply
      // the prescribed inviscid bank width below; psi is undefined on the cut.
      const t = base.tangentDerivative, length = Math.hypot(t.x, t.y);
      wakeSeed = { x: base.center.x + profile.epsilon * t.x / length,
        y: base.center.y + profile.epsilon * t.y / length };
      if (!fluid(wakeSeed)) throw new Error('Finite-base wake-center seed is inside a solid body.');
      diagnostics.finiteBaseWakes ??= [];
      diagnostics.finiteBaseWakes.push({ body: b, width: base.width,
        method: 'Velocity-traced center guess with constant normal bank width; not a single-valued panel-psi cut.' });
    } else wakeSeed = onPsi({ x: te.x + profile.epsilon, y: te.y }, profile.level, { x: 0, y: 1 });
    let delta;
    try {
      if (base) {
        const continuation = finiteBaseWakePotentialDifference({ trailingEdge: te, wakeSeed,
          direction: base.tangentDerivative, element: b, field });
        delta = continuation.value;
        if (continuation.detour) diagnostics.finiteBaseWakes.at(-1).potentialContinuation = continuation;
      } else delta = potentialDifference(te, wakeSeed, field);
    }
    catch (error) { throw new Error(`Body ${b + 1} trailing-edge-to-wake potential: ${error.message}`); }
    profile.wake = trace(wakeSeed, delta, delta + potentialBudget, { normal: { x: 1, y: 0 }, offset: xmax });
  });
  const fullCurve = level => {
    const p = inlet(level), phi = potentialDifference(reference, p, field);
    return trace(p, phi, phi + potentialBudget, { normal: { x: 1, y: 0 }, offset: xmax });
  };
  const outerPaths = [fullCurve(input.captureLevels[0]), fullCurve(input.captureLevels.at(-1))];
  let guides, outer;
  if (crosslinePlacement === 'potential') {
    let matched;
    for (let attempt = 0; attempt <= profiles.length; attempt++) {
      // A failed map can already have changed its local body station counts.
      // Commit only a complete map; previously sampled obsolete guide points
      // must not survive installation of a geometric incoming connector.
      const candidateInput = structuredClone(input);
      try {
        matched = matchPanelPotentialGuides({ input: candidateInput, profiles, surfaceFractions: surfaceDemands ?? system.fractions, outerPaths, sample, sampleCut, sampleCutX,
          growth: crosslineGrowth, maxIntervals: maxCrosslineIntervals, outerSpread: outerCrosslineSpread, resolveTurns: resolveSurfaceTurning,
          surfaceMetric: surfaceCrosslineMetric, cutStationSpacing, surfaceStationPlacement, passageIntervalCounts, passageCountPlanning });
        Object.assign(input, candidateInput); break;
      } catch (error) {
        if (error.code !== 'STAGNATION_GUIDE_REFERENCE_MISMATCH' || attempt === profiles.length) throw error;
        const b = error.body, profile = profiles[b], path = profile.upstream;
        if (profile.geometricStagnationConnector) throw error;
        const panel = field.panels.filter(p => p.element === b)[profile.stagnationSegment];
        const joinX = Math.min(panel.a.x, panel.b.x), anchor = profile.curve.evaluate(profile.stag);
        const join = sampleUniqueGuideX(path, joinX, { sample, velocity, tolerance });
        const connector = createStagnationGuideConnector({ join, anchor: anchor.point,
          anchorDirection: { x: -anchor.derivative.y, y: anchor.derivative.x },
          stagnationPotential: profile.phiStag, field, tolerance,
          streamfunctionLevel: profile.level, admissible: fluid, admissibleSegment: segment });
        const prefix = path.filter(p => p.potential < join.potential);
        const connected = [...prefix, ...connector.points];
        if (connected.some((p, k) => k && (p.x <= connected[k - 1].x || p.potential <= connected[k - 1].potential)))
          throw new Error('The geometric stagnation connector does not join an ordered incoming guide.');
        profile.upstream = connected; geometricCuts.set(connected, connector);
        profile.geometricStagnationConnector = connector;
        diagnostics.geometricStagnationConnectors ??= [];
        diagnostics.geometricStagnationConnectors.push({ body: b, trigger: error.code, rejectedRequestX: error.requestedX,
          originalFailure: error.cause?.message, originalSampledMaximumX: Math.max(...path.map(p => p.x)),
          sourcePanel: profile.stagnationSegment, sourcePanelEndpoints: [panel.a, panel.b],
          ...connector.diagnostics });
      }
    }
    ({ guides, outer } = matched); diagnostics.potentialCrosslines = matched.diagnostics;
    system = createStreamtubeBodySystem(input); ({ nx, elements } = system.layout);
    if (input.gridSpacing) {
      input.gridSpacing.preMatching = { inlet: input.gridSpacing.inlet, outlet: input.gridSpacing.outlet };
      input.gridSpacing.surfaceIntervalsByElement = input.bodies.map(b => ({ element: b.element, intervals: b.trailingIndex - b.leadingIndex }));
      const distribution = cutStationSpacing === 'physical-x' ? 'reconciled physical outline rank' : 'matched potential metric';
      input.gridSpacing.inlet = { intervals: Math.min(...input.bodies.map(b => b.leadingIndex)), distribution };
      input.gridSpacing.outlet = { intervals: nx - Math.max(...input.bodies.map(b => b.trailingIndex)), distribution };
      input.gridSpacing.coordinate = cutStationSpacing === 'physical-x' ? 'reconciled physical outline rank' : 'matched potential';
      diagnostics.gridSpacing = structuredClone(input.gridSpacing);
    }
    if (diagnostics.surfaceDistributions) diagnostics.surfaceDistributions = diagnostics.surfaceDistributions.map(d => {
      const body = input.bodies[d.body], curve = profiles[d.body].curve, arc = createContourArc(curve);
      const fractions = system.fractions[d.body][d.side], distances = fractions.map(f => arc.at(curve.branch(d.side, f, profiles[d.body].stag).parameter));
      const averageSpacing = .5 * arc.length / fractions.length;
      return { ...d, preMatching: d, count: fractions.length, averageSpacing,
        actualLeadingSpacingRatio: Math.abs(distances[1] - distances[0]) / averageSpacing,
        actualTrailingSpacingRatio: Math.abs(distances.at(-1) - distances.at(-2)) / averageSpacing,
        method: 'curvature demand with matched potential cross-lines' };
    });
  } else {
    guides = profiles.map((profile, b) => {
      const range = input.bodies[b], upper = [], lower = [], wake = profile.wake;
      for (let i = 0; i <= nx; i++) for (const [side, row] of [['upper', upper], ['lower', lower]]) {
        if (i < range.leadingIndex) row.push(atX(profile.upstream, input.cutPaths[b][i].x));
        else if (i > range.trailingIndex) {
          const p = atX(wake, input.cutPaths[b][i].x);
          row.push({ ...p, potential: p.potential + profile.phase(side === 'upper' ? 0 : profile.curve.length) });
        } else {
          const v = profile.curve.branch(side, system.fractions[b][side][i - range.leadingIndex], profile.stag);
          row.push({ ...v.point, potential: profile.phase(v.parameter) });
        }
      }
      return { upper, lower };
    });
    outer = outerPaths.map((path, side) => (side === 0 ? input.outerLower : input.outerUpper).map(p => atX(path, p.x)));
  }
  guides.forEach((rows, b) => {
    for (const row of Object.values(rows)) {
      if (row.some((p, i) => i && p.potential <= row[i - 1].potential)) throw new Error('Body guide potential must increase.');
      for (let i = input.bodies[b].leadingIndex; i <= input.bodies[b].trailingIndex; i++)
        diagnostics.maxSurfaceStreamfunctionDefect = Math.max(diagnostics.maxSurfaceStreamfunctionDefect, Math.abs(psi(row[i]) - profiles[b].level));
    }
  });
  if (normalSpacing === 'stagnation' || normalSpacing === 'automatic') {
    diagnostics.normalAllocation = matchStagnationNormalSpacing({ input, profiles, guides, outer, aspectRatio: stagnationAspectRatio,
      expandOnly: normalSpacing === 'automatic',
      maximumAdjacentRatio: resolveSurfaceTurning ? 3 : Infinity,
      pointAt: (g, i, eta, potential) => {
        const level = (1 - eta) * input.captureLevels[g] + eta * input.captureLevels[g + 1];
        return sample(fullCurve(level), potential);
      } });
    system = createStreamtubeBodySystem(input);
  }
  if (interiorInitialization !== 'traced') diagnostics.interiorInitialization = interiorInitialization;
  const potentialCoordinates = recordPotentialCoordinates ? [] : null;
  const arcCoordinates = recordPotentialCoordinates && interiorStationPlacement === 'boundary-arc' ? [] : null;
  if (interiorStationPlacement === 'boundary-arc') {
    diagnostics.interiorStationPlacement = interiorStationPlacement;
    diagnostics.boundaryArcTransfer = [];
  }
  let nodes = input.weights.map((weights, g) => {
    const lower = g === 0 ? outer[0] : guides[g - 1].upper, upper = g === elements ? outer[1] : guides[g].lower;
    const arcMap = interiorStationPlacement === 'boundary-arc' ? createBoundaryArcStationMap({ lower, upper }) : null;
    const lineArcs = arcMap ? [] : null;
    const n = weights.length, sum = weights.reduce((a, b) => a + b, 0), eta = [0];
    for (const w of weights) eta.push(eta.at(-1) + w / sum); eta[n] = 1;
    const lines = eta.map((fraction, j) => {
      let line;
      const level = (1 - fraction) * input.captureLevels[g] + fraction * input.captureLevels[g + 1];
      if (j === 0) line = lower;
      else if (j === n) line = upper;
      else if (interiorInitialization === 'linear') {
        line = lower.map((p, i) => ({ x: (1 - fraction) * p.x + fraction * upper[i].x,
          y: (1 - fraction) * p.y + fraction * upper[i].y }));
      } else {
        try {
          if (arcMap) {
            const p = inlet(level), phi = potentialDifference(reference, p, field);
            const path = traceArc(p, phi, 0, 4 * (xmax - xmin), { normal: { x: 1, y: 0 }, offset: xmax });
            line = arcMap.at(fraction).map(s => sampleArc(path, s * path.at(-1).arc));
          } else {
            const path = fullCurve(level);
            line = Array.from({ length: nx + 1 }, (_, i) => i === 0 ? path[0] : i === nx ? path.at(-1)
              : sample(path, (1 - fraction) * lower[i].potential + fraction * upper[i].potential));
          }
        } catch (error) {
          throw new Error(`Fluid region ${g + 1}, streamline ${j} of ${n}: ${error.message}`, { cause: error });
        }
      }
      if (lineArcs) lineArcs[j] = j === 0 || j === n
        ? arcMap.at(fraction).map(s => s * (j === 0 ? arcMap.lowerLength : arcMap.upperLength)) : line.map(p => p.arc);
      // Wall guides are only approximate panel streamlines; report their
      // defect separately from the actual fluid traces.
      for (let i = 0; i <= nx; i++) {
        const wall = (j === 0 && g > 0 && system.layout.active(g - 1, i))
          || (j === n && g < elements && system.layout.active(g, i));
        const centerCut = (j === 0 && g > 0 && bodyTopologies[g - 1] && i > input.bodies[g - 1].trailingIndex)
          || (j === n && g < elements && bodyTopologies[g] && i > input.bodies[g].trailingIndex);
        if (!wall && !centerCut) diagnostics.maxStreamfunctionDrift = Math.max(diagnostics.maxStreamfunctionDrift, Math.abs(psi(line[i]) - level));
      }
      return line;
    });
    onProgress?.({ phase: 'fluid-group', group: g });
    if (arcMap) {
      let maximumInteriorArcRatio = 1;
      for (let j = 1; j < n; j++) for (let i = 1; i < nx; i++) {
        const a = lineArcs[j][i] - lineArcs[j][i - 1], b = lineArcs[j][i + 1] - lineArcs[j][i];
        maximumInteriorArcRatio = Math.max(maximumInteriorArcRatio, a / b, b / a);
      }
      diagnostics.boundaryArcTransfer.push({ region: g, metric: arcMap.metric,
        boundaryMaximumAdjacentRatio: arcMap.maximumAdjacentRatio, maximumInteriorArcRatio,
        lineLengths: lineArcs.map(s => s.at(-1)) });
      if (arcCoordinates) arcCoordinates[g] = Array.from({ length: nx + 1 }, (_, i) => lineArcs.map(line => line[i]));
    }
    if (potentialCoordinates) potentialCoordinates[g] = Array.from({ length: nx + 1 }, (_, i) => lines.map((line, j) => ({
      x: line[i].potential ?? (1 - eta[j]) * lower[i].potential + eta[j] * upper[i].potential,
      y: (1 - eta[j]) * input.captureLevels[g] + eta[j] * input.captureLevels[g + 1],
    })));
    return Array.from({ length: nx + 1 }, (_, i) => lines.map(line => ({ x: line[i].x, y: line[i].y })));
  });
  if (bodyTopologies.some(Boolean)) {
    const physical = createStreamtubeDisplacement({ layout: system.layout, curves: system.curves,
      fractions: system.fractions, thicknesses: system.displacement }).apply(nodes, system.initialStagnation).nodes;
    // Traced interiors already include the finite-base panel source field.
    // Moving them again by the bank opening introduces a discontinuity at
    // the TE and can reverse fine first-wake cells. Only linear interiors
    // need the bank offsets propagated through their interpolation weights.
    // Explicitly prescribed BL displacement still moves either interior.
    // Both constructions retain the physical wall and prescribed bank gap.
    const retainTracedInterior = interiorInitialization === 'traced' && system.inviscidBaseWake;
    diagnostics.finiteBaseInteriorTreatment = retainTracedInterior
      ? 'Retain panel-traced interior streamlines; open centerline banks only'
      : 'Interpolate prescribed displacement offsets through interior weights';
    nodes = nodes.map((group, g) => {
      const weights = input.weights[g], total = weights.reduce((a, b) => a + b, 0), eta = [0];
      weights.forEach(w => eta.push(eta.at(-1) + w / total)); eta[weights.length] = 1;
      return group.map((row, i) => {
        const a = physical[g][i][0], z = physical[g][i].at(-1);
        const dl = { x: a.x - row[0].x, y: a.y - row[0].y };
        const du = { x: z.x - row.at(-1).x, y: z.y - row.at(-1).y };
        return row.map((p, j) => j === 0 ? a : j === row.length - 1 ? z
          : retainTracedInterior ? p : { x: p.x + (1 - eta[j]) * dl.x + eta[j] * du.x,
            y: p.y + (1 - eta[j]) * dl.y + eta[j] * du.y });
      });
    });
  }
  // The normal-width construction can reverse the first wake edge when a
  // finite solid base has a tangential TE separation. Only that typed
  // first-wake failure permits a declared material-TE pairing; healthy
  // original grids and all other failed topologies retain their old path.
  const wakeRecovery = recoverFiniteBaseWakeTangency({ input, system, nodes }, { admissibleNode: fluid });
  if (wakeRecovery) {
    input.bodies = wakeRecovery.input.bodies; system = wakeRecovery.system; nodes = wakeRecovery.nodes;
    diagnostics.finiteBaseWakeTangency = wakeRecovery.report;
  }
  if (interiorInitialization === 'traced') {
    const wakeInterior = recoverFiniteBaseWakeInterior({ input, system, nodes }, { admissibleNode: fluid });
    if (wakeInterior) {
      nodes = wakeInterior.nodes;
      diagnostics.finiteBaseWakeInterior = wakeInterior.report;
    }
  }
  const initial = system.adoptGeometry(system.initial, nodes);
  if (diagnostics.normalAllocation) diagnostics.normalAllocation.groups.forEach((group, g) => {
    for (const [end, target] of Object.entries(group.targets)) {
      const j = end === 'first' ? 0 : input.weights[g].length, adjacent = end === 'first' ? 1 : j - 1;
      const i = input.bodies[target.body].leadingIndex, p = nodes[g][i][j], q = nodes[g][i][adjacent];
      target.actualDistance = Math.hypot(p.x - q.x, p.y - q.y); target.actualAspect = target.actualDistance / target.along;
    }
  });
  initial[system.layout.globals.circulation] = -panel.elements.reduce((sum, e) => sum + e.circulation, 0) / system.conditions.lengthScale;
  // Geometry repair must obey the same body-intersection guards as tracing,
  // and report its changed relation to the panel guide field explicitly.
  const guideField = {
    // Analytic evaluation of the same simultaneous panel reference, for
    // checking arbitrary refined grids independently of their discretization.
    streamfunctionAt: psi,
    velocityAt: velocity,
    admissibleNode: p => fluid(p),
    admissibleMove: (a, b) => fluid(b) && segment(a, b),
    diagnosticsForNodes: candidate => {
      const result = structuredClone(diagnostics);
      result.tracedMaxStreamfunctionDrift = diagnostics.maxStreamfunctionDrift;
      result.maxStreamfunctionDrift = 0;
      result.interiorStreamfunctionError = [];
      input.weights.forEach((weights, g) => {
        const sum = weights.reduce((a, b) => a + b, 0), eta = [0];
        weights.forEach(w => eta.push(eta.at(-1) + w / sum)); eta[weights.length] = 1;
        const span = input.captureLevels[g + 1] - input.captureLevels[g];
        const error = { region: g, maximumAbsolute: 0, maximumLocalIntervals: 0, location: null };
        candidate[g].forEach((row, i) => row.forEach((p, j) => {
          if (system.layout.nodes[g][i][j].kind === 'wall') return;
          const level = (1 - eta[j]) * input.captureLevels[g] + eta[j] * input.captureLevels[g + 1];
          const drift = Math.abs(psi(p) - level);
          result.maxStreamfunctionDrift = Math.max(result.maxStreamfunctionDrift, drift);
          if (i > 0 && i < system.layout.nx && j > 0 && j < weights.length) {
            const relative = drift / (span * Math.min(eta[j] - eta[j - 1], eta[j + 1] - eta[j]));
            error.maximumAbsolute = Math.max(error.maximumAbsolute, drift);
            if (relative > error.maximumLocalIntervals) { error.maximumLocalIntervals = relative; error.location = { i, j }; }
          }
        }));
        result.interiorStreamfunctionError.push(error);
      });
      result.normalAllocation?.groups.forEach((group, g) => {
        for (const [end, target] of Object.entries(group.targets)) {
          const j = end === 'first' ? 0 : input.weights[g].length, adjacent = end === 'first' ? 1 : j - 1;
          const i = input.bodies[target.body].leadingIndex, p = candidate[g][i][j], q = candidate[g][i][adjacent];
          target.actualDistance = Math.hypot(p.x - q.x, p.y - q.y); target.actualAspect = target.actualDistance / target.along;
        }
      });
      return result;
    },
  };
  const decodedNodes = system.decode(initial).nodes;
  // Research reference refinement samples the prescribed physical boundary,
  // then evaluates its potential. Uniform potential subdivision is singular
  // in physical distance at stagnation and is not surface-curve refinement.
  const boundaryAt = potentialCoordinates ? (g, i, j) => {
    const nt = decodedNodes[g]?.[0]?.length - 1;
    if (!Number.isInteger(g) || g < 0 || g > elements || !Number.isFinite(i) || !Number.isFinite(j)
      || i < 0 || i > nx || j < 0 || j > nt || (i !== 0 && i !== nx && j !== 0 && j !== nt))
      throw new Error('Potential reference sampling requires a region boundary index.');
    if (Number.isInteger(i) && Number.isInteger(j))
      return { point: { ...decodedNodes[g][i][j] }, coordinate: { ...potentialCoordinates[g][i][j] } };
    if (i === 0 || i === nx) {
      const k = Math.min(nt - 1, Math.floor(j)), t = j - k;
      const level = potentialCoordinates[g][i][k].y + t * (potentialCoordinates[g][i][k + 1].y - potentialCoordinates[g][i][k].y);
      const a = decodedNodes[g][i][k], b = decodedNodes[g][i][k + 1];
      const p = onPsi({ x: a.x + t * (b.x - a.x), y: a.y + t * (b.y - a.y) }, level, { x: 0, y: 1 });
      return { point: { x: p.x, y: p.y }, coordinate: { x: potentialCoordinates[g][i][k].x + potentialDifference(a, p, field), y: level } };
    }
    const k = Math.min(nx - 1, Math.floor(i)), t = i - k, b = j === 0 ? g - 1 : g, side = j === 0 ? 'upper' : 'lower';
    const x = decodedNodes[g][k][j].x + t * (decodedNodes[g][k + 1][j].x - decodedNodes[g][k][j].x);
    const sampleX = path => sampleCutX(path, x);
    let point, potential;
    if (b < 0 || b >= elements) {
      point = sampleX(outerPaths[b < 0 ? 0 : 1]); potential = point.potential;
    } else {
      const body = input.bodies[b], profile = profiles[b];
      if (i < body.leadingIndex) { point = sampleX(profile.upstream); potential = point.potential; }
      else if (i > body.trailingIndex) {
        point = sampleX(profile.wake); potential = point.potential + profile.phase(side === 'upper' ? 0 : profile.curve.length);
      } else {
        const fractions = system.fractions[b][side], first = fractions[k - body.leadingIndex];
        const fraction = first + t * (fractions[k - body.leadingIndex + 1] - first);
        const value = profile.curve.branch(side, fraction, profile.stag);
        point = value.point; potential = profile.phase(value.parameter);
      }
    }
    return { point: { x: point.x, y: point.y }, coordinate: { x: potential, y: potentialCoordinates[g][k][j].y } };
  } : null;
  return { input, system, initial, nodes: decodedNodes, diagnostics, guideField,
    ...(potentialCoordinates ? { potentialCoordinates, potentialBoundaryAt: boundaryAt } : {}),
    ...(arcCoordinates ? { arcCoordinates } : {}),
    status: 'research grid; supplied topology; flow not initialized' };
}

export function initializePanelStreamtubeBody(supplied, controls) {
  const r = createPanelStreamtubeGrid(supplied, controls), { system, diagnostics } = r;
  const initial = system.conditions.flowModel === 'incompressible' ? r.initial : initializeStreamtubeDensities(system, r.initial);
  const evaluated = system.evaluate(initial); // Reject folds, nonpositive Pi and sonic flow before returning.
  diagnostics[system.conditions.flowModel === 'incompressible' ? 'initialIncompressible' : 'initialEuler'] = evaluated.diagnostics;
  return { ...r, initial, nodes: evaluated.nodes, diagnostics,
    status: 'research initializer; supplied topology; physical acceptance incomplete' };
}
