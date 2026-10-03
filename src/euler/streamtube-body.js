// SPDX-License-Identifier: GPL-2.0-or-later
// Intrinsic Euler body verification system. All element cuts, densities,
// grid positions, stagnation parameters, captured masses and farfield
// strengths are solved together. Whole-residual finite differences remain
// an independent oracle for the analytic backend. Shared momentum upwinding
// is an explicit research option; existing callers retain the old equations.
import { createStreamtubeBodyLayout, allocateStreamtubeMasses } from './streamtube-body-layout.js';
import { createContourCurve } from '../geometry/contour-curve.js';
import { createContourTopology, createSurfaceContourCurve } from '../geometry/contour-topology.js';
import { streamtubeBaseGeometry, initialStreamtubeDisplacement } from './streamtube-geometry.js';
import { prepareContour, validateAssembly } from '../geometry/airfoil.js';
import { evaluateStreamtubeCell } from './streamtube-cell.js';
import { evaluateIncompressibleStreamtubeCell } from './incompressible-streamtube-cell.js';
import { multipoleBasis } from '../potential/farfield.js';
import { solveNewton } from '../numerics/newton.js';
import { createStreamtubeBodyJacobian } from './streamtube-body-jacobian.js';
import { solveSparseDirect } from '../numerics/klu.js';
import { solveLinear } from '../numerics/linear.js';
import { takeDoglegStep } from '../numerics/dogleg.js';
import { createStreamtubeDisplacement } from './streamtube-displacement.js';
import { streamtubeCornerConstraints } from './streamtube-corner-constraints.js';
import { proposeDensityNewton } from './streamtube-density-newton.js';
import { streamtubeMotionDirections } from './streamtube-geometry.js';
import { streamtubeWakeGap } from './streamtube-wake-geometry.js';
import { prepareStreamtubeTransportChain } from './streamtube-transport-chain.js';
import { evaluateStreamtubeMomentumBlend } from './streamtube-momentum-blend.js';
import { normalizeStreamtubeEquationSelection, validateStreamtubeEquationRegionTopology, streamtubeEquationAt } from './streamtube-equation-selection.js';

const sub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });
const mean = (a, b) => ({ x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 });
const cross = (a, b) => a.x * b.y - a.y * b.x;
const dot = (a, b) => a.x * b.x + a.y * b.y;
const finitePoint = p => p && Number.isFinite(p.x) && Number.isFinite(p.y);

export function createStreamtubeBodySystem({ bodies, outerLower, outerUpper, cutPaths, weights,
  primaryBody = 0, captureLevels, alpha = 0, mach = .2, gamma = 1.4, pressureCorrectionFactor = .1,
  flowModel = 'compressible', streamwiseMode = flowModel === 'incompressible' ? 'bernoulli' : 'momentum', displacement,
  normalStencil = 'centered', stagnationMotion = 'interpolated', geometryDomain = 'convex', wakeGeometry = 'centerline', wakeOutlet = 'centerline',
  wakeDisplacementMotion = 'fixed', upwind, hybrid }) {
  if (!['centerline', 'banks'].includes(wakeOutlet) || wakeOutlet === 'banks' && displacement === undefined)
    throw new Error('Bank outlet tangency requires displacement boundaries.');
  if (!['centerline', 'independent-banks'].includes(wakeGeometry)) throw new Error('Unknown wake geometry formulation.');
  const independentWakeBanks = wakeGeometry === 'independent-banks';
  if (independentWakeBanks && displacement === undefined) throw new Error('Independent wake banks require displacement boundary data (zero is allowed).');
  if (!['fixed', 'te-center'].includes(wakeDisplacementMotion)
    || wakeDisplacementMotion === 'te-center' && !independentWakeBanks)
    throw new Error('TE-center displacement motion requires independent wake banks.');
  if (!['convex', 'positive-simple'].includes(geometryDomain) || flowModel === 'incompressible' && geometryDomain !== 'convex')
    throw new Error('Unsupported intrinsic body geometry domain.');
  if (!['centered', 'body-stations'].includes(normalStencil) || !['interpolated', 'walls-only'].includes(stagnationMotion))
    throw new Error('Unknown intrinsic body motion parameterization.');
  if (!['compressible', 'incompressible'].includes(flowModel)) throw new Error('Unknown intrinsic body flow model.');
  const incompressible = flowModel === 'incompressible';
  if (!(incompressible ? ['bernoulli'] : ['momentum', 'isentropic', 'hybrid']).includes(streamwiseMode)) throw new Error('Unknown intrinsic body streamwise mode.');
  if (streamwiseMode === 'hybrid') {
    if (!upwind || !hybrid || typeof hybrid !== 'object' || !Number.isFinite(hybrid.epsilonP) || hybrid.epsilonP <= 0)
      throw new Error('Hybrid body rows require explicit upwinding and positive epsilonP.');
    hybrid = normalizeStreamtubeEquationSelection(hybrid);
  } else if (hybrid !== undefined) throw new Error('Hybrid controls require hybrid streamwise rows.');
  if (upwind !== undefined) {
    if (incompressible || !['momentum', 'hybrid'].includes(streamwiseMode)) throw new Error('Body upwinding requires compressible conservative momentum or hybrid rows.');
    if (!upwind || typeof upwind !== 'object' || Array.isArray(upwind)
      || upwind.boundary?.kind !== 'unfiltered-first-two' || Object.keys(upwind.boundary).some(key => key !== 'kind'))
      throw new Error('Body upwinding requires an explicit unfiltered-first-two inlet closure.');
    const { mucon = 1, mcrit = .99 } = upwind;
    if (![mucon, mcrit].every(Number.isFinite) || mcrit < 0 || mcrit > 1) throw new Error('Invalid body upwind controls.');
    upwind = Object.freeze({ mucon, mcrit, boundary: Object.freeze({ kind: upwind.boundary.kind }) });
  }
  if (!Array.isArray(outerLower) || !Array.isArray(outerUpper) || outerLower.length !== outerUpper.length
    || ![...outerLower, ...outerUpper].every(finitePoint) || !Array.isArray(weights)) throw new Error('Invalid body grid boundary paths.');
  if (![alpha, mach, gamma, pressureCorrectionFactor].every(Number.isFinite)
    || (incompressible ? mach !== 0 : mach <= 0 || mach >= 1) || gamma <= 1 || pressureCorrectionFactor < 0)
    throw new Error('The intrinsic body verification system requires finite subcritical gas conditions.');
  const layout = createStreamtubeBodyLayout({ segments: outerLower.length - 1, tubes: weights.map(w => w.length), bodies, primaryBody, densityUnknowns: !incompressible, independentWakeBanks });
  if (wakeDisplacementMotion === 'te-center') layout.wakeDisplacementMotion = wakeDisplacementMotion;
  const { nx, elements, tubes, globals } = layout;
  if (hybrid?.entropyRegions) validateStreamtubeEquationRegionTopology(hybrid.entropyRegions, layout);
  bodies = bodies.map(b => ({ ...b, surfaceFractions: b.surfaceFractions === undefined ? undefined : structuredClone(b.surfaceFractions) }));
  if (!Array.isArray(cutPaths) || cutPaths.length !== elements
    || !cutPaths.every(path => Array.isArray(path) && path.length === nx + 1 && path.every(finitePoint)))
    throw new Error('Supply one complete initial dividing/wake path for every body.');
  const topologies = bodies.map(b => b.trailingEdge?.kind === 'finite-base' ? createContourTopology(b.points, { trailingEdge: b.trailingEdge }) : null);
  const contours = bodies.map((b, k) => topologies[k]?.points ?? prepareContour(b.points)); validateAssembly(contours);
  const curves = contours.map((p, k) => topologies[k] ? createSurfaceContourCurve(p, { trailingEdge: bodies[k].trailingEdge }) : createContourCurve(p));
  const trailingCenters = contours.map((p, k) => topologies[k]?.trailingEdge.center ?? p[0]);
  const baseGeometry = streamtubeBaseGeometry(bodies, curves), hasFiniteBase = baseGeometry.some(Boolean);
  baseGeometry.forEach((base, b) => { if (base) base.points = topologies[b].base.points; });
  const inviscidBaseWake = hasFiniteBase && displacement === undefined;
  if (inviscidBaseWake) displacement = initialStreamtubeDisplacement(layout, baseGeometry);
  const initialStagnation = curves.map((curve, b) => {
    if (bodies[b].stagnationParameter !== undefined) return bodies[b].stagnationParameter;
    const points = topologies[b]?.surface.points ?? contours[b], te = trailingCenters[b]; let index = 1;
    for (let j = 2; j < points.length - 1; j++)
      if (Math.hypot(points[j].x - te.x, points[j].y - te.y) > Math.hypot(points[index].x - te.x, points[index].y - te.y)) index = j;
    return curve.knots[index];
  });
  const fractions = bodies.map((b, k) => {
    const count = b.trailingIndex - b.leadingIndex;
    if (!(initialStagnation[k] > 0 && initialStagnation[k] < curves[k].length)) throw new Error('Invalid body stagnation parameterization.');
    return Object.fromEntries(['upper', 'lower'].map(side => {
      const source = Array.isArray(b.surfaceFractions) ? b.surfaceFractions : b.surfaceFractions?.[side];
      const row = source === undefined && b.surfaceFractions === undefined ? Array.from({ length: count + 1 }, (_, i) => i / count) : source;
      if (!Array.isArray(row) || row.length !== count + 1 || row[0] !== 0 || row.at(-1) !== 1 || !row.every(Number.isFinite)
        || row.some((v, i) => i && v <= row[i - 1])) throw new Error('Invalid body surface parameterization.');
      return [side, row.slice()];
    }));
  });
  let displaced = displacement === undefined ? null : createStreamtubeDisplacement({ layout, curves, fractions, thicknesses: displacement });
  layout.displacedBoundaries = Boolean(displaced);
  outerLower = outerLower.map(p => ({ ...p })); outerUpper = outerUpper.map(p => ({ ...p }));
  cutPaths = cutPaths.map(path => path.map(p => ({ ...p }))); weights = weights.map(w => w.slice());
  const angle = alpha * Math.PI / 180, freestream = { x: Math.cos(angle), y: Math.sin(angle) };
  const inlet = [outerLower[0], ...cutPaths.map(path => path[0]), outerUpper[0]];
  const levels = captureLevels?.slice() ?? [0];
  if (captureLevels === undefined) for (let i = 1; i < inlet.length; i++) levels.push(levels.at(-1) + cross(freestream, sub(inlet[i], inlet[i - 1])));
  const initialAllocation = allocateStreamtubeMasses(levels, weights, primaryBody), massScale = initialAllocation.totalMass;
  const primaryTE = trailingCenters[primaryBody], primaryLE = curves[primaryBody].evaluate(initialStagnation[primaryBody]).point;
  const lengthScale = Math.hypot(primaryTE.x - primaryLE.x, primaryTE.y - primaryLE.y);
  const center = { x: .75 * primaryLE.x + .25 * primaryTE.x, y: .75 * primaryLE.y + .25 * primaryTE.y };
  // Incompressible pressure is relative to the common total pressure, with
  // rho_inf = q_inf = 1. It has no perfect-gas h0 or stagnation density.
  const pInf = incompressible ? -.5 : 1 / (gamma * mach * mach);
  const pressureScale = incompressible ? 1 : pInf;
  const h0 = incompressible ? null : 1 / ((gamma - 1) * mach * mach) + .5;
  const rhoTotal = incompressible ? null : (1 + .5 * (gamma - 1) * mach * mach) ** (1 / (gamma - 1));
  if (!(lengthScale > 0) || !Number.isFinite(pInf) || (!incompressible && ![h0, rhoTotal].every(Number.isFinite))) throw new Error('Invalid intrinsic body normalization.');
  const eta = weights.map(row => {
    const result = [0], sum = row.reduce((a, b) => a + b, 0);
    for (const w of row) result.push(result.at(-1) + w / sum); result[result.length - 1] = 1; return result;
  });
  let wakeSeeds = null;
  const baseNodes = stagnation => {
    const cuts = curves.map((curve, body) => Array.from({ length: nx + 1 }, (_, i) => {
      if (layout.active(body, i)) return Object.fromEntries(['lower', 'upper'].map(side => [side,
        curve.branch(side, fractions[body][side][i - bodies[body].leadingIndex], stagnation[body]).point]));
      if (wakeSeeds && i > bodies[body].trailingIndex) return wakeSeeds[body][i];
      const point = { ...cutPaths[body][i] };
      if (i < bodies[body].leadingIndex) {
        const change = sub(curve.evaluate(stagnation[body]).point, curve.evaluate(initialStagnation[body]).point), weight = i / bodies[body].leadingIndex;
        point.x += weight * change.x; point.y += weight * change.y;
      }
      return { lower: point, upper: point };
    }));
    return tubes.map((count, group) => Array.from({ length: nx + 1 }, (_, i) => {
      const lower = group === 0 ? outerLower[i] : cuts[group - 1][i].upper;
      const upper = group === elements ? outerUpper[i] : cuts[group][i].lower;
      return Array.from({ length: count + 1 }, (_, j) => j === 0 ? lower : j === count ? upper
        : { x: (1 - eta[group][j]) * lower.x + eta[group][j] * upper.x, y: (1 - eta[group][j]) * lower.y + eta[group][j] * upper.y });
    }));
  };
  // Construct the independent coordinate bases once. The optional TE-center
  // displacement chart adds its absolute translation during decode; these
  // bases always represent the zero-translation coordinates.
  if (independentWakeBanks) {
    const seed = createStreamtubeDisplacement({ layout: { ...layout, independentWakeBanks: false }, curves, fractions, thicknesses: displacement })
      .apply(baseNodes(initialStagnation), initialStagnation).nodes;
    wakeSeeds = bodies.map((_, b) => Array.from({ length: nx + 1 }, (_, i) => ({ lower: seed[b][i].at(-1), upper: seed[b + 1][i][0] })));
  }
  const originalNodes = baseNodes(initialStagnation);
  const initialPhysicalNodes = displaced ? displaced.apply(originalNodes, initialStagnation).nodes : originalNodes;
  const directions = streamtubeMotionDirections(layout, initialPhysicalNodes, normalStencil), offsets = new Map();
  // The source-comparison chart keeps free-node tangential coordinates
  // independent of stagnation motion. Wall points still follow the contour;
  // their neighboring free nodes respond through the simultaneous equations.
  const freeBase = (nodes, g, i, j) => (stagnationMotion === 'walls-only' ? originalNodes : nodes)[g][i][j];
  const decode = state => {
    if (state.length !== layout.n || !state.every(Number.isFinite)) throw new Error('Invalid intrinsic body state.');
    const stagnation = initialStagnation.map((s, b) => s + (globals.stagnation[b] === null ? 0 : curves[b].length * state[globals.stagnation[b]]));
    const captured = levels.map((v, i) => i === 0 || i === levels.length - 1 || globals.capture[i - 1] === null ? v : v + massScale * state[globals.capture[i - 1]]);
    const allocation = allocateStreamtubeMasses(captured, weights, primaryBody), nodes = baseNodes(stagnation), moved = new Map();
    for (let g = 0; g <= elements; g++) for (let i = 0; i <= nx; i++) for (let j = 0; j <= tubes[g]; j++) {
      const col = layout.nodes[g][i][j].column; if (col === null) continue;
      if (!moved.has(col)) {
        const normal = directions.get(col), base = freeBase(nodes, g, i, j), distance = lengthScale * state[col], offset = offsets.get(col) ?? { x: 0, y: 0 };
        moved.set(col, { x: base.x + offset.x + distance * normal.x, y: base.y + offset.y + distance * normal.y });
      }
      nodes[g][i][j] = moved.get(col);
    }
    return { stagnation, captured, allocation, ...(displaced
      ? { nodes: displaced.apply(nodes, stagnation).nodes, undisplacedNodes: nodes } : { nodes }) };
  };
  // Differentiate the current normal-coordinate chart. Its offsets and
  // normals stay fixed within Newton; rebasing changes these directions
  // only after an accepted step. The interpolated chart also moves free-node
  // base positions with stagnation; the walls-only chart omits that motion.
  const geometryDerivatives = (state, { includeDisplacement = false } = {}) => {
    const { stagnation, undisplacedNodes } = decode(state);
    const boundary = (body, side, i) => {
      const col = globals.stagnation[body], row = new Map();
      if (col === null || i > bodies[body].trailingIndex) return row;
      const curve = curves[body], derivative = layout.active(body, i)
        ? curve.branch(side, fractions[body][side][i - bodies[body].leadingIndex], stagnation[body]).stagnationDerivative
        : curve.evaluate(stagnation[body]).derivative;
      const weight = curve.length * (i < bodies[body].leadingIndex ? i / bodies[body].leadingIndex : 1);
      row.set(col, { x: weight * derivative.x, y: weight * derivative.y }); return row;
    };
    const derivatives = tubes.map((count, g) => Array.from({ length: nx + 1 }, (_, i) => {
      const lower = g === 0 ? new Map() : boundary(g - 1, 'upper', i);
      const upper = g === elements ? new Map() : boundary(g, 'lower', i);
      return Array.from({ length: count + 1 }, (_, j) => {
        const row = new Map();
        for (const [source, weight] of [[lower, 1 - eta[g][j]], [upper, eta[g][j]]])
          for (const [col, p] of source) if (weight !== 0) {
            const old = row.get(col) ?? { x: 0, y: 0 }; row.set(col, { x: old.x + weight * p.x, y: old.y + weight * p.y });
          }
        const col = layout.nodes[g][i][j].column;
        if (col !== null) {
          if (stagnationMotion === 'walls-only') row.clear();
          const normal = directions.get(col); row.set(col, { x: lengthScale * normal.x, y: lengthScale * normal.y });
        }
        return row;
      });
    }));
    return displaced ? displaced.apply(undisplacedNodes, stagnation, derivatives, includeDisplacement).derivatives : derivatives;
  };
  const evaluate = state => {
    const decoded = decode(state), { nodes, allocation } = decoded;
    const cells = Array.from({ length: nx - 1 }, () => tubes.map(() => [])), sections = Array.from({ length: nx }, () => tubes.map(() => []));
    const hybridCells = hybrid ? Array.from({ length: nx - 1 }, () => tubes.map(() => [])) : null;
    const hybridDiagnostics = hybrid ? { ...hybrid, entropyCells: 0, momentumCells: 0, blendedCells: 0, maxMomentumDeparture: 0 } : null;
    let transportSpeeds;
    if (upwind) {
      transportSpeeds = Array.from({ length: nx }, () => tubes.map(() => []));
      for (let g = 0; g <= elements; g++) for (let j = 0; j < tubes[g]; j++) {
        let chain;
        try { chain = prepareStreamtubeTransportChain({ lower: nodes[g].map(row => row[j]), upper: nodes[g].map(row => row[j + 1]),
          densities: Array.from({ length: nx }, (_, i) => Math.exp(state[layout.densityIndex(i, g, j)])),
          massFlow: allocation.groups[g][j].massFlow, stagnationEnthalpy: h0, gamma, geometryDomain, upwind }); }
        catch (error) { throw Object.assign(new Error(`Body transport chain group=${g}, tube=${j}: ${error.message}`, { cause: error }),
          error.code ? { code: error.code, diagnostics: { ...error.diagnostics, group: g, tube: j } } : {}); }
        // The remote boundary equations retain their subsonic domain. Test
        // the actual gas here, independently of the multipole model's speed.
        for (const i of [0, nx - 1]) if (chain.sections[i].machSquared >= 1) {
          const section = chain.sections[i];
          throw Object.assign(new Error(`Body inlet/outlet section must remain subsonic: i=${i}, group=${g}, tube=${j}.`), {
            code: 'streamtube-boundary-subsonic', diagnostics: {
              i, group: g, tube: j, boundary: i === 0 ? 'inlet' : 'outlet',
              machSquared: section.machSquared, rho: section.rho, p: section.p, q: section.q,
              enthalpy: section.enthalpy, machSquaredUpperBound: 1, subsonicMargin: 1 - section.machSquared,
            },
          });
        }
        for (let i = 0; i < nx; i++) transportSpeeds[i][g][j] = chain.transportSpeeds[i];
      }
    }
    let maxMachSquared = 0, maxEntropyJump = 0, maxMomentumResidual = 0, maxStagnationPressureError = 0;
    const p0Inf = incompressible ? null : (gamma - 1) / gamma * rhoTotal * h0;
    for (let i = 1; i < nx; i++) for (let group = 0; group <= elements; group++) for (let tube = 0; tube < tubes[group]; tube++) {
      let cell;
      try { cell = (incompressible ? evaluateIncompressibleStreamtubeCell : evaluateStreamtubeCell)({ lower: [nodes[group][i - 1][tube], nodes[group][i][tube], nodes[group][i + 1][tube]],
        upper: [nodes[group][i - 1][tube + 1], nodes[group][i][tube + 1], nodes[group][i + 1][tube + 1]],
        ...(incompressible ? {} : { densities: [Math.exp(state[layout.densityIndex(i - 1, group, tube)]), Math.exp(state[layout.densityIndex(i, group, tube)])] }),
        massFlow: allocation.groups[group][tube].massFlow, stagnationEnthalpy: h0, gamma, pressureCorrectionFactor, geometryDomain,
        ...(upwind ? { transportSpeeds: [transportSpeeds[i - 1][group][tube], transportSpeeds[i][group][tube]] } : {}) }); }
      catch (error) { throw Object.assign(new Error(`Body cell i=${i}, group=${group}, tube=${tube}: ${error.message}`, { cause: error }),
        error.code ? { code: error.code, diagnostics: { ...error.diagnostics, cell: { i, group, tube } } } : {}); }
      if (!incompressible) for (const s of cell.states) {
        if (!upwind && s.machSquared >= 1) throw new Error('Local sonic flow is outside this body verification system.');
        maxMachSquared = Math.max(maxMachSquared, s.machSquared);
        const logTotalPressureRatio = Math.log(s.p / p0Inf) + gamma / (gamma - 1) * Math.log(h0 / s.enthalpy);
        maxStagnationPressureError = Math.max(maxStagnationPressureError, Math.abs(Math.expm1(logTotalPressureRatio)));
      }
      if (!incompressible) maxEntropyJump = Math.max(maxEntropyJump, Math.abs(cell.entropyJump));
      maxMomentumResidual = Math.max(maxMomentumResidual, Math.abs(cell.streamwiseResidual / pressureScale));
      cells[i - 1][group][tube] = cell;
      if (hybrid) {
        const equation = streamtubeEquationAt({ hybrid, bodies, tubes, nx, i, group, tube });
        const blend = equation === 'hybrid' ? evaluateStreamtubeMomentumBlend({ ...cell, epsilonP: hybrid.epsilonP })
          : { fraction: equation === 'momentum' ? 1 : 0,
            residual: equation === 'momentum' ? cell.streamwiseResidual : cell.isentropicResidual };
        hybridCells[i - 1][group][tube] = blend;
        hybridDiagnostics[blend.fraction === 0 ? 'entropyCells' : blend.fraction === 1 ? 'momentumCells' : 'blendedCells']++;
        hybridDiagnostics.maxMomentumDeparture = Math.max(hybridDiagnostics.maxMomentumDeparture,
          Math.abs((1 - blend.fraction) * (cell.streamwiseResidual - cell.isentropicResidual) / pressureScale));
      }
      if (i === 1) sections[0][group][tube] = cell.states[0]; sections[i][group][tube] = cell.states[1];
    }
    const strengths = { circulation: lengthScale * state[globals.circulation], source: lengthScale * state[globals.source],
      doubletX: lengthScale ** 2 * state[globals.doubletX], doubletY: lengthScale ** 2 * state[globals.doubletY] };
    // multipoleBasis uses positive CCW circulation. The body's global Gamma
    // follows MSES's positive-lift convention, hence the minus sign here.
    const coefficients = [-strengths.circulation, strengths.source, strengths.doubletX, strengths.doubletY, strengths.circulation ** 2];
    const farfield = point => {
      const basis = multipoleBasis(point, { center, alpha, mach, gamma }), velocity = { ...freestream };
      for (let k = 0; k < coefficients.length; k++) { velocity.x += coefficients[k] * basis.velocity[k][0]; velocity.y += coefficients[k] * basis.velocity[k][1]; }
      const q2 = dot(velocity, velocity), temperature = 1 + .5 * (gamma - 1) * mach * mach * (1 - q2);
      if (!(temperature > 0) || mach * mach * q2 / temperature >= 1) throw new Error('Non-subcritical body farfield state.');
      return { velocity, p: incompressible ? -.5 * q2 : pInf * temperature ** (gamma / (gamma - 1)), basis };
    };
    const bodyPressure = (body, i, side) => side === 'lower'
      ? cells[i - 1][body].at(-1).interfacePressure.upper : cells[i - 1][body + 1][0].interfacePressure.lower;
    // Drela (1986), Eq. 3.20: each physical outlet streamline follows the
    // local multipole velocity. In a coupled wake the mean belongs to Euler;
    // the bank difference replaces the extra terminal BL velocity closure.
    const outletBankTangency = bodies.map((_, b) => ['lower', 'upper'].map(side => {
      const g = side === 'lower' ? b : b + 1, j = side === 'lower' ? tubes[g] : 0;
      const a = nodes[g][nx - 1][j], z = nodes[g][nx][j], edge = sub(z, a);
      return cross(edge, mean(farfield(a).velocity, farfield(z).velocity)) / Math.hypot(edge.x, edge.y);
    }));
    const match = [0, 0, 0], matchScale = [0, 0, 0];
    for (const [group, j, tube] of [[0, 0, 0], [elements, tubes[elements], tubes[elements] - 1]]) for (let i = 0; i < nx; i++) {
      const first = nodes[group][i][j], second = nodes[group][i + 1][j], edge = sub(second, first), length = Math.hypot(edge.x, edge.y);
      const q = sections[i][group][tube].q, velocity = { x: q * edge.x / length, y: q * edge.y / length }, ff = farfield(mean(first, second));
      const defect = cross(ff.velocity, velocity);
      for (let k = 0; k < 3; k++) {
        const dv = ff.basis.velocity[k + 1], mode = dv[0] * velocity.y - dv[1] * velocity.x;
        match[k] += length * defect * mode; matchScale[k] += length * mode * mode;
      }
    }
    if (matchScale.some(s => !(s > 0))) throw new Error('Degenerate body farfield matching mode.');
    const residual = Float64Array.from(layout.rows, row => {
      if (row.kind === 'streamwise') {
        const cell = cells[row.i - 1][row.group][row.tube];
        return (hybrid ? hybridCells[row.i - 1][row.group][row.tube].residual
          : streamwiseMode === 'isentropic' ? cell.isentropicResidual : cell.streamwiseResidual) / pressureScale;
      }
      if (row.kind === 'inletDensity') {
        const s = sections[0][row.group][row.tube]; return Math.log(s.rho / rhoTotal) - Math.log(s.enthalpy / h0) / (gamma - 1);
      }
      if (row.kind === 'internalPressure') return (cells[row.i - 1][row.group][row.j - 1].interfacePressure.upper - cells[row.i - 1][row.group][row.j].interfacePressure.lower) / pressureScale;
      if (row.kind === 'farfieldPressure') {
        const group = row.side === 'lower' ? 0 : elements, j = row.side === 'lower' ? 0 : tubes[group];
        const cell = cells[row.i - 1][group][row.side === 'lower' ? 0 : tubes[group] - 1];
        return (cell.interfacePressure[row.side] - farfield(nodes[group][row.i][j]).p) / pressureScale;
      }
      if (['cutPressure', 'trailingKutta', 'leadingKutta'].includes(row.kind)) return (bodyPressure(row.body, row.i, 'upper') - bodyPressure(row.body, row.i, 'lower')) / pressureScale;
      if (row.kind === 'wakeGap') {
        const indices = [row.i - 1, row.i, Math.min(nx, row.i + 1)];
        const lower = indices.map(i => nodes[row.body][i].at(-1)), upper = indices.map(i => nodes[row.body + 1][i][0]);
        const { gap } = streamtubeWakeGap(lower, upper);
        // Roundoff in coordinate subtraction must not imply penetration of
        // coincident banks. Keep the signed residual; never clip the gap.
        const roundoff = 64 * Number.EPSILON * Math.max(lengthScale, ...[...lower, ...upper].flatMap(p => [Math.abs(p.x), Math.abs(p.y)]));
        if (gap < -roundoff) throw new Error(`Wake banks overlap at body=${row.body}, station=${row.i}.`);
        return (gap - displaced.values.wakes[row.body][row.i - bodies[row.body].trailingIndex - 1]) / lengthScale;
      }
      if (row.kind === 'endTangency') {
        if (wakeOutlet === 'banks' && row.node.kind === 'cut' && row.i === nx)
          return .5 * (outletBankTangency[row.node.body][0] + outletBankTangency[row.node.body][1]);
        const group = row.node.kind === 'cut' ? row.node.body : row.node.group, j = row.node.kind === 'cut' ? tubes[group] : row.node.j;
        const point = i => displaced && row.node.kind === 'cut' && row.i === nx
          ? mean(nodes[group][i].at(-1), nodes[group + 1][i][0]) : nodes[group][i][j];
        const i = row.i === 0 ? 0 : nx - 1, first = point(i), second = point(i + 1), edge = sub(second, first);
        const velocity = mean(farfield(first).velocity, farfield(second).velocity);
        return cross(edge, velocity) / Math.hypot(edge.x, edge.y);
      }
      if (row.kind === 'farfieldMatch') { const k = ['source', 'doubletX', 'doubletY'].indexOf(row.mode); return match[k] / Math.sqrt(matchScale[k]); }
      throw new Error(`Unknown intrinsic body row ${row.kind}.`);
    });
    if (!residual.every(Number.isFinite)) throw new Error('Nonfinite intrinsic body residual.');
    const residualByFamily = {};
    layout.rows.forEach((row, i) => { residualByFamily[row.kind] = Math.max(residualByFamily[row.kind] ?? 0, Math.abs(residual[i])); });
    return { ...decoded, cells, sections, residual, strengths, bodyPressure, outletBankTangency,
      ...(upwind ? { transportSpeeds } : {}),
      ...(hybrid ? { hybridCells } : {}),
      diagnostics: { residual: residual.reduce((peak, value) => Math.max(peak, Math.abs(value)), -Infinity), residualByFamily, maxMach: Math.sqrt(maxMachSquared),
        flowModel, streamwiseMode, maxEntropyJump: incompressible ? null : maxEntropyJump, maxMomentumResidual,
        maxStagnationPressureError: incompressible ? null : maxStagnationPressureError,
        ...(upwind ? { upwind } : {}), ...(hybrid ? { hybrid: hybridDiagnostics } : {}) } };
  };
  const rejectedTrials = new Set();
  const admissible = state => { try { evaluate(state); return true; } catch (error) { rejectedTrials.add(error.message); return false; } };
  const adoptGeometry = (state, nodes) => {
    const { stagnation } = decode(state), base = baseNodes(stagnation), visited = new Set(), next = state.slice();
    const physicalNodes = nodes;
    const nextOffsets = new Map();
    if (!Array.isArray(nodes) || nodes.length !== tubes.length || !nodes.every((group, g) => Array.isArray(group) && group.length === nx + 1
      && group.every(row => Array.isArray(row) && row.length === tubes[g] + 1 && row.every(finitePoint)))) throw new Error('Incompatible body restart grid dimensions.');
    if (displaced) nodes = displaced.restore(nodes, stagnation, base, 1e-10 * lengthScale);
    for (let g = 0; g <= elements; g++) for (let i = 0; i <= nx; i++) for (let j = 0; j <= tubes[g]; j++) {
      const node = layout.nodes[g][i][j];
      if (node.column === null && Math.hypot(nodes[g][i][j].x - base[g][i][j].x, nodes[g][i][j].y - base[g][i][j].y) > 1e-10 * lengthScale)
        throw new Error('Body restart wall does not match the current geometry and stagnation state.');
      if (node.kind === 'cut' && !node.side && g === node.body + 1 && Math.hypot(nodes[g][i][j].x - nodes[g - 1][i].at(-1).x, nodes[g][i][j].y - nodes[g - 1][i].at(-1).y) > 1e-10 * lengthScale)
        throw new Error('Body restart has a disconnected dividing/wake cut.');
    }
    for (let g = 0; g <= elements; g++) for (let i = 0; i <= nx; i++) for (let j = 0; j <= tubes[g]; j++) {
      const col = layout.nodes[g][i][j].column; if (col === null || visited.has(col)) continue; visited.add(col);
      nextOffsets.set(col, sub(nodes[g][i][j], freeBase(base, g, i, j))); next[col] = 0;
    }
    let nextDirections;
    // Offsets belong to the underlying chart, but NCALC motion directions
    // belong to the current physical flow grid. In particular the first
    // wake secant must start at the displaced TE center, not the bare foil.
    try { nextDirections = streamtubeMotionDirections(layout, physicalNodes, normalStencil); }
    catch (error) {
      if (error.message === 'Degenerate body grid movement direction.') throw new Error('Degenerate rebased body normal.');
      throw error;
    }
    // A rejected restart must leave the previous coordinate chart usable.
    for (const [col, offset] of nextOffsets) offsets.set(col, offset);
    for (const [col, normal] of nextDirections) directions.set(col, normal);
    return next;
  };
  // Normals stay fixed within each Newton linearization. Recenter them at
  // the accepted physical grid; this coordinate change preserves its state.
  const rebase = state => adoptGeometry(state, decode(state).nodes);
  // At a zero-increment chart origin, refreshing NCALC directions cannot
  // move a node. Keep the existing offsets instead of restoring displaced
  // coordinates and subtracting/adding the TE translation a second time.
  const refreshGeometryDirections = state => {
    if (layout.positions.some(({ column }) => state[column] !== 0))
      throw new Error('Refreshing body directions requires zero free-position increments.');
    const next = streamtubeMotionDirections(layout, decode(state).nodes, normalStencil);
    for (const [column, direction] of next) directions.set(column, direction);
  };
  const geometryChart = () => layout.positions.map(({ column }) => ({ column, offset: { ...(offsets.get(column) ?? { x: 0, y: 0 }) }, normal: { ...directions.get(column) } }));
  const system = { layout, curves, fractions, originalNodes, initialStagnation, initial: new Float64Array(layout.n), decode, evaluate,
    ...(hasFiniteBase ? { baseGeometry, inviscidBaseWake } : {}),
    get displacement() { return displaced?.values ?? null; },
    get displacementParameters() { return displaced?.parameters ?? []; },
    setDisplacement: thicknesses => {
      if (!displaced) throw new Error('Enable displacement boundaries when constructing the body system.');
      // Validate before replacing the current data; this never changes the
      // normal-coordinate chart or the Euler state during a residual call.
      displaced = createStreamtubeDisplacement({ layout, curves, fractions, thicknesses });
    },
    residual: state => evaluate(state).residual, admissible, rebase, adoptGeometry, refreshGeometryDirections, geometryChart, geometryDerivatives, rejectedTrials,
    conditions: { mach, alpha, gamma, pInf, pressureScale, h0, rhoTotal, lengthScale, massScale, center, pressureCorrectionFactor, streamwiseMode, flowModel, normalStencil, stagnationMotion, geometryDomain, wakeGeometry, wakeOutlet,
      ...(wakeDisplacementMotion === 'te-center' ? { wakeDisplacementMotion } : {}),
      ...(upwind ? { upwind } : {}), ...(hybrid ? { hybrid } : {}) } };
  system.jacobian = createStreamtubeBodyJacobian(system);
  system.evaluateJacobian = system.jacobian.evaluateJacobian;
  return system;
}

export function solveStreamtubeBody(system, { initial = system.initial, maxIterations = 40, tolerance = 1e-10, onIteration, onMesh,
  jacobianBackend = 'analytic', linearBackend = jacobianBackend === 'analytic' ? 'klu' : 'dense', linearTolerance = 1e-10,
  stepMethod = 'newton', initialTrustRadius = 1, projectedSteps = stepMethod === 'dogleg',
  secondOrderSteps = projectedSteps, ...controls } = {}) {
  if (!['analytic', 'finite-difference'].includes(jacobianBackend) || !['dense', 'klu'].includes(linearBackend)
    || (linearBackend === 'klu' && jacobianBackend !== 'analytic') || !(linearTolerance > 0) || !Number.isFinite(linearTolerance))
    throw new Error('Invalid intrinsic body solver backends or linear tolerance.');
  if (!['newton', 'dogleg', 'density-newton'].includes(stepMethod) || (stepMethod === 'dogleg' && jacobianBackend !== 'analytic')
    || !Number.isFinite(initialTrustRadius) || initialTrustRadius <= 0 || initialTrustRadius > 1e6)
    throw new Error('Invalid intrinsic body step method or trust radius.');
  if (stepMethod === 'density-newton' && (jacobianBackend !== 'analytic' || !system.layout.densityCount))
    throw new Error('Density Newton requires compressible analytic body equations.');
  if (typeof projectedSteps !== 'boolean' || (projectedSteps && stepMethod !== 'dogleg'))
    throw new Error('Projected Euler steps require the dogleg controller.');
  if (typeof secondOrderSteps !== 'boolean' || (secondOrderSteps && !projectedSteps))
    throw new Error('Second-order Euler steps require projected steps.');
  if (!Number.isInteger(maxIterations) || maxIterations < 0 || !(tolerance > 0) || !Number.isFinite(tolerance) || !system.admissible(initial))
    throw new Error('Invalid intrinsic body solve controls or initial state.');
  const jacobian = jacobianBackend === 'analytic' ? state => system.jacobian(state, { sparse: linearBackend === 'klu' }) : undefined;
  const linearDiagnostics = linearBackend === 'klu' ? { solves: 0, maxRelativeResidual: 0, maxFactorNonzeros: 0, refinements: 0 } : null;
  const linearSolve = linearBackend === 'klu' ? (matrix, rhs) => {
    const r = solveSparseDirect(matrix, rhs, { tolerance: linearTolerance });
    linearDiagnostics.solves++; linearDiagnostics.maxRelativeResidual = Math.max(linearDiagnostics.maxRelativeResidual, r.relativeResidual);
    linearDiagnostics.maxFactorNonzeros = Math.max(linearDiagnostics.maxFactorNonzeros, r.factorNonzeros);
    linearDiagnostics.refinements += r.refinements; return r.x;
  } : undefined;
  let state = Float64Array.from(initial), norm = system.evaluate(state).diagnostics.residual, reason = 'iteration limit', trustRadius = initialTrustRadius;
  let lastRejectedStep = null;
  // The Euler-only path needs the same geometric step constraints as the
  // coupled solve. Shrinking an outward Newton/gradient ray cannot pass a
  // corner limit. Projection changes the proposal, not the Euler equations;
  // nonlinear gas and signed-cell admissibility still decide acceptance.
  const stepConstraints = x => streamtubeCornerConstraints(system, x).map(({ value, gradient }) =>
    ({ value, gradient, lower: -.9 * value }));
  const constraintValues = x => streamtubeCornerConstraints(system, x, { derivatives: false }).map(c => c.value);
  const history = [{ iteration: 0, residual: norm, step: 0 }]; onIteration?.(history[0]);
  for (let iteration = 0; iteration < maxIterations && norm > tolerance; iteration++) {
    system.rejectedTrials.clear();
    let next, entry;
    if (stepMethod === 'dogleg') {
      const residual = system.residual(state), matrix = jacobian(state);
      let newtonDirection;
      try { newtonDirection = (linearSolve ?? solveLinear)(matrix, residual.map(v => -v)); }
      catch (error) { reason = error.message; break; }
      const step = takeDoglegStep({ initial: state, currentResidual: residual, matrix, newtonDirection,
        residual: system.residual, admissible: system.admissible, radius: trustRadius,
        linearizedConstraints: projectedSteps ? stepConstraints : undefined,
        constraintValues: secondOrderSteps ? constraintValues : undefined });
      trustRadius = step.radius;
      if (!step.accepted) { reason = step.reason; break; }
      next = step.x;
      // This ratio is the scaled length relative to the Newton direction,
      // not a scalar line-search multiple of that direction. The certified
      // KLU residual belongs to the full Newton solve; dogleg deliberately
      // retains a nonzero linearized residual until the Newton step fits.
      entry = { step: step.scaledStepNorm / step.scaledNewtonNorm, stepKind: step.kind, trustRadius,
        trialRadius: step.trialRadius, reductionRatio: step.ratio, actualReduction: step.actualReduction,
        predictedReduction: step.predictedReduction, linearizedResidualNorm: step.linearizedResidualNorm,
        trials: step.trials.length, ...(step.projection ? { projection: step.projection } : {}),
        ...(step.doglegSegment ? { doglegSegment: step.doglegSegment } : {}),
        ...(step.correction ? { correction: step.correction } : {}) };
    } else if (stepMethod === 'density-newton') {
      // A source-comparison path, not the full ISES grid/update algorithm.
      // No residual line search or directional fallback: expose the first
      // invalid candidate so missing grid maintenance can be investigated.
      let proposal;
      try {
        const residual = system.residual(state), matrix = jacobian(state);
        const direction = (linearSolve ?? solveLinear)(matrix, residual.map(v => -v));
        proposal = proposeDensityNewton(system, state, direction);
        ({ x: next, ...entry } = proposal);
        system.evaluate(next);
      } catch (error) {
        if (proposal) { const { x, ...details } = proposal; lastRejectedStep = details; }
        reason = `Density Newton proposal rejected: ${error.message}`; break;
      }
    } else {
      const step = solveNewton({ ...controls, initial: state, maxIterations: 1, tolerance, residual: system.residual, admissible: system.admissible, jacobian, linearSolve });
      if (step.history.length < 2) { reason = step.reason; break; }
      next = step.x; entry = step.history.at(-1);
    }
    state = system.rebase(next);
    const accepted = system.evaluate(state); norm = accepted.diagnostics.residual;
    history.push({ ...entry, iteration: iteration + 1, residual: norm }); onIteration?.(history.at(-1));
    // Publish the accepted physical grid, including the iteration it belongs
    // to. Rejected line-search/trust-region trial grids are never displayed.
    onMesh?.({ system, nodes: accepted.nodes, flow: accepted, iteration: { ...history.at(-1) } });
  }
  const r = { x: state, history, converged: norm <= tolerance, reason: norm <= tolerance ? 'residual' : reason,
    projectedSteps, secondOrderSteps, ...(lastRejectedStep ? { lastRejectedStep } : {}) };
  const value = system.evaluate(r.x);
  const { surfaces, diagnosticForces: forces } = streamtubeSurfaceDiagnostics(system, value);
  const { bodyPressure, ...serializable } = value;
  return { ...r, ...serializable, stepMethod, geometryChart: system.geometryChart(), rejectedTrials: [...system.rejectedTrials], surfaces, diagnosticForces: forces,
    jacobianBackend, linearBackend, linearDiagnostics, streamwiseMode: system.conditions.streamwiseMode, flowModel: system.conditions.flowModel,
    displacement: system.displacement,
    ...(system.inviscidBaseWake ? { inviscidBaseWake: true, baseGeometry: system.baseGeometry,
      baseForceModel: 'Mean upper/lower TE pressure on every retained base segment; zero pressure jump is imposed by Kutta.' } : {}),
    forceStatus: system.displacement ? 'Unvalidated displacement-surface traction diagnostics; not solid-wall aerodynamic forces.' : 'Unvalidated body discretization diagnostics; not released aerodynamic predictions.',
    formulation: `${system.conditions.flowModel === 'incompressible' ? 'Incompressible moving-streamtube grid with constant density, mass and Bernoulli eliminated' : `Intrinsic moving-streamtube flow with ${system.conditions.streamwiseMode} streamwise rows`}${system.displacement ? '; prescribed wall displacement and separated wake banks' : ''}; no BL coupling${system.conditions.upwind ? '; MSES shared-speed shock dissipation with locally transonic flow and subsonic remote endpoints' : ' or shock dissipation'}; ${jacobianBackend === 'analytic' ? 'analytic chain-rule' : 'whole-residual finite-difference'} Newton with ${stepMethod === 'dogleg' ? 'scaled dogleg trust region' : stepMethod === 'density-newton' ? 'additive physical density and common scalar damping; ISES grid redistribution absent' : 'residual line search'}, ${linearBackend === 'klu' ? 'KLU/WASM sparse LU' : 'dense LU'}.` };
}

// Integrate an already evaluated flow; shared by final and live reporting.
export function streamtubeSurfaceDiagnostics(system, value) {
  const { lengthScale, pInf, center, alpha } = system.conditions;
  const surfaces = [], forces = [];
  for (let body = 0; body < system.layout.elements; body++) {
    const range = system.layout.bodies[body]; let cx = 0, cy = 0, cm = 0;
    for (const side of ['lower', 'upper']) {
      const group = side === 'lower' ? body : body + 1, j = side === 'lower' ? system.layout.tubes[group] : 0, sign = side === 'lower' ? 1 : -1;
      const points = [];
      for (let i = range.leadingIndex; i <= range.trailingIndex; i++) points.push({ ...value.nodes[group][i][j], cp: 2 * (value.bodyPressure(body, i, side) - pInf) });
      for (let i = 0; i < points.length - 1; i++) {
        const a = points[i], b = points[i + 1], middle = mean(a, b);
        // Constant interface pressure on each half-edge matches the actual
        // bent dual-volume wall traction, including LE/TE half segments.
        for (const [first, second, cp] of [[a, middle, a.cp], [middle, b, b.cp]]) {
          const edge = sub(second, first), location = sub(mean(first, second), center), fx = -sign * cp * edge.y, fy = sign * cp * edge.x;
          cx += fx / lengthScale; cy += fy / lengthScale; cm -= (location.x * fy - location.y * fx) / lengthScale ** 2;
        }
      }
      surfaces.push({ body, side, points });
    }
    if (system.inviscidBaseWake && system.baseGeometry[body]) {
      const base = system.baseGeometry[body];
      const cp = value.bodyPressure(body, range.trailingIndex, 'upper')
        + value.bodyPressure(body, range.trailingIndex, 'lower') - 2 * pInf;
      for (let k = 1; k < base.points.length; k++) {
        const a = base.points[k - 1], b = base.points[k], edge = sub(b, a);
        const location = sub(mean(a, b), center), fx = -cp * edge.y, fy = cp * edge.x;
        cx += fx / lengthScale; cy += fy / lengthScale;
        cm -= (location.x * fy - location.y * fx) / lengthScale ** 2;
      }
    }
    const angle = alpha * Math.PI / 180;
    forces.push({ cx, cy, cm, cl: cy * Math.cos(angle) - cx * Math.sin(angle), cd: cx * Math.cos(angle) + cy * Math.sin(angle) });
  }
  return { surfaces, diagnosticForces: forces };
}
