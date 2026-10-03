// SPDX-License-Identifier: GPL-2.0-or-later
// Native integral BL equations on all intrinsic body/wake stations. This
// Material forced trips stay fixed on the contour; their active intervals
// can change when the stagnation/grid coordinates carry a node past a trip.
import { createIntegralKernel } from '../viscous/integral.js';
import { createContourArc } from '../geometry/contour-arc.js';
import { solveNewton } from '../numerics/newton.js';
import { boundedResidualPartial } from '../numerics/bounded-partial.js';
import { checkAutomaticTransition, evaluateTransitionInterval } from '../viscous/transition-interval.js';
import { prepareSurfaceTransition, selectSurfaceTransition, chooseMixedTransitionShear, minimizeMixedTransitionShear } from '../viscous/transition-selection.js';
import { evaluateLeadingTransitionInterval } from '../viscous/leading-transition-interval.js';
import { initializeXfoilSurfaces } from '../viscous/xfoil-surface-initializer.js';
import { createXfoilDeadAirGap } from '../viscous/xfoil-dead-air-gap.js';
import { projectXfoilDisplacement } from './streamtube-coupled-xfoil-update.js';
import { solveLinear, normInf } from '../numerics/linear.js';

const add = (row, col, value) => row.set(col, (row.get(col) ?? 0) + value);
const zero = { x: 0, y: 0 };
const mean = (a, b) => ({ x: .5 * (a.x + b.x), y: .5 * (a.y + b.y) });

// Local MRCHUE inverse update (xbl.f:737-774). It is used only by the
// boundary-layer initializer, so keeping it here makes that startup path
// inspectable as one unit.
function solveXfoilInverseWakeStation({ initial, residual, admissible, scale, mach, gamma,
  wakeGap = 0, maxIterations = 30, tolerance = 1e-10 }) {
  if (initial.length !== 4 || !initial.every(Number.isFinite) || !(scale > 0)
    || !Number.isFinite(scale) || !Number.isInteger(maxIterations) || maxIterations < 0
    || !(tolerance > 0) || !Number.isFinite(tolerance)) throw new Error('Invalid native inverse-wake controls.');
  let x = Float64Array.from(initial), r = Float64Array.from(residual(x));
  if (r.length !== 4 || !r.every(Number.isFinite) || !admissible(x) || x.some(v => v <= 0))
    throw new Error('Invalid native inverse-wake initial state.');
  const history = [{ iteration: 0, residual: normInf(r), step: 0 }];
  const finish = reason => ({ converged: normInf(r) <= tolerance, x, history,
    reason: normInf(r) <= tolerance ? 'residual' : reason });
  for (let iteration = 1; iteration <= maxIterations && normInf(r) > tolerance; iteration++) {
    const matrix = new Float64Array(16);
    for (let col = 0; col < 4; col++) {
      const h = Math.cbrt(Number.EPSILON) * Math.max(1, Math.abs(x[col]));
      const p = x.slice(), m = x.slice(); p[col] += h; m[col] -= h;
      const plus = admissible(p), minus = admissible(m);
      if (!plus && !minus) return finish('No admissible inverse-wake Jacobian perturbation.');
      const rp = plus ? residual(p) : r, rm = minus ? residual(m) : r, width = h * (plus && minus ? 2 : 1);
      for (let row = 0; row < 4; row++) matrix[4 * row + col] = (rp[row] - rm[row]) / width;
    }
    const direction = solveLinear(matrix, r.map(v => -v));
    const maximumRelativeIncrement = Math.max(...direction.map((v, col) => Math.abs(v / x[col])));
    let step = Math.min(1, .3 / maximumRelativeIncrement), accepted = false;
    for (let trial = 0; trial <= 20; trial++, step *= .5) {
      const candidate = x.map((v, col) => v + step * direction[col]);
      candidate[0] = Math.max(1e-7, Math.min(.30, candidate[0]));
      let projection;
      try {
        projection = projectXfoilDisplacement({ theta: scale * candidate[1], deltaStar: scale * candidate[2],
          ue: candidate[3], wakeGap, wake: true, mach, gamma });
      } catch (error) {
        if (error.message === 'DSLIM cannot repair a nonpositive thermal state.') continue;
        throw error;
      }
      if (projection.active) candidate[2] = projection.deltaStar / scale;
      if (!admissible(candidate)) continue;
      const next = Float64Array.from(residual(candidate));
      if (next.length !== 4 || !next.every(Number.isFinite)) throw new Error('Invalid native inverse-wake residual.');
      x = candidate; r = next; accepted = true;
      history.push({ iteration, residual: normInf(r), step, maximumRelativeIncrement,
        ...(trial ? { domainBacktracks: trial } : {}), ...(projection.active ? { projection } : {}) });
      break;
    }
    if (!accepted) return finish('No admissible native inverse-wake step.');
  }
  return finish('iteration limit');
}

export function evaluateFixedTripInterval(kernel, input) {
  const value = kernel.interval(input);
  if (input.regime !== 'transition' || value.transition?.forced) return value;
  // TRCHEK2's strict XT/XIFORC comparison can mislabel a root selected by
  // the trip. Use the same earlier-event check as automatic transition.
  // This leaves the forced mixed-interval equations unchanged and still
  // rejects a genuinely earlier natural transition in fixed-trip mode.
  const selected = evaluateTransitionInterval(kernel, input);
  if (!selected.transition?.forced) throw new Error('Natural transition precedes the fixed-trip research interval.');
  return selected;
}

// Local direct/inverse initialization only; the final solve still couples
// every surface/wake station and Euler unknown. `interval` supplies all native
// thermodynamic and transition-branch checks, including for FD perturbations.
export function initializeStreamtubeBLStation({ interval, properties, upstream, s, ue, regime, tripS, reynolds, initialState,
  allowSupersonicEdge = false, wakeGap = 0, mach = 0, gamma = 1.4 }) {
  if (typeof allowSupersonicEdge !== 'boolean') throw new Error('Invalid BL edge-domain option.');
  const scale = 1 / Math.sqrt(reynolds);
  const theta = initialState?.theta ?? upstream?.theta ?? .292 * Math.sqrt(s / (reynolds * ue));
  let deltaStar = initialState?.deltaStar ?? (upstream?.wakeGap !== undefined
    ? upstream.deltaStar - (upstream.wakeGap ?? 0) + wakeGap
    : upstream?.deltaStar ?? .647 * Math.sqrt(s / (reynolds * ue)));
  let seedProjection;
  // Copying an upstream wake thickness to a faster edge velocity can give
  // H>1 yet compressible Hk<=1. Initialize that otherwise-positive guess
  // with MRCHUE's fluid-thickness DSLIM (xbl.f:767-774), before its first
  // equation evaluation. Preserve every already-admissible seed exactly;
  // the local and global solves still require the original residuals.
  if (regime === 'wake' && mach > 0 && theta > 0 && deltaStar - wakeGap > theta) {
    const projected = projectXfoilDisplacement({ theta, deltaStar, ue, wakeGap, wake: true, mach, gamma });
    if (projected.rawHk <= 1) { seedProjection = projected; deltaStar = projected.deltaStar; }
  }
  const turbulent = ['transition', 'turbulent', 'wake'].includes(regime);
  // Amplification N and turbulent shear Ctau share a slot, not a physical
  // meaning. XFOIL xbl.f explicitly initializes a newly turbulent Ctau=.03.
  const aux = regime === 'transition' ? .03 : turbulent ? Math.max(.03, upstream?.aux ?? 0) : upstream?.aux ?? 0;
  const state = z => ({ s, ue: z.length === 4 ? z[3] : ue, aux: z[0], theta: scale * z[1], deltaStar: scale * z[2],
    ...(wakeGap ? { wakeGap } : {}) });
  const residual = z => interval({ upstream, downstream: state(z), regime, tripS }).residual;
  const admissible = z => {
    if (!z.every(Number.isFinite) || !(z[1] > 0) || !(z[2] - wakeGap / scale > z[1]) || (turbulent && !(z[0] > 0))) return false;
    try { residual(z); return true; } catch { return false; }
  };
  const initial = [aux, theta / scale, deltaStar / scale];
  // Preserve the actual physical failure reason for an invalid initial guess.
  residual(initial);
  const r = solveNewton({ initial, maxIterations: 30, tolerance: 1e-10, admissible, residual });
  if (r.converged) return { ...r, state: state(r.x), mode: 'direct', ...(seedProjection ? { seedProjection } : {}) };
  if (!properties || !upstream || regime === 'similarity')
    throw new Error(`BL station initialization failed (${regime}, s=${s}): ${r.reason}`);
  // MRCHUE in the supplied XFOIL source switches to a prescribed kinematic
  // shape parameter near direct-mode separation. This is only an initial
  // guess: the global system restores Euler pressure matching at this node.
  const hk = properties(upstream).rawHk, ds = s - upstream.s, hmax = regime === 'laminar' ? 3.8 : 2.5;
  if (!(hk > 1) || !Number.isFinite(hk)) throw new Error('Inverse BL initialization requires the raw kinematic shape parameter.');
  let targetHK;
  if (regime === 'wake') {
    const c = .03 * ds / upstream.theta; targetHK = hk;
    for (let k = 0; k < 3; k++) targetHK -= (targetHK + c * (targetHK - 1) ** 3 - hk) / (1 + 3 * c * (targetHK - 1) ** 2);
    targetHK = Math.max(targetHK, 1.01);
  } else {
    const transitionS = regime === 'transition' ? interval({ upstream, downstream: state(r.x), regime, tripS }).transition.s : tripS;
    const growth = regime === 'laminar' ? .03 * ds : regime === 'transition'
      ? .03 * (transitionS - upstream.s) - .15 * (s - transitionS) : -.15 * ds;
    targetHK = Math.max(hmax, hk + growth / upstream.theta);
  }
  // Prescribe the actual kinematic shape, not BLVAR's clipped closure Hk.
  // The latter has a flat derivative at a failed direct iterate and creates
  // a spurious zero row in the inverse Jacobian. At the inverse root both
  // values agree because every target lies above its closure lower bound.
  const inverseResidual = z => [...residual(z), properties(state(z)).rawHk - targetHK];
  const inverseAdmissible = z => {
    if (!(z[3] > 0) || !admissible(z)) return false;
    try { const p = properties(state(z)); return Number.isFinite(p.rawHk) && p.rawHk > 1
      && (allowSupersonicEdge ? Number.isFinite(p.machSquared) && p.machSquared >= 0 : p.machSquared < 1); } catch { return false; }
  };
  let inverse = solveNewton({ initial: [...r.x, ue], maxIterations: 30, tolerance: 1e-10,
    residual: inverseResidual, admissible: inverseAdmissible });
  let inverseRecovery;
  if (!inverse.converged && regime === 'wake') {
    // Native MRCHUE accepts bounded, admissible inverse updates even when
    // their intermediate residual increases. A monotone local line search
    // can stall before the same physical wake root on a decelerating seed.
    inverseRecovery = { method: 'xfoil-mrchue-bounded-inverse', reason: inverse.reason, history: inverse.history };
    inverse = solveXfoilInverseWakeStation({ initial: [...r.x, ue], residual: inverseResidual,
      admissible: inverseAdmissible, scale, mach, gamma, wakeGap });
  }
  if (!inverse.converged) throw Object.assign(new Error(`BL station inverse initialization failed (${regime}, s=${s}): ${inverse.reason}`), {
    code: 'BL_STATION_INVERSE_INITIALIZATION', diagnostics: { regime, reynolds, mach, gamma, tripS, targetHK,
      upstream: structuredClone(upstream), requestedUe: ue, originalSeed: state(initial), directState: state(r.x),
      retainedState: state(inverse.x), directReason: r.reason, history: inverse.history,
      ...(inverseRecovery ? { inverseRecovery } : {}) },
  });
  return { ...inverse, state: state(inverse.x), mode: 'inverse', targetHK, ...(seedProjection ? { seedProjection } : {}),
    ...(inverseRecovery ? { inverseRecovery } : {}),
    requestedUe: ue, directReason: r.reason, directHistory: r.history };
}

export function createStreamtubeBoundaryLayers(euler, initialEuler, { reynolds = 1e6, ncrit = 9, tripFractions, transitionMode = 'fixed-trip',
  allowSupersonicEdge = false, hkFloorLinearization = 'exact' } = {}) {
  const { layout, curves, fractions, initialStagnation } = euler, { lengthScale, mach, gamma } = euler.conditions, ne = layout.n;
  if (typeof allowSupersonicEdge !== 'boolean' || allowSupersonicEdge
    && (euler.conditions.flowModel !== 'compressible' || euler.conditions.streamwiseMode !== 'hybrid' || !euler.conditions.upwind))
    throw new Error('Supersonic BL initialization requires the historical hybrid Euler path.');
  if (!(reynolds > 0) || !Number.isFinite(reynolds)) throw new Error('Invalid coupled Reynolds number.');
  if (!['fixed-trip', 'automatic'].includes(transitionMode)) throw new Error('Unknown coupled transition mode.');
  const automatic = transitionMode === 'automatic';
  // Automatic transition is an eliminated nonlinear root in this Newton
  // system. Native TRCHEK's 5e-5 stopping error can jump when its iteration
  // count changes, corrupting local differences. Resolve the same root;
  // fixed-trip and original-Fortran reference behavior retain their settings.
  const kernel = createIntegralKernel({ reynolds, mach, gamma, ncrit, exactJacobian: true,
    hkFloorLinearization, ...(automatic ? { transitionTolerance: 1e-12 } : {}) });
  const scale = 1 / Math.sqrt(reynolds), arcs = curves.map(createContourArc), stations = [], surfaces = [], wakes = [];
  const hasFiniteBase = !!euler.baseGeometry?.some(Boolean);
  const deadAir = hasFiniteBase ? euler.baseGeometry.map(base => base && createXfoilDeadAirGap({ normalGap: base.width,
    upperDerivative: base.upperDerivative, lowerDerivative: base.lowerDerivative })) : [];
  const trips = tripFractions ?? layout.bodies.map(() => automatic ? [1, 1] : [.31, .31]);
  if (!Array.isArray(trips) || trips.length !== layout.elements || trips.some(row => !Array.isArray(row) || row.length !== 2
    || row.some(t => !Number.isFinite(t) || t <= 0 || (automatic ? t > 1 : t >= 1)))) throw new Error('Supply resolved upper/lower material trip fractions for each body.');
  const append = s => { const id = stations.length; stations.push({ ...s, id }); return id; };
  for (let body = 0; body < layout.elements; body++) {
    const b = layout.bodies[body];
    for (const side of ['upper', 'lower']) {
      const ids = []; for (let i = b.leadingIndex + 1; i <= b.trailingIndex; i++) ids.push(append({ kind: 'surface', body, side, i, k: i - b.leadingIndex }));
      const tripParameter = curves[body].branch(side, trips[body][side === 'upper' ? 0 : 1], initialStagnation[body]).parameter;
      surfaces.push({ body, side, ids, tripArc: arcs[body].at(tripParameter), tripParameter,
        ...(automatic ? { prescribedTrip: trips[body][side === 'upper' ? 0 : 1] < 1 } : {}) });
    }
    const ids = []; for (let i = b.trailingIndex; i <= layout.nx; i++) ids.push(append({ kind: 'wake', body, i, k: i - b.trailingIndex }));
    wakes.push({ body, ids });
  }
  const geometry = (eulerState, jacobian = false) => {
    const flow = euler.decode(eulerState), maps = jacobian
      ? euler.geometryDerivatives(eulerState, { includeDisplacement: true }) : null;
    const centerMap = (b, i) => {
      const result = new Map();
      if (maps) for (const col of new Set([...maps[b][i].at(-1).keys(), ...maps[b + 1][i][0].keys()])) {
        const delta = mean(maps[b][i].at(-1).get(col) ?? zero, maps[b + 1][i][0].get(col) ?? zero);
        if (col < ne) result.set(col, delta);
        else for (const [column, factor] of thicknessMap[col - ne]) {
          const prior = result.get(column) ?? zero;
          result.set(column, { x: prior.x + factor * delta.x, y: prior.y + factor * delta.y });
        }
      }
      return result;
    };
    const coordinates = [], surfaceData = surfaces.map(s => {
      const { body, side } = s, stag = flow.stagnation[body], arc = arcs[body], col = layout.globals.stagnation[body];
      const orientation = side === 'upper' ? -1 : 1, tripS = orientation * (s.tripArc - arc.at(stag)) / lengthScale;
      const tripDerivatives = new Map(); if (col !== null) tripDerivatives.set(col, -orientation * arc.speed(stag) * curves[body].length / lengthScale);
      for (const id of s.ids) {
        const f = fractions[body][side][stations[id].k], parameter = curves[body].branch(side, f, stag).parameter;
        const distance = orientation * (arc.at(parameter) - arc.at(stag)) / lengthScale, derivatives = new Map();
        if (col !== null) derivatives.set(col, orientation * ((1 - f) * arc.speed(parameter) - arc.speed(stag)) * curves[body].length / lengthScale);
        coordinates[id] = { s: distance, derivatives };
      }
      // With no prescribed trip, a zero-length turbulent portion at the TE
      // converts the auxiliary slot to shear before wake matching. The
      // surface thickness equations remain laminar up to that endpoint.
      return automatic && !s.prescribedTrip
        ? { tripS: coordinates[s.ids.at(-1)].s, tripDerivatives: new Map(coordinates[s.ids.at(-1)].derivatives), terminal: true }
        : { tripS, tripDerivatives };
    });
    for (const w of wakes) {
      const b = w.body, lower = surfaces.find(s => s.body === b && s.side === 'lower');
      const base = euler.baseGeometry?.[b], gapModel = deadAir[b];
      // Continue the physical wake center from the displaced TE. Starting
      // at the solid midpoint adds a fictitious first interval when the two
      // surface thicknesses differ, prematurely closing the dead-air gap.
      // Its thickness derivatives are needed even in a fixed wake chart.
      const te = layout.bodies[b].trailingIndex;
      let previous = mean(flow.nodes[b][te].at(-1), flow.nodes[b + 1][te][0]), previousMap = centerMap(b, te),
        distance = coordinates[lower.ids.at(-1)].s;
      let wakeDistance = 0; const wakeDistanceDerivatives = new Map();
      let derivatives = new Map(coordinates[lower.ids.at(-1)].derivatives);
      coordinates[w.ids[0]] = { s: distance, derivatives: new Map(derivatives),
        ...(base ? { wakeDistance: 0, wakeGap: base.width / lengthScale, wakeGapDerivatives: new Map() } : {}) };
      for (const id of w.ids.slice(1)) {
        const i = stations[id].i, p = mean(flow.nodes[b][i].at(-1), flow.nodes[b + 1][i][0]), currentMap = centerMap(b, i);
        const dx = p.x - previous.x, dy = p.y - previous.y, length = Math.hypot(dx, dy);
        if (!(length > 0)) throw new Error('Collapsed BL wake interval.');
        distance += length / lengthScale;
        if (base) wakeDistance += length;
        for (const col of new Set([...currentMap.keys(), ...previousMap.keys()])) {
          const a = currentMap.get(col) ?? zero, q = previousMap.get(col) ?? zero;
          add(derivatives, col, (dx * (a.x - q.x) + dy * (a.y - q.y)) / (length * lengthScale));
          if (base) add(wakeDistanceDerivatives, col, (dx * (a.x - q.x) + dy * (a.y - q.y)) / length);
        }
        const gap = base ? gapModel.at(wakeDistance) : null;
        coordinates[id] = { s: distance, derivatives: new Map(derivatives),
          ...(base ? { wakeDistance: wakeDistance / lengthScale, wakeGap: gap.gap / lengthScale,
            wakeGapDerivatives: new Map([...wakeDistanceDerivatives].map(([col, d]) => [col, gap.dDistance * d / lengthScale])) } : {}) };
        previous = p; previousMap = currentMap;
      }
    }
    return { coordinates, surfaceData };
  };
  const initialGeometry = geometry(initialEuler);
  const assignRegimes = (s, transition) => {
    s.transition = transition;
    s.ids.forEach((id, j) => { stations[id].regime = automatic && j === 0 && transition === 0 ? 'leading-transition'
      : j === 0 ? 'similarity' : j < transition ? 'laminar' : j === transition ? 'transition' : 'turbulent'; });
  };
  surfaces.forEach((s, k) => {
    const transition = s.ids.findIndex(id => initialGeometry.coordinates[id].s >= initialGeometry.surfaceData[k].tripS);
    if (automatic ? transition < 0 : transition < 1) throw new Error('Resolve at least one laminar station before each fixed transition trip.');
    assignRegimes(s, transition);
  });
  wakes.forEach(w => w.ids.forEach((id, k) => { stations[id].regime = k === 0 ? 'trailing-edge' : 'wake'; }));
  const thicknesses = x => ({ surfaces: layout.bodies.map((_, b) => {
    const s = Object.fromEntries(['upper', 'lower'].map(side => [side, surfaces.find(s => s.body === b && s.side === side)]));
    const leading = .5 * scale * lengthScale * (x[4 * s.upper.ids[0] + 2] + x[4 * s.lower.ids[0] + 2]);
    return Object.fromEntries(['upper', 'lower'].map(side => [side, [leading, ...s[side].ids.map(id => scale * lengthScale * x[4 * id + 2])]]));
  }), wakes: wakes.map(w => w.ids.slice(1).map(id => scale * lengthScale * x[4 * id + 2])) });
  const thicknessMap = euler.displacementParameters.map(p => {
    const row = new Map(), factor = scale * lengthScale;
    if (p.kind === 'wake') row.set(ne + 4 * wakes[p.body].ids[p.index + 1] + 2, factor);
    else for (const side of p.side === 'both' ? ['upper', 'lower'] : [p.side]) {
      const s = surfaces.find(s => s.body === p.body && s.side === side), id = s.ids[Math.max(0, p.index - 1)];
      row.set(ne + 4 * id + 2, factor * (p.side === 'both' ? .5 : 1));
    }
    return row;
  });
  const decode = (x, geo) => stations.map(({ id }) => ({ s: geo.coordinates[id].s, aux: x[4 * id], theta: scale * x[4 * id + 1], deltaStar: scale * x[4 * id + 2], ue: x[4 * id + 3],
    ...(geo.coordinates[id].wakeGap !== undefined ? { wakeGap: geo.coordinates[id].wakeGap } : {}) }));
  const snapshotActive = () => surfaces.map(s => s.transition);
  const restoreActive = saved => {
    if (!Array.isArray(saved) || saved.length !== surfaces.length || saved.some((j, k) => !Number.isInteger(j) || j < (automatic ? 0 : 1) || j >= surfaces[k].ids.length))
      throw new Error('Invalid streamtube transition active set.');
    surfaces.forEach((s, k) => assignRegimes(s, saved[k]));
  };
  const withSurfaceDiagnostic = (surface, operation) => {
    try { return operation(); }
    catch (error) {
      if (Number.isInteger(error.stationIndex) && surface.ids[error.stationIndex] !== undefined) {
        const id = surface.ids[error.stationIndex], station = stations[id];
        const element = layout.bodies[surface.body].element ?? surface.body;
        error.diagnostics = { ...error.diagnostics, stationLocation: { element, body: surface.body,
          side: surface.side, surfaceStationIndex: error.stationIndex, id, gridStation: station.i,
          regime: station.regime, indexBase: 0 } };
        error.message += ` Element ${element + 1}, ${surface.side}, station ${error.stationIndex + 1} (BL id ${id}, grid station ${station.i}).`;
      }
      throw error;
    }
  };
  const activeTargets = (eulerState, x) => {
    const geo = geometry(eulerState);
    if (automatic) {
      if (!x) throw new Error('Automatic transition selection requires the complete BL state.');
      const states = decode(x, geo);
      return surfaces.map((s, k) => {
        const target = withSurfaceDiagnostic(s, () => selectSurfaceTransition(kernel, s.ids.map(id => states[id]), { tripS: geo.surfaceData[k].tripS }));
        if (target.index === null) throw new Error('Missing terminal transition-to-wake interval.');
        const up = target.index ? states[s.ids[target.index - 1]].s : 0, down = states[s.ids[target.index]].s;
        // Selection integrates N from stagnation, whereas the mixed row uses
        // the actual packed upstream N. Away from a root those can disagree
        // even when every packed N is below Ncrit and the index is unchanged.
        // Signal only an unusable mixed interval; valid Newton N remains free.
        // The caller repairs the actual N prefix before evaluating the same
        // simultaneous equations, rather than substituting N inside a row.
        const amplificationReconciliation = target.index === s.transition && target.index > 0
          && states[s.ids[target.index - 1]].aux < ncrit
          && !withSurfaceDiagnostic(s, () => checkAutomaticTransition(kernel, {
            upstream: states[s.ids[target.index - 1]], downstream: states[s.ids[target.index]],
            tripS: geo.surfaceData[k].tripS,
          })).transition;
        return { body: s.body, side: s.side, from: s.transition, to: target.index, tripS: geo.surfaceData[k].tripS,
          transitionS: target.s, kind: target.kind === 'forced' && geo.surfaceData[k].terminal ? 'trailing-edge' : target.kind,
          fraction: (target.s - up) / (down - up), ...(amplificationReconciliation ? { amplificationReconciliation: true } : {}) };
      });
    }
    return surfaces.map((s, k) => {
      const tripS = geo.surfaceData[k].tripS;
      const target = s.ids.findIndex(id => geo.coordinates[id].s >= tripS);
      if (target < 1) throw new Error('A fixed trip requires at least one resolved upstream laminar station.');
      const up = geo.coordinates[s.ids[target - 1]].s, down = geo.coordinates[s.ids[target]].s;
      return { body: s.body, side: s.side, from: s.transition, to: target, tripS,
        fraction: (tripS - up) / (down - up) };
    });
  };
  const updateActive = (x, eulerState, { reinitializeAmplification = false, reconcileAmplificationSurfaces = [] } = {}) => {
    if (!Array.isArray(reconcileAmplificationSurfaces) || reconcileAmplificationSurfaces.some(k => !Number.isInteger(k) || k < 0 || k >= surfaces.length))
      throw new Error('Invalid amplification-reconciliation surfaces.');
    if ((reinitializeAmplification || reconcileAmplificationSurfaces.length) && !automatic) throw new Error('Amplification reinitialization requires automatic transition.');
    const reconcile = new Set(reconcileAmplificationSurfaces);
    const targets = activeTargets(eulerState, x), changes = targets.filter(t => t.from !== t.to);
    if (!changes.length && !reinitializeAmplification && !reconcile.size) return { changed: false, changes: [] };
    // Prepare every conversion before committing any phase or state. An
    // invalid native transition calculation cannot leave a half-updated BL.
    const candidate = x.slice(), geo = geometry(eulerState), states = decode(candidate, geo);
    if (automatic) {
      const prepared = surfaces.map((s, k) => withSurfaceDiagnostic(s, () => prepareSurfaceTransition(kernel, s.ids.map(id => states[id]),
        { previousIndex: s.transition, tripS: geo.surfaceData[k].tripS, reinitializeAmplification: reinitializeAmplification || reconcile.has(k) })));
      // Convert N/shear where the regime changes, using both affected BL
      // intervals. Downstream onset can preserve an existing Ctau; upstream
      // onset needs a new shear guess in place of laminar amplification.
      prepared.forEach((p, k) => {
        const s = surfaces[k], j = p.index;
        if (p.kind !== 'natural' || !(j !== s.transition && j > 0 && j + 1 < s.ids.length)) return;
        const conversion = p.converted.find(c => c.index === j);
        const select = j < s.transition ? minimizeMixedTransitionShear : chooseMixedTransitionShear;
        const selected = select(kernel, {
          upstream: { ...states[s.ids[j - 1]], aux: p.auxiliary[j - 1] },
          downstream: states[s.ids[j]], following: { ...states[s.ids[j + 1]], aux: p.auxiliary[j + 1] },
          tripS: geo.surfaceData[k].tripS, initializedAux: p.auxiliary[j], shearWeight: 20,
        });
        p.auxiliary[j] = selected.aux; conversion.aux = selected.aux;
        conversion.shearSelection = selected.diagnostics;
        if (selected.diagnostics.selected !== 'initialized') {
          if (selected.diagnostics.selected === 'existing') conversion.shearPreserved = true;
          delete conversion.transitionShearInitialization;
          delete conversion.shearInitialization;
        }
      });
      // A terminal fallback initializes the junction below, but does not
      // change the meaning of any existing downstream wake unknown.
      prepared.forEach((p, k) => {
        surfaces[k].ids.forEach((id, j) => { candidate[4 * id] = p.auxiliary[j]; });
        targets[k].converted = p.converted.map(c => ({ ...c, id: surfaces[k].ids[c.index] }));
        if (p.changed) targets[k].downstreamShearPreserved = true;
      });
      // Keep the merged TE shear consistent only for a terminal fallback event.
      // This is the existing algebraic theta-weighted junction,
      // using fluid momentum thickness; a finite-base gap is not a weight.
      // Preserve the downstream wake solution. Re-marching shear alone at
      // unchanged thickness/velocity can mask a bad transition trial by
      // closing just its shear rows, then drive subsequent thickness updates
      // onto the Hk floor. The full Newton solve closes the wake transport.
      // All values remain private until every junction has succeeded.
      const affectedBodies = new Set(prepared.flatMap((p, k) =>
        p.converted.some(c => c.terminalShearInitialization) ? [surfaces[k].body] : []));
      for (const body of affectedBodies) {
        const upper = surfaces.find(s => s.body === body && s.side === 'upper').ids.at(-1);
        const lower = surfaces.find(s => s.body === body && s.side === 'lower').ids.at(-1);
        const wake = wakes.find(w => w.body === body), id = wake.ids[0], a = candidate[4 * upper + 1], b = candidate[4 * lower + 1];
        const aux = (a * candidate[4 * upper] + b * candidate[4 * lower]) / (a + b);
        if (!(aux > 0) || !Number.isFinite(aux)) throw new Error('Invalid terminal-transition wake shear initialization.');
        const initialization = { id, oldAux: candidate[4 * id], aux, method: 'theta-weighted-TE-shear' };
        candidate[4 * id] = aux;
        targets.forEach(t => { if (t.body === body && t.converted.some(c => c.transitionShearInitialization)) {
          t.terminalWakeShearInitialization = initialization;
          t.downstreamWakeShearPreserved = true;
        } });
      }
      const amplificationReconciliations = [...reconcile].map(k => ({ body: surfaces[k].body, side: surfaces[k].side,
        from: surfaces[k].transition, to: prepared[k].index, criterion: ncrit,
        stations: surfaces[k].ids.filter(id => candidate[4 * id] !== x[4 * id])
          .map(id => ({ id, oldAux: x[4 * id], aux: candidate[4 * id] })) })).filter(r => r.stations.length);
      x.set(candidate); surfaces.forEach((s, k) => assignRegimes(s, prepared[k].index));
      return { changed: changes.length > 0 || amplificationReconciliations.length > 0,
        changes: [...changes, ...amplificationReconciliations.map(r => ({ ...r, kind: 'amplification-reconciliation', auxiliaryOnly: true }))],
        ...(amplificationReconciliations.length ? { amplificationReconciliations } : {}) };
    }
    targets.forEach((target, k) => {
      if (target.from === target.to) return;
      const surface = surfaces[k]; target.converted = [];
      for (let j = Math.min(target.from, target.to); j < Math.max(target.from, target.to); j++) {
        const id = surface.ids[j], oldAux = candidate[4 * id]; let aux;
        if (j < target.to) {
          const upstream = states[surface.ids[j - 1]], downstream = states[id];
          // Only the newly laminar amplification is initialized. Velocity
          // and both thickness equations remain in the global coupled solve.
          const check = kernel.transitionCheck({ upstream, downstream, tripS: target.tripS });
          if (check.transition || !(check.amplification < ncrit)) throw new Error('Natural transition precedes the fixed material trip during an interval change.');
          aux = check.amplification;
        } else aux = kernel.station({ ...states[id], aux: .03 }, 'turbulent').transitionShear;
        if (!Number.isFinite(aux) || (j >= target.to && !(aux > 0))) throw new Error('Invalid transition auxiliary-state transfer.');
        candidate[4 * id] = aux; states[id] = { ...states[id], aux };
        target.converted.push({ id, oldAux, aux, to: j < target.to ? 'laminar' : 'turbulent' });
      }
    });
    x.set(candidate); surfaces.forEach((s, k) => assignRegimes(s, targets[k].to));
    return { changed: true, changes };
  };
  const interval = input => {
    if (automatic && input.regime === 'transition') return evaluateTransitionInterval(kernel, input);
    return evaluateFixedTripInterval(kernel, input);
  };
  const blocks = (states, geo) => {
    const result = [];
    surfaces.forEach((s, k) => {
      const tripS = geo.surfaceData[k].tripS, up = s.ids[s.transition - 1], down = s.ids[s.transition];
      if (!automatic && !(states[up].s < tripS && tripS <= states[down].s)) throw new Error('A material trip left its resolved research interval.');
      s.ids.forEach((id, j) => {
        const previous = j ? s.ids[j - 1] : id, regime = stations[id].regime;
        if ((regime === 'laminar' || regime === 'similarity') && states[id].aux >= ncrit) throw new Error(automatic
          ? 'Amplification left the active laminar interval.' : 'Natural transition is outside the fixed-trip research model.');
        result.push({ id, previous, regime, trip: geo.surfaceData[k], input: { upstream: states[previous], downstream: states[id], regime, tripS } });
      });
    });
    for (const w of wakes) for (const id of w.ids.slice(1)) result.push({ id, previous: id - 1, regime: 'wake', input: { upstream: states[id - 1], downstream: states[id], regime: 'wake' } });
    return result;
  };
  const evaluate = (x, eulerState, { jacobian = false } = {}) => {
    const geo = geometry(eulerState, jacobian), states = decode(x, geo), residual = new Float64Array(4 * stations.length);
    const rows = jacobian ? Array.from({ length: residual.length }, () => new Map()) : null;
    const transitions = [];
    for (const block of blocks(states, geo)) {
      const { id, previous, regime, input, trip } = block;
      const value = regime === 'leading-transition' ? evaluateLeadingTransitionInterval(kernel, input, { jacobian })
        : automatic && regime === 'transition' ? evaluateTransitionInterval(kernel, input, { jacobian }) : interval(input);
      const turbulent = ['leading-transition', 'transition', 'turbulent', 'wake'].includes(regime);
      if (automatic && value.transition) transitions.push({ body: stations[id].body, side: stations[id].side, id,
        ...value.transition, kind: value.transition.forced ? (trip.terminal ? 'trailing-edge' : 'forced') : 'natural' });
      const rowScale = [turbulent ? 20 : 1, 1, 1];
      for (let r = 0; r < 3; r++) residual[4 * id + r] = rowScale[r] * value.residual[r];
      if (!jacobian) continue;
      // Native Hk-floor linearization retains the XFOIL quasi-Newton shape
      // sensitivity below the closure floor; other kernel partials use
      // exact mode. Transition retains its resolved-root local partials.
      // Transition still needs local numerical partials; these do not
      // finite-difference the Euler system.
      const numerical = regime === 'transition';
      const derivative = (side, key, k) => {
        if (value.partials) return value.partials[side].map(row => row[k]);
        if (!numerical) return value[side].map(row => row[k]);
        const v = input[side][key], h = Math.cbrt(Number.EPSILON) * Math.max(Math.abs(v), k === 0 ? .01 : k < 3 ? 1e-7 : 1e-6);
        if (regime === 'transition' && key === 's') return boundedResidualPartial(
          s => interval({ ...input, [side]: { ...input[side], s } }).residual, v, {
            step: h, base: value.residual,
            lower: side === 'upstream' ? 0 : input.tripS,
            upper: side === 'upstream' ? input.tripS : Infinity });
        const p = interval({ ...input, [side]: { ...input[side], [key]: v + h } }).residual;
        const m = interval({ ...input, [side]: { ...input[side], [key]: v - h } }).residual;
        return p.map((v, r) => (v - m[r]) / (2 * h));
      };
      for (const side of ['similarity', 'leading-transition'].includes(regime) ? ['downstream'] : ['upstream', 'downstream']) {
        const node = side === 'upstream' ? previous : id;
        ['aux', 'theta', 'deltaStar', 'ue', 's'].forEach((key, k) => {
          let d;
          try { d = derivative(side, key, k); }
          catch (error) { throw new Error(`BL ${regime} derivative at station ${id}, ${side}.${key}: ${error.message}`, { cause: error }); }
          for (let r = 0; r < 3; r++) {
            if (k < 4) add(rows[4 * id + r], ne + 4 * node + k, rowScale[r] * d[r] * (k === 1 || k === 2 ? scale : 1));
            else for (const [col, factor] of geo.coordinates[node].derivatives) add(rows[4 * id + r], col, rowScale[r] * d[r] * factor);
          }
        });
        if (regime === 'wake' && geo.coordinates[node].wakeGapDerivatives?.size) {
          // BLPRV's fluid displacement is total deltaStar minus DW. The
          // only additional explicit DW dependence is HWA in momentum and
          // energy (XFOIL BLDIF). No finite differences of the Euler rows.
          const dDelta = derivative(side, 'deltaStar', 2), theta = input[side].theta;
          const explicit = Math.log(input.downstream.ue / input.upstream.ue) / (2 * theta);
          const dGap = [-dDelta[0], -dDelta[1] + explicit, -dDelta[2] - explicit];
          for (let r = 0; r < 3; r++) for (const [col, factor] of geo.coordinates[node].wakeGapDerivatives)
            add(rows[4 * id + r], col, rowScale[r] * dGap[r] * factor);
        }
      }
      if (['transition', 'leading-transition'].includes(regime)) {
        const h = Math.cbrt(Number.EPSILON) * input.tripS;
        const d = value.partials?.trip ?? boundedResidualPartial(tripS => interval({ ...input, tripS }).residual, input.tripS,
          { step: h, base: value.residual, lower: input.upstream.s, upper: input.downstream.s });
        for (let r = 0; r < 3; r++) for (const [col, factor] of trip.tripDerivatives) add(rows[4 * id + r], col, rowScale[r] * d[r] * factor);
      }
    }
    for (const w of wakes) {
      const id = w.ids[0], upper = surfaces.find(s => s.body === w.body && s.side === 'upper').ids.at(-1), lower = surfaces.find(s => s.body === w.body && s.side === 'lower').ids.at(-1);
      const a = states[upper], b = states[lower], sum = a.theta + b.theta, aux = (a.theta * a.aux + b.theta * b.aux) / sum;
      const r = kernel.trailingEdge(a, b, states[id], (euler.baseGeometry?.[w.body]?.width ?? 0) / lengthScale).residual;
      for (let k = 0; k < 3; k++) residual[4 * id + k] = r[k] * (k === 0 ? 20 : 1 / scale);
      if (rows) {
        for (const node of [upper, lower]) {
          add(rows[4 * id], ne + 4 * node, -20 * states[node].theta / sum);
          add(rows[4 * id], ne + 4 * node + 1, -20 * (states[node].aux - aux) * scale / sum);
          for (const k of [1, 2]) add(rows[4 * id + k], ne + 4 * node + k, -1);
        }
        for (let k = 0; k < 3; k++) add(rows[4 * id + k], ne + 4 * id + k, k === 0 ? 20 : 1);
      }
    }
    return { states, geometry: geo, residual, rows, ...(automatic ? { transitions } : {}) };
  };
  const initialization = [];
  const initialize = (edgeVelocity, { method = 'direct-inverse' } = {}) => {
    if (!['direct-inverse', 'mrchue'].includes(method)) throw new Error('Unknown BL initialization method.');
    initialization.length = 0;
    const x = new Float64Array(4 * stations.length), geo = initialGeometry;
    const put = (id, a) => x.set([a.aux, a.theta / scale, a.deltaStar / scale, a.ue], 4 * id);
    const get = id => ({ ...decode(x, geo)[id] });
    const march = (id, previous, regime, tripS, initialState) => {
      const ue = edgeVelocity(stations[id]), s = geo.coordinates[id].s, upstream = previous === null ? null : get(previous);
      try {
        const r = initializeStreamtubeBLStation({ interval: regime === 'leading-transition'
          ? input => evaluateLeadingTransitionInterval(kernel, input) : interval,
          properties: state => kernel.station(state, regime === 'wake' ? 'wake' : turbulentRegime(regime)),
          upstream: upstream ?? undefined, s, ue, regime: regime === 'leading-transition' ? 'transition' : regime, tripS, reynolds, mach, gamma, initialState,
          ...(geo.coordinates[id].wakeGap !== undefined ? { wakeGap: geo.coordinates[id].wakeGap } : {}),
          ...(allowSupersonicEdge ? { allowSupersonicEdge } : {}) });
        put(id, r.state);
        initialization.push({ id, body: stations[id].body, regime, mode: r.mode, targetHK: r.targetHK,
          requestedUe: ue, initialUe: r.state.ue, residual: r.history.at(-1).residual,
          iterations: r.history.length - 1, directReason: r.directReason,
          ...(r.seedProjection ? { seedProjection: r.seedProjection } : {}),
          ...(r.inverseRecovery ? { inverseRecovery: r.inverseRecovery } : {}) });
        return r.state;
      } catch (error) {
        throw new Error(`BL initialization failed at body ${stations[id].body}, ${regime} station ${id}: ${error.message}`, { cause: error });
      }
    };
    if (!automatic) surfaces.forEach((s, k) => s.ids.forEach((id, j) => march(id, j ? s.ids[j - 1] : null, stations[id].regime, geo.surfaceData[k].tripS)));
    else {
      const native = new Map();
      if (method === 'mrchue') for (let body = 0; body < layout.elements; body++) {
        // MRCHUE's geometry-free adapter uses the terminal material trip.
        // Explicit earlier trips retain the resolved-interval initializer,
        // including its virtual leading station; final equations are shared.
        const free = surfaces.filter(s => s.body === body && !s.prescribedTrip);
        if (!free.length) continue;
        const profiles = free.map(s => s.ids.map(id => ({ s: geo.coordinates[id].s, ue: edgeVelocity(stations[id]) })));
        const r = initializeXfoilSurfaces(profiles.length === 2 ? profiles : [profiles[0], profiles[0]], { reynolds, mach, gamma, ncrit });
        free.forEach((s, k) => native.set(s, { ...r.surfaces[k], messages: r.messages, localConvergenceWarnings: r.localConvergenceWarnings }));
      }
      surfaces.forEach((s, k) => {
      if (native.has(s)) {
        const r = native.get(s); assignRegimes(s, r.transition);
        s.ids.forEach((id, j) => {
          put(id, r.states[j]);
          initialization.push({ id, body: s.body, regime: stations[id].regime, method: 'mrchue',
            mode: r.targetHK[j] > 0 ? 'inverse' : 'direct', targetHK: r.targetHK[j],
            requestedUe: edgeVelocity(stations[id]), initialUe: r.states[j].ue,
            ...(j === 0 ? { messages: r.messages, localConvergenceWarnings: r.localConvergenceWarnings } : {}) });
        });
        return;
      }
      const tripS = geo.surfaceData[k].tripS;
      let transition = s.transition === 0 ? 0 : null;
      s.ids.forEach((id, j) => {
        const previous = j ? s.ids[j - 1] : null;
        if (j === 0) { march(id, null, transition === 0 ? 'leading-transition' : 'similarity', tripS); return; }
        if (transition !== null) { march(id, previous, 'turbulent', tripS); return; }
        // Direct/inverse laminar initialization only. The final coupled
        // equations retain N, thickness, edge velocity and all grid unknowns.
        const guess = march(id, previous, 'laminar', tripS);
        const check = checkAutomaticTransition(kernel, { upstream: get(previous), downstream: guess, tripS });
        if (check.transition) {
          transition = j; assignRegimes(s, j);
          march(id, previous, 'transition', tripS, guess);
        }
      });
      if (transition === null) throw new Error('Automatic initialization missed the terminal transition interval.');
      assignRegimes(s, transition);
      });
    }
    for (const w of wakes) {
      const upper = get(surfaces.find(s => s.body === w.body && s.side === 'upper').ids.at(-1)), lower = get(surfaces.find(s => s.body === w.body && s.side === 'lower').ids.at(-1));
      const theta = upper.theta + lower.theta;
      put(w.ids[0], { ue: edgeVelocity(stations[w.ids[0]]), theta,
        deltaStar: upper.deltaStar + lower.deltaStar + (euler.baseGeometry?.[w.body]?.width ?? 0) / lengthScale,
        aux: (upper.aux * upper.theta + lower.aux * lower.theta) / theta });
      w.ids.slice(1).forEach(id => march(id, id - 1, 'wake'));
    }
    return x;
  };
  const admissible = x => x.length === 4 * stations.length && x.every(Number.isFinite) && stations.every(({ id, regime }) =>
    x[4 * id + 1] > 0 && x[4 * id + 2] > x[4 * id + 1] && x[4 * id + 3] > 0
    && (!['leading-transition', 'transition', 'turbulent', 'trailing-edge', 'wake'].includes(regime) || x[4 * id] > 0));
  return { stations, surfaces, wakes, kernel, scale, thicknessMap, thicknesses, geometry, evaluate, initialize, initialization, admissible, trips,
    ...(hasFiniteBase ? { hasFiniteBase, initialWakeGaps: Float64Array.from(initialGeometry.coordinates, p => p.wakeGap ?? 0) } : {}),
    snapshotActive, restoreActive, activeTargets, updateActive, transitionMode };
}

function turbulentRegime(regime) {
  return ['leading-transition', 'transition', 'turbulent'].includes(regime) ? 'turbulent' : 'laminar';
}
