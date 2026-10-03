// SPDX-License-Identifier: GPL-2.0-or-later
// Research simultaneous Euler/grid/BL/wake equations on every element.
// Fixed resolved material trips; smooth subcritical pressure/edge matching.
import { createStreamtubeBodySystem, solveStreamtubeBody } from './streamtube-body.js';
import { coupledConvergenceSatisfied } from './streamtube-coupled-convergence.js';
import { createStreamtubeBoundaryLayers } from './streamtube-boundary-layers.js';
import { streamtubeEdgeState, streamtubeEdgeVelocity, streamtubeEdgePressure } from './streamtube-edge-velocity.js';
import { streamtubeMeshSnapshot, streamtubeMeshQuality, STREAMTUBE_MINIMUM_CORNER_SINE } from './streamtube-mesh-preview.js';
import { streamtubeCornerConstraints } from './streamtube-corner-constraints.js';
import { solveNewton } from '../numerics/newton.js';
import { takeDoglegStep } from '../numerics/dogleg.js';
import { solveSparseDirect } from '../numerics/klu.js';
import { sparseMatrixFromRows } from '../numerics/sparse-rows.js';
import { createPhysicalBLDomain } from '../viscous/physical-domain.js';
import { transferStreamtubeGeometry } from './streamtube-geometry.js';
import { seedStreamtubeWakeBanks } from './streamtube-displacement.js';
import { streamtubeBLEdgeThermodynamics } from './streamtube-bl-edge-thermodynamics.js';
import { initialStreamtubeDisplacement } from './streamtube-geometry.js';
import { XFOIL_BL_EDGE_SPEED_LIMIT } from './streamtube-coupled-xfoil-update.js';

const maximum = a => a.reduce((peak, value) => Math.max(peak, Math.abs(value)), -Infinity);
const add = (row, col, value) => row.set(col, (row.get(col) ?? 0) + value);
export function createCoupledStreamtubeBody(input, { reynolds = 1e6, ncrit = 9, tripFractions, initialEuler, initialBL,
  edgeMatching = 'pressure', transitionMode = 'fixed-trip', transitionState, blInitialization = 'direct-inverse', blThermodynamics,
  hkFloorLinearization, geometryReplay, initialEdgeVelocity } = {}) {
  if (initialEdgeVelocity !== undefined && typeof initialEdgeVelocity !== 'function')
    throw new Error('Initial BL edge velocity must be a function.');
  if (input.displacement !== undefined) throw new Error('The coupled system computes displacement from its BL unknowns.');
  if (hkFloorLinearization !== undefined && !['exact', 'native'].includes(hkFloorLinearization))
    throw new Error('Unknown coupled Hk-floor linearization policy.');
  // A numerical direction policy only. Omission and explicit exact retain
  // legacy kernels/checkpoint schemas; native is recorded when requested.
  const hkPolicy = hkFloorLinearization === 'native' ? { hkFloorLinearization } : {};
  if (geometryReplay !== undefined && !['legacy', 'preserve-undisplaced'].includes(geometryReplay))
    throw new Error('Unknown coupled geometry replay policy.');
  const preserveUndisplaced = geometryReplay === 'preserve-undisplaced';
  if (preserveUndisplaced && (!initialEuler?.undisplacedNodes || initialBL == null))
    throw new Error('Preserving coupled geometry requires a supplied undisplaced chart and BL state.');
  const geometryPolicy = preserveUndisplaced ? { geometryReplay } : {};
  if (!['pressure', 'section-velocity', 'section-velocity-distance'].includes(edgeMatching)) throw new Error('Unknown streamtube edge matching formulation.');
  const historicalTransonic = blThermodynamics === 'historical-common-isentrope';
  if (blThermodynamics !== undefined && !historicalTransonic) throw new Error('Unknown coupled BL thermodynamic model.');
  if (historicalTransonic && (input.streamwiseMode !== 'hybrid' || !input.upwind
    || input.flowModel !== undefined && input.flowModel !== 'compressible' || edgeMatching !== 'section-velocity'))
    throw new Error('Historical transonic BL requires compressible hybrid Euler, upwinding and arithmetic section-velocity matching.');
  const edgeInterpolation = edgeMatching === 'section-velocity-distance' ? 'distance-weighted' : 'arithmetic';
  let displacement = { surfaces: input.bodies.map(b => ({ upper: Array(b.trailingIndex - b.leadingIndex + 1).fill(0), lower: Array(b.trailingIndex - b.leadingIndex + 1).fill(0) })),
    wakes: input.bodies.map(b => Array(input.outerLower.length - 1 - b.trailingIndex).fill(0)) };
  const euler = createStreamtubeBodySystem({ ...input, displacement }), ne = euler.layout.n;
  if (euler.baseGeometry) {
    displacement = initialStreamtubeDisplacement(euler.layout, euler.baseGeometry);
    euler.setDisplacement(displacement);
  }
  let start;
  if (initialEuler) start = euler.adoptGeometry(initialEuler.x, initialEuler.undisplacedNodes ?? initialEuler.nodes);
  else {
    // The inviscid precursor has coincident banks. Solve its reduced,
    // zero-gap system and transfer the complete physical flow to the new
    // coordinates before activating the simultaneous viscous equations.
    const precursor = euler.layout.independentWakeBanks ? createStreamtubeBodySystem({ ...input, displacement, wakeGeometry: 'centerline', wakeDisplacementMotion: 'fixed' }) : euler;
    const inviscid = solveStreamtubeBody(precursor, { maxIterations: 25, tolerance: 1e-11 });
    if (!inviscid.converged) throw new Error(`Coupled Euler initialization failed: ${inviscid.reason}`);
    start = precursor === euler ? inviscid.x : transferStreamtubeGeometry(precursor, inviscid.x, euler);
  }
  const bl = createStreamtubeBoundaryLayers(euler, start, { reynolds, ncrit, tripFractions, transitionMode,
    ...hkPolicy,
    ...(historicalTransonic ? { allowSupersonicEdge: true } : {}) });
  if (transitionMode === 'automatic' && initialBL != null && transitionState === undefined)
    throw new Error('A supplied automatic BL state requires its transition-interval map.');
  if (transitionState !== undefined) {
    if (transitionMode !== 'automatic') throw new Error('An explicit transition-interval map requires automatic mode.');
    bl.restoreActive(transitionState);
  }
  const blDomain = createPhysicalBLDomain(bl.kernel.parameters);
  const bankOutlet = s => input.wakeOutlet === 'banks' && s.kind === 'wake' && s.i === euler.layout.nx;
  const pressure = (flow, s) => s.kind === 'surface' ? flow.bodyPressure(s.body, s.i, s.side)
    : .5 * (flow.bodyPressure(s.body, Math.min(s.i, euler.layout.nx - 1), 'upper') + flow.bodyPressure(s.body, Math.min(s.i, euler.layout.nx - 1), 'lower'));
  const edge = (flow, s) => {
    if (edgeMatching === 'pressure') return { pressure: pressure(flow, s), matching: edgeMatching };
    const sides = s.kind === 'surface' ? [s.side] : ['upper', 'lower'];
    const values = sides.map(side => {
      const g = side === 'lower' ? s.body : s.body + 1, j = side === 'lower' ? euler.layout.tubes[g] - 1 : 0;
      return streamtubeEdgeVelocity(flow.cells[Math.min(s.i, euler.layout.nx - 1) - 1][g][j], .2, edgeInterpolation);
    });
    // This assembly currently has one merged wake BL: average its two
    // bank velocities. It is not a model for unequal-entropy wake confluence.
    const ue = values.reduce((sum, v) => sum + v.ue, 0) / values.length;
    if (euler.conditions.flowModel !== 'incompressible') {
      const { h0, gamma } = euler.conditions, h = h0 - .5 * ue * ue;
      if (historicalTransonic) {
        // Giles SETBL pp.226–227 reconstructs BL gas from one RSTOUT and
        // physical Ue even when the Euler flow is supersonic. Keep that
        // approximation explicit; never inject Euler shock entropy into
        // the unchanged BL equations. Each wake bank must remain thermal.
        if (!(h > 0) || !Number.isFinite(h) || values.some(v =>
          !(v.ue > 0) || !Number.isFinite(v.ue) || !(h0 - .5 * v.ue * v.ue > 0)))
          throw new Error('Corrected historical BL edge velocity exceeds its physical thermal domain.');
      } else if (!(h > 0) || ue * ue >= (gamma - 1) * h) throw new Error('Corrected BL edge velocity is sonic or inadmissible.');
    }
    return { ue, sides: values, matching: edgeMatching };
  };
  // A saved coupled state need not have an admissible zero-displacement
  // Euler grid. Evaluate that auxiliary grid only when it is actually
  // needed to initialize a missing BL profile.
  let blStart = initialBL;
  if (blStart == null) {
    const inviscid = euler.evaluate(start);
    const geometry = initialEdgeVelocity ? bl.geometry(start) : null;
    blStart = bl.initialize(s => initialEdgeVelocity ? initialEdgeVelocity(s, bl, geometry)
      : edgeMatching === 'pressure' ? streamtubeEdgeState(pressure(inviscid, s), euler.conditions).ue
        : edge(inviscid, s).ue, { method: blInitialization });
  }
  if (!bl.admissible(blStart)) throw new Error('Invalid simultaneous BL initial state.');
  const n = ne + blStart.length, initial = new Float64Array(n); initial.set(start); initial.set(blStart, ne);
  // Restore the normal chart on the actual displaced flow grid too. The
  // zero-displacement adoption above recovers coordinates, not the final
  // wake motion directions for this supplied/initialized BL state.
  euler.setDisplacement(bl.thicknesses(blStart));
  if (euler.layout.independentWakeBanks && initialBL == null)
    initial.set(euler.adoptGeometry(initial.subarray(0, ne), seedStreamtubeWakeBanks(euler, initial.subarray(0, ne))));
  if (preserveUndisplaced) euler.refreshGeometryDirections(initial.subarray(0, ne));
  else initial.set(euler.rebase(initial.subarray(0, ne)));
  const evaluateInternal = (x, jacobian, linearizeEuler) => {
    if (x.length !== n || !x.every(Number.isFinite) || !bl.admissible(x.subarray(ne))) throw new Error('Inadmissible coupled streamtube state.');
    euler.setDisplacement(bl.thicknesses(x.subarray(ne)));
    // Euler value and derivatives share one fresh evaluation. Nothing is
    // retained across residual calls, trial steps, charts or transition events.
    const linearized = linearizeEuler ? euler.evaluateJacobian(x.subarray(0, ne), { sparse: true, includeDisplacement: true }) : null;
    const outer = linearized ? linearized.value : euler.evaluate(x.subarray(0, ne));
    const layers = bl.evaluate(x.subarray(ne), x.subarray(0, ne), { jacobian });
    const residual = new Float64Array(n); residual.set(outer.residual); residual.set(layers.residual, ne);
    const edges = bl.stations.map(s => bankOutlet(s) ? { ue: x[ne + 4 * s.id + 3], matching: 'outlet-bank-tangency' } : edge(outer, s));
    bl.stations.forEach(s => {
      const [lower, upper] = outer.outletBankTangency[s.body];
      residual[ne + 4 * s.id + 3] = bankOutlet(s) ? .5 * (lower - upper)
        : edgeMatching === 'pressure' ? streamtubeEdgePressure(x[ne + 4 * s.id + 3], euler.conditions).pressure - edges[s.id].pressure
          : x[ne + 4 * s.id + 3] - edges[s.id].ue;
    });
    const families = { euler: maximum(outer.residual), boundaryLayer: maximum(layers.residual), edgeMatching: bl.stations.reduce((peak, { id }) => Math.max(peak, Math.abs(residual[ne + 4 * id + 3])), -Infinity) };
    return { residual, outer, layers, edges, families, ...(linearized ? { eulerJacobian: linearized.jacobian } : {}) };
  };
  const evaluate = (x, { jacobian = false } = {}) => evaluateInternal(x, jacobian, false);
  const jacobian = (x, { sparse = true } = {}) => {
    const value = evaluateInternal(x, true, true), j = value.eulerJacobian;
    const rows = Array.from({ length: n }, () => new Map());
    const chain = (row, col, d) => {
      if (col < ne) add(row, col, d);
      else for (const [blCol, scale] of bl.thicknessMap[col - ne]) add(row, blCol, d * scale);
    };
    for (let row = 0; row < ne; row++) {
      for (let p = j.state.rowPtr[row]; p < j.state.rowPtr[row + 1]; p++) add(rows[row], j.state.colIndex[p], j.state.values[p]);
      for (const [col, d] of j.displacement[row]) chain(rows[row], ne + col, d);
    }
    value.layers.rows.forEach((row, k) => { rows[ne + k] = new Map(row); });
    for (const s of bl.stations) {
      if (bankOutlet(s)) {
        const row = rows[ne + 4 * s.id + 3];
        for (const [side, weight] of [['lower', .5], ['upper', -.5]])
          for (const [col, d] of j.boundaryOutletTangency(s.body, side).derivatives) chain(row, col, weight * d);
        continue;
      }
      const row = rows[ne + 4 * s.id + 3];
      add(row, ne + 4 * s.id + 3, edgeMatching === 'pressure'
        ? streamtubeEdgePressure(x[ne + 4 * s.id + 3], euler.conditions).derivative : 1);
      for (const side of s.kind === 'surface' ? [s.side] : ['upper', 'lower']) {
        const p = edgeMatching === 'pressure' ? j.boundaryPressure(s.body, s.i, side) : j.boundaryEdgeVelocity(s.body, s.i, side, edgeInterpolation);
        const weight = s.kind === 'surface' ? 1 : .5;
        for (const [col, d] of p.derivatives) chain(row, col, -weight * d);
      }
    }
    if (sparse) return sparseMatrixFromRows(rows);
    const matrix = new Float64Array(n * n);
    rows.forEach((row, i) => { for (const [col, d] of row) {
      if (!Number.isFinite(d)) throw new Error('Nonfinite coupled streamtube Jacobian.');
      matrix[i * n + col] = d;
    } });
    return matrix;
  };
  const stationDomain = (x, id, geometry) => {
    const k = ne + 4 * id;
    return blDomain({ theta: x[k + 1], deltaStar: x[k + 2], ue: x[k + 3],
      ...(geometry ? { wakeGap: (geometry.coordinates[id].wakeGap ?? 0) / bl.scale } : {}) });
  };
  const admissibleValue = (x, { requireConvex = true, onFailure } = {}) => {
    // Diagnostics observe the existing first failing check. They never
    // evaluate another residual or let an observer alter domain acceptance.
    const reject = details => {
      if (typeof onFailure === 'function') { try { onFailure(details); } catch { /* observer only */ } }
      return false;
    };
    let stage = 'bl-geometry';
    try {
      if (bl.hasFiniteBase) euler.setDisplacement(bl.thicknesses(x.subarray(ne)));
      const geometry = bl.hasFiniteBase ? bl.geometry(x.subarray(0, ne)) : null;
      // Algebraically equivalent raw-Hk checks can disagree by roundoff
      // on the boundary. Accepted states must also be strictly inside the
      // exact polynomial halfspaces used to construct the next step.
      stage = 'bl-domain';
      for (const station of bl.stations) {
        const { id } = station;
        const domain = stationDomain(x, id, geometry);
        if (!(domain.shape > 0) || !(domain.enthalpy > 0)) {
          if (typeof onFailure !== 'function') return false;
          const k = ne + 4 * id, theta = x[k + 1], deltaStar = x[k + 2], ue = x[k + 3],
            wakeGap = geometry ? (geometry.coordinates[id].wakeGap ?? 0) / bl.scale : 0;
          const { mach, gamma } = bl.kernel.parameters, machFactor = mach ** 2 / (1 + .5 * (gamma - 1) * mach ** 2);
          const edgeMachSquared = domain.enthalpy > 0 ? machFactor * ue ** 2 / domain.enthalpy : null;
          const rawHk = theta > 0 && edgeMachSquared !== null
            ? ((deltaStar - wakeGap) / theta - .29 * edgeMachSquared) / (1 + .113 * edgeMachSquared) : null;
          return reject({ kind: 'bl-domain', stage, station: { ...station },
            failedConstraints: [...(!(domain.shape > 0) ? ['kinematic-shape'] : []), ...(!(domain.enthalpy > 0) ? ['static-enthalpy'] : [])],
            shape: domain.shape, enthalpy: domain.enthalpy, threshold: 0,
            packed: { aux: x[k], theta, deltaStar, ue, wakeGap },
            physicalLengthPerPackedThickness: bl.scale * euler.conditions.lengthScale,
            edgeMachSquared, rawHk, interpretation: 'Exact polynomial halfspaces control acceptance; rawHk is diagnostic only.' });
        }
      }
      stage = 'evaluation';
      const value = evaluate(x);
      // Explicit research drivers may use the Euler positive-simple domain
      // during iteration. Final coupled acceptance still requires convexity.
      if (!requireConvex) return value;
      stage = 'final-grid';
      const quality = streamtubeMeshQuality({ system: euler, nodes: value.outer.nodes });
      if (quality.valid) return value;
      if (typeof onFailure !== 'function') return false;
      const mesh = streamtubeMeshSnapshot({ system: euler, nodes: value.outer.nodes });
      const cell = mesh.quality.invalidCells[0], { nx, tubes } = euler.layout;
      let local = cell, group = 0;
      while (group < tubes.length - 1 && local >= nx * tubes[group]) local -= nx * tubes[group++];
      return reject({ kind: 'grid-convexity', stage, requireConvex,
        invalidCellCount: mesh.quality.invalidCells.length,
        firstCell: { cell, group, interval: Math.floor(local / tubes[group]), tube: local % tubes[group],
          points: mesh.cells[cell].map(id => ({ ...mesh.vertices[id] })) },
        minArea: mesh.quality.minArea, minCornerSine: mesh.quality.minCornerSine,
        minimumCornerSine: STREAMTUBE_MINIMUM_CORNER_SINE });
    } catch (error) {
      return typeof onFailure === 'function' ? reject({ kind: 'evaluation-error', stage,
        message: error?.message ?? String(error), ...(error?.code === undefined ? {} : { code: error.code }),
        ...(error?.diagnostics === undefined ? {} : { diagnostics: structuredClone(error.diagnostics) }),
        interpretation: 'Original exception from the existing evaluation; no unreported station or cell is inferred.' }) : false;
    }
  };
  // Return the freshly checked evaluation; never cache across states or phases.
  const admissible = (x, options) => !!admissibleValue(x, options);
  const rebase = x => {
    evaluate(x); const next = x.slice(); next.set(euler.rebase(x.subarray(0, ne))); return next;
  };
  const constraintsAt = (x, derivatives = true) => {
    euler.setDisplacement(bl.thicknesses(x.subarray(ne)));
    const corners = streamtubeCornerConstraints(euler, x.subarray(0, ne), {
      displacementMap: bl.thicknessMap, minimumSine: STREAMTUBE_MINIMUM_CORNER_SINE, derivatives });
    // Fraction-to-boundary applies to the current positive quantities. It
    // introduces no absolute cell-size floor or new governing equation. The
    // corner margin matches the existing published-mesh resolution check.
    const constraints = corners.map(({ value, gradient }) => ({ value, gradient, lower: -.9 * value }));
    const geometry = bl.hasFiniteBase ? bl.geometry(x.subarray(0, ne), derivatives) : null;
    for (const { id, regime } of bl.stations) {
      const k = ne + 4 * id;
      const domain = stationDomain(x, id, geometry);
      const makeGradient = entries => new Map(entries.map((v, i) => [k + 1 + i, v]).filter(([, v]) => v !== 0));
      const shapeGradient = makeGradient(domain.shapeGradient);
      if (geometry) for (const [col, derivative] of geometry.coordinates[id].wakeGapDerivatives ?? [])
        add(shapeGradient, col, -domain.enthalpy * derivative / bl.scale);
      for (const [kind, gradient, value] of [['momentum-thickness', new Map([[k + 1, 1]]), x[k + 1]],
        ['kinematic-shape', shapeGradient, domain.shape],
        ['static-enthalpy', makeGradient(domain.enthalpyGradient), domain.enthalpy],
        ['edge-speed', new Map([[k + 3, 1]]), x[k + 3]],
        ...(['leading-transition', 'transition', 'turbulent', 'trailing-edge', 'wake'].includes(regime) ? [['positive-shear', new Map([[k, 1]]), x[k]]] : [])])
        constraints.push({ kind, station: id, value, gradient, lower: -.9 * value });
    }
    return constraints;
  };
  const stepConstraints = x => constraintsAt(x);
  const constraintValues = x => constraintsAt(x, false).map(c => c.value);
  return { n, ne, initial, euler, bl, evaluate, residual: x => evaluate(x).residual, jacobian, admissible, admissibleValue, rebase, stepConstraints, constraintValues,
    initialization: { suppliedBL: initialBL !== undefined, boundaryLayer: structuredClone(bl.initialization) },
    conditions: { reynolds, ncrit, mach: euler.conditions.mach, trips: bl.trips, referenceChord: euler.conditions.lengthScale, edgeMatching,
      wakeOutlet: euler.conditions.wakeOutlet, ...(transitionMode === 'automatic' ? { transitionMode } : {}),
      ...(historicalTransonic ? { blThermodynamics } : {}), ...hkPolicy, ...geometryPolicy } };
}

export function solveCoupledStreamtubeBody(system, { initial = system.initial, maxIterations = 40, tolerance = 1e-8,
  stepMethod = 'newton', initialTrustRadius = 1, transitionEvents = true, projectedSteps = stepMethod === 'dogleg',
  secondOrderSteps = projectedSteps, onIteration, onMesh } = {}) {
  if (!Number.isInteger(maxIterations) || maxIterations < 0 || !Number.isFinite(tolerance) || tolerance <= 0 || !system.admissible(initial))
    throw new Error('Invalid coupled streamtube solve controls or initial state.');
  if (!['newton', 'dogleg'].includes(stepMethod) || !(initialTrustRadius > 0) || initialTrustRadius > 1e6 || !Number.isFinite(initialTrustRadius))
    throw new Error('Invalid coupled streamtube step method or trust radius.');
  if (typeof transitionEvents !== 'boolean')
    throw new Error('Invalid coupled material-trip event control.');
  if (typeof projectedSteps !== 'boolean' || (projectedSteps && stepMethod !== 'dogleg'))
    throw new Error('Projected steps require the coupled dogleg controller.');
  if (typeof secondOrderSteps !== 'boolean' || (secondOrderSteps && !projectedSteps))
    throw new Error('Second-order steps require projected coupled steps.');
  let x = Float64Array.from(initial), value = system.evaluate(x), reason = 'iteration limit', trustRadius = initialTrustRadius;
  const history = [], linearDiagnostics = { solves: 0, maxRelativeResidual: 0, orderingFallbacks: 0 };
  let preferredOrdering = 'amd';
  const report = (iteration, step, details = {}) => {
    const h = { ...details, iteration, step, residual: maximum(value.residual), ...value.families };
    history.push(h); onIteration?.(h);
  };
  report(0, 0);
  const linearSolve = (a, b) => {
    const r = solveSparseDirect(a, b, { tolerance: 1e-10, preferredOrdering });
    preferredOrdering = r.ordering ?? preferredOrdering;
    if (r.attempts?.length > 1) linearDiagnostics.orderingFallbacks++;
    linearDiagnostics.solves++; linearDiagnostics.maxRelativeResidual = Math.max(linearDiagnostics.maxRelativeResidual, r.relativeResidual); return r.x;
  };
  const activeSet = transitionEvents ? coupledStreamtubeTripEvents(system) : undefined;
  for (let k = 0; k < maxIterations && maximum(value.residual) > tolerance; k++) {
    let next, entry;
    if (stepMethod === 'dogleg') {
      let matrix, newtonDirection;
      try { matrix = system.jacobian(x); }
      catch (error) { reason = `Coupled Jacobian assembly failed: ${error.message}`; break; }
      try { newtonDirection = linearSolve(matrix, value.residual.map(v => -v)); }
      catch (error) { reason = error.message; break; }
      const step = takeDoglegStep({ initial: x, currentResidual: value.residual, matrix, newtonDirection,
        residual: system.residual, admissible: system.admissible, radius: trustRadius, activeSet,
        linearizedConstraints: projectedSteps ? system.stepConstraints : undefined,
        constraintValues: secondOrderSteps ? system.constraintValues : undefined });
      trustRadius = step.radius;
      if (!step.accepted) { reason = step.reason; break; }
      next = step.x;
      entry = { step: step.scaledStepNorm / step.scaledNewtonNorm, stepKind: step.kind, trustRadius,
        trialRadius: step.trialRadius, reductionRatio: step.ratio, actualReduction: step.actualReduction,
        predictedReduction: step.predictedReduction, linearizedResidualNorm: step.linearizedResidualNorm, trials: step.trials.length,
        ...(step.projection ? { projection: step.projection } : {}),
        ...(step.doglegSegment ? { doglegSegment: step.doglegSegment } : {}),
        ...(step.correction ? { correction: step.correction } : {}),
        ...(step.activeChange ? { activeChange: true, transitionChanges: step.changes, meritComparable: false } : {}) };
    } else {
      let step;
      try { step = solveNewton({ initial: x, maxIterations: 1, tolerance, residual: system.residual, jacobian: system.jacobian, admissible: system.admissible, linearSolve, activeSet }); }
      catch (error) { reason = `Coupled Newton step failed: ${error.message}`; break; }
      if (step.history.length < 2) { reason = step.reason; break; }
      next = step.x; entry = step.history.at(-1);
    }
    x = system.rebase(next); value = system.evaluate(x); report(k + 1, entry.step, entry);
    onMesh?.(streamtubeMeshSnapshot({ system: system.euler, nodes: value.outer.nodes, iteration: history.at(-1) }));
  }
  return coupledStreamtubeResult(system, x, { tolerance, reason, history, linearDiagnostics, stepMethod,
    transitionEvents, projectedSteps, secondOrderSteps });
}

// Shared auxiliary-state transfer for simultaneous update drivers. Fixed
// material trips retain the original event controls. XFOIL UPDATE already
// limits global BL increments; native SETBL/MRCHDU then scans the full surface
// for natural onset or its terminal fallback, without an old-index +/-1 or
// new-interval fraction cap. Disappearance of natural onset does not move a trip.
export function coupledStreamtubeTripEvents(system, { blUpdate = 'giles' } = {}) {
  if (!['giles', 'xfoil'].includes(blUpdate)) throw new Error('Unknown BL update policy for transition events.');
  const nativeAutomatic = blUpdate === 'xfoil' && system.bl.transitionMode === 'automatic';
  const nativeSurfaceSelection = target => nativeAutomatic && (target.kind === 'natural' || target.kind === 'trailing-edge');
  return {
    snapshot: system.bl.snapshotActive, restore: system.bl.restoreActive,
    prepare: (candidate, current) => {
      const targets = system.bl.activeTargets(candidate.subarray(0, system.ne), candidate.subarray(system.ne)), changes = targets.filter(t => t.from !== t.to);
      // The native selector marches N from the stagnation boundary, while N
      // is also a simultaneous unknown. Away from a root the two can differ:
      // a packed laminar N can cross Ncrit without moving the selected index,
      // or stay below Ncrit while making the selected mixed interval unusable.
      // Reconcile only that surface with the same native amplification law.
      // This is an auxiliary-state event, not a relaxed laminar-domain guard.
      // XFOIL SETBL/MRCHDU likewise establishes N and transition together;
      // no full MRCHDU thickness/velocity sweep is introduced here.
      const reconcileAmplificationSurfaces = system.bl.transitionMode === 'automatic'
        ? system.bl.surfaces.flatMap((s, k) => targets[k].from === targets[k].to
          && (targets[k].amplificationReconciliation
            || s.ids.slice(0, s.transition).some(id => candidate[system.ne + 4 * id] >= system.bl.kernel.parameters.ncrit)) ? [k] : []) : [];
      if (!changes.length && !reconcileAmplificationSurfaces.length) return { changed: false };
      // Retain the workbench's bounded migration for material trips and
      // legacy drivers. Natural onset and its fixed terminal fallback under
      // XFOIL UPDATE are reselected over the whole candidate profile. A TE
      // target has fraction 1; applying the migration cap would prohibit it.
      // TRCHEK's separate inner N-root limiter,
      // global variable limits and subsequent physical/grid gates stay active.
      const blockedTargets = changes.filter(t => !nativeSurfaceSelection(t)
        && (Math.abs(t.to - t.from) > 1 || (t.to > t.from ? t.fraction : 1 - t.fraction) > .25));
      if (blockedTargets.length) throw Object.assign(new Error('Material-trip event step crosses too much of the new interval.'), {
        code: 'COUPLED_TRANSITION_MIGRATION_LIMIT',
        diagnostics: { stage: 'material-trip transfer', blUpdate, transitionMode: system.bl.transitionMode,
          targets: blockedTargets.map(target => ({ ...target })) },
      });
      // Natural/terminal onset and amplification reconciliation change only auxiliary
      // unknowns. Use the already-global XFOIL speed bound for these events,
      // rather than imposing a smaller bound only when conversion is needed.
      // Classify PRE-transfer targets: reconciliation records lose their
      // underlying natural/forced kind. Other event policies remain unchanged.
      const nativeSpeedLimit = nativeAutomatic && changes.every(nativeSurfaceSelection)
        && reconcileAmplificationSurfaces.every(k => nativeSurfaceSelection(targets[k]));
      for (let i = system.ne; i < candidate.length; i++) {
        const k = (i - system.ne) % 4;
        if (k === 1 || k === 2) {
          const minimum = .5 * current[i];
          // XFOIL UPDATE already limits the relative thickness decrease to
          // -.5. Forming old + step*direction can round just below old/2;
          // this repeated event check must not halve that limited step again.
          // Scale to these operands, not an absolute unit thickness, and
          // retain the legacy Giles event gate exactly.
          const roundoff = blUpdate === 'xfoil' && Number.isFinite(candidate[i]) && Number.isFinite(current[i])
            ? 8 * Number.EPSILON * Math.max(Math.abs(candidate[i]), Math.abs(current[i])) : 0;
          if (minimum - candidate[i] > roundoff) throw Object.assign(
            new Error('Material-trip event reduces a thickness too far.'), {
              code: 'COUPLED_TRANSITION_THICKNESS_LIMIT',
              diagnostics: { stage: 'material-trip transfer', blUpdate, station: Math.floor((i - system.ne) / 4),
                variable: k === 1 ? 'theta' : 'delta-star', current: current[i], candidate: candidate[i],
                minimum, roundoff, normalizedIncrement: (candidate[i] - current[i]) / current[i] },
            });
        }
        if (k === 3) {
          const limit = nativeSpeedLimit ? XFOIL_BL_EDGE_SPEED_LIMIT : .2;
          // Adding the limited increment and subtracting the old speed can
          // exceed the exact .375 by roundoff. Preserve the legacy .2 test.
          const roundoff = nativeSpeedLimit
            ? 8 * Number.EPSILON * Math.max(1, Math.abs(candidate[i]), Math.abs(current[i])) : 0;
          if (Math.abs(candidate[i] - current[i]) > limit + roundoff)
            throw new Error('Material-trip event changes edge velocity too far.');
        }
      }
      return system.bl.updateActive(candidate.subarray(system.ne), candidate.subarray(0, system.ne), { reconcileAmplificationSurfaces });
    },
  };
}

// Report a physically admissible research endpoint even when its final
// convex-grid gate fails. This function never promotes that endpoint.
export function coupledStreamtubeResult(system, x, { tolerance = 1e-8, reason = 'iteration limit', ...details } = {}) {
  const value = system.evaluate(x);
  const flow = solveStreamtubeBody(system.euler, { initial: x.subarray(0, system.ne), maxIterations: 0, tolerance });
  const mesh = streamtubeMeshSnapshot({ system: system.euler, nodes: flow.nodes });
  const residualConverged = maximum(value.residual) <= tolerance;
  const changesConverged = coupledConvergenceSatisfied({ families: value.families, convergence: details.convergence }, tolerance);
  const converged = (residualConverged || changesConverged) && mesh.quality.valid;
  mesh.initialization.flowSolved = converged;
  return { ...details, x, converged, residualConverged, reason: converged ? (residualConverged ? 'residual' : 'solution changes') : mesh.quality.valid ? reason : 'Invalid final displacement grid',
    // Preserve the numerical stop when final convexity overrides its display reason.
    ...(!mesh.quality.valid ? { solverStopReason: residualConverged ? 'residual' : reason } : {}),
    residual: value.residual, families: value.families, flow, mesh,
    ...(system.conditions.blThermodynamics === 'historical-common-isentrope'
      ? { edgeThermodynamics: streamtubeBLEdgeThermodynamics({ flow: value.outer, euler: system.euler, bl: system.bl, states: value.layers.states }) } : {}),
    boundaryLayer: { stations: system.bl.stations.map((s, id) => ({ ...s, ...value.layers.states[id] })), surfaces: system.bl.surfaces, wakes: system.bl.wakes,
      ...(value.layers.transitions ? { transitions: value.layers.transitions, transitionState: system.bl.snapshotActive() } : {}) },
    initialization: system.initialization,
    conditions: system.conditions, cl: null, cd: null, cm: null,
    status: converged ? (residualConverged ? 'research-coupled-equations-converged' : 'research-coupled-changes-converged') : 'unconverged',
    limitations: system.conditions.blThermodynamics === 'historical-common-isentrope'
      ? `Research ${system.bl.transitionMode} simultaneous hybrid Euler/BL/displacement/wake equations with historical common-isentrope BL thermodynamics. Physical Euler entropy is retained separately; one merged BL wake per body approximates unequal-entropy wake banks. Transonic coupling, shock/BL interaction and physical/refinement force acceptance remain unvalidated.`
      : `Research ${system.bl.transitionMode} simultaneous Euler/BL/displacement/wake equations; ${system.conditions.edgeMatching} edge matching. Smooth subcritical BL thermodynamics and one merged wake per body. Natural-transition movement, shocks, general grid robustness and physical/refinement force acceptance remain unvalidated.` };
}
