// SPDX-License-Identifier: GPL-2.0-or-later
// Construct an admissible initial guess for the simultaneous system. A local
// BL profile need not be a geometrically admissible displacement of its
// inviscid precursor. Backtrack ALL initial thicknesses together; Reynolds
// number, closure equations and final displacement coupling remain unchanged.
import { createCoupledStreamtubeBody } from './streamtube-coupled.js';
import { extendStreamtubeDisplacement } from './streamtube-displacement.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';
import { initializeStreamtubeWakeCorrespondence } from './streamtube-wake-correspondence.js';
import { relaxStreamtubeInitialGrid } from './streamtube-elliptic-initializer.js';
import { solveViscousAssembly } from '../viscous/result.js';

// A converged incompressible panel/BL profile supplies a cold surface/wake
// guess. The subsequent Euler/BL solve retains the requested
// Mach, Reynolds number, transition model, geometry and grid resolution.
export function panelCoupledEdgeGuess(input, bodies, lengthScale) {
  if (!(lengthScale > 0) || !Number.isFinite(lengthScale)) throw new Error('Panel edge guess needs the Euler length scale.');
  const result = solveViscousAssembly({ ...input, mach: 0, maxIterations: 60,
    elements: input.elements.map((element, i) => ({ ...element,
      ...(input.materialTrips?.[i] ? { trips: input.materialTrips[i] } : {}) })) });
  if (result.status !== 'solved') return { accepted: false, reason: result.diagnostics?.reason };
  const profiles = result.boundaryLayer;
  const interpolate = (rows, distance, key = 'ue') => {
    if (!rows.length) throw new Error(`Panel BL profile has no ${key} samples.`);
    const i = rows.findIndex(row => row.s >= distance);
    if (i < 0) return rows.at(-1)[key];
    if (i === 0) return rows[0][key];
    const a = rows[i - 1], b = rows[i], weight = (distance - a.s) / (b.s - a.s);
    return a[key] * (1 - weight) + b[key] * weight;
  };
  return { accepted: true, residual: result.diagnostics.equationResidual,
    initialEdgeVelocity(station, bl, geometry) {
      const element = bodies[station.body].element;
      if (station.kind === 'surface') {
        const surface = bl.surfaces.find(s => s.body === station.body && s.side === station.side);
        const rows = profiles.surfaces.find(s => s.element === element && s.side === station.side).stations;
        const end = geometry.coordinates[surface.ids.at(-1)].s;
        return interpolate([{ s: 0, ue: 0 }, ...rows], geometry.coordinates[station.id].s / end * rows.at(-1).s);
      }
      const rows = profiles.wakes[element].stations, wake = bl.wakes[station.body];
      const distance = (geometry.coordinates[station.id].s - geometry.coordinates[wake.ids[0]].s) * lengthScale;
      return interpolate(rows, distance + rows[0].s);
    },
    initialBoundaryLayer(system, { thicknessFactor = 1 } = {}) {
      if (!(thicknessFactor > 0 && thicknessFactor <= 1)) throw new Error('Invalid retained BL thickness factor.');
      const { bl, ne, initial } = system, eulerState = initial.slice(0, ne);
      const geometry = bl.geometry(eulerState), initialBL = initial.slice(ne), transitionState = [];
      const savedTransition = bl.snapshotActive();
      for (const surface of bl.surfaces) {
        const element = bodies[surface.body].element;
        const rows = profiles.surfaces.find(s => s.element === element && s.side === surface.side).stations;
        const laminar = rows.filter(r => r.regime === 'laminar'), turbulent = rows.filter(r => r.regime === 'turbulent');
        const end = geometry.coordinates[surface.ids.at(-1)].s;
        const distance = id => geometry.coordinates[id].s / end * rows.at(-1).s;
        const firstTurbulent = turbulent.length ? surface.ids.findIndex(id => distance(id) >= turbulent[0].s) : -1;
        const transition = firstTurbulent < 0 ? surface.ids.length : firstTurbulent;
        transitionState.push(transition);
        surface.ids.forEach((id, j) => {
          const s = distance(id);
          // Keep the local similarity seed ahead of the first panel station;
          // clamping its thickness to a downstream sample breaks the LE row.
          if (s < rows[0].s) {
            // A geometrically admitted fallback may have thinned this local
            // seed. Restore its physical thickness when transferring the
            // full panel profile; otherwise the leading rows start off-root.
            initialBL[4 * id + 1] /= thicknessFactor;
            initialBL[4 * id + 2] /= thicknessFactor;
            return;
          }
          initialBL.set([
            j < transition ? interpolate(laminar, s, 'amplification') : interpolate(turbulent, s, 'ctau'),
            interpolate(rows, s, 'theta') / (lengthScale * bl.scale),
            interpolate(rows, s, 'deltaStar') / (lengthScale * bl.scale), interpolate(rows, s),
          ], 4 * id);
        });
      }
      for (const wake of bl.wakes) {
        const rows = profiles.wakes[bodies[wake.body].element].stations;
        for (const id of wake.ids) {
          const coordinate = geometry.coordinates[id];
          const s = (coordinate.s - geometry.coordinates[wake.ids[0]].s) * lengthScale + rows[0].s;
          initialBL.set([interpolate(rows, s, 'ctau'), interpolate(rows, s, 'theta') / (lengthScale * bl.scale),
            // Transfer viscous thickness only; the Euler wake owns its
            // prescribed finite-base gap and physical wake coordinates.
            (interpolate(rows, s, 'viscousDeltaStar') / lengthScale + (coordinate.wakeGap ?? 0)) / bl.scale,
            interpolate(rows, s)], 4 * id);
        }
      }
      try {
        bl.restoreActive(transitionState);
        // Interpolated N does not satisfy amplification on the new station
        // spacing. Reconcile it before testing transition or thinning a seed.
        bl.updateActive(initialBL, eulerState, { reinitializeAmplification: true });
        // A coarse target mesh can miss the laminar amplification peak.
        // Reintegrating N through the interpolated turbulent tail then sends
        // transition to the TE and converts that tail back to laminar, while
        // retaining its turbulent thicknesses and the transferred wake. Such
        // a profile is admissible but is not a useful transferred solution.
        // Retain the locally marched profile for that entire element,
        // including both surfaces and its merged wake. A transfer failure
        // on one element must not discard usable profiles on the others.
        const lostTransitions = bl.surfaces.filter((surface, k) =>
          transitionState[k] < surface.ids.length - 1 && surface.transition === surface.ids.length - 1);
        const fallbackBodies = new Set(lostTransitions.map(surface => surface.body));
        if (fallbackBodies.size === bodies.length) {
          const first = lostTransitions[0];
          throw new Error(`Panel profile transfer lost resolved transition on body ${first.body + 1}, ${first.side}; retained the locally marched panel-speed seed.`);
        }
        if (fallbackBodies.size) {
          const phases = bl.snapshotActive();
          bl.surfaces.forEach((surface, k) => {
            if (fallbackBodies.has(surface.body)) phases[k] = savedTransition[k];
          });
          for (const station of bl.stations) if (fallbackBodies.has(station.body)) {
            const k = 4 * station.id;
            initialBL.set(initial.subarray(ne + k, ne + k + 4), k);
            if (thicknessFactor !== 1) {
              // Undo only viscous thinning; the finite trailing-edge gap
              // belongs to the geometry and must never be scaled.
              const gap = (geometry.coordinates[station.id].wakeGap ?? 0) / bl.scale;
              initialBL[k + 1] /= thicknessFactor;
              initialBL[k + 2] = (initialBL[k + 2] - gap) / thicknessFactor + gap;
            }
          }
          bl.restoreActive(phases);
        }
        return { initialBL, transitionState: bl.snapshotActive(),
          profileTransfer: { mappedBodies: bodies.flatMap((_, b) => fallbackBodies.has(b) ? [] : [b]),
            fallbackBodies: [...fallbackBodies],
            rejectedSurfaces: lostTransitions.map(({ body, side }) => ({ body, side, reason: 'unresolved transferred transition' })) } };
      } finally { bl.restoreActive(savedTransition); }
    } };
}

// A displaced trailing edge can pinch the adjacent passage interiors even
// when all cells remain convex. Relax that part of a fresh transonic seed,
// after its initial tangential redistribution. Whole-grid relaxation can
// destroy a good Euler seed elsewhere, so leave the inlet and forward body
// alone. The cosine window overlaps the aft body to avoid a kink at the TE.
// Its width (20% of body arc) and half correction are initialization controls,
// not changes to the Euler/BL equations or a claim of exact MSES parity.
export function relaxCoupledWakeSeed(input, options, system, state, value) {
  if (!(value.outer.diagnostics.maxMach > 1) || !system.euler.layout.independentWakeBanks
    || value.families.euler <= Math.max(value.families.boundaryLayer, value.families.edgeMatching)) return null;
  const squared = residual => residual.reduce((sum, v) => sum + v * v, 0);
  const beforeSquaredNorm = squared(value.residual), nodes = value.outer.nodes;
  const report = { method: 'elliptic-wake-seed', accepted: false, beforeSquaredNorm,
    amplitude: .5, arcWidthFraction: .2, boundaryNodesFixed: true, equationsChanged: false };
  try {
    const relaxed = relaxStreamtubeInitialGrid({ system: system.euler, initial: state.subarray(0, system.ne), nodes },
      { seed: 'supplied', maxSweeps: 1000, tolerance: 1e-8, omega: 1, requireConvex: true });
    report.regions = relaxed.regions.map(r => ({ converged: r.converged, reason: r.reason,
      iterations: r.history.at(-1)?.iteration, residual: r.history.at(-1)?.residual }));
    if (!relaxed.converged) return { report: { ...report, reason: 'Elliptic seed correction did not converge.' } };
    const masks = system.euler.layout.bodies.map((body, b) => {
      const arc = [0], bank = nodes[b + 1].map(row => row[0]);
      for (let i = 1; i < bank.length; i++)
        arc.push(arc.at(-1) + Math.hypot(bank[i].x - bank[i - 1].x, bank[i].y - bank[i - 1].y));
      const width = .2 * (arc[body.trailingIndex] - arc[body.leadingIndex]);
      if (!(width > 0)) throw new Error('Wake smoothing needs a positive body outline length.');
      return arc.map(s => {
        const distance = (s - arc[body.trailingIndex]) / width;
        return distance <= -1 ? 0 : distance >= 1 ? 1 : (1 + Math.sin(Math.PI * distance / 2)) / 2;
      });
    });
    const target = nodes.map((group, g) => group.map((row, i) => {
      const adjacent = [masks[g - 1], masks[g]].filter(Boolean);
      const weight = adjacent.length === 1 ? adjacent[0][i]
        : 1 - adjacent.reduce((remaining, mask) => remaining * (1 - mask[i]), 1);
      return row.map((p, j) => ({
        x: p.x + .5 * weight * (relaxed.nodes[g][i][j].x - p.x),
        y: p.y + .5 * weight * (relaxed.nodes[g][i][j].y - p.y),
      }));
    }));
    const bl = state.slice(system.ne);
    const candidate = createCoupledStreamtubeBody(input, { ...options, transitionState: system.bl.snapshotActive(),
      initialBL: bl, initialEuler: { x: state.slice(0, system.ne), nodes,
        undisplacedNodes: value.outer.undisplacedNodes } });
    const next = candidate.initial.slice();
    candidate.euler.setDisplacement(candidate.bl.thicknesses(bl));
    next.set(candidate.euler.adoptGeometry(next.subarray(0, candidate.ne), target));
    const nextValue = candidate.admissibleValue(next, { requireConvex: true,
      onFailure: failure => { report.admissibility = failure; } });
    if (!nextValue) return { report: { ...report, reason: 'Corrected seed is not physically admissible.' } };
    report.afterSquaredNorm = squared(nextValue.residual);
    report.afterFamilies = nextValue.families;
    if (!(report.afterSquaredNorm < beforeSquaredNorm))
      return { report: { ...report, reason: 'Corrected seed does not reduce the coupled residual.' } };
    report.accepted = true;
    return { system: candidate, state: next, value: nextValue, report };
  } catch (error) {
    return { report: { ...report, reason: error.message, ...(error.code ? { code: error.code } : {}) } };
  }
}

// The coordinate initializer only changes downstream wake banks and their
// passage interiors. Do not invoke it for a bad cell outside that support.
// This is an eligibility check, never a substitute for the final grid test.
export function coldStreamtubeWakeCorrespondenceEligible({ layout, quality } = {}) {
  const { nx, elements, tubes, bodies, independentWakeBanks } = layout ?? {};
  if (independentWakeBanks !== true || !Number.isInteger(nx) || nx < 2
    || !Number.isInteger(elements) || elements < 1 || !Array.isArray(tubes)
    || tubes.length !== elements + 1 || tubes.some(n => !Number.isInteger(n) || n < 1)
    || !Array.isArray(bodies) || bodies.length !== elements || bodies.some(b =>
      !Number.isInteger(b?.trailingIndex) || b.trailingIndex < 1 || b.trailingIndex >= nx)
    || quality?.valid !== false || !Array.isArray(quality.invalidCells) || !quality.invalidCells.length) return false;
  const starts = [0]; for (const n of tubes) starts.push(starts.at(-1) + nx * n);
  return quality.invalidCells.every(id => {
    if (!Number.isInteger(id) || id < 0 || id >= starts.at(-1)) return false;
    const g = starts.findIndex((start, j) => j < tubes.length && id >= start && id < starts[j + 1]);
    const i = Math.floor((id - starts[g]) / tubes[g]);
    return bodies.some((body, b) => (g === b || g === b + 1) && i >= body.trailingIndex);
  });
}

// Giles ISET/TDCALC pp.183–184 starts its wake with constant momentum
// thickness and index-linear H approaching 1.1. Use the solved surface
// guess's merged TE values here, rather than ISET's uniform surface H=2.4.
// Only an initial guess: every wake equation remains in the Newton solve.
function isetWakeGuess(bl, original) {
  const result = original.slice();
  for (const wake of bl.wakes) {
    const first = 4 * wake.ids[0], theta = original[first + 1];
    const h = (original[first + 2] - (bl.initialWakeGaps?.[wake.ids[0]] ?? 0) / bl.scale) / theta;
    const count = wake.ids.length - 1;
    wake.ids.slice(1).forEach((id, k) => {
      result[4 * id + 1] = theta;
      result[4 * id + 2] = theta * (h + (1.1 - h) * (k + 1) / count) + (bl.initialWakeGaps?.[id] ?? 0) / bl.scale;
    });
  }
  return result;
}

// The dead-air contribution is prescribed geometry, not a viscous unknown
// that can be thinned. Native BL variables contain TOTAL displacement.
export function scaleStreamtubeBLThicknesses(bl, original, factor) {
  if (!Number.isFinite(factor) || factor <= 0 || factor > 1) throw new Error('Invalid BL thickness scale.');
  if (!bl.hasFiniteBase) return original.map((v, k) => k % 4 === 1 || k % 4 === 2 ? factor * v : v);
  return original.map((v, k) => {
    if (k % 4 === 1) return factor * v;
    if (k % 4 !== 2 || factor === 1) return v;
    const gap = bl.initialWakeGaps[Math.floor(k / 4)] / bl.scale;
    return factor * (v - gap) + gap;
  });
}

export function initializeCoupledStreamtubeBody(input, options, { maximumBacktracks = 10, initialThicknessFactor = 1, onAttempt } = {}) {
  if (!options?.initialEuler) throw new Error('Supply an Euler state and its physical grid for coupled initialization.');
  if (!Number.isInteger(maximumBacktracks) || maximumBacktracks < 0 || maximumBacktracks > 30)
    throw new Error('Invalid coupled initialization backtracking limit.');
  if (!Number.isFinite(initialThicknessFactor) || initialThicknessFactor <= 0 || initialThicknessFactor > 1)
    throw new Error('Invalid initial coupled thickness factor.');
  // MRCHUE checks transition during each local update and switches to
  // inverse Hk before accepting a strongly separated direct initializer.
  // The low-level constructor retains its controlled direct/inverse seed;
  // this public admissible-grid initializer uses MRCHUE for natural surfaces.
  let seeded = createCoupledStreamtubeBody(input, { ...options,
    ...(options.transitionMode === 'automatic' ? { blInitialization: 'mrchue' } : {}) });
  let surfaceRecovery;
  const warnings = seeded.initialization?.boundaryLayer;
  const incompleteSurfaceMarch = Array.isArray(warnings) && warnings.some(row => row?.method === 'mrchue'
    && Array.isArray(row.localConvergenceWarnings) && row.localConvergenceWarnings.length > 0);
  const automatic = seeded.bl.transitionMode === 'automatic';
  let initialOptions = { ...options, ...(automatic ? { transitionState: seeded.bl.snapshotActive() } : {}) };
  let originalBL = seeded.initial.slice(seeded.ne);
  const history = [];
  const guesses = Array.from({ length: maximumBacktracks + 1 }, (_, attempt) => ({ factor: initialThicknessFactor * 2 ** -attempt }));
  // Preserve a completed MRCHUE surface march's existing first guess.
  // An unresolved local surface march may leave a poor merged wake seed;
  // try the existing ISET wake guess first only for that recorded failure.
  // Supplied BL states and explicit thinner retries keep their old order.
  if (automatic && options.initialBL == null && initialThicknessFactor === 1) {
    guesses.splice(incompleteSurfaceMarch ? 0 : 1, 0, { factor: 1, wakeInitialization: 'iset-linear-shape' });
  }
  for (let index = 0; index < guesses.length; index++) {
    let guess = guesses[index];
    // Keep every admissible full-thickness MRCHUE/ISET seed. If both fail,
    // thinning an unresolved surface profile can amplify its equation error.
    // Try one resolved march before resorting to those thinner guesses.
    if (guess.factor < initialThicknessFactor && incompleteSurfaceMarch && !surfaceRecovery
      && automatic && options.initialBL == null && initialThicknessFactor === 1) {
      surfaceRecovery = { method: 'direct-inverse', trigger: 'incomplete-mrchue', accepted: false,
        warnings: warnings.filter(row => row?.method === 'mrchue'
          && Array.isArray(row.localConvergenceWarnings) && row.localConvergenceWarnings.length).map(row => ({
          id: row.id, body: row.body, messages: [...row.localConvergenceWarnings],
        })) };
      try {
        // Same operating condition and physical grid. Only the initial BL
        // guess changes; final coupling, transition and closure stay intact.
        seeded = createCoupledStreamtubeBody(input, { ...options, blInitialization: 'direct-inverse' });
        originalBL = seeded.initial.slice(seeded.ne);
        initialOptions = { ...options, transitionState: seeded.bl.snapshotActive() };
        surfaceRecovery.accepted = true;
        guess = { factor: initialThicknessFactor };
        guesses.splice(index, 0, guess);
      } catch (error) {
        surfaceRecovery.reason = error.message;
      }
    }
    const factor = guess.factor, baseBL = guess.wakeInitialization ? isetWakeGuess(seeded.bl, originalBL) : originalBL;
    const initialBL = scaleStreamtubeBLThicknesses(seeded.bl, baseBL, factor);
    const row = { attempt: history.length, thicknessFactor: factor, accepted: false,
      ...(guess.wakeInitialization ? { wakeInitialization: guess.wakeInitialization } : {}) };
    history.push(row);
    let accepted;
    try {
      // Each attempt owns its chart: failed geometry/gas checks cannot leave
      // partial offsets in a subsequent candidate or mutate the input state.
      const trial = createCoupledStreamtubeBody(input, { ...initialOptions, initialBL });
      const euler = trial.euler, state = trial.initial.slice(0, trial.ne);
      if (automatic) {
        // Thickness backtracking changes amplification growth. Convert the
        // initial guess's phases before evaluating its coupled equations.
        const event = trial.bl.updateActive(initialBL, state);
        row.transitionChanges = event.changes;
      }
      const trialOptions = { ...options, ...(automatic ? { transitionState: trial.bl.snapshotActive() } : {}) };
      euler.setDisplacement(trial.bl.thicknesses(initialBL));
      let nodes = extendStreamtubeDisplacement(euler, state);
      row.quality = streamtubeMeshSnapshot({ system: euler, nodes }).quality;
      if (!row.quality.valid && row.quality.invalidCells.some(cell => {
        let group = 0;
        while (group < euler.layout.tubes.length - 1 && cell >= euler.layout.nx * euler.layout.tubes[group])
          cell -= euler.layout.nx * euler.layout.tubes[group++];
        const i = Math.floor(cell / euler.layout.tubes[group]);
        return [group - 1, group].some(b => euler.layout.bodies[b]?.leadingIndex === i + 1);
      })) {
        row.inletDisplacement = { attempted: true, beforeQuality: row.quality };
        nodes = extendStreamtubeDisplacement(euler, state, { inletDisplacement: true });
        row.quality = streamtubeMeshSnapshot({ system: euler, nodes }).quality;
        row.inletDisplacement.afterQuality = row.quality;
      }
      // Mapped profiles and finite-base cold profiles can offset the two TE
      // banks tangentially. Pair their wake coordinates before discarding a
      // resolved full-thickness BL guess. The correction preserves the wake
      // center and prescribed normal gap; surface defects remain ineligible.
      if ((options.initialBL != null || seeded.bl.hasFiniteBase)
        && coldStreamtubeWakeCorrespondenceEligible({ layout: euler.layout, quality: row.quality })) {
        row.wakeCorrespondence = { attempted: true, accepted: false, beforeQuality: row.quality,
          initialGuessOnly: true, equationsChanged: false, thicknessChanged: false };
        const massFractions = euler.decode(state).allocation.groups.map(group => {
          const total = group.reduce((sum, tube) => sum + tube.massFlow, 0), fractions = [0];
          for (const tube of group) fractions.push(fractions.at(-1) + tube.massFlow / total);
          fractions[fractions.length - 1] = 1; return fractions;
        });
        const candidate = initializeStreamtubeWakeCorrespondence({ nodes, layout: euler.layout, massFractions });
        nodes = candidate.nodes;
        row.wakeCorrespondence.geometry = candidate.diagnostics;
        row.quality = streamtubeMeshSnapshot({ system: euler, nodes }).quality;
        row.wakeCorrespondence.afterQuality = row.quality;
      }
      if (!row.quality.valid) throw new Error('Initial BL displacement extension has invalid cells.');
      // Preserve the supplied density and global flow unknowns. A fresh
      // isentropic density inversion would discard a progressed Euler state
      // even as factor tends to zero, defeating continuous backtracking.
      const gas = euler.adoptGeometry(state, nodes), flow = euler.evaluate(gas);
      const system = createCoupledStreamtubeBody(input, { ...trialOptions, initialBL, initialEuler: { ...flow, x: gas } });
      const value = system.evaluate(system.initial), mesh = streamtubeMeshSnapshot({ system: system.euler, nodes: value.outer.nodes });
      if (!mesh.quality.valid) throw new Error('Restored coupled grid has invalid cells.');
      row.accepted = true; row.families = value.families; row.maxMach = value.outer.diagnostics.maxMach;
      if (row.wakeCorrespondence) row.wakeCorrespondence.accepted = true;
      accepted = { system, mesh, initialization: { method: surfaceRecovery?.accepted ? 'Resolved direct/inverse BL initialization'
        : guess.wakeInitialization ? 'MRCHUE surfaces and ISET-style wake guess' : 'uniform initial thickness backtracking',
        thicknessFactor: factor, history, originalBLInitialization: seeded.initialization,
        ...(surfaceRecovery ? { surfaceRecovery } : {}),
        ...(guess.wakeInitialization ? { wakeInitialization: guess.wakeInitialization, surfaceThicknessFactor: 1 } : {}),
        equationsChanged: false, flowSolved: false } };
    } catch (error) {
      row.reason = error.message;
      if (row.wakeCorrespondence) row.wakeCorrespondence.reason = error.message;
    }
    onAttempt?.(structuredClone(row));
    if (accepted) return accepted;
  }
  const error = new Error(`Coupled initialization failed after ${history.length} thickness guesses: ${history.at(-1).reason}`);
  error.initialization = { method: 'uniform initial thickness backtracking', history,
    ...(surfaceRecovery ? { surfaceRecovery } : {}), equationsChanged: false, flowSolved: false };
  throw error;
}
