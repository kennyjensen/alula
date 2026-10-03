// SPDX-License-Identifier: GPL-2.0-or-later
import { initializeCoupledStreamtubeBody } from './streamtube-coupled-initializer.js';
import { solveCoupledStreamtubeIses } from './streamtube-coupled-ises.js';

// One alternative initial guess on the SAME grid and operating point.
// Panel speeds omit shocks; a failed panel-seeded solve can therefore benefit
// from a fresh BL march using the already computed Euler precursor. A failed
// alternative never replaces the retained complete coupled state.
export function retryCoupledEulerSeed({ result, input, options, initialEuler, initialization,
  maxIterations, tolerance, onIteration, onMesh, onCheckpoint }) {
  if (result.converged || !(input.mach >= .7) || options.transitionMode !== 'automatic'
    || options.edgeMatching !== 'section-velocity' || !result.checkpoint || !result.mesh.quality.valid
    || initialization.panelEdgeGuess?.accepted !== true || maxIterations < 1) return null;
  const diagnostics = { method: 'euler-edge-seed-retry', accepted: false, equationsChanged: false,
    operatingConditionsChanged: false, precursorReused: true, originalFamilies: { ...result.families },
    originalIterations: result.history.length - 1 };
  let prepared;
  try { prepared = initializeCoupledStreamtubeBody(input, { ...options, initialEuler }); }
  catch (error) { return { diagnostics: { ...diagnostics, reason: error.message, iterations: 0 } }; }
  const s = prepared.system, v = s.evaluate(s.initial), controls = result.checkpoint.continuation;
  const next = solveCoupledStreamtubeIses(input, { ...options, transitionState: s.bl.snapshotActive(),
    initialEuler: { x: s.initial.slice(0, s.ne), nodes: v.outer.nodes, undisplacedNodes: v.outer.undisplacedNodes },
    initialBL: s.initial.slice(s.ne), maxIterations, tolerance, convergence: controls.convergence,
    iterationGeometry: controls.iterationGeometry, stepAcceptance: controls.stepAcceptance,
    stagnationLimiter: controls.stagnationLimiter, blUpdate: controls.blUpdate,
    projectionGeometry: controls.projectionGeometry, linearOrdering: controls.linearOrdering,
    shearCoordinate: controls.shearCoordinate, hkProjectionRecovery: controls.hkProjectionRecovery ?? false,
    directionRecovery: controls.directionRecovery ?? false, wakeGridInitialization: true,
    eliminateAmplification: controls.eliminateAmplification ?? false,
    onIteration, onMesh, onCheckpoint });
  Object.assign(diagnostics, { accepted: next.converged, reason: next.reason, families: { ...next.families },
    iterations: next.history.length - 1, history: next.history, linearDiagnostics: next.linearDiagnostics, initialization: prepared.initialization });
  if (!next.converged) return { diagnostics };
  return { result: next, prepared, diagnostics };
}
