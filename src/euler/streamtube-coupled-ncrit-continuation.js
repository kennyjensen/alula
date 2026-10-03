// SPDX-License-Identifier: GPL-2.0-or-later
import { observableFlow, observableBL } from './streamtube-flow-preview.js';
// Bounded initial-guess continuation. Every accepted source is a complete
// coupled root; only the final requested Ncrit can certify the user's case.
import { prepareCoupledNcritRestart } from './streamtube-coupled-ncrit-restart.js';
import { solveCoupledStreamtubeIses } from './streamtube-coupled-ises.js';
import { createCoupledStreamtubeBody } from './streamtube-coupled.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';

export function continueCoupledNcrit(source, { targetNcrit, phase, maxIterations = 40, tolerance = 1e-10, maximumStep = 1,
  normalization, startupAttempt = 1, onStage, onIteration, onMesh, onFlow, onIterationCheckpoint } = {}) {
  if (!['coarse', 'fine'].includes(phase) || source?.converged !== true || !source.checkpoint?.restart?.options
    || !Number.isFinite(targetNcrit) || !Number.isFinite(source.checkpoint.restart.options.ncrit)
    || !(targetNcrit >= source.checkpoint.restart.options.ncrit)
    || !Number.isFinite(maximumStep) || maximumStep > 1 || maximumStep < (phase === 'coarse' ? .5 : .125)
    || !Number.isInteger(maxIterations) || maxIterations < 0 || !Number.isFinite(tolerance) || tolerance <= 0)
    throw new Error('Invalid coupled Ncrit continuation source or controls.');
  let retained = source, step = maximumStep, observerFailed = false, observerError;
  // Work selection only: limiting an Ncrit increment never changes the target
  // criterion, equations or convergence tolerance. Existing callers keep 1.
  const gridProvenance = source.gridSequence ?? source.automaticRefinement ?? source.mesh?.initialization?.gridRefinement;
  const attempts = [], stage = phase === 'coarse' ? 'coupled-coarse-initialization' : 'coupled';
  const observe = fn => (...args) => { try { return fn?.(...args); }
    catch (error) { observerFailed = true; observerError = error; throw error; } };
  const emitStage = observe(onStage), emitIteration = observe(onIteration), emitMesh = observe(onMesh);
  const emitFlow = observe(onFlow), emitCheckpoint = observe(onIterationCheckpoint);
  const labels = (actualNcrit, stateConverged = false) => ({ actualNcrit, targetNcrit,
    ncritContinuation: { actualNcrit, targetNcrit, reachedTarget: actualNcrit === targetNcrit,
      stateConverged, phase } });
  // At most one coarse halving. A difficult coarse transition should be
  // resolved on the requested finer grid rather than by repeated coarse solves.
  // Fine continuation has a bounded adaptive step and retains the last root.
  const minimumStep = phase === 'coarse' ? .5 : .125;
  for (let attempt = 0; attempt < 32 && retained.checkpoint.restart.options.ncrit < targetNcrit; attempt++) {
    const actual = retained.checkpoint.restart.options.ncrit;
    const nextInteger = Math.floor(actual + 1e-12) + 1;
    const selected = Math.min(targetNcrit, actual + step, nextInteger);
    const record = { sourceNcrit: actual, selectedNcrit: selected, phase, converged: false };
    attempts.push(record);
    emitStage({ stage, startupAttempt, ...labels(selected), ncritStartup: true });
    let candidate, checkpoint, system;
    const publish = state => {
      const mesh = streamtubeMeshSnapshot(state);
      if (source.mesh?.initialization?.gridSmoothing) mesh.initialization.gridSmoothing = source.mesh.initialization.gridSmoothing;
      if (gridProvenance) mesh.initialization.gridRefinement = gridProvenance;
      emitMesh(mesh, state.iteration.iteration ? 'solving' : 'initial', stage);
      if (onFlow && checkpoint) emitFlow(structuredClone({ checkpoint,
        flow: observableFlow(state.flow),
        bl: observableBL(system.bl),
        bodies: state.system.layout.bodies, normalization, iteration: state.iteration,
        mach: checkpoint.restart.input.mach, stage, startupAttempt,
        ...labels(checkpoint.restart.options.ncrit, state.iteration.ncritContinuation?.stateConverged === true) }));
    };
    try {
      const prepared = prepareCoupledNcritRestart(selected, retained.checkpoint, { requestedNcrit: targetNcrit,
        tolerance, blPredictor: phase === 'coarse' ? 'xfoil-mrchdu' : 'preserve' });
      record.preparation = prepared.diagnostics;
      if (prepared.diagnostics.warnings.length) record.reason = 'Native BL predictor did not converge locally.';
      else {
        const f = prepared.checkpoint.restart;
        const { iterationGeometry, stepAcceptance, stagnationLimiter } = prepared.checkpoint.continuation;
        system = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
        candidate = solveCoupledStreamtubeIses(undefined, { resume: prepared.checkpoint, maxIterations, tolerance,
          iterationGeometry, stepAcceptance, stagnationLimiter,
          onIteration: h => emitIteration({ ...h, stage, startupAttempt, ...labels(selected) }),
          onCheckpoint: (value, details) => {
            checkpoint = value;
            emitCheckpoint(structuredClone(value), { ...details, stage, startupAttempt, ...labels(selected) });
          }, onMesh: publish });
        record.converged = candidate.converged; record.reason = candidate.reason;
        record.families = candidate.families; record.iterations = candidate.history.length - 1;
        if (candidate.converged) {
          retained = { ...retained, ...candidate };
          // Keep the published grid provenance across complete-state restarts.
          if (source.mesh?.initialization?.gridSmoothing) retained.mesh.initialization.gridSmoothing = source.mesh.initialization.gridSmoothing;
          if (gridProvenance) retained.mesh.initialization.gridRefinement = gridProvenance;
          continue;
        }
      }
    } catch (error) {
      if (observerFailed) throw observerError;
      record.reason = error.message; record.failureCode = error.code;
      record.diagnostics = error.diagnostics;
    }
    // A failed candidate is never a source for a subsequent change of Ncrit.
    // Restore the matching retained grid and pressure checkpoint together.
    const f = retained.checkpoint.restart;
    checkpoint = retained.checkpoint;
    system = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
    const value = system.evaluate(system.initial);
    emitCheckpoint(structuredClone(checkpoint), { stage, startupAttempt, ...labels(actual, true), retained: true,
      history: structuredClone(retained.history), retryReason: record.reason });
    publish({ system: system.euler, nodes: value.outer.nodes, flow: value.outer,
      iteration: { ...retained.history.at(-1), ...labels(actual, true) } });
    emitStage({ stage, startupAttempt, ...labels(actual, true), ncritStartup: true, retryReason: record.reason });
    const attemptedIncrement = selected - actual;
    if (attemptedIncrement <= minimumStep) break;
    // Integer/target clipping can make the attempted increment smaller than
    // step; halve the increment actually tried, never repeat the same seed.
    step = Math.max(minimumStep, attemptedIncrement / 2);
  }
  const actualNcrit = retained.checkpoint.restart.options.ncrit;
  return { result: retained, diagnostics: { attempted: true, phase, actualNcrit, targetNcrit, maximumStep,
    reachedTarget: actualNcrit === targetNcrit, stateConverged: retained.converged, attempts } };
}
