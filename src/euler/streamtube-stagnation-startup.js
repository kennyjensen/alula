// SPDX-License-Identifier: GPL-2.0-or-later
import { createStreamtubeBodySystem } from './streamtube-body.js';
import { solveStreamtubeIses } from './streamtube-ises-update.js';

// A stalled wall-only chart can be changed without running initial SMOVE a
// second time. Re-express and re-evaluate the same admitted physical state;
// the receiving checkpoint still has to pass the ordinary exact replay gate.
export function recoverStreamtubeStagnationMotion(before, { maxIterations, tolerance = 1e-10,
  onIteration, onMesh } = {}) {
  const source = before?.checkpoint, used = before?.history?.at(-1)?.iteration;
  if (!Number.isInteger(maxIterations) || maxIterations < 0) throw new Error('Invalid stagnation recovery budget.');
  if (before?.converged || !source || !Number.isInteger(used) || used >= maxIterations
    || before.lastRejectedStep?.code !== 'EULER_GRID_STAGNATION'
    || source.input.stagnationMotion !== 'walls-only' || source.input.bodies.length < 2
    || !(source.input.mach <= .3) || !(before.diagnostics.maxMach < 1)) return before;
  const resume = structuredClone(source);
  resume.input.stagnationMotion = 'interpolated';
  const system = createStreamtubeBodySystem(resume.input);
  const state = system.adoptGeometry(Float64Array.from(resume.initialEuler.x), resume.initialEuler.nodes);
  const value = system.evaluate(state);
  let maximumCoordinateChange = 0;
  value.nodes.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => {
    const old = source.initialEuler.nodes[g][i][j];
    const change = Math.hypot(p.x - old.x, p.y - old.y);
    if (change > 64 * Number.EPSILON * Math.max(1, Math.abs(old.x), Math.abs(old.y)))
      throw new Error('Stagnation coordinate recovery changed the physical grid.');
    maximumCoordinateChange = Math.max(maximumCoordinateChange, change);
  })));
  if (state.some((x, i) => x !== source.initialEuler.x[i]))
    throw new Error('Stagnation coordinate recovery changed the packed flow state.');
  resume.initialEuler = { x: state, nodes: value.nodes };
  resume.residual = value.residual;
  const result = solveStreamtubeIses(undefined, { ...resume.continuation, resume,
    maxIterations: maxIterations - used, tolerance, retainCheckpoint: true,
    onIteration: h => onIteration?.({ ...h, iteration: used + h.iteration,
      startupStrategy: 'interpolated-stagnation-recovery' }),
    onMesh: frame => onMesh?.({ ...frame, iteration: { ...frame.iteration,
      iteration: used + frame.iteration.iteration, startupStrategy: 'interpolated-stagnation-recovery' } }),
  });
  result.history = [...before.history.slice(0, -1), ...result.history.map(h => ({ ...h,
    iteration: used + h.iteration, startupStrategy: 'interpolated-stagnation-recovery' }))];
  const first = before.linearDiagnostics, last = result.linearDiagnostics;
  result.linearDiagnostics = { ...last,
    solves: first.solves + last.solves, refinements: first.refinements + last.refinements,
    orderingFallbacks: first.orderingFallbacks + last.orderingFallbacks,
    maxRelativeResidual: Math.max(first.maxRelativeResidual, last.maxRelativeResidual),
    maxFactorNonzeros: Math.max(first.maxFactorNonzeros, last.maxFactorNonzeros) };
  result.stagnationRecovery = { method: 'interpolated-stagnation-recovery',
    initialIterations: used, remainingIterations: maxIterations - used,
    maximumCoordinateChange, repeatedLimit: before.lastRejectedStep.diagnostics,
    converged: result.converged };
  return result;
}
