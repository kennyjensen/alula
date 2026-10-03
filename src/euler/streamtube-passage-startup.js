// SPDX-License-Identifier: GPL-2.0-or-later
import { createStreamtubeBodySystem } from './streamtube-body.js';
import { initializeStreamtubeStartup } from './streamtube-startup.js';
import { solveStreamtubeIses } from './streamtube-ises-update.js';
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid, gilesStreamwiseCoordinates } from '../geometry/elliptic-streamtube-grid.js';
import { assertConvexStreamtubeGrid } from '../geometry/streamtube-convex-step.js';

// The main element's station coordinate can be a poor label for a passage
// containing another element's nose. Both banks already prescribe the common
// block intersections. Average their normalized Giles coordinates to smooth
// the interior without moving a wall, cut, farfield or endpoint.
export function smoothStreamtubePassages(prepared) {
  const { system, initial, nodes } = prepared;
  const { allocation } = system.decode(initial), reports = [];
  const smoothed = nodes.map((group, g) => {
    const lower = gilesStreamwiseCoordinates(group.map(row => row[0]));
    const upper = gilesStreamwiseCoordinates(group.map(row => row.at(-1)));
    const xi = lower.map((value, i) => .5 * (value + upper[i]));
    const grid = createEllipticStreamtubeGrid({ nodes: group,
      massFlows: allocation.groups[g].map(tube => tube.massFlow),
      streamwiseCoordinates: xi, discretization: 'giles-1985' });
    const result = smoothEllipticStreamtubeGrid(grid, { maxSweeps: 1000,
      tolerance: 1e-8, omega: 1, requireConvex: true, detectCoordinateRoundoff: true });
    reports.push({ group: g, converged: result.converged, reason: result.reason,
      iterations: result.history.at(-1)?.iteration, residual: result.history.at(-1)?.residual,
      minCornerSine: result.history.at(-1)?.minCornerSine });
    return result.nodes;
  });
  if (!reports.every(r => r.converged)) return { converged: false, reports };
  assertConvexStreamtubeGrid(smoothed);
  smoothed.forEach((group, g) => group.forEach((row, i) => row.forEach((point, j) => {
    if (i && i < group.length - 1 && j && j < row.length - 1) return;
    const original = nodes[g][i][j];
    if (point.x !== original.x || point.y !== original.y)
      throw new Error('Passage smoothing moved a prescribed grid boundary.');
  })));
  return { converged: true, nodes: smoothed, reports };
}

export function recoverStreamtubePassageGrid(before, prepared, { maxIterations, tolerance = 1e-10,
  onIteration, onMesh } = {}) {
  if (!Number.isInteger(maxIterations) || maxIterations < 0) throw new Error('Invalid passage recovery budget.');
  const used = before.history.at(-1).iteration;
  if (before.converged || !['EULER_GRID_STAGNATION', 'EULER_PASSAGE_STAGNATION'].includes(before.lastRejectedStep?.code)
    || used >= maxIterations || prepared.input.bodies.length < 2
    || !(prepared.input.mach <= .3) || !(before.diagnostics.maxMach < 1)) return before;
  const smoothing = smoothStreamtubePassages(prepared);
  if (!smoothing.converged) return { ...before, passageRecovery: { attempted: true,
    converged: false, stage: 'smoothing', smoothing: smoothing.reports } };
  const input = { ...prepared.input, stagnationMotion: 'interpolated' };
  const system = createStreamtubeBodySystem(input);
  const state = system.adoptGeometry(prepared.initial, smoothing.nodes);
  const gas = initializeStreamtubeStartup(system, state);
  const result = solveStreamtubeIses(input, { initialEuler: { x: gas.initial, nodes: gas.flow.nodes },
    initialRedistributionMode: 'supplied', tangentialGridRecovery: true,
    maxIterations: maxIterations - used, tolerance, retainCheckpoint: true, retainBestCheckpoint: true,
    iterationGeometry: 'ises-sampled', stepAcceptance: 'armijo', adaptiveMcrit: before.checkpoint?.continuation.adaptiveMcrit ?? true,
    ...(before.checkpoint?.continuation.targetMcrit === undefined ? {} : { targetMcrit: before.checkpoint.continuation.targetMcrit }),
    stopOnGridStagnation: true, gridCorrectionBacktracking: 'before-damping',
    onIteration: h => onIteration?.({ ...h, iteration: used + h.iteration, startupStrategy: 'passage-grid-recovery' }),
    onMesh: frame => onMesh?.({ ...frame, iteration: { ...frame.iteration,
      iteration: used + frame.iteration.iteration, startupStrategy: 'passage-grid-recovery' } }) });
  result.history = [...before.history.slice(0, -1), ...result.history.map(h => ({ ...h,
    iteration: used + h.iteration, startupStrategy: 'passage-grid-recovery' }))];
  if (result.bestCheckpoint) result.bestCheckpoint.iteration += used;
  const a = before.linearDiagnostics, b = result.linearDiagnostics;
  result.linearDiagnostics = { ...b, solves: a.solves + b.solves, refinements: a.refinements + b.refinements,
    orderingFallbacks: a.orderingFallbacks + b.orderingFallbacks,
    maxRelativeResidual: Math.max(a.maxRelativeResidual, b.maxRelativeResidual),
    maxFactorNonzeros: Math.max(a.maxFactorNonzeros, b.maxFactorNonzeros) };
  result.stagnationRecovery = before.stagnationRecovery;
  result.passageRecovery = { attempted: true, converged: result.converged,
    initialIterations: used, remainingIterations: maxIterations - used,
    smoothing: smoothing.reports, initialResidual: gas.flow.diagnostics.residual,
    method: 'fixed-boundary passage smoothing and coordinated tangential Newton recovery' };
  return result;
}
