// SPDX-License-Identifier: GPL-2.0-or-later
// Match equations to unknowns before AMD sees A+A^T. A streamwise sweep
// leaves a wide front on many-tube grids; AMD can eliminate across tubes too.
import { solveSparseDirectAligned } from '../numerics/klu.js';
import { createStreamtubeStationOrdering } from './streamtube-station-ordering.js';
import { solveCoupledLinearSystem } from './streamtube-coupled-linear-solve.js';

export function solveCoupledAlignedSystem(matrix, rhs, controls = {}) {
  const started = performance.now();
  let matchingMilliseconds = 0;
  try {
    const matching = createStreamtubeStationOrdering({ matrix, layout: controls.layout, stations: controls.stations });
    // Matching returns position-to-original row AND column indices. Only
    // permute rows here: variables retain their original packed ordering.
    const rows = new Int32Array(matrix.n);
    for (let k = 0; k < matrix.n; k++) rows[matching.Quser[k]] = matching.Puser[k];
    matchingMilliseconds = performance.now() - started;
    const linear = solveSparseDirectAligned(matrix, rhs, rows, {
      preferredOrdering: controls.preferredOrdering, pivotTolerance: controls.pivotTolerance });
    return { linear, ordering: null, stationPolicy: {
      version: 1, mode: 'aligned-auto', attempted: true, accepted: linear.equationOrdering === 'aligned',
      fallback: linear.equationOrdering !== 'aligned', recommendation: 'aligned-auto',
      selected: { ordering: linear.ordering, equationOrdering: linear.equationOrdering, pivotTolerance: linear.pivotTolerance },
      timings: { matchingMilliseconds, alignedSolveMilliseconds: performance.now() - started - matchingMilliseconds,
        totalMilliseconds: performance.now() - started },
    } };
  } catch (error) {
    // Preserve cancellation/programming errors. Only a numerical rejection
    // can try the established station/automatic policy on the SAME system.
    if (['AbortError', 'CanceledError', 'CancelledError'].includes(error?.name)
      || ['ABORT_ERR', 'ERR_CANCELED', 'ERR_CANCELLED'].includes(error?.code)
      || !(error?.code === 'STREAMTUBE_STATION_MATCHING_FAILED' || error?.code === 'KLU_RESIDUAL_LIMIT'
        || error?.attempts?.length && error.attempts.every(a => [1, -2].includes(a.status)))) throw error;
    const alignedFailure = { code: error.code, message: error.message };
    const result = solveCoupledLinearSystem(matrix, rhs, { ...controls, mode: 'station-auto' });
    result.stationPolicy = { ...result.stationPolicy, mode: 'aligned-auto', fallback: true,
      recommendation: 'aligned-auto', alignedFailure,
      timings: { ...result.stationPolicy.timings, totalMilliseconds: performance.now() - started } };
    return result;
  }
}
