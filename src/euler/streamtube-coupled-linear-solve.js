// SPDX-License-Identifier: GPL-2.0-or-later
// Ordering policy only: the complete original Jacobian/RHS and KLU accuracy
// gate are unchanged. This helper does not construct or evaluate a flow.
import { solveSparseDirect } from '../numerics/klu.js';
import { createStreamtubeStationOrdering } from './streamtube-station-ordering.js';

const cancelled = error => ['AbortError', 'CanceledError', 'CancelledError'].includes(error?.name)
  || ['ABORT_ERR', 'ERR_CANCELED', 'ERR_CANCELLED'].includes(error?.code);
// Upstream klu.h status values. The typed factor failure is returned only
// after native/bridge/JS cleanup, so a lower-fill ordering may recover OOM.
const KLU_SINGULAR = 1, KLU_OUT_OF_MEMORY = -2;
const numericalRejection = error => !cancelled(error) && (error?.code === 'KLU_RESIDUAL_LIMIT'
  || error?.code === 'KLU_GIVEN_FACTORIZATION' && [KLU_SINGULAR, KLU_OUT_OF_MEMORY].includes(error.status));
const matchingRejection = error => !cancelled(error) && error?.code === 'STREAMTUBE_STATION_MATCHING_FAILED';
const failure = (error, stage) => ({ stage, name: error.name, message: error.message,
  ...Object.fromEntries(['code', 'status', 'relativeResidual', 'attempts', 'diagnostics']
    .filter(key => error[key] !== undefined).map(key => [key, structuredClone(error[key])])) });

export function solveCoupledLinearSystem(matrix, rhs, { layout, stations, preferredOrdering = 'amd',
  pivotTolerance = .001, mode = 'auto', stationEnabled = true,
  stationSkipReason = 'previous-station-fallback' } = {}) {
  if (!['auto', 'station', 'station-auto'].includes(mode))
    throw new Error('Unknown coupled linear ordering policy.');
  const automatic = () => solveSparseDirect(matrix, rhs, { preferredOrdering, pivotTolerance });
  const match = () => createStreamtubeStationOrdering({ matrix, layout, stations });
  const station = ordering => solveSparseDirect(matrix, rhs, { ordering: 'given',
    rowPermutation: ordering.Puser, columnPermutation: ordering.Quser,
    btf: false, pivotTolerance, pivotFallback: false });
  // Preserve legacy options, errors and complete linear payloads exactly.
  if (mode === 'auto') return { linear: automatic(), ordering: null };
  if (mode === 'station') {
    const ordering = match();
    return { linear: station(ordering), ordering };
  }
  if (typeof stationEnabled !== 'boolean' || !['amd', 'colamd'].includes(preferredOrdering)
    || !Number.isFinite(pivotTolerance) || pivotTolerance <= 0 || pivotTolerance > 1
    || !stationEnabled && (typeof stationSkipReason !== 'string' || !stationSkipReason.trim()))
    throw new Error('Invalid automatic station ordering controls.');
  const started = performance.now(), timings = { matchingMilliseconds: 0, stationSolveMilliseconds: 0,
    autoSolveMilliseconds: 0, totalMilliseconds: 0 };
  const policy = { version: 1, mode, attempted: stationEnabled, accepted: false, fallback: false,
    recommendation: stationEnabled ? 'station-auto' : 'auto', timings,
    ...(!stationEnabled ? { skippedReason: stationSkipReason } : {}) };
  const finish = linear => {
    timings.totalMilliseconds = performance.now() - started;
    if (linear) policy.selected = { ordering: linear.ordering ?? 'none',
      ...(linear.pivotTolerance === undefined ? {} : { pivotTolerance: linear.pivotTolerance }) };
    return structuredClone(policy);
  };
  let ordering;
  if (stationEnabled) {
    let stage = 'matching', clock = performance.now();
    try {
      ordering = match(); timings.matchingMilliseconds = performance.now() - clock;
      stage = 'station-solve'; clock = performance.now();
      const linear = station(ordering); timings.stationSolveMilliseconds = performance.now() - clock;
      policy.accepted = true;
      return { linear, ordering, stationPolicy: finish(linear) };
    } catch (error) {
      if (stage === 'matching') timings.matchingMilliseconds = performance.now() - clock;
      else timings.stationSolveMilliseconds = performance.now() - clock;
      if (!(stage === 'matching' ? matchingRejection(error) : numericalRejection(error))) throw error;
      policy.fallback = true; policy.recommendation = 'auto';
      policy.stationFailure = failure(error, stage);
    }
  }
  const clock = performance.now();
  try {
    const linear = automatic(); timings.autoSolveMilliseconds = performance.now() - clock;
    return { linear, ordering: null, stationPolicy: finish(linear) };
  } catch (error) {
    // Never retry this failure. Known factor/certificate failures retain both
    // policies; cancellation and unexpected runtime errors pass through intact.
    const automaticFactorFailure = error?.attempts?.length > 0
      && error.attempts.every(attempt => ['amd', 'colamd'].includes(attempt.ordering))
      && [KLU_SINGULAR, KLU_OUT_OF_MEMORY].includes(error.attempts.at(-1).status);
    if (!cancelled(error) && (numericalRejection(error) || automaticFactorFailure)) {
      timings.autoSolveMilliseconds = performance.now() - clock;
      policy.autoFailure = failure(error, 'auto-solve');
      const stationPolicy = finish();
      error.stationPolicy = stationPolicy;
      error.diagnostics = { ...error.diagnostics, stationPolicy };
    }
    throw error;
  }
}
