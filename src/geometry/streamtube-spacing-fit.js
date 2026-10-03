// SPDX-License-Identifier: GPL-2.0-or-later
// Fit an actual streamline distance. Every bracket endpoint is measured;
// the wall is never substituted for an unevaluated zero-mass streamline.
export function fitNormalStreamtubeFraction({ initialFraction, targetDistance, distanceAtFraction,
  relativeTolerance = 1e-4, maxIterations = 30 }) {
  if (![initialFraction, targetDistance, relativeTolerance].every(Number.isFinite)
    || initialFraction <= 0 || targetDistance <= 0 || relativeTolerance <= 0
    || !Number.isInteger(maxIterations) || maxIterations < 1 || typeof distanceAtFraction !== 'function')
    throw new Error('Invalid normal-spacing fit controls.');
  // The local strain approximation may overestimate even an attainable
  // tube mass. An out-of-range estimate chooses an interior trial only;
  // every accepted mass still comes from the unchanged measured-distance
  // equation. Do not clip a fitted endpoint mass or alter its target.
  let fraction = initialFraction < 1 ? initialFraction : .5, lower, upper; const history = [];
  for (let i = 0; i < maxIterations; i++) {
    const distance = distanceAtFraction(fraction);
    if (!(distance > 0) || !Number.isFinite(distance)) throw new Error('Invalid measured normal streamline distance.');
    const relativeError = distance / targetDistance - 1;
    history.push({ fraction, distance, relativeError });
    if (Math.abs(relativeError) <= relativeTolerance) return { fraction, distance, relativeError, history };
    if ((lower && (fraction <= lower.fraction || distance <= lower.distance))
      || (upper && (fraction >= upper.fraction || distance >= upper.distance)))
      throw new Error('Normal cross-line distance is not monotone in streamtube mass.');
    const value = { fraction, distance };
    if (distance < targetDistance) lower = value; else upper = value;
    if (lower && upper) {
      const t = Math.max(.05, Math.min(.95, Math.log(targetDistance / lower.distance) / Math.log(upper.distance / lower.distance)));
      fraction = Math.exp((1 - t) * Math.log(lower.fraction) + t * Math.log(upper.fraction));
    } else {
      // Linear stagnation scaling supplies the first bracket. Safeguards
      // limit extrapolation; they never replace a measured distance.
      const ratio = (targetDistance / distance) ** 2;
      const factor = upper ? Math.max(.1, Math.min(.9, ratio)) : Math.max(1.1, Math.min(10, ratio));
      fraction = upper ? fraction * factor : Math.min(.5 * (fraction + 1), fraction * factor);
    }
    if (!(fraction > 0 && fraction < 1)) throw new Error('Normal-spacing fit exhausted its mass interval.');
  }
  throw new Error('Normal-spacing distance fit did not converge.');
}
