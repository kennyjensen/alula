// SPDX-License-Identifier: GPL-2.0-or-later
// Vinokur NASA-CR-3313, Eqs. (44), (47), (50), printed pp. 13-14.
// Independent discrete calibration: fit the FIRST AND LAST INTERVALS,
// rather than interpreting an endpoint derivative as a finite cell length.
// This is not a recovered MSET endpoint-reconciliation rule.

function logSinhc(x) {
  if (x < .1) {
    const z = x * x;
    return z * (1 / 6 - z * (1 / 180 - z * (1 / 2835 - z / 37800)));
  }
  return x < 20 ? Math.log(Math.sinh(x) / x)
    : x + Math.log(-Math.expm1(-2 * x)) - Math.log(2 * x);
}
function logCosh(x) {
  return x < 20 ? Math.log1p(2 * Math.sinh(x / 2) ** 2)
    : x + Math.log1p(Math.exp(-2 * x)) - Math.LN2;
}
function logSinc(x) {
  if (x < .1) {
    const z = x * x;
    return -z * (1 / 6 + z * (1 / 180 + z * (1 / 2835 + z / 37800)));
  }
  return Math.log(Math.sin(x) / x);
}
function logSymmetric(t, beta, hyperbolic) {
  // 0 < t <= 1/2. The equivalent sin/sinh identity avoids subtracting
  // nearly equal tan/tanh values close to a clustered endpoint.
  const a = beta * t, b = beta / 2, c = beta * (.5 - t);
  return Math.log(t) + (hyperbolic
    ? logSinhc(a) - logSinhc(b) - logCosh(c)
    : logSinc(a) - logSinc(b) - Math.log(Math.cos(c)));
}
function logistic(x) {
  return x > 0 ? 1 / (1 + Math.exp(-x)) : Math.exp(x) / (1 + Math.exp(x));
}

export function createVinokurStations({ length, intervals, firstSpacing, lastSpacing } = {}) {
  if (!Number.isSafeInteger(intervals) || intervals < 3 || intervals > 1000000
    || ![length, firstSpacing, lastSpacing].every(v => Number.isFinite(v) && v > 0)
    || !(firstSpacing < length - lastSpacing))
    throw new Error('Vinokur stations require at least three intervals and positive end lengths whose sum is below the block length.');
  // Odds at t=1/N and 1-1/N eliminate A analytically. Only beta remains:
  // v^2 = a*b/((L-a)*(L-b)), A^2 = b*(L-a)/(a*(L-b)).
  const logFirstOdds = Math.log(firstSpacing) - Math.log(length - firstSpacing);
  const logLastOdds = Math.log(lastSpacing) - Math.log(length - lastSpacing);
  const logV = .5 * (logFirstOdds + logLastOdds), logA = .5 * (logLastOdds - logFirstOdds);
  const logTarget = logV - Math.log1p(Math.exp(logV));
  const t1 = 1 / intervals, logUniform = Math.log(t1);
  const branch = logTarget === logUniform ? 'rational' : logTarget < logUniform ? 'hyperbolic' : 'trigonometric';
  const hyperbolic = branch === 'hyperbolic';
  let parameter = 0;
  if (branch !== 'rational') {
    let lo = 0, hi = hyperbolic ? 1 : Math.PI - 2 * Number.EPSILON;
    if (hyperbolic) while (logSymmetric(t1, hi, true) > logTarget && hi < 65536) hi *= 2;
    if (hyperbolic ? logSymmetric(t1, hi, true) > logTarget : logSymmetric(t1, hi, false) < logTarget)
      throw new Error('Vinokur interval fit has no representable parameter.');
    for (let k = 0; k < 128; k++) {
      const mid = lo + .5 * (hi - lo), residual = logSymmetric(t1, mid, hyperbolic) - logTarget;
      parameter = mid;
      if (residual === 0 || mid === lo || mid === hi) break;
      if (hyperbolic ? residual > 0 : residual < 0) lo = mid; else hi = mid;
    }
    if (Math.abs(logSymmetric(t1, parameter, hyperbolic) - logTarget) > 128 * Number.EPSILON * Math.max(1, Math.abs(logTarget)))
      throw new Error('Vinokur interval fit failed its endpoint residual check.');
  }
  const value = t => {
    if (!Number.isFinite(t) || t < 0 || t > 1) throw new Error('Vinokur station argument must lie in [0,1].');
    if (t === 0) return 0; if (t === 1) return length;
    const lowT = Math.min(t, 1 - t);
    const logU = branch === 'rational' ? Math.log(lowT) : logSymmetric(lowT, parameter, hyperbolic);
    const logOneMinusU = Math.log1p(-Math.exp(logU));
    const logOdds = (t <= .5 ? logU - logOneMinusU : logOneMinusU - logU) - logA;
    return length * logistic(logOdds);
  };
  const positions = Array.from({ length: intervals + 1 }, (_, i) => value(i / intervals));
  const tolerance = 256 * Number.EPSILON * length;
  if (Math.abs(positions[1] - firstSpacing) > tolerance
    || Math.abs((length - positions.at(-2)) - lastSpacing) > tolerance)
    throw new Error('Vinokur sampled endpoint intervals do not match the requested lengths.');
  if (positions.some((s, i) => !Number.isFinite(s) || i && !(s > positions[i - 1])))
    throw new Error('Vinokur sampled intervals are not resolved and increasing.');
  return { positions, value, parameter, branch, logAsymmetry: logA,
    firstSpacing: positions[1], lastSpacing: length - positions.at(-2) };
}
