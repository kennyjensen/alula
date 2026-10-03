// SPDX-License-Identifier: GPL-2.0-or-later
// Steihaug CG on the Gauss-Newton quadratic in column-normalized variables.
// Products with A and A^T implement H=A^T A; H is never formed or factored.
// This returns a trust-subproblem proposal, not a Newton linear certificate.
const dot = (a, b) => a.reduce((sum, v, i) => sum + v * b[i], 0);
const norm = a => Math.sqrt(dot(a, a));

export function scaledLeastSquaresTrustCG(matrix, residual, { radius, gradientTolerance = .01, maxIterations = 100 } = {}) {
  const { n, rowPtr, colIndex, values } = matrix;
  if (residual.length !== n || !residual.every(Number.isFinite) || !(radius > 0) || !Number.isFinite(radius)
    || !(gradientTolerance > 0 && gradientTolerance < 1) || !Number.isInteger(maxIterations) || maxIterations < 1)
    throw new Error('Invalid scaled least-squares trust CG controls.');
  const scales = new Float64Array(n);
  for (let k = 0; k < values.length; k++) scales[colIndex[k]] += values[k] ** 2;
  scales.forEach((v, i) => { if (!(v > 0 && Number.isFinite(v))) throw new Error('Unresolved trust CG column.'); scales[i] = Math.sqrt(v); });
  const a = Float64Array.from(values, (v, k) => v / scales[colIndex[k]]);
  const product = z => {
    const y = new Float64Array(n);
    for (let row = 0; row < n; row++) for (let k = rowPtr[row]; k < rowPtr[row + 1]; k++) y[row] += a[k] * z[colIndex[k]];
    return y;
  };
  const transpose = y => {
    const z = new Float64Array(n);
    for (let row = 0; row < n; row++) for (let k = rowPtr[row]; k < rowPtr[row + 1]; k++) z[colIndex[k]] += a[k] * y[row];
    return z;
  };
  const gradient = transpose(residual), gradientNorm = norm(gradient);
  if (!(gradientNorm > 0)) throw new Error('Stationary least-squares model; residual convergence is unproven.');
  let z = new Float64Array(n), r = gradient.map(v => -v), p = r.slice(), rr = dot(r, r), reason = 'iteration limit';
  const history = [];
  const boundaryLength = () => {
    const aa = dot(p, p), b = dot(z, p), c = (norm(z) - radius) * (norm(z) + radius);
    const root = Math.sqrt(b * b - aa * c);
    return b >= 0 ? -c / (b + root) : (-b + root) / aa;
  };
  for (let iteration = 1; iteration <= maxIterations; iteration++) {
    const ap = product(p), hp = transpose(ap), curvature = dot(ap, ap);
    const alpha = curvature > 0 ? rr / curvature : Infinity;
    const next = z.map((v, i) => v + alpha * p[i]);
    if (!Number.isFinite(alpha) || norm(next) >= radius) {
      const tau = boundaryLength(); z = z.map((v, i) => v + tau * p[i]);
      history.push({ iteration, scaledNorm: norm(z), boundary: true }); reason = 'trust boundary'; break;
    }
    z = next;
    const nextR = r.map((v, i) => v - alpha * hp[i]), nextRR = dot(nextR, nextR);
    history.push({ iteration, scaledNorm: norm(z), relativeRecursiveGradient: Math.sqrt(nextRR) / gradientNorm });
    if (Math.sqrt(nextRR) <= gradientTolerance * gradientNorm) { reason = 'gradient tolerance'; break; }
    const beta = nextRR / rr; p = nextR.map((v, i) => v + beta * p[i]); r = nextR; rr = nextRR;
  }
  const az = product(z), linearized = residual.map((v, i) => v + az[i]);
  const relativeGradientResidual = norm(transpose(linearized)) / gradientNorm;
  const predictedReduction = -dot(residual, az) - .5 * dot(az, az);
  if (!(predictedReduction > 64 * Number.EPSILON * .5 * dot(residual, residual)) || !z.every(Number.isFinite))
    throw new Error('Trust CG has no resolved model decrease.');
  return { direction: z.map((v, i) => v / scales[i]), scales, scaledStepNorm: norm(z),
    predictedReduction, relativeGradientResidual, gradientConverged: relativeGradientResidual <= gradientTolerance,
    linearRelativeResidual: norm(linearized) / norm(residual), reason, iterations: history.length, history };
}
