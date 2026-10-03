// SPDX-License-Identifier: GPL-2.0-or-later
// Standalone XFOIL XICALC dead-air wake-gap cubic. This is a viscous
// displacement contribution, not an Euler constant-width base wake.
// ANTE and endpoint parameter derivatives are supplied by the caller.
export function createXfoilDeadAirGap({ normalGap, upperDerivative, lowerDerivative, sharp = false }) {
  if (!Number.isFinite(normalGap) || normalGap < 0 || typeof sharp !== 'boolean'
    || ![upperDerivative, lowerDerivative].every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y)))
    throw new Error('Invalid XFOIL dead-air gap geometry.');
  const a = Math.hypot(upperDerivative.x, upperDerivative.y), b = Math.hypot(lowerDerivative.x, lowerDerivative.y);
  if (!(a > 0 && b > 0) || !Number.isFinite(a) || !Number.isFinite(b))
    throw new Error('XFOIL dead-air gap requires nonzero finite endpoint derivatives.');
  const cross = (upperDerivative.x / a) * (lowerDerivative.y / b)
    - (upperDerivative.y / a) * (lowerDerivative.x / b);
  if (!Number.isFinite(cross) || Math.abs(cross) > 1 + 8 * Number.EPSILON)
    throw new Error('Invalid normalized trailing-edge tangent cross product.');
  const lengthToGapRatio = 2.5, limit = 3 / lengthToGapRatio;
  const threshold = limit / Math.sqrt(1 + limit * limit);
  // Evaluate the capped branch before sqrt(1-cross^2). This is the same
  // XICALC clamp, also regular at perpendicular tangents; it is not a clip
  // to force an out-of-domain or negative width into a valid one.
  const rawSlope = Math.abs(cross) < 1 ? cross / Math.sqrt(1 - cross * cross) : null;
  const slope = Math.abs(cross) >= threshold ? Math.sign(cross) * limit : rawSlope;
  const width = sharp ? 0 : normalGap, closureDistance = lengthToGapRatio * width;
  if (!Number.isFinite(closureDistance)) throw new Error('Nonfinite dead-air gap closure distance.');
  const aa = 3 + lengthToGapRatio * slope, bb = -2 - lengthToGapRatio * slope;
  const at = distance => {
    if (!Number.isFinite(distance) || distance < 0) throw new Error('Dead-air gap distance must be finite and nonnegative.');
    if (width === 0 || distance >= closureDistance) return { gap: 0, dDistance: 0 };
    const z = 1 - distance / closureDistance;
    return { gap: width * (aa + bb * z) * z * z,
      dDistance: -(2 * aa * z + 3 * bb * z * z) / lengthToGapRatio };
  };
  return { normalGap, closureDistance, slope, rawSlope, rawSlopeUnbounded: rawSlope === null,
    normalizedTangentCross: cross, slopeLimited: Math.abs(cross) >= threshold,
    lengthToGapRatio, sharp, at };
}
