// SPDX-License-Identifier: GPL-2.0-or-later
// Numerical update metadata only. Never place this policy in BL physics options.
export function coupledShearCoordinate(value) {
  if (value === undefined) return 'linear';
  if (value !== 'linear' && value !== 'logarithmic')
    throw new Error('Invalid coupled shear-coordinate policy.');
  return value;
}

export function coupledCheckpointShearCoordinate(checkpoint) {
  const c = checkpoint?.continuation;
  // Legacy transfer sources may lack continuation metadata; their original
  // policy is linear. Actual checkpoint resumes retain the driver's guards.
  if (c === undefined) return 'linear';
  if (!c || typeof c !== 'object' || Array.isArray(c))
    throw new Error('Invalid coupled checkpoint continuation controls.');
  const value = coupledShearCoordinate(c.shearCoordinate);
  if (value === 'logarithmic' && c.blUpdate !== 'xfoil')
    throw new Error('Logarithmic shear requires the XFOIL BL update policy.');
  return value;
}

export function coupledFreshShearCoordinate(prepared, recoveryPlan) {
  const c = prepared?.iterationControls;
  if (c !== undefined && (!c || typeof c !== 'object' || Array.isArray(c)))
    throw new Error('Invalid prepared coupled iteration controls.');
  const inherited = coupledShearCoordinate(c?.shearCoordinate);
  if (recoveryPlan === undefined || recoveryPlan === null) return inherited;
  const selected = coupledShearCoordinate(recoveryPlan.shearCoordinate);
  if (c !== undefined && selected !== inherited)
    throw new Error('Prepared grid and recovery disagree on shear-coordinate policy.');
  return selected;
}

export function coupledResultShearCoordinate(result) {
  return result?.checkpoint === undefined
    ? coupledShearCoordinate(result?.shearCoordinate)
    : coupledCheckpointShearCoordinate(result.checkpoint);
}
