// Two adjacent doubles at the native wake DSLIM/BLVAR contact.
export function wakeHkContactCase() {
  const parameters = { reynolds: 992889.5740457802, mach: .2, ncrit: 9, gamma: 1.4,
    velocityConvention: 'physical', exactJacobian: true, transitionTolerance: 1e-12 };
  const minimumHk = 1.00005, theta = .0044, ue = .99, gm1 = parameters.gamma - 1;
  const hstinv = gm1 * parameters.mach ** 2 / (1 + .5 * gm1 * parameters.mach ** 2);
  const machSquared = ue * ue * hstinv / (gm1 * (1 - .5 * ue * ue * hstinv));
  const deltaStar = theta * (minimumHk * (1 + .113 * machSquared) + .29 * machSquared);
  const bits = new DataView(new ArrayBuffer(8)); bits.setFloat64(0, deltaStar);
  bits.setBigUint64(0, bits.getBigUint64(0) - 1n); const belowDeltaStar = bits.getFloat64(0);
  const downstream = { s: 2.5, theta, deltaStar, ue, aux: .02, wakeGap: 0 };
  const upstream = { s: 2, theta: .0046, deltaStar: .005, ue: .98, aux: .021, wakeGap: 0 };
  return { parameters, minimumHk, below: { upstream, downstream: { ...downstream, deltaStar: belowDeltaStar }, regime: 'wake' },
    contact: { upstream, downstream, regime: 'wake' }, scope: 'Manufactured physically admissible ordinary wake interval at the same Mach/Re as the saved RAE case. It is a local branch discriminant, not a captured full RAE state.' };
}
