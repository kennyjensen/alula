// Independent perfect-gas formulas. No runtime gas, BL or cell imports.
// Units: freestream density and speed are one; R is absorbed into p/rho.
export function referenceGas({ gamma = 1.4, freestreamMach = .3 } = {}) {
  if (!(gamma > 1) || !(freestreamMach > 0 && freestreamMach < 1)) throw Error('Invalid reference gas.');
  const pInf = 1 / (gamma * freestreamMach ** 2);
  const h0 = gamma * pInf / (gamma - 1) + .5;
  const totalRatio = 1 + .5 * (gamma - 1) * freestreamMach ** 2;
  return { gamma, freestreamMach, pInf, h0,
    rhoTotal: totalRatio ** (1 / (gamma - 1)),
    pTotal: pInf * totalRatio ** (gamma / (gamma - 1)) };
}

export function isentropeAtSpeed(q, gas) {
  const { gamma, h0, pTotal, rhoTotal } = gas, enthalpy = h0 - .5 * q ** 2;
  if (!(q > 0) || !(enthalpy > 0)) throw Error('Invalid isentropic speed.');
  return { q, enthalpy, rho: rhoTotal * (enthalpy / h0) ** (1 / (gamma - 1)),
    p: pTotal * (enthalpy / h0) ** (gamma / (gamma - 1)),
    machSquared: q ** 2 / ((gamma - 1) * enthalpy) };
}

export function isentropeAtMach(mach, gas) {
  const h = gas.h0 / (1 + .5 * (gas.gamma - 1) * mach ** 2);
  return isentropeAtSpeed(mach * Math.sqrt((gas.gamma - 1) * h), gas);
}

export function normalShock(mach, gas) {
  if (!Number.isFinite(mach) || mach < 1) throw Error('A normal shock requires upstream Mach >= 1.');
  const upstream = isentropeAtMach(mach, gas), { gamma } = gas;
  const densityRatio = (gamma + 1) * mach ** 2 / (2 + (gamma - 1) * mach ** 2);
  const pressureRatio = (2 * gamma * mach ** 2 - (gamma - 1)) / (gamma + 1);
  const downstream = { q: upstream.q / densityRatio, rho: upstream.rho * densityRatio,
    p: upstream.p * pressureRatio };
  downstream.enthalpy = gamma / (gamma - 1) * downstream.p / downstream.rho;
  downstream.machSquared = (1 + .5 * (gamma - 1) * mach ** 2) / (gamma * mach ** 2 - .5 * (gamma - 1));
  const entropyOverR = (Math.log(pressureRatio) - gamma * Math.log(densityRatio)) / (gamma - 1);
  return { upstream, downstream, densityRatio, pressureRatio, entropyOverR };
}

// p(q)=K(h0-q²/2)^n, n=gamma/(gamma-1). The fourth-derivative
// bound is the triangle inequality over the entire symmetric speed interval.
export function pressureTaylorBound(q, fullSpeedDifference, gas) {
  const { gamma, h0, pTotal } = gas, n = gamma / (gamma - 1), k = pTotal / h0 ** n;
  const half = .5 * fullSpeedDifference, qMax = q + half, qMin = q - half;
  if (!(qMin > 0 && fullSpeedDifference > 0)) throw Error('Invalid symmetric speed interval.');
  const h = h0 - .5 * q ** 2, hMin = h0 - .5 * qMax ** 2, hMax = h0 - .5 * qMin ** 2;
  if (!(hMin > 0)) throw Error('Nonpositive static enthalpy in Taylor interval.');
  const powerBound = exponent => (exponent < 0 ? hMin : hMax) ** exponent;
  const pSecond = k * (-n * h ** (n - 1) + n * (n - 1) * q ** 2 * h ** (n - 2));
  const fourthMagnitudeBound = k * (
    3 * Math.abs(n * (n - 1)) * powerBound(n - 2)
    + 6 * Math.abs(n * (n - 1) * (n - 2)) * qMax ** 2 * powerBound(n - 3)
    + Math.abs(n * (n - 1) * (n - 2) * (n - 3)) * qMax ** 4 * powerBound(n - 4));
  return { pSecond, leadingDifference: pSecond * fullSpeedDifference ** 2 / 8,
    fourthMagnitudeBound, remainderBound: fourthMagnitudeBound * fullSpeedDifference ** 4 / 384 };
}
