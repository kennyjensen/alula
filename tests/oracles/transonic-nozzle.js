// Independent perfect-gas, quasi-one-dimensional nozzle reference.
// No numerical solver, cell residual, speed filter, or runtime imports.
// NASA Glenn perfect-gas area-ratio and normal-shock relations:
// https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
// https://www.grc.nasa.gov/www/k-12/airplane/normal.html
// Reservoir conditions and downstream static pressure locate an internal
// shock; both ends must stay subsonic. This slender-flow reference is not
// an exact two-dimensional solution on the moving streamline grid.

function gas(gamma) {
  if (!Number.isFinite(gamma) || gamma <= 1) throw new Error('Require finite gamma > 1.');
}

export function areaMachRatio(mach, gamma = 1.4) {
  gas(gamma);
  if (!Number.isFinite(mach) || mach <= 0) throw new Error('Require positive finite Mach number.');
  // log1p makes the sonic value exactly one without subtracting two logs.
  const exponent = (gamma + 1) / (2 * (gamma - 1));
  return Math.exp(exponent * Math.log1p((gamma - 1) / (gamma + 1) * (mach * mach - 1)) - Math.log(mach));
}

export function machFromAreaRatio(ratio, { branch, gamma = 1.4 } = {}) {
  gas(gamma);
  if (!Number.isFinite(ratio) || ratio < 1 || !['subsonic', 'supersonic'].includes(branch))
    throw new Error('Require area ratio >= 1 and an explicit subsonic/supersonic branch.');
  if (ratio === 1) return 1;
  let lo = branch === 'subsonic' ? 0 : 1, hi = branch === 'subsonic' ? 1 / ratio : 2;
  while (branch === 'supersonic' && areaMachRatio(hi, gamma) < ratio) hi *= 2;
  for (let k = 0; k < 100; k++) {
    const mid = .5 * (lo + hi);
    if (mid === lo || mid === hi) break;
    const below = areaMachRatio(mid, gamma) < ratio;
    if (below === (branch === 'supersonic')) lo = mid;
    else hi = mid;
  }
  return .5 * (lo + hi);
}

export function normalShock(machUpstream, { gamma = 1.4 } = {}) {
  gas(gamma);
  if (!Number.isFinite(machUpstream) || machUpstream < 1)
    throw new Error('A normal compression shock requires upstream Mach >= 1.');
  const m2 = machUpstream * machUpstream;
  if (!Number.isFinite(m2)) throw new Error('Nonfinite normal-shock Mach squared.');
  const densityRatio = (gamma + 1) * m2 / ((gamma - 1) * m2 + 2);
  const pressureRatio = (2 * gamma * m2 - (gamma - 1)) / (gamma + 1);
  const logTotalPressureRatio = gamma / (gamma - 1) * Math.log(densityRatio)
    - Math.log(pressureRatio) / (gamma - 1);
  return {
    machUpstream,
    machDownstream: Math.sqrt((1 + .5 * (gamma - 1) * m2) / (gamma * m2 - .5 * (gamma - 1))),
    densityRatio, pressureRatio, temperatureRatio: pressureRatio / densityRatio,
    totalPressureRatio: Math.exp(logTotalPressureRatio),
    entropyRiseOverR: -logTotalPressureRatio,
  };
}

export function transonicNozzleReference({ gamma = 1.4, stagnationEnthalpy = 3.5,
  stagnationDensity = 1, throatArea = .02, areaRatio = 2.5, shockMach = 2,
  extensionLength = .25 } = {}) {
  gas(gamma);
  if (![stagnationEnthalpy, stagnationDensity, throatArea].every(v => Number.isFinite(v) && v > 0)
    || !Number.isFinite(areaRatio) || areaRatio <= 1 || !Number.isFinite(extensionLength) || extensionLength < 0)
    throw new Error('Invalid nozzle reservoir, area, or extension parameters.');
  const pTotal = (gamma - 1) / gamma * stagnationDensity * stagnationEnthalpy;
  const chokedMass = throatArea * stagnationDensity * Math.sqrt((gamma - 1) * stagnationEnthalpy)
    * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)));
  const area = x => {
    if (!Number.isFinite(x)) throw new Error('Require a finite nozzle coordinate.');
    const t = Math.min(1, Math.abs(x));
    return throatArea * (1 + (areaRatio - 1) * Math.sin(Math.PI * t / 2) ** 2);
  };
  const areaSlope = x => {
    if (!Number.isFinite(x)) throw new Error('Require a finite nozzle coordinate.');
    return Math.abs(x) >= 1 ? 0 : throatArea * (areaRatio - 1) * Math.PI / 2 * Math.sin(Math.PI * x);
  };
  const root = (ratio, branch) => machFromAreaRatio(ratio, { gamma, branch });
  const machBranches = x => ({ subsonic: root(area(x) / throatArea, 'subsonic'),
    supersonic: root(area(x) / throatArea, 'supersonic') });
  const stateForMach = (mach, localTotalPressure = pTotal) => {
    if (!Number.isFinite(mach) || mach < 0 || !Number.isFinite(localTotalPressure) || localTotalPressure <= 0)
      throw new Error('Invalid physical Mach number or total pressure.');
    const enthalpy = stagnationEnthalpy / (1 + .5 * (gamma - 1) * mach * mach);
    const localTotalDensity = stagnationDensity * localTotalPressure / pTotal;
    const rho = localTotalDensity * (enthalpy / stagnationEnthalpy) ** (1 / (gamma - 1));
    return { mach, machSquared: mach * mach, q: mach * Math.sqrt((gamma - 1) * enthalpy),
      rho, p: (gamma - 1) / gamma * rho * enthalpy, enthalpy, stagnationEnthalpy,
      stagnationDensity: localTotalDensity, totalPressure: localTotalPressure,
      entropyRiseOverR: -Math.log(localTotalPressure / pTotal) };
  };
  const shockAt = x => {
    if (!Number.isFinite(x) || x < 0 || x > 1) throw new Error('Shock position must lie on the divergent segment [0,1].');
    const shock = normalShock(root(area(x) / throatArea, 'supersonic'), { gamma });
    const downstreamCriticalArea = throatArea / shock.totalPressureRatio;
    const exitMach = root(areaRatio * shock.totalPressureRatio, 'subsonic');
    const exit = stateForMach(exitMach, pTotal * shock.totalPressureRatio);
    return { ...shock, x, area: area(x), downstreamCriticalArea, exitMach, backPressure: exit.p };
  };
  if (!Number.isFinite(shockMach) || shockMach <= 1)
    throw new Error('Choose a shock Mach placing a nonzero shock strictly inside the divergent nozzle.');
  const shockAreaRatio = areaMachRatio(shockMach, gamma);
  if (shockAreaRatio <= 1 || shockAreaRatio >= areaRatio)
    throw new Error('Choose a shock Mach placing a nonzero shock strictly inside the divergent nozzle.');
  const shockX = 2 / Math.PI * Math.asin(Math.sqrt((shockAreaRatio - 1) / (areaRatio - 1)));
  const shock = shockAt(shockX);
  const backPressureRange = { minimum: shockAt(1).backPressure, maximum: shockAt(0).backPressure };
  const shockLocationForBackPressure = pressure => {
    if (!Number.isFinite(pressure) || pressure < backPressureRange.minimum || pressure > backPressureRange.maximum)
      throw new Error('Back pressure does not put a normal shock inside this nozzle.');
    if (pressure === backPressureRange.maximum) return 0;
    if (pressure === backPressureRange.minimum) return 1;
    let lo = 0, hi = 1;
    for (let k = 0; k < 80; k++) {
      const mid = .5 * (lo + hi);
      if (mid === lo || mid === hi) break;
      if (shockAt(mid).backPressure > pressure) lo = mid;
      else hi = mid;
    }
    return .5 * (lo + hi);
  };
  const stateAt = (x, { shockSide = 'downstream' } = {}) => {
    if (!['upstream', 'downstream'].includes(shockSide)) throw new Error('Invalid side of shock.');
    const a = area(x), downstream = x > shockX || (x === shockX && shockSide === 'downstream');
    const branch = downstream || x <= 0 ? 'subsonic' : 'supersonic';
    const mach = root(a / (downstream ? shock.downstreamCriticalArea : throatArea), branch);
    const state = stateForMach(mach, pTotal * (downstream ? shock.totalPressureRatio : 1));
    return { ...state, x, area: a, branch: x === 0 ? 'sonic' : branch, massFlow: a * state.rho * state.q };
  };
  return { gamma, stagnationEnthalpy, stagnationDensity, totalPressure: pTotal,
    throatArea, areaRatio, extensionLength, domain: [-1 - extensionLength, 1 + extensionLength],
    throatX: 0, shockX, shock, chokedMass, backPressure: shock.backPressure,
    backPressureRange, area, areaSlope, machBranches, stateForMach, stateAt, shockAt, shockLocationForBackPressure };
}
