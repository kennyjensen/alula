// SPDX-License-Identifier: GPL-2.0-or-later
// Diagnostics only: keep physical Euler gas separate from the historical
// BL common isentrope. No interpolation of shock entropy into a BL equation.
export function streamtubePhysicalSectionGas(section, { gamma, h0, rhoTotal }) {
  const { rho, q, p, enthalpy } = section;
  if (![gamma, h0, rhoTotal, rho, q, p, enthalpy].every(Number.isFinite)
    || gamma <= 1 || h0 <= 0 || rhoTotal <= 0 || rho <= 0 || q <= 0 || p <= 0 || enthalpy <= 0)
    throw new Error('Physical Euler edge diagnostics require positive finite compressible gas.');
  const logRhoTotal = Math.log(rho) + Math.log(h0 / enthalpy) / (gamma - 1);
  const entropyOverR = Math.log(rhoTotal) - logRhoTotal;
  return { rho, q, p, enthalpy, machSquared: q * q / ((gamma - 1) * enthalpy), logRhoTotal, entropyOverR,
    stagnationPressureRatio: Math.exp(-entropyOverR),
    // Compare at THIS section's physical speed, not an interpolated BL speed.
    commonIsentropeRelativeDensityError: Math.expm1(entropyOverR) };
}

export function streamtubeBLEdgeThermodynamics({ flow, euler, bl, states }) {
  const { gamma, h0, rhoTotal, pInf, flowModel } = euler.conditions, { layout } = euler;
  if (flowModel !== 'compressible' || states?.length !== bl.stations.length)
    throw new Error('BL edge gas diagnostics require a complete compressible coupled state.');
  let maxDensityMismatch = 0, maxWakeBankEntropyDifference = 0;
  const stations = bl.stations.map(s => {
    const ue = states[s.id].ue, h = h0 - .5 * ue * ue;
    if (!Number.isFinite(ue) || ue <= 0 || !(h > 0)) throw new Error('BL gas has nonpositive speed or static enthalpy.');
    const i = Math.min(s.i, layout.nx - 1);
    const banks = (s.kind === 'surface' ? [s.side] : ['upper', 'lower']).map(side => {
      const group = side === 'lower' ? s.body : s.body + 1, tube = side === 'lower' ? layout.tubes[group] - 1 : 0;
      const cell = flow.cells[i - 1][group][tube];
      const sections = cell.states.map((physical, k) => {
        const gas = streamtubePhysicalSectionGas(physical, euler.conditions);
        maxDensityMismatch = Math.max(maxDensityMismatch, Math.abs(gas.commonIsentropeRelativeDensityError));
        return { index: i - 1 + k, ...gas };
      });
      const interfacePressure = side === 'lower' ? cell.interfacePressure.upper : cell.interfacePressure.lower;
      if (!Number.isFinite(interfacePressure) || interfacePressure <= 0) throw new Error('Invalid Euler interface pressure.');
      return { side, group, tube, interfaceIndex: i, interfacePressure, cp: 2 * (interfacePressure - pInf), sections };
    });
    if (banks.length === 2) for (let k = 0; k < 2; k++)
      maxWakeBankEntropyDifference = Math.max(maxWakeBankEntropyDifference,
        Math.abs(banks[0].sections[k].entropyOverR - banks[1].sections[k].entropyOverR));
    return { id: s.id, kind: s.kind, body: s.body, ...(s.kind === 'surface' ? { side: s.side } : {}),
      boundaryLayer: { ue, rho: rhoTotal * (h / h0) ** (1 / (gamma - 1)), enthalpy: h,
        machSquared: ue * ue / ((gamma - 1) * h) }, banks };
  });
  return { model: 'historical-common-isentrope', stations, maxDensityMismatch, maxWakeBankEntropyDifference,
    interpretation: 'BL density uses the common reference isentrope. Euler section gas and interface pressure are retained independently. Density errors compare gas at the same section speed; adjacent sections are not claimed to be an exact wall gas state. Wake banks remain separate.' };
}
