// SPDX-License-Identifier: GPL-2.0-or-later
import { streamtubeCellGeometry } from '../streamtube-cell.js';

// Giles's ISET listing (thesis p. 172) initializes uniform stagnation density.
// Here the minimum signed section area also supplies a sonic-capacity scale
// for each FREE mass unknown. This is only a Newton guess: neither choking,
// a sonic branch, shock position nor downstream entropy is prescribed.
export function initializeUpwindStreamtubeChannel(system) {
  const { nx, nt, initial: original, conditions, massIndex } = system;
  if (!Number.isInteger(nx) || nx < 3 || !Number.isInteger(nt) || nt < 1
    || typeof massIndex !== 'function' || !original || original.length !== system.n)
    throw new Error('Supply a complete upwind streamtube channel for initialization.');
  const { stagnationDensity: rhoTotal, stagnationEnthalpy: h0, gamma,
    referenceDensity, massFlows: referenceMass } = conditions;
  // decode constructs geometry without evaluating the possibly inadmissible
  // arbitrary reference-mass guess. No back-pressure data enter this seed.
  const { nodes } = system.decode(original);
  const normalAreas = Array.from({ length: nx }, () => Array(nt));
  for (let i = 1; i < nx; i++) for (let j = 0; j < nt; j++) {
    const lower = [nodes[i - 1][j], nodes[i][j], nodes[i + 1][j]];
    const upper = [nodes[i - 1][j + 1], nodes[i][j + 1], nodes[i + 1][j + 1]];
    const geometry = streamtubeCellGeometry(lower, upper);
    if (i === 1) normalAreas[0][j] = geometry.normalAreas[0];
    normalAreas[i][j] = geometry.normalAreas[1];
  }
  const initial = Float64Array.from(original);
  const tubes = Array.from({ length: nt }, (_, j) => {
    let minimumNormalArea = Infinity, minimumAreaSection = -1;
    for (let i = 0; i < nx; i++) if (normalAreas[i][j] < minimumNormalArea) {
      minimumNormalArea = normalAreas[i][j]; minimumAreaSection = i;
    }
    const sonicDensity = rhoTotal[j] * (2 / (gamma + 1)) ** (1 / (gamma - 1));
    const sonicSpeed = Math.sqrt(2 * (gamma - 1) / (gamma + 1) * h0[j]);
    const sonicMassFlux = sonicDensity * sonicSpeed;
    const capacityMassFlow = minimumNormalArea * sonicMassFlux;
    if (![minimumNormalArea, sonicDensity, sonicSpeed, capacityMassFlow].every(v => Number.isFinite(v) && v > 0))
      throw new Error('Invalid signed section capacity in upwind-channel initialization.');
    initial[massIndex(j)] = Math.log(capacityMassFlow / referenceMass[j]);
    const density = Math.log(rhoTotal[j] / referenceDensity);
    for (let i = 0; i < nx; i++) initial[i * nt + j] = density;
    return { tube: j, minimumNormalArea, minimumAreaSection, sonicDensity, sonicSpeed,
      sonicMassFlux, capacityMassFlow };
  });
  // The complete target residual enforces positive geometry, thermodynamics
  // and subsonic inlet/outlet. Its nonzero residual is retained, never hidden.
  const flow = system.evaluate(initial);
  let initialMaxMach = 0;
  for (const row of flow.sections) for (const state of row)
    initialMaxMach = Math.max(initialMaxMach, Math.sqrt(state.machSquared));
  const residual = flow.residual.reduce((largest, r) => Math.max(largest, Math.abs(r)), 0);
  if (![initialMaxMach, residual].every(Number.isFinite)) throw new Error('Nonfinite upwind-channel startup diagnostic.');
  return { initial, flow, diagnostics: {
    method: 'stagnation-density-sonic-capacity', initialGuessOnly: true, targetEquationsUnchanged: true,
    massFlowsAreUnknowns: true, minimumNormalArea: Math.min(...tubes.map(t => t.minimumNormalArea)),
    capacityMassFlows: tubes.map(t => t.capacityMassFlow), initialMaxMach, residual, tubes,
  } };
}
