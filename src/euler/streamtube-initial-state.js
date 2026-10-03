// SPDX-License-Identifier: GPL-2.0-or-later
import { streamtubeCellGeometry } from './streamtube-cell.js';

function sonicCapacityError(massFlux, mach, gamma) {
  const criticalMassFlux = isentropicSonicMassFlux({ mach, gamma });
  return Object.assign(new Error('Initial section mass/area reaches or exceeds the sonic limit.'), {
    code: 'streamtube-sonic-capacity',
    diagnostics: { mach, massFlux, criticalMassFlux, capacityRatio: massFlux / criticalMassFlux },
  });
}
// Maximum rho*q on the uniform-entropy subsonic branch. Density and
// freestream velocity are nondimensionalized to one by the body solver.
export function isentropicSonicMassFlux({ mach, gamma = 1.4 }) {
  if (!Number.isFinite(mach) || mach <= 0 || mach >= 1 || !Number.isFinite(gamma) || gamma <= 1)
    throw new Error('Invalid sonic mass-flux conditions.');
  const gm = gamma - 1, speed = Math.sqrt((2 / gm + mach * mach) / ((1 + 2 / gm) * mach * mach));
  const temperature = 1 + .5 * gm * mach * mach * (1 - speed * speed);
  return speed * temperature ** (1 / gm);
}
// Recover the subsonic isentropic branch from section mass/normal area.
// Freestream density and speed are one, as in streamtube-body.js.
export function isentropicSectionDensity({ massFlux, mach, gamma = 1.4 }) {
  if (![massFlux, mach, gamma].every(Number.isFinite) || massFlux <= 0 || mach <= 0 || mach >= 1 || gamma <= 1)
    throw new Error('Invalid isentropic section initialization.');
  const gm = gamma - 1, sonicSpeed = Math.sqrt((2 / gm + mach * mach) / ((1 + 2 / gm) * mach * mach));
  const density = q => (1 + .5 * gm * mach * mach * (1 - q * q)) ** (1 / gm);
  const sonicMassFlux = isentropicSonicMassFlux({ mach, gamma });
  if (massFlux >= sonicMassFlux) throw sonicCapacityError(massFlux, mach, gamma);
  let lower = 0, upper = sonicSpeed;
  for (let iteration = 0; iteration < 60; iteration++) {
    const q = .5 * (lower + upper);
    if (density(q) * q < massFlux) lower = q; else upper = q;
  }
  return density(.5 * (lower + upper));
}

// Geometry is adopted into a fresh research system before initializing its
// densities. The returned state still must pass system.evaluate/admissible.
export function initializeStreamtubeDensities(system, state) {
  if (system.conditions.flowModel === 'incompressible') throw new Error('Incompressible density is prescribed; gas initialization is not applicable.');
  const initial = Float64Array.from(state), { nodes, allocation } = system.decode(initial);
  const { gamma, mach } = system.conditions, { nx, tubes } = system.layout;
  // Check the entire grid before any gas inversion. Positive mean area alone
  // does not exclude a concave/folded quadrilateral. Otherwise an earlier
  // sonic section can conceal a downstream geometry defect, or a nominally
  // initialized state is rejected only when the full Euler solve starts.
  const areas = tubes.map(count => Array.from({ length: nx }, () => new Float64Array(count)));
  for (let g = 0; g < nodes.length; g++) for (let i = 1; i < nx; i++) for (let j = 0; j < tubes[g]; j++) {
    try {
      const geometry = streamtubeCellGeometry([nodes[g][i - 1][j], nodes[g][i][j], nodes[g][i + 1][j]],
        [nodes[g][i - 1][j + 1], nodes[g][i][j + 1], nodes[g][i + 1][j + 1]]);
      if (i === 1) areas[g][0][j] = geometry.normalAreas[0];
      areas[g][i][j] = geometry.normalAreas[1];
    } catch (error) { throw new Error(`Initial grid i=${i}, group=${g}, tube=${j}: ${error.message}`); }
  }
  // Report the largest capacity excess, not whichever sonic section the
  // density loop happens to encounter first. Geometry is checked first;
  // this is a constraint on the initial grid, not a physical critical Mach.
  let bottleneck = null;
  for (let g = 0; g < nodes.length; g++) for (let i = 0; i < nx; i++) for (let j = 0; j < tubes[g]; j++) {
    const massFlux = allocation.groups[g][j].massFlow / areas[g][i][j];
    if (!bottleneck || massFlux > bottleneck.massFlux) bottleneck = { i, group: g, tube: j, massFlux };
  }
  if (bottleneck && bottleneck.massFlux >= isentropicSonicMassFlux({ mach, gamma })) {
    const { i, group, tube, massFlux } = bottleneck, error = sonicCapacityError(massFlux, mach, gamma);
    error.message = `Initial section i=${i}, group=${group}, tube=${tube}: ${error.message}`;
    Object.assign(error.diagnostics, { section: { i, group, tube }, stage: 'gas-initialization' });
    throw error;
  }
  for (let g = 0; g < nodes.length; g++) for (let i = 0; i < nx; i++) for (let j = 0; j < tubes[g]; j++) {
    const area = areas[g][i][j];
    try {
      if (!(area > 0) || !Number.isFinite(area)) throw new Error('Nonpositive initial normal section area.');
      initial[system.layout.densityIndex(i, g, j)] = Math.log(isentropicSectionDensity({
        massFlux: allocation.groups[g][j].massFlow / area, mach, gamma }));
    } catch (error) { throw new Error(`Initial section i=${i}, group=${g}, tube=${j}: ${error.message}`); }
  }
  return initial;
}
