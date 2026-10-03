// SPDX-License-Identifier: GPL-2.0-or-later
// Transfer physical geometry, captured masses, stagnation locations and
// multipoles into a fresh system. Density is recomputed for the target gas;
// an incompressible grid is never relabeled as a compressible solution.
import { createStreamtubeBodySystem } from '../streamtube-body.js';
import { initializeStreamtubeDensities } from '../streamtube-initial-state.js';

export function initializeStreamtubeBodyFromGrid(input, source) {
  const system = createStreamtubeBodySystem(input), { globals, elements } = system.layout;
  const state = system.initial.slice(), base = system.decode(state), { lengthScale, massScale, flowModel } = system.conditions;
  if (!source || !Array.isArray(source.captured) || source.captured.length !== elements + 2 || !source.captured.every(Number.isFinite)
    || !Array.isArray(source.stagnation) || source.stagnation.length !== elements || !source.stagnation.every(Number.isFinite)
    || !source.strengths || !['circulation', 'source', 'doubletX', 'doubletY'].every(k => Number.isFinite(source.strengths[k])))
    throw new Error('Supply physical body-grid capture, stagnation and farfield data for restart.');
  for (let k = 0; k < source.captured.length; k++) {
    const col = k === 0 || k === source.captured.length - 1 ? null : globals.capture[k - 1];
    const difference = source.captured[k] - base.captured[k];
    if (col === null && Math.abs(difference) > 1e-12 * massScale)
      throw new Error('Restart changes a prescribed outer or primary captured-mass level.');
    if (col !== null) state[col] = difference / massScale;
  }
  for (let b = 0; b < elements; b++) {
    const col = globals.stagnation[b], difference = source.stagnation[b] - base.stagnation[b];
    if (col === null && Math.abs(difference) > 1e-12 * system.curves[b].length)
      throw new Error('Restart changes a prescribed stagnation location.');
    if (col !== null) state[col] = difference / system.curves[b].length;
  }
  for (const key of ['circulation', 'source', 'doubletX', 'doubletY'])
    state[globals[key]] = source.strengths[key] / (key.startsWith('doublet') ? lengthScale ** 2 : lengthScale);
  const geometry = system.adoptGeometry(state, source.nodes);
  const initial = flowModel === 'incompressible' ? geometry : initializeStreamtubeDensities(system, geometry);
  const evaluated = system.evaluate(initial);
  return { input: structuredClone(input), system, initial, nodes: evaluated.nodes,
    diagnostics: { initialFlow: evaluated.diagnostics, densityInitialization: flowModel === 'incompressible' ? 'constant density' : 'target-Mach subsonic isentropic mass/area inversion' },
    status: 'research grid restart; target flow still requires convergence' };
}
