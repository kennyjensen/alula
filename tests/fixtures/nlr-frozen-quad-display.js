// SPDX-License-Identifier: GPL-2.0-or-later
// Presentation fixture reconstructed only from the saved coupled seed.
// No flow constructor, initializer, residual, Jacobian or Newton operation.
import fs from 'node:fs';
import { createSurfaceContourCurve } from '../../src/geometry/contour-topology.js';

export const seedPath = new URL('../../docs/nlr-finite-base/coupled-seed-state.json', import.meta.url);
export function nlrFrozenQuadDisplay() {
  const saved = JSON.parse(fs.readFileSync(seedPath, 'utf8'));
  const { input, options, initialBL, initialEuler } = saved.restart;
  const stations = [], surfaces = [], wakes = [], nx = input.outerLower.length - 1;
  const append = metadata => {
    const id = stations.length;
    stations.push({ ...metadata, id, ...saved.layers.states[id] }); return id;
  };
  input.bodies.forEach((body, b) => {
    const curve = createSurfaceContourCurve(body.points, body);
    for (const side of ['upper', 'lower']) {
      const index = surfaces.length, transition = options.transitionState[index], ids = [];
      for (let i = body.leadingIndex + 1; i <= body.trailingIndex; i++) {
        const j = i - body.leadingIndex - 1;
        const regime = j === 0 && transition === 0 ? 'leading-transition' : j === 0 ? 'similarity'
          : j < transition ? 'laminar' : j === transition ? 'transition' : 'turbulent';
        ids.push(append({ kind: 'surface', body: b, side, i, k: j + 1, regime }));
      }
      const trip = options.tripFractions[b][side === 'upper' ? 0 : 1];
      surfaces.push({ body: b, side, ids, transition, prescribedTrip: trip < 1,
        tripParameter: curve.branch(side, trip, saved.outer.stagnation[b]).parameter });
    }
    const ids = [];
    for (let i = body.trailingIndex; i <= nx; i++) ids.push(append({ kind: 'wake', body: b, i,
      k: i - body.trailingIndex, regime: i === body.trailingIndex ? 'trailing-edge' : 'wake' }));
    wakes.push({ body: b, ids });
  });
  if (stations.length !== saved.layers.states.length) throw new Error('Saved NLR station topology changed.');
  return { input: saved.caseInput, saved, raw: {
    model: 'research-streamtube-euler-bl', physicalAcceptance: false, converged: false,
    reason: 'frozen initialization; no coupled Newton updates', alpha: input.alpha, mach: input.mach,
    ...saved.settings.normalization, solverInput: input, restart: saved.restart,
    checkpoint: { version: 1, restart: saved.restart }, x: [...initialEuler.x, ...initialBL],
    flow: saved.outer, mesh: saved.mesh, families: saved.families, history: [{ iteration: 0 }],
    boundaryLayer: { stations, surfaces, wakes, transitions: saved.layers.transitions,
      transitionState: options.transitionState },
    conditions: { ...options, mach: input.mach },
    initialization: { ...saved.initialization, euler: { gridSmoothing: undefined } },
  } };
}
