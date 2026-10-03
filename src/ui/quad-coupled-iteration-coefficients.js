// SPDX-License-Identifier: GPL-2.0-or-later
// Read the accepted iterate; no flow/BL evaluation or change to its state.
import { quadCoupledCoefficients } from './quad-coupled-coefficients.js';
import { isentropicState } from '../potential/isentropic.js';
import { quadCoupledWakeGap } from './quad-coupled-wake-gap.js';

export function quadCoupledIterationCoefficients({ checkpoint, flow, bl, bodies,
  solverLength, referenceChord, momentReference, alpha, mach, gamma = 1.4 }) {
  const packed = checkpoint?.restart?.initialBL;
  if (checkpoint?.version !== 1 || packed?.length !== 4 * bl.stations.length
    || !(solverLength > 0) || !(bl.scale > 0)) throw new Error('Missing coupled iterate for coefficient reporting.');
  const surfaces = bl.surfaces.map(branch => {
    const group = branch.side === 'upper' ? branch.body + 1 : branch.body;
    return { element: bodies[branch.body].element, side: branch.side,
      stations: branch.ids.map(id => {
        const row = flow.undisplacedNodes[group][bl.stations[id].i];
        const point = branch.side === 'upper' ? row[0] : row.at(-1);
        return { ...point, cp: isentropicState(packed[4 * id + 3], 0, { mach, gamma }).cp };
      }) };
  });
  const wakes = bl.wakes.map(wake => {
    const id = wake.ids.at(-1);
    const finite = bodies[wake.body].trailingEdge?.kind === 'finite-base';
    return { element: bodies[wake.body].element, stations: [{ index: id,
      theta: packed[4 * id + 1] * bl.scale * solverLength,
      deltaStar: packed[4 * id + 2] * bl.scale * solverLength, ue: packed[4 * id + 3],
      ...(finite ? { wakeGap: quadCoupledWakeGap({ body: bodies[wake.body], bodyIndex: wake.body,
        flow, station: bl.stations[id], solverLength }) } : {}) }] };
  }).sort((a, b) => a.element - b.element);
  const stagnation = [];
  bodies.forEach((body, b) => { stagnation[body.element] = flow.undisplacedNodes[b + 1][body.leadingIndex][0]; });
  return quadCoupledCoefficients({ surfaces, wakes, stagnation, bodies, referenceChord, momentReference, alpha, mach, gamma });
}
