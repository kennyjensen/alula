// SPDX-License-Identifier: GPL-2.0-or-later
// Sample the current observed state without evaluating flow or BL equations.
// Match the terminal Cp plot: wetted surface stations only, in element order.
import { isentropicState } from '../potential/isentropic.js';

export function quadCoupledIterationPressure({ checkpoint, flow, bl, bodies, elements, referenceChord, mach }) {
  const restart = checkpoint?.restart, input = restart?.input, packed = restart?.initialBL;
  const actualMach = input?.mach, gamma = input?.gamma ?? 1.4;
  if (checkpoint?.version !== 1 || !restart?.options || !packed
    || !Number.isFinite(referenceChord) || referenceChord <= 0
    || !Number.isFinite(actualMach) || actualMach <= 0 || actualMach >= 1
    || !Number.isFinite(gamma) || gamma <= 1 || mach !== undefined && mach !== actualMach
    || !Array.isArray(bl?.stations) || packed.length !== 4 * bl.stations.length
    || !Array.isArray(bl?.surfaces) || !Array.isArray(bodies) || !Array.isArray(elements)
    || bodies.length !== elements.length || bodies.length !== input.bodies?.length
    || bl.surfaces.length !== 2 * bodies.length)
    throw new Error('Missing or mismatched coupled pressure frame.');
  const historical = restart.options.blThermodynamics === 'historical-common-isentrope';
  const pInf = 1 / (gamma * actualMach * actualMach), used = new Set();
  const samePoint = (a, b) => a && b && Number.isFinite(a.x) && Number.isFinite(a.y) && a.x === b.x && a.y === b.y;
  const output = Array(elements.length);
  for (let body = 0; body < bodies.length; body++) {
    const b = bodies[body], saved = input.bodies[body], element = b.element;
    if (!Number.isInteger(element) || element < 0 || element >= elements.length || used.has(element)
      || b.leadingIndex !== saved.leadingIndex || b.trailingIndex !== saved.trailingIndex || element !== saved.element)
      throw new Error('Pressure element mapping differs from the current checkpoint.');
    used.add(element);
    const sides = {};
    for (const side of ['upper', 'lower']) {
      const branches = bl.surfaces.filter(s => s.body === body && s.side === side), group = side === 'upper' ? body + 1 : body;
      if (branches.length !== 1 || branches[0].ids.length !== b.trailingIndex - b.leadingIndex)
        throw new Error('Pressure surface stations differ from the current grid.');
      const bank = row => side === 'upper' ? row?.[0] : row?.at(-1);
      sides[side] = branches[0].ids.map((id, j) => {
        const station = bl.stations[id], i = b.leadingIndex + 1 + j;
        if (station?.id !== id || station.body !== body || station.side !== side || station.i !== i)
          throw new Error('Pressure station metadata belongs to another grid.');
        const point = bank(flow?.undisplacedNodes?.[group]?.[i]);
        if (!samePoint(point, bank(restart.initialEuler?.undisplacedNodes?.[group]?.[i]))
          || !samePoint(bank(flow?.nodes?.[group]?.[i]), bank(restart.initialEuler?.nodes?.[group]?.[i])))
          throw new Error('Pressure flow coordinates differ from the current checkpoint.');
        let cp;
        if (historical) {
          const cells = flow.cells?.[i - 1]?.[group];
          const pressure = side === 'upper' ? cells?.[0]?.interfacePressure?.lower : cells?.at(-1)?.interfacePressure?.upper;
          if (!Number.isFinite(pressure) || pressure <= 0) throw new Error('Current solid pressure is unavailable.');
          cp = 2 * (pressure - pInf);
        } else {
          const ue = packed[4 * id + 3];
          if (!Number.isFinite(ue) || ue <= 0) throw new Error('Current BL edge speed is unavailable.');
          cp = isentropicState(ue, 0, { mach: actualMach, gamma }).cp;
        }
        if (!Number.isFinite(cp)) throw new Error('Nonfinite current pressure coefficient.');
        return { x: point.x, y: point.y, cp };
      });
    }
    output[element] = { name: elements[element].name ?? `Element ${element + 1}`, cp: [...sides.upper.reverse(), ...sides.lower] };
  }
  return { elements: output, referenceChord, actualMach,
    pressureKind: historical ? 'physical Euler interface pressure on the undisplaced solid contour'
      : 'isentropic BL-edge pressure from solved edge speed', provisional: true, physicalAcceptance: false };
}
