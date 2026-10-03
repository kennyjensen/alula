// SPDX-License-Identifier: GPL-2.0-or-later
// Exact M -> 0 limit of the smooth, constant-total-pressure intrinsic cell.
// Pressure is measured relative to that common total pressure. It may be
// negative; this is not a perfect-gas state or a relaxation of Euler guards.
import { streamtubeCellGeometry } from './streamtube-cell.js';

const dot = (a, b) => a.x * b.x + a.y * b.y;
const cross = (a, b) => a.x * b.y - a.y * b.x;

export function evaluateIncompressibleStreamtubeCell({ lower, upper, massFlow, density = 1, pressureCorrectionFactor = .1 }) {
  if (![massFlow, density, pressureCorrectionFactor].every(Number.isFinite)
    || massFlow <= 0 || density <= 0 || pressureCorrectionFactor < 0)
    throw new Error('Invalid incompressible streamtube controls.');
  const geometry = streamtubeCellGeometry(lower, upper);
  const states = geometry.normalAreas.map(area => {
    const q = massFlow / (density * area);
    return { rho: density, q, p: -.5 * density * q * q };
  });
  const [a, b] = states, { directions, transverse, streamwise, sections, sides, area } = geometry;
  // gamma*pMean*MMean^2 -> rho*(q_a^2+q_b^2)/2. The curvature
  // correction remains finite in the incompressible limit.
  const pressureCorrection = pressureCorrectionFactor * .5 * density * (a.q * a.q + b.q * b.q) * geometry.pressureCurvature;
  const pressureDifference = (massFlow * (a.q * dot(directions[0], transverse) - b.q * dot(directions[1], transverse))
    + pressureCorrection * cross(sections[0], sections[1])) / area;
  const pressureSum = a.p + b.p + 2 * pressureCorrection;
  const interfacePressure = { lower: .5 * (pressureSum - pressureDifference), upper: .5 * (pressureSum + pressureDifference) };
  const streamwiseResidual = massFlow * (b.q * dot(directions[1], streamwise) - a.q * dot(directions[0], streamwise)) / area
    + b.p - a.p - pressureCorrection * cross(sides.lower, sides.upper) / area;
  if (![pressureCorrection, streamwiseResidual, ...Object.values(interfacePressure), ...states.flatMap(Object.values)].every(Number.isFinite))
    throw new Error('Nonfinite incompressible streamtube cell.');
  return { geometry, states, pressureCorrection, interfacePressure, streamwiseResidual,
    pressureReference: 'relative to the common incompressible total pressure' };
}
