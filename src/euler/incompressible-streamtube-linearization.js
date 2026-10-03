// SPDX-License-Identifier: GPL-2.0-or-later
import { evaluateIncompressibleStreamtubeCell } from './incompressible-streamtube-cell.js';
import { linearizeStreamtubeGeometry } from './streamtube-cell.js';
const dot = (a, b) => a.x * b.x + a.y * b.y;
const cross = (a, b) => a.x * b.y - a.y * b.x;
const ddot = (a, b, da, db) => dot(da, b) + dot(a, db);
const dcross = (a, b, da, db) => cross(da, b) + cross(a, db);

export function linearizeIncompressibleStreamtubeCell(parameters) {
  const value = evaluateIncompressibleStreamtubeCell(parameters), g = value.geometry;
  const geometry = linearizeStreamtubeGeometry(parameters.lower, parameters.upper, g);
  const { massFlow: m, density: rho = 1, pressureCorrectionFactor: factor = .1 } = parameters;
  const [a, b] = value.states, { pressureCorrection: pc } = value;
  const normal = a.q * dot(g.directions[0], g.transverse) - b.q * dot(g.directions[1], g.transverse);
  const streamwise = b.q * dot(g.directions[1], g.streamwise) - a.q * dot(g.directions[0], g.streamwise);
  const crossA = cross(...g.sections), crossB = cross(g.sides.lower, g.sides.upper);
  const apply = ({ lower, upper, massFlow: dm = 0, density: drho = 0, pressureCorrectionFactor: dfactor = 0 } = {}) => {
    if (![dm, drho, dfactor].every(Number.isFinite)) throw new Error('Invalid incompressible-cell tangent.');
    const dg = geometry.apply({ lower, upper });
    const states = value.states.map((s, i) => {
      const q = s.q * (dm / m - drho / rho - dg.normalAreas[i] / g.normalAreas[i]);
      return { rho: drho, q, p: -.5 * drho * s.q ** 2 - rho * s.q * q };
    });
    const [da, db] = states, speedSum = a.q ** 2 + b.q ** 2;
    const dpc = .5 * (dfactor * rho * speedSum * g.pressureCurvature
      + factor * ((drho * speedSum + 2 * rho * (a.q * da.q + b.q * db.q)) * g.pressureCurvature + rho * speedSum * dg.pressureCurvature));
    const dn = da.q * dot(g.directions[0], g.transverse) + a.q * ddot(g.directions[0], g.transverse, dg.directions[0], dg.transverse)
      - db.q * dot(g.directions[1], g.transverse) - b.q * ddot(g.directions[1], g.transverse, dg.directions[1], dg.transverse);
    const ds = db.q * dot(g.directions[1], g.streamwise) + b.q * ddot(g.directions[1], g.streamwise, dg.directions[1], dg.streamwise)
      - da.q * dot(g.directions[0], g.streamwise) - a.q * ddot(g.directions[0], g.streamwise, dg.directions[0], dg.streamwise);
    const quotient = (numerator, derivative) => (derivative - numerator * dg.area / g.area) / g.area;
    const difference = quotient(m * normal, dm * normal + m * dn)
      + quotient(pc * crossA, dpc * crossA + pc * dcross(...g.sections, ...dg.sections));
    const sum = da.p + db.p + 2 * dpc;
    const interfacePressure = { lower: .5 * (sum - difference), upper: .5 * (sum + difference) };
    const streamwiseResidual = quotient(m * streamwise, dm * streamwise + m * ds) + db.p - da.p
      - quotient(pc * crossB, dpc * crossB + pc * dcross(g.sides.lower, g.sides.upper, dg.sides.lower, dg.sides.upper));
    if (![streamwiseResidual, dpc, ...Object.values(interfacePressure), ...states.flatMap(Object.values)].every(Number.isFinite))
      throw new Error('Nonfinite incompressible-cell derivative.');
    return { states, interfacePressure, pressureCorrection: dpc, streamwiseResidual, geometry: dg };
  };
  return { value, apply };
}
