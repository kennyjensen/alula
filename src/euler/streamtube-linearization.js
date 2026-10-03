// SPDX-License-Identifier: GPL-2.0-or-later
// Exact directional chain rule for the intrinsic Euler cell, including
// optional shared-section transport speeds supplied by an upstream stencil.
// Both coordinates of every node participate. Mass flow and total enthalpy
// are differentiated here even though the verification channel prescribes
// them. The reference value kernel remains independently implemented.
import { evaluateStreamtubeCell, linearizeStreamtubeGeometry } from './streamtube-cell.js';

const cross = (a, b) => a.x * b.y - a.y * b.x;
const dot = (a, b) => a.x * b.x + a.y * b.y;
const dcross = (a, b, da, db) => cross(da, b) + cross(a, db);
const ddot = (a, b, da, db) => dot(da, b) + dot(a, db);
const zeroRow = Object.freeze(Array.from({ length: 3 }, () => Object.freeze({ x: 0, y: 0 })));
const zeroPair = Object.freeze([0, 0]);
function finiteRow(row) {
  if (!Array.isArray(row) || row.length !== 3) return false;
  for (let i = 0; i < 3; i++) if (!Number.isFinite(row[i].x) || !Number.isFinite(row[i].y)) return false;
  return true;
}
const finitePair = a => Array.isArray(a) && a.length === 2 && Number.isFinite(a[0]) && Number.isFinite(a[1]);
const finiteState = a => Number.isFinite(a.rho) && Number.isFinite(a.q) && Number.isFinite(a.p)
  && Number.isFinite(a.enthalpy) && Number.isFinite(a.machSquared);
const sectionDerivative = (partial, rho, mass, area, enthalpy) => partial.density * rho + partial.massFlow * mass
  + partial.normalArea * area + partial.stagnationEnthalpy * enthalpy;

export function linearizeStreamtubeCell(parameters) {
  const value = evaluateStreamtubeCell(parameters);
  const lower = parameters.lower.map(p => ({ ...p })), upper = parameters.upper.map(p => ({ ...p }));
  const { massFlow: m, gamma = 1.4, pressureCorrectionFactor: factor = .1 } = parameters;
  const { geometry: g, states: [a, b], pressureCorrection: pc } = value;
  const geometry = linearizeStreamtubeGeometry(lower, upper, g);
  const pMean = .5 * (a.p + b.p), machMean = .5 * (a.machSquared + b.machSquared);
  const { streamwise: s, transverse: n, area, directions: t, sections: section, sides } = g;
  const crossA = cross(...section), crossB = cross(sides.lower, sides.upper);
  const biased = parameters.transportSpeeds !== undefined;
  const [qa, qb] = biased ? value.transportSpeeds : [a.q, b.q];
  const normalTerm = qa * dot(t[0], n) - qb * dot(t[1], n);
  const streamwiseTerm = qb * dot(t[1], s) - qa * dot(t[0], s);

  const apply = ({ lower: dl = zeroRow, upper: du = zeroRow, densities: drho = zeroPair,
    massFlow: dm = 0, stagnationEnthalpy: dh0 = 0, pressureCorrectionFactor: dfactor = 0,
    transportSpeeds: dtransport = zeroPair } = {}) => {
    if (!finiteRow(dl) || !finiteRow(du) || !finitePair(drho) || !finitePair(dtransport)
      || !Number.isFinite(dm) || !Number.isFinite(dh0) || !Number.isFinite(dfactor)) throw new Error('Invalid intrinsic-cell tangent.');
    if (!biased && (dtransport[0] !== 0 || dtransport[1] !== 0)) throw new Error('Transport-speed tangent requires explicit transport speeds.');
    const { directions: dt, sections: dsection, sides: dsides, streamwise: ds, transverse: dn,
      area: darea, normalAreas: dnormalAreas, pressureCurvature: dcurvature, streamwiseLengths: dlengths } = geometry.apply({ lower: dl, upper: du });
    const states = new Array(2);
    for (let i = 0; i < 2; i++) {
      const state = i === 0 ? a : b;
      const q = sectionDerivative(state.derivatives.q, drho[i], dm, dnormalAreas[i], dh0);
      const p = sectionDerivative(state.derivatives.p, drho[i], dm, dnormalAreas[i], dh0), enthalpy = dh0 - state.q * q;
      const machSquared = (2 * state.q * q - state.machSquared * (gamma - 1) * enthalpy) / ((gamma - 1) * state.enthalpy);
      states[i] = { rho: drho[i], q, p, enthalpy, machSquared };
    }
    const [da, db] = states, dpMean = .5 * (da.p + db.p), dmachMean = .5 * (da.machSquared + db.machSquared);
    const [dqa, dqb] = biased ? dtransport : [da.q, db.q];
    // Differentiate the active smooth branch. At mean M²=1 this switch has
    // no unique two-sided derivative; the channel excludes sonic states.
    const dpc = machMean < 1 ? gamma * (
      dfactor * pMean * machMean * (1 - machMean) * g.pressureCurvature
      + factor * ((dpMean * machMean * (1 - machMean) + pMean * (1 - 2 * machMean) * dmachMean) * g.pressureCurvature
        + pMean * machMean * (1 - machMean) * dcurvature)) : 0;
    const dnormalTerm = dqa * dot(t[0], n) + qa * ddot(t[0], n, dt[0], dn)
      - dqb * dot(t[1], n) - qb * ddot(t[1], n, dt[1], dn);
    const dstreamwiseTerm = dqb * dot(t[1], s) + qb * ddot(t[1], s, dt[1], ds)
      - dqa * dot(t[0], s) - qa * ddot(t[0], s, dt[0], ds);
    const dquotient = (numerator, dnumerator) => (dnumerator - numerator * darea / area) / area;
    const dnormal = dquotient(m * normalTerm, dm * normalTerm + m * dnormalTerm);
    const dpressureDifference = dnormal + dquotient(pc * crossA, dpc * crossA + pc * dcross(...section, ...dsection));
    const dpressureSum = da.p + db.p + 2 * dpc;
    const interfacePressure = { lower: .5 * (dpressureSum - dpressureDifference), upper: .5 * (dpressureSum + dpressureDifference) };
    const streamwiseResidual = dquotient(m * streamwiseTerm, dm * streamwiseTerm + m * dstreamwiseTerm)
      + db.p - da.p - dquotient(pc * crossB, dpc * crossB + pc * dcross(sides.lower, sides.upper, dsides.lower, dsides.upper));
    const entropyJump = (db.enthalpy / b.enthalpy - da.enthalpy / a.enthalpy) / (gamma - 1) - db.rho / b.rho + da.rho / a.rho;
    const artificialEnthalpies = biased ? [dh0 - qa * dqa, dh0 - qb * dqb] : undefined;
    const artificialEntropyJump = biased ? gamma / (gamma - 1) *
      (artificialEnthalpies[1] / value.artificialEnthalpies[1] - artificialEnthalpies[0] / value.artificialEnthalpies[0])
      - db.p / b.p + da.p / a.p : entropyJump;
    const isentropicResidual = -dpMean * (biased ? value.artificialEntropyJump : value.entropyJump) - pMean * artificialEntropyJump;
    if (!Number.isFinite(streamwiseResidual) || !Number.isFinite(entropyJump) || !Number.isFinite(isentropicResidual)
      || !Number.isFinite(dpc) || !Number.isFinite(interfacePressure.lower) || !Number.isFinite(interfacePressure.upper)
      || !finiteState(da) || !finiteState(db))
      throw new Error('Nonfinite intrinsic-cell derivative.');
    return { streamwiseResidual, entropyJump, isentropicResidual, interfacePressure, pressureCorrection: dpc, states,
      ...(biased ? { transportSpeeds: [...dtransport], artificialEnthalpies, artificialEntropyJump } : {}),
      geometry: { sections: dsection, directions: dt, sides: dsides, streamwise: ds, transverse: dn,
        area: darea, normalAreas: dnormalAreas, pressureCurvature: dcurvature, streamwiseLengths: dlengths } };
  };
  return { value, apply };
}
