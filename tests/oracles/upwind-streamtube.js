// SPDX-License-Identifier: GPL-2.0-or-later
// Independent midpoint-polygon flux integration; no solver/kernel imports.
export function directUpwindChannelConservation(result, stagnationEnthalpy, gamma = 1.4) {
  const midpoint = (a, b) => ({ x: .5 * (a.x + b.x), y: .5 * (a.y + b.y) });
  const nx = result.nodes.length - 1, nt = result.nodes[0].length - 1;
  const h0 = Array.isArray(stagnationEnthalpy) ? stagnationEnthalpy : Array(nt).fill(stagnationEnthalpy);
  const summary = () => ({ total: [0, 0, 0, 0], external: [0, 0, 0, 0], maxLocal: [0, 0, 0, 0] });
  const biased = summary(), physical = summary(), local = [];
  let maximumMassMismatch = 0, maximumEnthalpyMismatch = 0, maximumTransportMomentumDifference = 0;
  for (let i = 1; i < nx; i++) for (let j = 0; j < nt; j++) {
    const l = [result.nodes[i - 1][j], result.nodes[i][j], result.nodes[i + 1][j]];
    const u = [result.nodes[i - 1][j + 1], result.nodes[i][j + 1], result.nodes[i + 1][j + 1]];
    const polygon = [midpoint(l[0], l[1]), l[1], midpoint(l[1], l[2]),
      midpoint(u[1], u[2]), u[1], midpoint(u[0], u[1])];
    const cell = result.cells[i - 1][j], a = result.sections[i - 1][j], b = result.sections[i][j];
    const pressure = [cell.interfacePressure.lower, cell.interfacePressure.lower, b.p,
      cell.interfacePressure.upper, cell.interfacePressure.upper, a.p];
    const fluxes = { biased: [0, 0, 0, 0], physical: [0, 0, 0, 0] };
    for (let f = 0; f < 6; f++) {
      const start = polygon[f], end = polygon[(f + 1) % 6];
      const normal = { x: end.y - start.y, y: start.x - end.x };
      const face = [0, pressure[f] * normal.x, pressure[f] * normal.y, 0], real = face.slice();
      if (f === 2 || f === 5) {
        const k = f === 5 ? 0 : 1, section = k === 0 ? a : b, sign = k === 0 ? -1 : 1;
        const c0 = midpoint(l[k], u[k]), c1 = midpoint(l[k + 1], u[k + 1]);
        const length = Math.hypot(c1.x - c0.x, c1.y - c0.y);
        const direction = { x: (c1.x - c0.x) / length, y: (c1.y - c0.y) / length };
        const mass = sign * result.massFlows[j];
        const geometricMass = section.rho * section.q * (direction.x * normal.x + direction.y * normal.y);
        maximumMassMismatch = Math.max(maximumMassMismatch, Math.abs(geometricMass - mass));
        const physicalEnthalpy = gamma / (gamma - 1) * section.p / section.rho + .5 * section.q ** 2;
        maximumEnthalpyMismatch = Math.max(maximumEnthalpyMismatch, Math.abs(physicalEnthalpy - h0[j]));
        face[0] = real[0] = mass;
        face[3] = real[3] = mass * h0[j];
        for (const [index, coordinate] of [[1, 'x'], [2, 'y']]) {
          face[index] += mass * cell.transportSpeeds[k] * direction[coordinate];
          real[index] += mass * section.q * direction[coordinate];
          maximumTransportMomentumDifference = Math.max(maximumTransportMomentumDifference, Math.abs(face[index] - real[index]));
        }
      }
      const boundary = (j === 0 && f < 2) || (j === nt - 1 && (f === 3 || f === 4))
        || (i === 1 && f === 5) || (i === nx - 1 && f === 2);
      for (const [key, values, sums] of [['biased', face, biased], ['physical', real, physical]])
        values.forEach((v, k) => { fluxes[key][k] += v; if (boundary) sums.external[k] += v; });
    }
    for (const [key, sums] of [['biased', biased], ['physical', physical]])
      fluxes[key].forEach((v, k) => { sums.total[k] += v; sums.maxLocal[k] = Math.max(sums.maxLocal[k], Math.abs(v)); });
    const transverse = { x: .5 * (polygon[5].x - polygon[0].x + polygon[3].x - polygon[2].x),
      y: .5 * (polygon[5].y - polygon[0].y + polygon[3].y - polygon[2].y) };
    local.push({ i, j, transverse, ...fluxes });
  }
  for (const value of [biased, physical]) value.internalCancellation = value.total.map((v, k) => v - value.external[k]);
  return { order: ['mass', 'xMomentum', 'yMomentum', 'totalEnthalpy'], biased, physical, local,
    maximumMassMismatch, maximumEnthalpyMismatch, maximumTransportMomentumDifference,
    physicalMinusBiased: Object.fromEntries(['total', 'external'].map(key => [key,
      physical[key].map((v, k) => v - biased[key][k])])) };
}
