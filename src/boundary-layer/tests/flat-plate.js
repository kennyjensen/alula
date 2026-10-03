// SPDX-License-Identifier: GPL-2.0-or-later
// Research harness, NOT an airfoil viscous solver. Incompressible, laminar,
// zero pressure gradient, no transition. Log integral equations specialize
// XFOIL BLDIF's momentum/shape equations; endpoint trapezoidal quadrature.
import { cfl, dil, hsl } from './closures.js';
import { solveNewton } from '../../numerics/newton.js';

export function solveFlatPlate({ reynolds = 1e6, stations = 41, initialFactor = 1.4 } = {}) {
  if (!(reynolds > 0) || !Number.isFinite(reynolds) || !Number.isInteger(stations) || stations < 3 || stations > 101
    || !Number.isFinite(initialFactor) || initialFactor <= 0) throw new Error('Invalid flat-plate test parameters.');
  const s = Float64Array.from({ length: stations }, (_, i) => 10 ** (-4 + 4 * i / (stations - 1)));
  const scale = Math.sqrt(reynolds);
  const initial = new Float64Array(2 * stations);
  for (let i = 0; i < stations; i++) {
    initial[2 * i] = Math.log(initialFactor * 0.664 * Math.sqrt(s[i]));
    initial[2 * i + 1] = 2.8;
  }
  const stateAt = (state, i) => {
    const theta = Math.exp(state[2 * i]) / scale;
    const h = state[2 * i + 1];
    const rt = reynolds * theta;
    return { theta, h, hs: hsl(h).hs, cf: cfl(h, rt).cf, di: dil(h, rt).di };
  };
  const residual = state => {
    const r = new Float64Array(2 * stations);
    let previous = stateAt(state, 0);
    // Leading similarity station: d(log theta)/d(log x)=1/2, dH*/dx=0.
    r[0] = 0.5 - 0.5 * previous.cf * s[0] / previous.theta;
    r[1] = (0.5 * previous.cf - previous.di) * s[0] / previous.theta;
    for (let i = 1; i < stations; i++) {
      const next = stateAt(state, i);
      const xlog = Math.log(s[i] / s[i - 1]);
      const friction = 0.5 * (previous.cf * s[i - 1] / previous.theta + next.cf * s[i] / next.theta);
      const dissipation = 0.5 * (previous.di * s[i - 1] / previous.theta + next.di * s[i] / next.theta);
      r[2 * i] = Math.log(next.theta / previous.theta) - 0.5 * xlog * friction;
      r[2 * i + 1] = Math.log(next.hs / previous.hs) + xlog * (0.5 * friction - dissipation);
      previous = next;
    }
    return r;
  };
  const solved = solveNewton({ residual, initial,
    admissible: state => state.every((v, i) => i % 2 ? v > 1.1 && v < 4.3 : v > -30 && v < 10) });
  return { ...solved, model: 'laminar-flat-plate-global-newton', reynolds,
    stations: Array.from(s, (x, i) => ({ x, ...stateAt(solved.x, i) })) };
}
