// Continuum channel cases and direct flux integration. No intrinsic-cell or
// channel-solver imports: these oracles do not reuse the reduced equations.
import { isentropicNozzle } from './euler.js';

export function channelGas(mach = .3, gamma = 1.4) {
  const pressure = 1 / (gamma * mach ** 2);
  return { gamma, referencePressure: pressure,
    stagnationEnthalpy: gamma / (gamma - 1) * pressure + .5,
    stagnationDensity: (1 + .5 * (gamma - 1) * mach ** 2) ** (1 / (gamma - 1)) };
}

export function nozzleChannel(nx, nt = 3) {
  const parameters = { mach: .3, height: .08, contraction: .12, length: 2 };
  const x = Array.from({ length: nx + 1 }, (_, i) => parameters.length * i / nx);
  return { ...channelGas(), x, lower: x.map(() => 0), upper: x.map(x => isentropicNozzle(x, parameters).area),
    massFlows: Array(nt).fill(parameters.height / nt),
    exact: x => isentropicNozzle(x, parameters), parameters };
}

export function vortexChannel(nx, nt) {
  const mach = .3, gas = channelGas(mach), { gamma, referencePressure } = gas;
  const exact = radius => {
    const q = 1 / radius, temperature = 1 + .5 * (gamma - 1) * mach ** 2 * (1 - q * q);
    return { q, rho: temperature ** (1 / (gamma - 1)), p: referencePressure * temperature ** (gamma / (gamma - 1)) };
  };
  const x = Array.from({ length: nx + 1 }, (_, i) => -.4 + .8 * i / nx);
  const radii = Array.from({ length: nt + 1 }, (_, j) => 1.25 - .25 * j / nt);
  const ordinate = (x, radius) => -Math.sqrt(radius * radius - x * x);
  // Simpson quadrature of the continuum rho*q dr, independent of the
  // discrete normal area and midpoint state. All tubes share inlet entropy.
  const massFlows = radii.slice(1).map((radius, j) => {
    const n = 200, dr = (radii[j] - radius) / n;
    let sum = 0;
    for (let i = 0; i <= n; i++) {
      const s = exact(radius + i * dr);
      sum += (i === 0 || i === n ? 1 : i % 2 ? 4 : 2) * s.rho * s.q;
    }
    return sum * dr / 3;
  });
  const slopes = i => radii.slice(1, -1).map(r => (ordinate(x[i + 1], r) - ordinate(x[i], r)) / (x[i + 1] - x[i]));
  return { ...gas, x, lower: x.map(x => ordinate(x, radii[0])), upper: x.map(x => ordinate(x, radii.at(-1))),
    massFlows, inletSlopes: slopes(0), outletSlopes: slopes(nx - 1), exact, radii, ordinate };
}

export function vortexErrors(result, oracle) {
  let pressureMax = 0, positionMax = 0;
  for (let i = 0; i < result.sections.length; i++) for (let j = 0; j < result.sections[i].length; j++) {
    const x = (oracle.x[i] + oracle.x[i + 1]) / 2;
    const y = (result.nodes[i][j].y + result.nodes[i + 1][j].y + result.nodes[i][j + 1].y + result.nodes[i + 1][j + 1].y) / 4;
    pressureMax = Math.max(pressureMax, Math.abs(result.sections[i][j].p - oracle.exact(Math.hypot(x, y)).p));
  }
  for (let i = 0; i < result.nodes.length; i++) for (let j = 1; j < result.nodes[i].length - 1; j++)
    positionMax = Math.max(positionMax, Math.abs(result.nodes[i][j].y - oracle.ordinate(oracle.x[i], oracle.radii[j])));
  return { pressureMax, positionMax };
}

export function directChannelConservation(result, gamma = 1.4, { flowModel = 'compressible' } = {}) {
  if (!['compressible', 'incompressible'].includes(flowModel)) throw new Error('Unknown conservation-oracle flow model.');
  const enthalpyFactor = flowModel === 'incompressible' ? 1 : gamma / (gamma - 1);
  const midpoint = (a, b) => ({ x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 });
  const maxLocal = [0, 0, 0, 0], total = [0, 0, 0, 0], external = [0, 0, 0, 0];
  const nx = result.nodes.length - 1, nt = result.nodes[0].length - 1;
  for (let i = 1; i < nx; i++) for (let j = 0; j < nt; j++) {
    const l = [result.nodes[i - 1][j], result.nodes[i][j], result.nodes[i + 1][j]];
    const u = [result.nodes[i - 1][j + 1], result.nodes[i][j + 1], result.nodes[i + 1][j + 1]];
    const polygon = [midpoint(l[0], l[1]), l[1], midpoint(l[1], l[2]), midpoint(u[1], u[2]), u[1], midpoint(u[0], u[1])];
    const pi = result.cells[i - 1][j].interfacePressure, a = result.sections[i - 1][j], b = result.sections[i][j];
    const pressures = [pi.lower, pi.lower, b.p, pi.upper, pi.upper, a.p], flux = [0, 0, 0, 0];
    for (let f = 0; f < 6; f++) {
      const first = polygon[f], next = polygon[(f + 1) % 6], normal = { x: next.y - first.y, y: first.x - next.x };
      const face = [0, pressures[f] * normal.x, pressures[f] * normal.y, 0];
      if (f === 2 || f === 5) {
        const k = f === 5 ? 0 : 1, s = f === 5 ? a : b;
        const c0 = midpoint(l[k], u[k]), c1 = midpoint(l[k + 1], u[k + 1]);
        const length = Math.hypot(c1.x - c0.x, c1.y - c0.y);
        const vx = s.q * (c1.x - c0.x) / length, vy = s.q * (c1.y - c0.y) / length;
        face[0] = s.rho * (vx * normal.x + vy * normal.y);
        face[1] += face[0] * vx; face[2] += face[0] * vy;
        face[3] = face[0] * (enthalpyFactor * s.p / s.rho + .5 * (vx * vx + vy * vy));
      }
      const boundary = (j === 0 && f < 2) || (j === nt - 1 && (f === 3 || f === 4))
        || (i === 1 && f === 5) || (i === nx - 1 && f === 2);
      face.forEach((v, k) => { flux[k] += v; if (boundary) external[k] += v; });
    }
    flux.forEach((v, k) => { maxLocal[k] = Math.max(maxLocal[k], Math.abs(v)); total[k] += v; });
  }
  return { order: ['mass', 'xMomentum', 'yMomentum', 'totalEnthalpy'], maxLocal, total, external,
    internalCancellation: total.map((v, k) => Math.abs(v - external[k])) };
}
