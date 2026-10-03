// Independent continuum solutions; no numerical flux or solver imports.
import { buildMesh } from '../../src/euler/mesh.js';

export function isentropicNozzle(x, { mach = 0.3, height = 0.08, contraction = 0.12, length = 2, gamma = 1.4 } = {}) {
  const area = height * (1 - contraction * Math.sin(Math.PI * x / length) ** 2);
  const areaMach = m => (2 / (gamma + 1) * (1 + (gamma - 1) / 2 * m * m)) ** ((gamma + 1) / (2 * (gamma - 1))) / m;
  const target = area / height * areaMach(mach);
  if (target < 1) throw new Error('Choked nozzle oracle.');
  let lo = mach; let hi = 1;
  for (let k = 0; k < 60; k++) { const m = (lo + hi) / 2; if (areaMach(m) > target) lo = m; else hi = m; }
  const m = (lo + hi) / 2;
  const temperature = (1 + (gamma - 1) / 2 * mach ** 2) / (1 + (gamma - 1) / 2 * m ** 2);
  const rho = temperature ** (1 / (gamma - 1));
  const p = temperature ** (gamma / (gamma - 1)) / (gamma * mach ** 2);
  return { area, mach: m, rho, p, u: m / mach * Math.sqrt(temperature), v: 0 };
}

// Exact steady 2-D compressible irrotational vortex: q=K/r, uniform entropy
// and stagnation enthalpy. Radial momentum gives dp/dr=rho*q^2/r.
export function annularVortex(nx = 12, ny = 3, { mach = 0.3, gamma = 1.4 } = {}) {
  const k = 1.5; const pRef = 1 / (gamma * mach * mach); const h0 = gamma / (gamma - 1) * pRef + 0.5;
  const exact = ({ x, y }) => {
    const radius = Math.hypot(x, y); const speed = k / radius;
    const a2 = (gamma - 1) * (h0 - 0.5 * speed * speed);
    const rho = (a2 / (gamma * pRef)) ** (1 / (gamma - 1));
    return { rho, p: rho * a2 / gamma, u: -speed * y / radius, v: speed * x / radius };
  };
  const vertices = []; const cells = []; const boundaries = []; const id = (i, j) => i * (ny + 1) + j;
  for (let i = 0; i <= nx; i++) for (let j = 0; j <= ny; j++) {
    const angle = Math.PI / 3 * i / nx; const radius = 1 + j / ny;
    vertices.push({ x: radius * Math.cos(angle), y: radius * Math.sin(angle) });
  }
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) cells.push([id(i, j), id(i, j + 1), id(i + 1, j + 1), id(i + 1, j)]);
  for (let i = 0; i < nx; i++) for (const j of [0, ny]) boundaries.push({ a: id(i, j), b: id(i + 1, j), type: 'wall', side: j === 0 ? 'inner' : 'outer' });
  for (const i of [0, nx]) for (let j = 0; j < ny; j++) {
    const a = id(i, j); const b = id(i, j + 1);
    boundaries.push({ a, b, type: 'farfield', state: exact({ x: (vertices[a].x + vertices[b].x) / 2, y: (vertices[a].y + vertices[b].y) / 2 }) });
  }
  return { mesh: buildMesh(vertices, cells, boundaries), exact };
}

export function weightedError(result, exact, variable) {
  const { mesh, states } = result;
  return mesh.cells.reduce((sum, c, i) => sum + c.area * Math.abs(states[i][variable] - exact(c)[variable]), 0) / mesh.cells.reduce((sum, c) => sum + c.area, 0);
}
