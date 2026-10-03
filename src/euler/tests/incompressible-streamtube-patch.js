// SPDX-License-Identifier: GPL-2.0-or-later
// Small verification patch for incompressible normal-pressure relaxation.
// All wall/outer/end nodes and tube masses are prescribed. Interior nodes
// move normally to the supplied streamline chart. This has no body cuts,
// free stagnation/capture conditions, compressible gas or BL coupling.
import { evaluateIncompressibleStreamtubeCell } from '../incompressible-streamtube-cell.js';
import { solveNewton } from '../../numerics/newton.js';

export function createIncompressibleStreamtubePatch({ nodes, massFlows, density = 1, pressureCorrectionFactor = .1 }) {
  const nx = nodes?.length - 1, nt = massFlows?.length;
  if (!Array.isArray(nodes) || nx < 2 || !Array.isArray(massFlows) || nt < 2
    || !massFlows.every(m => Number.isFinite(m) && m > 0)
    || !nodes.every(row => Array.isArray(row) && row.length === nt + 1 && row.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y)))
    || ![density, pressureCorrectionFactor].every(Number.isFinite) || density <= 0 || pressureCorrectionFactor < 0)
    throw new Error('Invalid incompressible streamtube patch.');
  nodes = nodes.map(row => row.map(p => ({ ...p }))); massFlows = massFlows.slice();
  const lengthScale = Math.hypot(nodes[0][nt].x - nodes[0][0].x, nodes[0][nt].y - nodes[0][0].y);
  if (!(lengthScale > 0)) throw new Error('Degenerate incompressible patch width.');
  const n = (nx - 1) * (nt - 1), index = (i, j) => (i - 1) * (nt - 1) + j - 1;
  // Dense numerical derivatives are an independent small-system oracle,
  // deliberately bounded rather than a production body-grid backend.
  if (n > 200) throw new Error('The dense verification patch is limited to 200 free nodes.');
  const directions = nodes.map((row, i) => row.map((p, j) => {
    const a = nodes[Math.max(0, i - 1)][j], b = nodes[Math.min(nx, i + 1)][j], length = Math.hypot(b.x - a.x, b.y - a.y);
    if (!(length > 0)) throw new Error('Degenerate incompressible patch tangent.');
    return { x: -(b.y - a.y) / length, y: (b.x - a.x) / length };
  }));
  const pressureScale = massFlows.reduce((s, m) => s + m, 0) ** 2 / (density * lengthScale ** 2);
  if (!(pressureScale > 0) || !Number.isFinite(pressureScale)) throw new Error('Unresolved incompressible patch pressure scale.');
  const decode = state => {
    if (state.length !== n || !state.every(Number.isFinite)) throw new Error('Invalid incompressible patch state.');
    return nodes.map((row, i) => row.map((p, j) => {
      const shift = i === 0 || i === nx || j === 0 || j === nt ? 0 : lengthScale * state[index(i, j)];
      return { x: p.x + shift * directions[i][j].x, y: p.y + shift * directions[i][j].y };
    }));
  };
  const evaluate = state => {
    const points = decode(state), residual = new Float64Array(n), cells = [];
    for (let i = 1; i < nx; i++) {
      const strip = massFlows.map((massFlow, j) => evaluateIncompressibleStreamtubeCell({
        lower: [points[i - 1][j], points[i][j], points[i + 1][j]],
        upper: [points[i - 1][j + 1], points[i][j + 1], points[i + 1][j + 1]], massFlow, density, pressureCorrectionFactor }));
      for (let j = 1; j < nt; j++) residual[index(i, j)] = (strip[j - 1].interfacePressure.upper - strip[j].interfacePressure.lower) / pressureScale;
      cells.push(strip);
    }
    return { nodes: points, cells, residual };
  };
  const admissible = state => { try { evaluate(state); return true; } catch { return false; } };
  return { n, initial: new Float64Array(n), decode, evaluate, residual: state => evaluate(state).residual, admissible,
    conditions: { density, massFlows, lengthScale, pressureScale, pressureCorrectionFactor } };
}

export function solveIncompressibleStreamtubePatch(system, { initial = system.initial, tolerance = 1e-10, maxIterations = 20 } = {}) {
  const r = solveNewton({ initial, tolerance, maxIterations, residual: system.residual, admissible: system.admissible });
  return { ...r, ...system.evaluate(r.x), formulation: 'Incompressible fixed-boundary verification patch; no body capture or Euler/BL coupling.' };
}
