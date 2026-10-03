// SPDX-License-Identifier: GPL-2.0-or-later
// Intrinsic Euler verification channel: all densities and internal streamline
// positions are solved simultaneously. Fixed x cross-lines and prescribed
// tube masses are channel constraints, not an airfoil mass-capture model.
import { evaluateStreamtubeCell } from '../streamtube-cell.js';
import { createStreamtubeChannelJacobian } from './streamtube-channel-jacobian.js';
import { solveNewton } from '../../numerics/newton.js';
import { solveSparseDirect } from '../../numerics/klu.js';

export function createStreamtubeChannel({ x, lower, upper, massFlows,
  stagnationEnthalpy, stagnationDensity, gamma = 1.4,
  referenceDensity = 1, referencePressure = 1, pressureCorrectionFactor = .1,
  inletSlopes, outletSlopes, streamwiseMode = 'momentum' }) {
  if (!['momentum', 'isentropic'].includes(streamwiseMode)) throw new Error('Unknown intrinsic channel streamwise mode.');
  const finiteArray = a => Array.isArray(a) && a.every(Number.isFinite);
  if (!finiteArray(x) || x.length < 3 || x.some((v, i) => i && v <= x[i - 1])
    || ![lower, upper].every(a => finiteArray(a) && a.length === x.length)
    || lower.some((v, i) => v >= upper[i])) throw new Error('Invalid channel cross-lines or walls.');
  if (!finiteArray(massFlows) || !massFlows.length || massFlows.some(v => v <= 0))
    throw new Error('Prescribe positive mass flow in every streamtube.');
  if (![referenceDensity, referencePressure, gamma, pressureCorrectionFactor].every(Number.isFinite)
    || Math.min(referenceDensity, referencePressure) <= 0 || gamma <= 1 || pressureCorrectionFactor < 0)
    throw new Error('Invalid channel scales or gas controls.');
  const nx = x.length - 1, nt = massFlows.length, densityCount = nx * nt;
  const positionCount = (nx + 1) * (nt - 1), n = densityCount + positionCount;
  const perTube = (v, name) => {
    const a = Array.isArray(v) ? v.slice() : Array(nt).fill(v);
    if (a.length !== nt || !a.every(v => Number.isFinite(v) && v > 0)) throw new Error(`Invalid channel ${name}.`);
    return a;
  };
  const h0 = perTube(stagnationEnthalpy, 'stagnation enthalpy');
  const rhoTotal = perTube(stagnationDensity, 'stagnation density');
  const slopeArray = v => {
    if (v !== undefined && !finiteArray(v)) throw new Error('Invalid internal-streamline end slopes.');
    const a = v === undefined ? Array(nt - 1).fill(0) : v.slice();
    if (!finiteArray(a) || a.length !== nt - 1) throw new Error('Invalid internal-streamline end slopes.');
    return a;
  };
  const slopesIn = slopeArray(inletSlopes), slopesOut = slopeArray(outletSlopes);
  // Copy prescribed data so callers cannot change the equations mid-solve.
  x = x.slice(); lower = lower.slice(); upper = upper.slice(); massFlows = massFlows.slice();
  const heightScale = upper[0] - lower[0], totalMass = massFlows.reduce((s, m) => s + m, 0);
  const fractions = [0];
  for (const m of massFlows) fractions.push(fractions.at(-1) + m / totalMass);
  const positionIndex = (i, j) => densityCount + i * (nt - 1) + j - 1;
  const initial = new Float64Array(n);
  const decode = state => {
    if (state.length !== n || !state.every(Number.isFinite)) throw new Error('Invalid channel state vector.');
    const nodes = x.map((xi, i) => Array.from({ length: nt + 1 }, (_, j) => ({ x: xi,
      y: j === 0 ? lower[i] : j === nt ? upper[i]
        : lower[i] + fractions[j] * (upper[i] - lower[i]) + heightScale * state[positionIndex(i, j)] })));
    const densities = Array.from({ length: nx }, (_, i) => Array.from({ length: nt }, (_, j) => referenceDensity * Math.exp(state[i * nt + j])));
    return { nodes, densities };
  };
  const evaluate = state => {
    const { nodes, densities } = decode(state), cells = [], residual = new Float64Array(n);
    let row = 0;
    for (let i = 1; i < nx; i++) {
      const strip = [];
      for (let j = 0; j < nt; j++) {
        const cell = evaluateStreamtubeCell({ lower: [nodes[i - 1][j], nodes[i][j], nodes[i + 1][j]],
          upper: [nodes[i - 1][j + 1], nodes[i][j + 1], nodes[i + 1][j + 1]],
          densities: [densities[i - 1][j], densities[i][j]], massFlow: massFlows[j],
          stagnationEnthalpy: h0[j], gamma, pressureCorrectionFactor });
        if (cell.states.some(s => s.machSquared >= 1)) throw new Error('Intrinsic channel currently requires subsonic flow.');
        residual[row++] = (streamwiseMode === 'isentropic' ? cell.isentropicResidual : cell.streamwiseResidual) / referencePressure;
        strip.push(cell);
      }
      cells.push(strip);
    }
    for (let j = 0; j < nt; j++) {
      const s = cells[0][j].states[0];
      residual[row++] = Math.log(s.rho / rhoTotal[j]) - Math.log(s.enthalpy / h0[j]) / (gamma - 1);
    }
    for (const strip of cells) for (let j = 0; j < nt - 1; j++)
      residual[row++] = (strip[j].interfacePressure.upper - strip[j + 1].interfacePressure.lower) / referencePressure;
    for (const [i, slopes] of [[0, slopesIn], [nx - 1, slopesOut]]) for (let j = 1; j < nt; j++)
      residual[row++] = (nodes[i + 1][j].y - nodes[i][j].y - slopes[j - 1] * (x[i + 1] - x[i])) / heightScale;
    if (row !== n || !residual.every(Number.isFinite)) throw new Error('Invalid intrinsic channel equation count or residual.');
    const sections = [cells[0].map(c => c.states[0]), ...cells.map(strip => strip.map(c => c.states[1]))];
    return { residual, nodes, densities, cells, sections };
  };
  const admissible = state => { try { evaluate(state); return true; } catch { return false; } };
  const system = { n, nx, nt, densityCount, positionCount, positionIndex, initial, decode, evaluate, admissible,
    residual: state => evaluate(state).residual,
    conditions: { massFlows, stagnationEnthalpy: h0, stagnationDensity: rhoTotal, gamma,
      referenceDensity, referencePressure, heightScale, pressureCorrectionFactor, streamwiseMode } };
  system.jacobian = createStreamtubeChannelJacobian(system);
  return system;
}

// The dense whole-residual finite-difference backend remains an independent
// oracle. KLU/WASM verifies each solve against the original sparse matrix.
export function solveStreamtubeChannel(system, { initial = system.initial,
  jacobianBackend = 'analytic', linearBackend = jacobianBackend === 'analytic' ? 'klu' : 'dense', ...controls } = {}) {
  if (!['analytic', 'finite-difference'].includes(jacobianBackend) || !['dense', 'klu'].includes(linearBackend)
    || (linearBackend === 'klu' && jacobianBackend !== 'analytic')) throw new Error('Invalid intrinsic channel solver backends.');
  const jacobian = jacobianBackend === 'analytic' ? state => system.jacobian(state, { sparse: linearBackend === 'klu' }) : undefined;
  const linearDiagnostics = linearBackend === 'klu' ? { solves: 0, maxRelativeResidual: 0, maxFactorNonzeros: 0, refinements: 0 } : null;
  const linearSolve = linearBackend === 'klu' ? (matrix, rhs) => {
    const r = solveSparseDirect(matrix, rhs);
    linearDiagnostics.solves++;
    linearDiagnostics.maxRelativeResidual = Math.max(linearDiagnostics.maxRelativeResidual, r.relativeResidual);
    linearDiagnostics.maxFactorNonzeros = Math.max(linearDiagnostics.maxFactorNonzeros, r.factorNonzeros);
    linearDiagnostics.refinements += r.refinements;
    return r.x;
  } : undefined;
  const result = solveNewton({ ...controls, initial, residual: system.residual, admissible: system.admissible, jacobian, linearSolve });
  return { ...result, ...system.evaluate(result.x), jacobianBackend, linearBackend, linearDiagnostics, streamwiseMode: system.conditions.streamwiseMode };
}
