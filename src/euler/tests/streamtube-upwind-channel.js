// SPDX-License-Identifier: GPL-2.0-or-later
// Transonic verification channel using the same intrinsic moving-grid cells
// as the airfoil solver. Reservoir inflow: h0, entropy, direction. Subsonic
// outflow: static pressure, with outlet direction free. Tube masses are
// unknown; fixed inlet streamline labels remove their repartition gauge.
// This assembly is not yet connected to the body/BL solver or GUI.
import { createStreamtubeChannel } from './streamtube-channel.js';
import { streamtubeCellGeometry, streamtubeSection, evaluateStreamtubeCell } from '../streamtube-cell.js';
import { linearizeStreamtubeCell } from '../streamtube-linearization.js';
import { evaluateStreamtubeSpeedStencil } from './streamtube-speed-stencil.js';
import { linearizeStreamtubeSpeedUpwind } from '../streamtube-speed-upwind.js';
import { evaluateStreamtubeMomentumBlend, linearizeStreamtubeMomentumBlend } from '../streamtube-momentum-blend.js';
import { sparseMatrix, sparseAdd } from '../../numerics/sparse.js';

const cross = (a, b) => a.x * b.y - a.y * b.x;

export function createUpwindStreamtubeChannel(input) {
  if (input.streamwiseMode !== undefined && !['momentum', 'hybrid'].includes(input.streamwiseMode))
    throw new Error('The transonic channel requires momentum or hybrid streamwise equations.');
  const hybridControls = input.streamwiseMode === 'hybrid' ? { epsilonP: input.hybrid?.epsilonP } : null;
  if (hybridControls && !(Number.isFinite(hybridControls.epsilonP) && hybridControls.epsilonP > 0))
    throw new Error('Hybrid streamwise equations require a positive finite epsilonP.');
  if (input.outletSlopes !== undefined) throw new Error('Outlet direction is free at the prescribed-pressure boundary.');
  if (input.upwind?.boundary?.kind !== 'unfiltered-first-two')
    throw new Error('Explicit unfiltered-first-two inlet closure required for this channel.');
  const upwindControls = structuredClone(input.upwind);
  const base = createStreamtubeChannel({ ...input, streamwiseMode: 'momentum' });
  const { nx, nt, densityCount, positionCount, positionIndex } = base;
  if (nx < 3) throw new Error('The upwind channel requires at least three sections.');
  const { stagnationEnthalpy: h0, stagnationDensity: rhoTotal, gamma,
    referenceDensity, referencePressure, heightScale, pressureCorrectionFactor } = base.conditions;
  const referenceMass = base.conditions.massFlows.slice();
  const outletPressure = Array.isArray(input.outletPressure) ? input.outletPressure.slice() : Array(nt).fill(input.outletPressure);
  if (outletPressure.length !== nt || !outletPressure.every(p => Number.isFinite(p) && p > 0))
    throw new Error('Prescribe positive outlet static pressure in every tube.');
  const inletSlopes = input.inletSlopes?.slice() ?? Array(nt - 1).fill(0);
  const x = input.x.slice(), n = base.n + nt, massIndex = j => base.n + j;
  const initial = new Float64Array(n); initial.set(base.initial);
  const conditions = { ...base.conditions, upwind: upwindControls, outletPressure,
    ...(hybridControls ? { streamwiseMode: 'hybrid', hybrid: { ...hybridControls } } : {}),
    boundary: 'inlet stagnation state/direction and fixed labels; outlet static pressure with free direction' };
  const decode = state => {
    if (state.length !== n || !state.every(Number.isFinite)) throw new Error('Invalid upwind-channel state.');
    const decoded = base.decode(state.slice(0, base.n));
    const massFlows = referenceMass.map((m, j) => m * Math.exp(state[massIndex(j)]));
    if (!massFlows.every(m => Number.isFinite(m) && m > 0)) throw new Error('Invalid upwind-channel mass flow.');
    return { ...decoded, massFlows };
  };

  function prepare(state, derivatives = false) {
    const { nodes, densities, massFlows } = decode(state);
    const sections = Array.from({ length: nx }, () => Array(nt));
    const sectionGeometry = Array.from({ length: nx }, () => Array(nt));
    const parameters = Array.from({ length: nx - 1 }, () => Array(nt));
    // Shared physical sections are built before any speed filtering. An
    // undissipated interface pressure is not an admissibility prerequisite.
    for (let i = 1; i < nx; i++) for (let j = 0; j < nt; j++) {
      const lower = [nodes[i - 1][j], nodes[i][j], nodes[i + 1][j]];
      const upper = [nodes[i - 1][j + 1], nodes[i][j + 1], nodes[i + 1][j + 1]];
      const geometry = streamtubeCellGeometry(lower, upper);
      parameters[i - 1][j] = { lower, upper, densities: [densities[i - 1][j], densities[i][j]],
        massFlow: massFlows[j], stagnationEnthalpy: h0[j], gamma, pressureCorrectionFactor };
      for (const slot of i === 1 ? [0, 1] : [1]) {
        const k = i - 1 + slot;
        sections[k][j] = streamtubeSection({ density: densities[k][j], massFlow: massFlows[j],
          normalArea: geometry.normalAreas[slot], stagnationEnthalpy: h0[j], gamma });
        sectionGeometry[k][j] = { normalArea: geometry.normalAreas[slot],
          length: geometry.streamwiseLengths[slot], direction: geometry.directions[slot], gap: geometry.sections[slot] };
      }
    }
    for (const end of [0, nx - 1]) if (sections[end].some(s => s.machSquared >= 1))
      throw new Error('Reservoir inlet and prescribed-pressure outlet must remain subsonic.');
    const sectionArcs = Array.from({ length: nt }, (_, j) => {
      let arc = 0;
      return sectionGeometry.map(row => { const length = row[j].length, s = arc + .5 * length; arc += length; return s; });
    });
    const filterParameters = sectionArcs.map((sectionArc, j) => ({ ...upwindControls, gamma, sectionArc,
      speeds: sections.map(row => row[j].q), machSquared: sections.map(row => row[j].machSquared) }));
    const upwind = filterParameters.map(evaluateStreamtubeSpeedStencil);
    const transportSpeeds = sections.map((_, i) => upwind.map(f => f.speeds[i]));
    const cells = [], linears = [], residual = new Float64Array(n);
    const hybrid = hybridControls ? [] : null, hybridLinears = hybridControls ? [] : null;
    let row = 0;
    for (let i = 1; i < nx; i++) {
      const strip = [], linearStrip = [];
      const hybridStrip = [], hybridLinearStrip = [];
      for (let j = 0; j < nt; j++) {
        const p = { ...parameters[i - 1][j], transportSpeeds: [transportSpeeds[i - 1][j], transportSpeeds[i][j]] };
        const linear = derivatives ? linearizeStreamtubeCell(p) : null;
        const cell = linear ? linear.value : evaluateStreamtubeCell(p);
        if (hybridControls) {
          const blendInput = { states: cell.states, transportSpeeds: cell.transportSpeeds,
            streamwiseResidual: cell.streamwiseResidual, isentropicResidual: cell.isentropicResidual, ...hybridControls };
          const blendLinear = derivatives ? linearizeStreamtubeMomentumBlend(blendInput) : null;
          const blend = blendLinear ? blendLinear.value : evaluateStreamtubeMomentumBlend(blendInput);
          residual[row++] = blend.residual / referencePressure;
          hybridStrip.push(blend); hybridLinearStrip.push(blendLinear);
        } else residual[row++] = cell.streamwiseResidual / referencePressure;
        strip.push(cell); linearStrip.push(linear);
      }
      cells.push(strip); linears.push(linearStrip);
      if (hybridControls) { hybrid.push(hybridStrip); hybridLinears.push(hybridLinearStrip); }
    }
    for (let j = 0; j < nt; j++) {
      const s = sections[0][j];
      residual[row++] = Math.log(s.rho / rhoTotal[j]) - Math.log(s.enthalpy / h0[j]) / (gamma - 1);
    }
    for (const strip of cells) for (let j = 0; j < nt - 1; j++)
      residual[row++] = (strip[j].interfacePressure.upper - strip[j + 1].interfacePressure.lower) / referencePressure;
    for (let j = 1; j < nt; j++)
      residual[row++] = (nodes[1][j].y - nodes[0][j].y - inletSlopes[j - 1] * (x[1] - x[0])) / heightScale;
    for (let j = 1; j < nt; j++) residual[row++] = state[positionIndex(0, j)];
    for (let j = 0; j < nt; j++) residual[row++] = (sections.at(-1)[j].p - outletPressure[j]) / referencePressure;
    if (row !== n || !residual.every(Number.isFinite)) throw new Error('Invalid upwind-channel residual.');
    const value = { residual, nodes, densities, massFlows, cells, sections, sectionArcs, transportSpeeds, upwind,
      ...(hybridControls ? { hybrid } : {}) };
    if (!derivatives) return { value };
    // Only local three-section filters are needed by a cell's Jacobian.
    // Stencil assembly below therefore remains O(nx*nt), not dense AD.
    const filterLinears = filterParameters.map(p => p.speeds.map((_, k) => k < 2 ? null : linearizeStreamtubeSpeedUpwind({
      speeds: p.speeds.slice(k - 2, k + 1), machSquared: p.machSquared.slice(k - 1, k + 1),
      spacing: [p.sectionArc[k - 1] - p.sectionArc[k - 2], p.sectionArc[k] - p.sectionArc[k - 1]],
      mucon: p.mucon, mcrit: p.mcrit, gamma })));
    return { value, linears, sectionGeometry, filterLinears, hybridLinears };
  }

  const columns = Array.from({ length: n }, () => new Set()), stencils = [];
  for (let i = 1; i < nx; i++) for (let j = 0; j < nt; j++) {
    const streamwiseRow = (i - 1) * nt + j, inletRow = i === 1 ? (nx - 1) * nt + j : null;
    const upperRow = j < nt - 1 ? densityCount + (i - 1) * (nt - 1) + j : null;
    const lowerRow = j > 0 ? densityCount + (i - 1) * (nt - 1) + j - 1 : null;
    const outletRow = i === nx - 1 ? base.n + j : null;
    const variables = [{ col: massIndex(j), kind: 'mass' }];
    for (let k = Math.max(0, i - 3); k <= i; k++) variables.push({ col: k * nt + j, kind: 'density', index: k });
    for (const streamline of [j, j + 1]) if (streamline > 0 && streamline < nt)
      for (let k = Math.max(0, i - 3); k <= i + 1; k++)
        variables.push({ col: positionIndex(k, streamline), kind: 'node', index: k, streamline });
    const rows = [streamwiseRow, inletRow, upperRow, lowerRow, outletRow].filter(r => r !== null);
    for (const row of rows) for (const { col } of variables) columns[row].add(col);
    stencils.push({ i, j, streamwiseRow, inletRow, upperRow, lowerRow, outletRow, variables });
  }
  const endRows = [];
  let row = densityCount + (nx - 1) * (nt - 1);
  for (let j = 1; j < nt; j++) {
    const left = positionIndex(0, j), right = positionIndex(1, j);
    columns[row].add(left); columns[row].add(right); endRows.push({ row: row++, left, right });
  }
  for (let j = 1; j < nt; j++) { const col = positionIndex(0, j); columns[row].add(col); endRows.push({ row: row++, anchor: col }); }
  if (row !== base.n) throw new Error('Upwind-channel boundary row count mismatch.');
  const pattern = sparseMatrix(columns);

  function jacobian(state, { sparse = false } = {}) {
    const { value, linears, sectionGeometry, filterLinears, hybridLinears } = prepare(state, true);
    const { sections, massFlows } = value;
    const matrix = sparse ? { ...pattern, values: new Float64Array(pattern.values.length) } : new Float64Array(n * n);
    const add = sparse ? (r, c, v) => sparseAdd(matrix, r, c, v) : (r, c, v) => { matrix[r * n + c] += v; };
    for (const { i, j, streamwiseRow, inletRow, upperRow, lowerRow, outletRow, variables } of stencils) {
      const cell = linears[i - 1][j];
      for (const variable of variables) {
        const { col, kind, index, streamline } = variable, dm = kind === 'mass' ? massFlows[j] : 0;
        const dy = (k, line) => kind === 'node' && index === k && streamline === line ? heightScale : 0;
        const sectionDerivative = k => {
          const g = sectionGeometry[k][j], s = sections[k][j];
          const segmentY = .5 * (dy(k + 1, j) - dy(k, j) + dy(k + 1, j + 1) - dy(k, j + 1));
          const gapY = .5 * (dy(k, j + 1) - dy(k, j) + dy(k + 1, j + 1) - dy(k + 1, j));
          const length = g.direction.y * segmentY;
          const direction = { x: -g.direction.x * length / g.length, y: (segmentY - g.direction.y * length) / g.length };
          const normalArea = cross(direction, g.gap) + g.direction.x * gapY;
          const rho = kind === 'density' && index === k ? s.rho : 0;
          const chain = d => d.density * rho + d.massFlow * dm + d.normalArea * normalArea;
          const q = chain(s.derivatives.q), p = chain(s.derivatives.p), enthalpy = -s.q * q;
          const machSquared = (2 * s.q * q - s.machSquared * (gamma - 1) * enthalpy) / ((gamma - 1) * s.enthalpy);
          return { rho, q, p, enthalpy, machSquared, length };
        };
        const ds = new Map();
        for (let k = Math.max(0, i - 3); k <= i; k++) ds.set(k, sectionDerivative(k));
        const biasedDerivative = k => {
          if (k < 2) return ds.get(k).q;
          const a = ds.get(k - 2), b = ds.get(k - 1), c = ds.get(k);
          return filterLinears[j][k].apply({ speeds: [a.q, b.q, c.q], machSquared: [b.machSquared, c.machSquared],
            spacing: [.5 * (a.length + b.length), .5 * (b.length + c.length)] }).speed;
        };
        const tangent = { massFlow: dm, densities: [ds.get(i - 1).rho, ds.get(i).rho],
          transportSpeeds: [biasedDerivative(i - 1), biasedDerivative(i)] };
        if (kind === 'node' && index >= i - 1 && index <= i + 1) {
          const side = streamline === j ? 'lower' : 'upper';
          tangent[side] = Array.from({ length: 3 }, (_, k) => ({ x: 0, y: index === i - 1 + k ? heightScale : 0 }));
        }
        const d = cell.apply(tangent);
        const streamwiseDerivative = hybridControls ? hybridLinears[i - 1][j].apply({ states: d.states,
          transportSpeeds: d.transportSpeeds, streamwiseResidual: d.streamwiseResidual,
          isentropicResidual: d.isentropicResidual }).residual : d.streamwiseResidual;
        add(streamwiseRow, col, streamwiseDerivative / referencePressure);
        if (upperRow !== null) add(upperRow, col, d.interfacePressure.upper / referencePressure);
        if (lowerRow !== null) add(lowerRow, col, -d.interfacePressure.lower / referencePressure);
        if (inletRow !== null) {
          const s = sections[0][j], v = ds.get(0);
          add(inletRow, col, v.rho / s.rho - v.enthalpy / ((gamma - 1) * s.enthalpy));
        }
        if (outletRow !== null) add(outletRow, col, ds.get(nx - 1).p / referencePressure);
      }
    }
    for (const { row, left, right, anchor } of endRows) {
      if (anchor !== undefined) add(row, anchor, 1);
      else { add(row, left, -1); add(row, right, 1); }
    }
    if (sparse ? !matrix.values.every(Number.isFinite) : !matrix.every(Number.isFinite))
      throw new Error('Nonfinite upwind-channel Jacobian.');
    return matrix;
  }
  const evaluate = state => prepare(state).value;
  const admissible = state => { try { evaluate(state); return true; } catch { return false; } };
  return { n, nx, nt, densityCount, positionCount, positionIndex, massIndex, initial, decode, evaluate, admissible,
    residual: state => evaluate(state).residual, jacobian, conditions };
}
