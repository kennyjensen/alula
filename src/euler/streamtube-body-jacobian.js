// SPDX-License-Identifier: GPL-2.0-or-later
// Full intrinsic body Jacobian: conservative cells, moving wall/cut charts,
// captured mass, LE/TE Kutta, end tangency and normalized farfield matching.
// The independently implemented residual remains the verification oracle.
import { linearizeStreamtubeCell } from './streamtube-linearization.js';
import { linearizeIncompressibleStreamtubeCell } from './incompressible-streamtube-linearization.js';
import { streamtubeEdgeVelocity, streamtubeEdgeVelocityTangent } from './streamtube-edge-velocity.js';
import { multipoleGeometryDerivatives } from '../potential/farfield-geometry.js';
import { sparseMatrixFromRows } from '../numerics/sparse-rows.js';
import { streamtubeWakeGap } from './streamtube-wake-geometry.js';
import { prepareStreamtubeTransportChain } from './streamtube-transport-chain.js';
import { createStreamtubeJacobianTransport } from './streamtube-jacobian-transport.js';
import { linearizeStreamtubeMomentumBlend } from './streamtube-momentum-blend.js';
import { streamtubeEquationAt } from './streamtube-equation-selection.js';

const zero = Object.freeze({ x: 0, y: 0 });
const sub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });
const mean = (a, b) => ({ x: .5 * (a.x + b.x), y: .5 * (a.y + b.y) });
const dot = (a, b) => a.x * b.x + a.y * b.y;
const cross = (a, b) => a.x * b.y - a.y * b.x;
const accumulate = (row, col, value) => row.set(col, (row.get(col) ?? 0) + value);
const sumPoints = (a, b, wa = .5, wb = .5) => new Map([...new Set([...a.keys(), ...b.keys()])].map(col => {
  const p = a.get(col) ?? zero, q = b.get(col) ?? zero;
  return [col, { x: wa * p.x + wb * q.x, y: wa * p.y + wb * q.y }];
}));

export function createStreamtubeBodyJacobian(system) {
  const { layout, conditions } = system, { nx, tubes, elements, globals, n } = layout;
  const { gamma, mach, alpha, pInf, pressureScale, h0, lengthScale: length, massScale, pressureCorrectionFactor, center, streamwiseMode, flowModel, geometryDomain, upwind, hybrid } = conditions;
  const incompressible = flowModel === 'incompressible';
  const angle = alpha * Math.PI / 180, freestream = { x: Math.cos(angle), y: Math.sin(angle) };
  const farfieldColumns = [globals.circulation, globals.source, globals.doubletX, globals.doubletY];
  // `raw` is private to one synchronous evaluation/assembly call. Callers
  // cannot supply an earlier evaluation after changing state or geometry.
  const assemble = (state, raw, { sparse = false, includeDisplacement = false } = {}) => {
    const { nodes, allocation, strengths } = raw;
    const geometry = system.geometryDerivatives(state, { includeDisplacement });
    // Each passage/tube is one continuous physical chain through the body
    // cuts and wake. Prepare its local speed-filter derivatives once; a
    // cell-column query below touches at most four physical sections.
    const chains = upwind ? tubes.map((count, g) => Array.from({ length: count }, (_, j) =>
      prepareStreamtubeTransportChain({ lower: nodes[g].map(line => line[j]), upper: nodes[g].map(line => line[j + 1]),
        densities: Array.from({ length: nx }, (_, k) => Math.exp(state[layout.densityIndex(k, g, j)])),
        massFlow: allocation.groups[g][j].massFlow, stagnationEnthalpy: h0, gamma, geometryDomain, upwind },
      { linearize: true }))) : null;
    const massDerivatives = allocation.groups.map(group => group.map(mass =>
      new Map([...mass.derivatives].map(([body, derivative]) => [globals.capture[body], derivative * massScale]))));
    const transport = chains ? createStreamtubeJacobianTransport({ chains, geometry, layout, massDerivatives }) : null;
    const cells = Array.from({ length: nx - 1 }, () => tubes.map(() => []));
    for (let i = 1; i < nx; i++) for (let g = 0; g <= elements; g++) for (let j = 0; j < tubes[g]; j++) {
      const lowerMaps = [geometry[g][i - 1][j], geometry[g][i][j], geometry[g][i + 1][j]];
      const upperMaps = [geometry[g][i - 1][j + 1], geometry[g][i][j + 1], geometry[g][i + 1][j + 1]];
      const densityColumns = incompressible ? [] : [layout.densityIndex(i - 1, g, j), layout.densityIndex(i, g, j)];
      const densities = densityColumns.map(col => Math.exp(state[col])), mass = allocation.groups[g][j];
      const massColumns = massDerivatives[g][j];
      const columns = new Set([...densityColumns, ...massColumns.keys(), ...[...lowerMaps, ...upperMaps].flatMap(row => [...row.keys()])]);
      if (upwind) {
        // Filtering section i reaches back to section i-2; filtering i-1
        // reaches i-3. Geometry maps also contain moving stagnation and
        // displacement dependencies which are not plain node columns.
        for (let k = Math.max(0, i - 3); k <= i; k++) columns.add(layout.densityIndex(k, g, j));
        for (let k = Math.max(0, i - 3); k <= i + 1; k++)
          for (const side of [j, j + 1]) for (const col of geometry[g][k][side].keys()) columns.add(col);
      }
      const cell = (incompressible ? linearizeIncompressibleStreamtubeCell : linearizeStreamtubeCell)({ lower: [nodes[g][i - 1][j], nodes[g][i][j], nodes[g][i + 1][j]],
        upper: [nodes[g][i - 1][j + 1], nodes[g][i][j + 1], nodes[g][i + 1][j + 1]],
        densities, massFlow: mass.massFlow, stagnationEnthalpy: h0, gamma, pressureCorrectionFactor, geometryDomain,
        ...(upwind ? { transportSpeeds: [raw.transportSpeeds[i - 1][g][j], raw.transportSpeeds[i][g][j]] } : {}) });
      const equation = hybrid ? streamtubeEquationAt({ hybrid, bodies: layout.bodies, tubes, nx: layout.nx, i, group: g, tube: j }) : null;
      const blend = equation === 'hybrid' ? linearizeStreamtubeMomentumBlend({ ...cell.value, epsilonP: hybrid.epsilonP }) : null;
      const derivatives = new Map();
      // The local kernel consumes inputs synchronously and returns owned
      // derivatives. Reuse only the input buffers across column queries.
      const tangent = { lower: new Array(3), upper: new Array(3),
        ...(incompressible ? {} : { densities: new Array(2) }), massFlow: 0 };
      for (const col of columns) {
        for (let k = 0; k < 3; k++) {
          tangent.lower[k] = lowerMaps[k].get(col) ?? zero; tangent.upper[k] = upperMaps[k].get(col) ?? zero;
        }
        if (!incompressible) for (let k = 0; k < 2; k++) tangent.densities[k] = densityColumns[k] === col ? densities[k] : 0;
        tangent.massFlow = massColumns.get(col) ?? 0;
        if (transport) tangent.transportSpeeds = transport.forCell(i, g, j, col);
        const d = cell.apply(tangent);
        if (hybrid) d.hybridResidual = blend ? blend.apply(d).residual
          : equation === 'isentropic' ? d.isentropicResidual : d.streamwiseResidual;
        // Retain only fields consumed by row assembly and boundary queries.
        // The full kernel and all finite/domain checks still run above;
        // hybrid blending must precede dropping transport-speed derivatives.
        derivatives.set(col, { streamwiseResidual: d.streamwiseResidual, isentropicResidual: d.isentropicResidual,
          ...(hybrid ? { hybridResidual: d.hybridResidual } : {}), interfacePressure: d.interfacePressure, states: d.states,
          geometry: { pressureCurvature: d.geometry.pressureCurvature, streamwiseLengths: d.geometry.streamwiseLengths } });
      }
      cells[i - 1][g][j] = { value: cell.value, derivatives };
    }
    const section = (i, g, j) => {
      const cell = cells[Math.max(0, i - 1)][g][j], slot = i === 0 ? 0 : 1;
      return { value: cell.value.states[slot], derivatives: new Map([...cell.derivatives].map(([col, d]) => [col, d.states[slot]])) };
    };
    const coefficients = [-strengths.circulation, strengths.source, strengths.doubletX, strengths.doubletY, strengths.circulation ** 2];
    const farfield = (point, pointDerivatives) => {
      const basis = multipoleGeometryDerivatives(point, { center, alpha, mach, gamma }), velocity = { ...freestream };
      for (let k = 0; k < coefficients.length; k++) { velocity.x += coefficients[k] * basis.velocity[k][0]; velocity.y += coefficients[k] * basis.velocity[k][1]; }
      const temperature = 1 + .5 * (gamma - 1) * mach * mach * (1 - dot(velocity, velocity));
      const pressure = incompressible ? -.5 * dot(velocity, velocity) : pInf * temperature ** (gamma / (gamma - 1)), derivatives = new Map();
      for (const col of new Set([...pointDerivatives.keys(), ...farfieldColumns])) {
        const dp = pointDerivatives.get(col) ?? zero;
        const modes = basis.velocityDerivatives.map(d => ({ x: d[0][0] * dp.x + d[0][1] * dp.y, y: d[1][0] * dp.x + d[1][1] * dp.y }));
        const dv = { x: 0, y: 0 };
        for (let k = 0; k < coefficients.length; k++) { dv.x += coefficients[k] * modes[k].x; dv.y += coefficients[k] * modes[k].y; }
        const changes = col === globals.circulation ? [-length, 0, 0, 0, 2 * strengths.circulation * length]
          : col === globals.source ? [0, length, 0, 0, 0] : col === globals.doubletX ? [0, 0, length ** 2, 0, 0]
          : col === globals.doubletY ? [0, 0, 0, length ** 2, 0] : [0, 0, 0, 0, 0];
        for (let k = 0; k < changes.length; k++) { dv.x += changes[k] * basis.velocity[k][0]; dv.y += changes[k] * basis.velocity[k][1]; }
        derivatives.set(col, { velocity: dv, p: incompressible ? -dot(velocity, dv) : -gamma * pressure * mach * mach / temperature * dot(velocity, dv), modes });
      }
      return { velocity, pressure, basis, derivatives };
    };
    const boundaryOutletTangency = (body, side) => {
      if (!Number.isInteger(body) || body < 0 || body >= elements || !['lower', 'upper'].includes(side))
        throw new Error('Invalid wake outlet bank.');
      const g = side === 'lower' ? body : body + 1, j = side === 'lower' ? tubes[g] : 0;
      const a = nodes[g][nx - 1][j], b = nodes[g][nx][j], da = geometry[g][nx - 1][j], db = geometry[g][nx][j];
      const edge = sub(b, a), norm = Math.hypot(edge.x, edge.y), fa = farfield(a, da), fb = farfield(b, db);
      const velocity = mean(fa.velocity, fb.velocity), value = cross(edge, velocity) / norm, derivatives = new Map();
      for (const col of new Set([...fa.derivatives.keys(), ...fb.derivatives.keys()])) {
        const de = sub(db.get(col) ?? zero, da.get(col) ?? zero), dn = dot(edge, de) / norm;
        const dv = mean(fa.derivatives.get(col)?.velocity ?? zero, fb.derivatives.get(col)?.velocity ?? zero);
        derivatives.set(col, (cross(de, velocity) + cross(edge, dv) - value * dn) / norm);
      }
      return { value, derivatives };
    };
    const rows = layout.rows.map(() => new Map());
    const cellTerm = (row, cell, component, sign = 1) => {
      for (const [col, d] of cell.derivatives) accumulate(row, col, sign * component(d));
    };
    // Local and cut/Kutta equations keep the residual's independent ordering.
    for (const row of layout.rows) {
      const out = rows[row.index];
      if (row.kind === 'streamwise') cellTerm(out, cells[row.i - 1][row.group][row.tube], d => (hybrid ? d.hybridResidual : streamwiseMode === 'isentropic' ? d.isentropicResidual : d.streamwiseResidual) / pressureScale);
      else if (row.kind === 'inletDensity') {
        const s = section(0, row.group, row.tube);
        for (const [col, d] of s.derivatives) accumulate(out, col, d.rho / s.value.rho - d.enthalpy / ((gamma - 1) * s.value.enthalpy));
      } else if (row.kind === 'internalPressure') {
        cellTerm(out, cells[row.i - 1][row.group][row.j - 1], d => d.interfacePressure.upper / pressureScale);
        cellTerm(out, cells[row.i - 1][row.group][row.j], d => d.interfacePressure.lower / pressureScale, -1);
      } else if (['cutPressure', 'leadingKutta', 'trailingKutta'].includes(row.kind)) {
        cellTerm(out, cells[row.i - 1][row.body + 1][0], d => d.interfacePressure.lower / pressureScale);
        cellTerm(out, cells[row.i - 1][row.body].at(-1), d => d.interfacePressure.upper / pressureScale, -1);
      } else if (row.kind === 'farfieldPressure') {
        const g = row.side === 'lower' ? 0 : elements, j = row.side === 'lower' ? 0 : tubes[g];
        cellTerm(out, cells[row.i - 1][g][row.side === 'lower' ? 0 : tubes[g] - 1], d => d.interfacePressure[row.side] / pressureScale);
        for (const [col, d] of farfield(nodes[g][row.i][j], geometry[g][row.i][j]).derivatives) accumulate(out, col, -d.p / pressureScale);
      } else if (row.kind === 'wakeGap') {
        const indices = [row.i - 1, row.i, Math.min(nx, row.i + 1)];
        const lowerMaps = indices.map(i => geometry[row.body][i].at(-1)), upperMaps = indices.map(i => geometry[row.body + 1][i][0]);
        const gap = streamtubeWakeGap(indices.map(i => nodes[row.body][i].at(-1)), indices.map(i => nodes[row.body + 1][i][0]));
        for (const col of new Set([...lowerMaps, ...upperMaps].flatMap(m => [...m.keys()])))
          accumulate(out, col, gap.apply(lowerMaps.map(m => m.get(col) ?? zero), upperMaps.map(m => m.get(col) ?? zero)) / length);
        if (includeDisplacement) {
          const index = system.displacementParameters.findIndex(p => p.kind === 'wake' && p.body === row.body && p.index === row.i - layout.bodies[row.body].trailingIndex - 1);
          if (index < 0) throw new Error('Missing independent wake displacement parameter.');
          accumulate(out, n + index, -1 / length);
        }
      } else if (row.kind === 'endTangency') {
        if (conditions.wakeOutlet === 'banks' && row.node.kind === 'cut' && row.i === nx) {
          for (const side of ['lower', 'upper']) for (const [col, d] of boundaryOutletTangency(row.node.body, side).derivatives)
            accumulate(out, col, .5 * d);
          continue;
        }
        const g = row.node.kind === 'cut' ? row.node.body : row.node.group, j = row.node.kind === 'cut' ? tubes[g] : row.node.j;
        const center = layout.displacedBoundaries && row.node.kind === 'cut' && row.i === nx;
        const point = i => center ? mean(nodes[g][i].at(-1), nodes[g + 1][i][0]) : nodes[g][i][j];
        const derivative = i => center ? sumPoints(geometry[g][i].at(-1), geometry[g + 1][i][0]) : geometry[g][i][j];
        const i = row.i === 0 ? 0 : nx - 1, a = point(i), b = point(i + 1), edge = sub(b, a), norm = Math.hypot(edge.x, edge.y);
        const da = derivative(i), db = derivative(i + 1), fa = farfield(a, da), fb = farfield(b, db), velocity = mean(fa.velocity, fb.velocity);
        for (const col of new Set([...fa.derivatives.keys(), ...fb.derivatives.keys()])) {
          const de = sub(db.get(col) ?? zero, da.get(col) ?? zero), dn = dot(edge, de) / norm;
          const dv = mean(fa.derivatives.get(col)?.velocity ?? zero, fb.derivatives.get(col)?.velocity ?? zero);
          accumulate(out, col, (cross(de, velocity) + cross(edge, dv)) / norm - cross(edge, velocity) * dn / norm ** 2);
        }
      } else if (row.kind !== 'farfieldMatch') throw new Error(`Unknown intrinsic body derivative row ${row.kind}.`);
    }
    // The three farfield rows include derivatives of both the least-squares
    // numerator and its normalization; freezing either changes Newton.
    const match = [0, 0, 0], scales = [0, 0, 0], dm = [new Map(), new Map(), new Map()], ds = [new Map(), new Map(), new Map()];
    for (const [g, j, tube] of [[0, 0, 0], [elements, tubes[elements], tubes[elements] - 1]]) for (let i = 0; i < nx; i++) {
      const a = nodes[g][i][j], b = nodes[g][i + 1][j], edge = sub(b, a), norm = Math.hypot(edge.x, edge.y), t = { x: edge.x / norm, y: edge.y / norm };
      const da = geometry[g][i][j], db = geometry[g][i + 1][j], state = section(i, g, tube), q = state.value.q;
      const velocity = { x: q * t.x, y: q * t.y }, ff = farfield(mean(a, b), sumPoints(da, db)), defect = cross(ff.velocity, velocity);
      const modes = ff.basis.velocity.slice(1, 4).map(([x, y]) => ({ x, y })), projection = modes.map(v => cross(v, velocity));
      for (let k = 0; k < 3; k++) { match[k] += norm * defect * projection[k]; scales[k] += norm * projection[k] ** 2; }
      for (const col of new Set([...state.derivatives.keys(), ...ff.derivatives.keys()])) {
        const de = sub(db.get(col) ?? zero, da.get(col) ?? zero), dn = dot(t, de), dq = state.derivatives.get(col)?.q ?? 0;
        const dv = { x: dq * t.x + q * (de.x - t.x * dn) / norm, y: dq * t.y + q * (de.y - t.y * dn) / norm };
        const dff = ff.derivatives.get(col), dDefect = cross(dff?.velocity ?? zero, velocity) + cross(ff.velocity, dv);
        for (let k = 0; k < 3; k++) {
          const mode = projection[k], dMode = cross(dff?.modes[k + 1] ?? zero, velocity) + cross(modes[k], dv);
          accumulate(dm[k], col, dn * defect * mode + norm * (dDefect * mode + defect * dMode));
          accumulate(ds[k], col, dn * mode ** 2 + 2 * norm * mode * dMode);
        }
      }
    }
    for (const row of layout.rows.filter(row => row.kind === 'farfieldMatch')) {
      const k = ['source', 'doubletX', 'doubletY'].indexOf(row.mode), scale = scales[k];
      if (!(scale > 0)) throw new Error('Degenerate intrinsic-body farfield Jacobian mode.');
      for (const [col, d] of dm[k]) accumulate(rows[row.index], col, d / Math.sqrt(scale) - .5 * match[k] * (ds[k].get(col) ?? 0) / scale ** 1.5);
    }
    let matrix = sparse ? null : new Float64Array(n * n);
    const displacementRows = includeDisplacement ? rows.map(() => new Map()) : null;
    for (let row = 0; row < n; row++) for (const [col, value] of rows[row]) {
      if (!Number.isFinite(value)) throw new Error('Nonfinite intrinsic-body Jacobian.');
      if (col >= n) {
        displacementRows[row].set(col - n, value);
        // These row Maps are private to this assembly. Preserve the separate
        // thickness block, then serialize the square state block directly.
        rows[row].delete(col);
      } else if (!sparse) matrix[row * n + col] = value;
    }
    if (sparse) matrix = sparseMatrixFromRows(rows);
    // Rectangular Euler-to-thickness block for a future simultaneous BL
    // assembly. Columns are physical lengths, ordered by the parameter list.
    if (!includeDisplacement) return matrix;
    const boundaryPressure = (body, i, side) => {
      if (!Number.isInteger(i) || i < 1 || i > nx || !['lower', 'upper'].includes(side)) throw new Error('Invalid boundary pressure station.');
      const g = side === 'lower' ? body : body + 1, j = side === 'lower' ? tubes[body] - 1 : 0, key = side === 'lower' ? 'upper' : 'lower';
      // Zero pressure gradient beyond the last interior wake station.
      const cell = cells[Math.min(i, nx - 1) - 1][g][j];
      return { value: cell.value.interfacePressure[key], derivatives: new Map([...cell.derivatives].map(([col, d]) => [col, d.interfacePressure[key]])) };
    };
    const boundaryEdgeVelocity = (body, i, side, interpolation = 'arithmetic') => {
      if (!Number.isInteger(body) || body < 0 || body >= elements || !Number.isInteger(i) || i < 1 || i > nx
        || !['lower', 'upper'].includes(side)) throw new Error('Invalid boundary velocity station.');
      const g = side === 'lower' ? body : body + 1, j = side === 'lower' ? tubes[body] - 1 : 0;
      // Retain the existing zero-gradient terminal wake closure.
      const cell = cells[Math.min(i, nx - 1) - 1][g][j];
      return { value: streamtubeEdgeVelocity(cell.value, .2, interpolation).ue,
        derivatives: new Map([...cell.derivatives].map(([col, d]) => [col, streamtubeEdgeVelocityTangent(cell.value, d, .2, interpolation)])) };
    };
    return { state: matrix, displacement: displacementRows, parameters: system.displacementParameters, boundaryPressure, boundaryEdgeVelocity, boundaryOutletTangency };
  };
  const jacobian = (state, options) => assemble(state, system.evaluate(state), options);
  jacobian.evaluateJacobian = (state, options) => {
    const value = system.evaluate(state);
    return { value, jacobian: assemble(state, value, options) };
  };
  return jacobian;
}
