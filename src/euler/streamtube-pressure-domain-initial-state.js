// SPDX-License-Identifier: GPL-2.0-or-later
import { streamtubeCellGeometry } from './streamtube-cell.js';
import { prepareStreamtubeTransportChain } from './streamtube-transport-chain.js';
import { assertConvexStreamtubeGrid } from '../geometry/streamtube-convex-step.js';

const dot = (a, b) => a.x * b.x + a.y * b.y;
const cross = (a, b) => a.x * b.y - a.y * b.x;

// Construct one fixed-geometry Newton guess, not a new equation of state or
// an isentropic solution. At uniform density r, put u_i = mass/normalArea_i,
// U = (u_0²+u_1²)/2 and k=(gamma-1)/gamma. The actual cell equations give
//   pressureSum = A*r-B/r+2*Pc, pressureDifference = C/r+K*Pc,
// where A=2*k*H0, B=k*U and C,K depend only on the fixed geometry and mass.
// Locally subsonic flow has h_i >= 2*H0/(gamma+1), which bounds the FULL
// auxiliary correction by |Pc| <= D/r, including its signed curvature:
//   D = correctionFactor*(gamma+1)*U*abs(curvature)/2.
// Both side pressures are therefore positive if
//   r² > (B+abs(C)+D*(2+abs(K)))/A.
// A separate sonic bound establishes the premise. The 0.95 fraction leaves
// an interior margin in reciprocal r², as in the grid fraction-to-boundary.
// Actual transport speeds must remain physical for the C/r proof to apply;
// active upwinding is rejected, never silently disabled. The caller MUST
// evaluate the complete selected Euler equations before adopting this seed.
export function initializeStreamtubePressureDomain(system, state) {
  const { gamma, h0, rhoTotal, geometryDomain, upwind, pressureCorrectionFactor = .1 } = system.conditions;
  if (system.conditions.flowModel === 'incompressible'
    || ![gamma, h0, rhoTotal, pressureCorrectionFactor].every(Number.isFinite)
    || gamma <= 1 || h0 <= 0 || rhoTotal <= 0 || pressureCorrectionFactor < 0)
    throw new Error('Invalid pressure-domain density initialization conditions.');
  const initial = Float64Array.from(state), { nodes, allocation } = system.decode(initial);
  const { nx, tubes, densityIndex } = system.layout;
  if (!Number.isInteger(nx) || nx < 2 || !Array.isArray(tubes) || !tubes.length
    || tubes.some(count => !Number.isInteger(count) || count < 1) || !initial.every(Number.isFinite))
    throw new Error('Invalid pressure-domain density initialization layout or state.');
  const k = (gamma - 1) / gamma, A = 2 * k * h0;
  const sonicSpeedSquared = 2 * (gamma - 1) * h0 / (gamma + 1);
  let pressure = { squaredDensity: 0 }, sonic = { squaredDensity: 0 }, thermal = { squaredDensity: 0 };
  for (let g = 0; g < tubes.length; g++) for (let i = 1; i < nx; i++) for (let j = 0; j < tubes[g]; j++) {
    const geometry = streamtubeCellGeometry([i - 1, i, i + 1].map(t => nodes[g][t][j]),
      [i - 1, i, i + 1].map(t => nodes[g][t][j + 1]), { geometryDomain });
    const massFlow = allocation.groups[g][j].massFlow;
    if (!Number.isFinite(massFlow) || massFlow <= 0) throw new Error('Invalid pressure-domain streamtube mass flow.');
    const u = geometry.normalAreas.map(area => massFlow / area), U = .5 * (u[0] ** 2 + u[1] ** 2);
    const B = k * U, C = massFlow * (u[0] * dot(geometry.directions[0], geometry.transverse)
      - u[1] * dot(geometry.directions[1], geometry.transverse)) / geometry.area;
    const K = cross(geometry.sections[0], geometry.sections[1]) / geometry.area;
    const D = pressureCorrectionFactor * (gamma + 1) / 2 * U * Math.abs(geometry.pressureCurvature);
    const squaredDensity = (B + Math.abs(C) + D * (2 + Math.abs(K))) / A;
    if (![squaredDensity, ...u].every(Number.isFinite) || squaredDensity < 0)
      throw new Error('Nonfinite pressure-domain density bound.');
    if (squaredDensity > pressure.squaredDensity)
      pressure = { squaredDensity, cell: { i, group: g, tube: j }, A, B, C, D, K, U,
        pressureCurvature: geometry.pressureCurvature };
    for (const [side, value] of u.entries()) {
      const section = { i: i - 1 + side, group: g, tube: j };
      if (value * value / sonicSpeedSquared > sonic.squaredDensity)
        sonic = { squaredDensity: value * value / sonicSpeedSquared, section, massFlux: value };
      if (value * value / (2 * h0) > thermal.squaredDensity)
        thermal = { squaredDensity: value * value / (2 * h0), section, massFlux: value };
    }
  }
  const safetyFraction = .95;
  const boundarySquaredDensity = Math.max(rhoTotal ** 2, pressure.squaredDensity, sonic.squaredDensity);
  const logDensity = Math.log(Math.sqrt(boundarySquaredDensity / safetyFraction));
  const density = Math.exp(logDensity);
  const minimumCertifiedPressure = .5 * A * (density - pressure.squaredDensity / density);
  const maximumSpeedSquared = sonic.squaredDensity * sonicSpeedSquared / density ** 2;
  const maximumMachSquared = maximumSpeedSquared / ((gamma - 1) * (h0 - .5 * maximumSpeedSquared));
  if (![density, minimumCertifiedPressure, maximumMachSquared].every(Number.isFinite)
    || density <= 0 || minimumCertifiedPressure <= 0 || maximumMachSquared >= 1)
    throw new Error('No finite interior pressure-domain density seed.');
  let densityUnknowns = 0;
  for (let g = 0; g < tubes.length; g++) for (let i = 0; i < nx; i++) for (let j = 0; j < tubes[g]; j++) {
    const index = densityIndex(i, g, j);
    if (!Number.isInteger(index) || index < 0 || index >= initial.length)
      throw new Error('Invalid pressure-domain density state index.');
    initial[index] = logDensity; densityUnknowns++;
  }
  const diagnostics = { density, densityUnknowns, safetyFraction, boundarySquaredDensity,
    bounds: { pressure, sonic, thermal, stagnationSquaredDensity: rhoTotal ** 2 }, minimumCertifiedPressure,
    transportCheck: { biasedSections: 0, maxTransportDeparture: 0, maximumMachSquared } };
  if (upwind) {
    for (let g = 0; g < tubes.length; g++) for (let j = 0; j < tubes[g]; j++) {
      const chain = prepareStreamtubeTransportChain({ lower: nodes[g].map(row => row[j]),
        upper: nodes[g].map(row => row[j + 1]), densities: Array(nx).fill(density),
        massFlow: allocation.groups[g][j].massFlow, stagnationEnthalpy: h0, gamma, geometryDomain, upwind });
      chain.sections.forEach((section, i) => {
        const departure = Math.abs(chain.transportSpeeds[i] - section.q);
        if (departure !== 0) diagnostics.transportCheck.biasedSections++;
        diagnostics.transportCheck.maxTransportDeparture = Math.max(diagnostics.transportCheck.maxTransportDeparture, departure);
      });
    }
    if (diagnostics.transportCheck.biasedSections)
      throw new StreamtubePressureDomainTransportError(initial, diagnostics, system, nodes);
  }
  return { initial, diagnostics };
}

const origins = new WeakMap();
const sameNodes = (a, b) => Array.isArray(a) && Array.isArray(b) && a.length === b.length
  && a.every((group, g) => Array.isArray(group) && group.length === b[g].length
    && group.every((row, i) => Array.isArray(row) && row.length === b[g][i].length
      && row.every((p, j) => Object.is(p.x, b[g][i][j].x) && Object.is(p.y, b[g][i][j].y))));
const sameChart = (a, b) => Array.isArray(a) && a.length === b.length && a.every((p, i) =>
  p.column === b[i].column && ['offset', 'normal'].every(key =>
    Object.is(p[key].x, b[i][key].x) && Object.is(p[key].y, b[i][key].y)));

// This exception says a sufficient unfiltered pressure proof is inapplicable.
// It does not assert failure of the actual selected equations. Retain the
// already generated candidate; do not search for another density or disable
// the transport operator that invalidated the sufficient proof.
export class StreamtubePressureDomainTransportError extends Error {
  constructor(initial, diagnostics, system, nodes) {
    super('Pressure-domain density certificate requires unchanged physical transport speeds.');
    this.name = 'StreamtubePressureDomainTransportError';
    this.code = 'streamtube-pressure-domain-transport';
    this.diagnostics = structuredClone(diagnostics);
    this.candidate = { initial: Float64Array.from(initial), conditions: structuredClone(system.conditions) };
    // The existing initializer already has these physical nodes. Keep a
    // private origin record so a branded error cannot be transplanted to
    // another system or reused after geometry adoption/direction refresh.
    origins.set(this, { system, nodes: structuredClone(nodes),
      chart: typeof system.geometryChart === 'function' ? structuredClone(system.geometryChart()) : null,
      initial: Float64Array.from(initial), diagnostics: structuredClone(diagnostics),
      conditions: structuredClone(system.conditions) });
  }
}

export function evaluateStreamtubePressureDomainCandidate(system, state, failure) {
  if (!(failure instanceof StreamtubePressureDomainTransportError)) throw failure;
  const require = (ok, message) => { if (!ok) throw new Error(message); };
  const original = failure.candidate, d = failure.diagnostics, origin = origins.get(failure);
  require(origin?.system === system && Array.isArray(origin.chart)
    && typeof system.geometryChart === 'function' && sameChart(system.geometryChart(), origin.chart),
    'Pressure-domain candidate changed its originating system or geometry chart.');
  require(original.initial.length === origin.initial.length
    && original.initial.every((v, i) => Object.is(v, origin.initial[i]))
    && JSON.stringify(d) === JSON.stringify(origin.diagnostics)
    && JSON.stringify(original.conditions) === JSON.stringify(origin.conditions),
  'Pressure-domain generated candidate or provenance was modified.');
  require(system.conditions.flowModel === 'compressible' && system.conditions.upwind
    && JSON.stringify(system.conditions) === JSON.stringify(original.conditions),
  'Pressure-domain candidate changed its selected gas/transport conditions.');
  require(d.transportCheck?.biasedSections > 0 && Number.isInteger(d.transportCheck.biasedSections)
    && d.transportCheck.maxTransportDeparture > 0 && Number.isFinite(d.transportCheck.maxTransportDeparture)
    && d.density > 0 && Number.isFinite(d.density), 'Invalid transport-inapplicable density candidate.');
  const initial = original.initial.slice(), { nx, tubes, densityIndex } = system.layout;
  require(state.length === initial.length && initial.every(Number.isFinite) && state.every(Number.isFinite)
    && Number.isInteger(nx) && nx >= 2 && Array.isArray(tubes) && tubes.length > 0
    && tubes.every(n => Number.isInteger(n) && n > 0), 'Invalid pressure-domain candidate state/layout.');
  const densities = new Set();
  for (let g = 0; g < tubes.length; g++) for (let i = 0; i < nx; i++) for (let j = 0; j < tubes[g]; j++) {
    const column = densityIndex(i, g, j);
    require(Number.isInteger(column) && column >= 0 && column < initial.length && !densities.has(column),
      'Invalid pressure-domain candidate density index.');
    densities.add(column);
    require(Math.exp(initial[column]) === d.density, 'Pressure-domain candidate density changed.');
  }
  require(densities.size === d.densityUnknowns, 'Pressure-domain candidate density count changed.');
  for (let i = 0; i < initial.length; i++) if (!densities.has(i))
    require(Object.is(initial[i], state[i]), 'Pressure-domain candidate changed a non-density unknown.');
  require(sameNodes(system.decode(initial).nodes, origin.nodes),
    'Pressure-domain candidate changed its originating physical geometry.');

  // One actual selected-equation evaluation, with all its physical/remote
  // boundary guards. No observer is called here and no exception is caught.
  const flow = system.evaluate(initial), grid = assertConvexStreamtubeGrid(flow.nodes);
  let minStaticPressure = Infinity, minStaticEnthalpy = Infinity, minInterfacePressure = Infinity,
    maximumMachSquared = 0, sections = 0, cells = 0;
  for (const row of flow.sections) for (const group of row) for (const s of group) {
    require([s.p, s.enthalpy, s.rho, s.machSquared].every(Number.isFinite)
      && s.p > 0 && s.enthalpy > 0 && s.rho > 0 && s.machSquared >= 0,
    'Actual pressure-domain candidate has an invalid physical gas section.');
    minStaticPressure = Math.min(minStaticPressure, s.p);
    minStaticEnthalpy = Math.min(minStaticEnthalpy, s.enthalpy);
    maximumMachSquared = Math.max(maximumMachSquared, s.machSquared); sections++;
  }
  for (const row of flow.cells) for (const group of row) for (const cell of group) {
    const values = [cell.interfacePressure.lower, cell.interfacePressure.upper];
    require(values.every(p => Number.isFinite(p) && p > 0),
      'Actual pressure-domain candidate has nonpositive/nonfinite interface pressure.');
    minInterfacePressure = Math.min(minInterfacePressure, ...values); cells++;
  }
  require(sections === densities.size && cells === (nx - 1) * tubes.reduce((a, b) => a + b, 0)
    && flow.residual.length === initial.length && flow.residual.every(Number.isFinite)
    && maximumMachSquared < 1, 'Incomplete/nonfinite/non-subsonic pressure-domain candidate evaluation.');
  return { initial, flow, diagnostics: { method: 'pressure-domain-density-evaluated',
    initialGuessOnly: true, targetEquationsUnchanged: true, density: d.density, densityUnknowns: d.densityUnknowns,
    pressureCertificateApplicable: false, certificateInapplicable: { code: failure.code, message: failure.message,
      diagnostics: structuredClone(d) },
    actualEvaluation: { accepted: true, completeSelectedEquations: true, convexGrid: true,
      minStaticPressure, minStaticEnthalpy, minInterfacePressure, maximumMachSquared,
      sections, cells, minCornerSine: grid.minCornerSine, evaluations: 1,
      transportOperatorChanged: false, nondensityUnknownsChanged: false,
      originSystemAndChartExact: true, originPhysicalNodesExact: true } } };
}
