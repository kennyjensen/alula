// SPDX-License-Identifier: GPL-2.0-or-later
import { initializeStreamtubeDensities } from './streamtube-initial-state.js';
import { initializeStreamtubePressureDomain } from './streamtube-pressure-domain-initial-state.js';
import { StreamtubePressureDomainTransportError, evaluateStreamtubePressureDomainCandidate } from './streamtube-pressure-domain-initial-state.js';

function startupDiagnostics(system, flow, method) {
  const { gamma, pInf } = system.conditions;
  let sectionCount = 0, maximumAbsoluteEntropyDeparture = 0;
  for (const row of flow.sections) for (const group of row) for (const section of group) {
    // Dimensionless entropy relative to the target freestream, scaled by cv.
    // It is generally nonzero for the stagnation-density Newton seed.
    const departure = Math.log(section.p / pInf) - gamma * Math.log(section.rho);
    if (!Number.isFinite(departure)) throw new Error('Nonfinite startup entropy diagnostic.');
    maximumAbsoluteEntropyDeparture = Math.max(maximumAbsoluteEntropyDeparture, Math.abs(departure));
    sectionCount++;
  }
  return { method, fallbackAttempted: method !== 'isentropic', initialGuessOnly: true, targetEquationsUnchanged: true,
    residual: flow.diagnostics.residual,
    residualByFamily: { ...flow.diagnostics.residualByFamily },
    maxMach: flow.diagnostics.maxMach,
    maxEntropyJump: flow.diagnostics.maxEntropyJump,
    maxStagnationPressureError: flow.diagnostics.maxStagnationPressureError,
    entropy: { definition: 'log(p/p_inf) - gamma*log(rho/rho_inf), with rho_inf=1',
      sectionCount, maximumAbsoluteDeparture: maximumAbsoluteEntropyDeparture } };
}

const failureSnapshot = error => ({ code: error.code, message: error.message,
  diagnostics: structuredClone(error.diagnostics) });
const rejectedCandidate = (method, error) => ({ attempted: true, method, admissible: false,
  reason: error.message, ...(error.code ? { code: error.code } : {}),
  ...(error.diagnostics === undefined ? {} : { diagnostics: structuredClone(error.diagnostics) }) });

function pressureDomainStartup(system, state, originalError, previous = {}) {
  const originalFailure = failureSnapshot(originalError);
  let candidate, certificateInapplicable;
  try {
    try { candidate = initializeStreamtubePressureDomain(system, state); }
    catch (error) {
      if (!(error instanceof StreamtubePressureDomainTransportError)) throw error;
      certificateInapplicable = failureSnapshot(error);
      candidate = evaluateStreamtubePressureDomainCandidate(system, state, error);
      return { initial: candidate.initial, flow: candidate.flow, diagnostics: {
        ...startupDiagnostics(system, candidate.flow, 'pressure-domain-density-evaluated'), ...candidate.diagnostics,
        strict: originalFailure.diagnostics, isentropicFailure: originalFailure, ...previous } };
    }
    // The certificate concerns local side pressures only. All selected Euler
    // rows, remote boundary conditions and physical/grid gates still run.
    const flow = system.evaluate(candidate.initial);
    return { initial: candidate.initial, flow, diagnostics: {
      ...startupDiagnostics(system, flow, 'pressure-domain-density'), ...candidate.diagnostics,
      strict: originalFailure.diagnostics, isentropicFailure: originalFailure, ...previous } };
  } catch (error) {
    originalError.diagnostics = { ...originalError.diagnostics, fallbackAttempted: true, ...previous,
      pressureDomainDensityFallback: { ...rejectedCandidate('pressure-domain-density', error),
        ...(certificateInapplicable ? { certificateInapplicable } : {}),
        ...(candidate ? { candidate: structuredClone(candidate.diagnostics) } : {}) } };
    throw originalError;
  }
}

// Exact freestream-entropy section inversion is the preferred seed. A valid
// initial grid can exceed that branch's fixed-geometry mass-flow capacity
// while still admitting a subsonic, nonzero-residual Newton iterate. Giles's
// ISET listing (thesis p. 172) uses stagnation density as its initial estimate.
// This alternate estimate changes densities alone; entropy and inlet rows,
// captured masses, geometry and total enthalpy retain their target equations.
export function initializeStreamtubeStartup(system, state) {
  let initial, isentropicFailure;
  try { initial = initializeStreamtubeDensities(system, state); }
  catch (error) {
    // Geometry errors and other invalid inputs retain the strict rejection.
    if (error.code !== 'streamtube-sonic-capacity') throw error;
    isentropicFailure = error;
  }
  if (!isentropicFailure) {
    try {
      const flow = system.evaluate(initial);
      return { initial, flow, diagnostics: startupDiagnostics(system, flow, 'isentropic') };
    } catch (error) {
      // A convex, locally subsonic isentropic seed can still require negative
      // side pressure to balance its grid curvature. Only that typed failure
      // enables the fixed-grid, non-isentropic pressure-domain initial guess.
      if (error.code !== 'streamtube-interface-pressure') throw error;
      return pressureDomainStartup(system, state, error);
    }
  }

  const originalFailure = { code: isentropicFailure.code, message: isentropicFailure.message,
    diagnostics: structuredClone(isentropicFailure.diagnostics) };
  try {
    const density = system.conditions.rhoTotal;
    if (!(density > 0) || !Number.isFinite(density)) throw new Error('Invalid physical stagnation density.');
    initial = Float64Array.from(state);
    const { nx, tubes, densityIndex } = system.layout;
    let densityUnknowns = 0;
    for (let g = 0; g < tubes.length; g++) for (let i = 0; i < nx; i++) for (let j = 0; j < tubes[g]; j++) {
      initial[densityIndex(i, g, j)] = Math.log(density); densityUnknowns++;
    }
    // The complete Euler evaluation enforces positive pressure/enthalpy,
    // subsonic flow, geometry and all other existing admissibility checks.
    const flow = system.evaluate(initial);
    const diagnostics = { ...startupDiagnostics(system, flow, 'stagnation-density'),
      density, densityUnknowns, strict: originalFailure.diagnostics, isentropicFailure: originalFailure };
    return { initial, flow, diagnostics };
  } catch (error) {
    const stagnationDensityFallback = rejectedCandidate('stagnation-density', error);
    // Severe fixed-grid mass/area compression can exhaust H0 before side
    // pressures are evaluated. The same density certificate bounds physical
    // section speed below sonic speed, hence also guarantees positive h.
    // Only these typed gas-domain failures authorize another initial guess.
    if (error.code === 'streamtube-interface-pressure' || error.code === 'streamtube-static-enthalpy')
      return pressureDomainStartup(system, state, isentropicFailure, { stagnationDensityFallback });
    // Preserve the original typed error and every original diagnostic. The
    // rejected alternative remains an explicit diagnostic, never a success.
    isentropicFailure.diagnostics = { ...isentropicFailure.diagnostics, fallbackAttempted: true,
      stagnationDensityFallback };
    throw isentropicFailure;
  }
}
