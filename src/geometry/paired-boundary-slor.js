// SPDX-License-Identifier: GPL-2.0-or-later
// Paired boundary-row SLOR with direct source elimination and Armijo
// backtracking inside each Newton line. The undamped whole-grid residual
// determines convergence; geometry stays admissible at every accepted row.
import { createWasmSmoother } from './wasm-smoothing.js';
import { slorErrorTermination, tagSlorTermination } from './slor-termination.js';
export function smoothPairedBoundaryGrid(system, { initial = system?.initial, maxSweeps = 600, tolerance = 1e-9, onSweep, backend = 'wasm' } = {}) {
  if (system?.coordinateEquations?.lineGrouping !== 'boundary-pairs' || system.coordinateEquations.lineSearch !== 'armijo'
    || system.coordinateEquations.controlUpdate || !Number.isInteger(maxSweeps) || maxSweeps < 0
    || !(tolerance > 0) || !Number.isFinite(tolerance))
    throw new Error('Paired SLOR requires direct boundary controls, Newton line backtracking and valid settings.');
  if (!['wasm', 'javascript'].includes(backend)) throw new Error('Unknown smoothing backend.');
  const kernel = backend === 'wasm' ? createWasmSmoother(system) : null;
  let nodes = structuredClone(initial), reason = 'sweep limit', converged = false, residual;
  let termination = { origin: 'solver', termination: 'sweep-limit' };
  const history = [];
  for (let iteration = 0; iteration <= maxSweeps; iteration++) {
    const quality = system.quality(nodes), row = { iteration, minCornerSine: quality.minCornerSine, invalidCells: quality.invalidCells.length };
    history.push(row);
    if (!quality.valid) { reason = 'Invalid starting state: folded grid'; termination = { origin: 'solver', termination: 'invalid-grid' }; break; }
    try {
      const state = kernel ? kernel.evaluate(nodes) : system.residuals(nodes); residual = state.residual;
      row.residual = residual; row.merit = kernel ? state.merit : .5 * state.rows.reduce((s, p) => s + p.x ** 2 + p.y ** 2, 0);
      if (!Number.isFinite(residual) || !Number.isFinite(row.merit)) throw new Error('Nonfinite paired SLOR residual or merit.');
      onSweep?.({ ...row }, structuredClone(nodes));
      if (residual <= tolerance) { converged = true; reason = 'converged'; break; }
      if (iteration === maxSweeps) break;
      // A line already below one tenth of the global target need not
      // decrease sub-roundoff corrections while other lines catch up.
      const next = (kernel ?? system).sweep(nodes, 1, { lineTolerance: .1 * tolerance });
      nodes = next.nodes; row.maxUpdate = next.maxUpdate; row.rowSteps = next.rowSteps;
    } catch (error) {
      if (error?.code === 'slor-observer-failed') throw error;
      reason = error.message; termination = slorErrorTermination(error); break;
    }
  }
  const result = tagSlorTermination({ nodes, history, converged, reason, residual, quality: system.quality(nodes), controls: { maxSweeps, tolerance, lineTolerance: .1 * tolerance },
    coordinateEquations: system.coordinateEquations,
    formulation: 'Paired boundary rows, full line metrics, direct boundary controls and row-wise Armijo backtracking; fixed boundaries and harmonic mass' }, termination);
  Object.defineProperty(result, 'backend', { value: kernel ? 'wasm' : 'javascript' });
  Object.defineProperty(result, 'referenceReplays', { value: kernel?.replays ?? 0 });
  return result;
}
