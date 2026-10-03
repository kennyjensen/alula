// SPDX-License-Identifier: GPL-2.0-or-later
import { solveSparseDirect } from '../numerics/klu.js';
import { sparseProduct } from '../numerics/sparse.js';
import { createDoglegModel } from '../numerics/dogleg.js';

// A feasible descent direction when damping the Newton ray cannot leave a
// domain boundary. Reuse the solver's physical and corner constraints in
// the column-scaled residual metric. This only proposes a direction: the
// ordinary nonlinear admission and merit tests still decide acceptance.
export function coupledConstrainedDirection(system, state, matrix, residual, { weights } = {}) {
  let modelMatrix = matrix, modelResidual = residual;
  if (weights) {
    modelMatrix = { ...matrix, values: matrix.values.slice() };
    modelResidual = residual.map((v, i) => v * weights[i]);
    for (let i = 0; i < matrix.n; i++) for (let p = matrix.rowPtr[i]; p < matrix.rowPtr[i + 1]; p++)
      modelMatrix.values[p] *= weights[i];
  }
  const model = createDoglegModel(modelMatrix, modelResidual, new Float64Array(matrix.n));
  let radius = 0;
  for (const v of modelResidual) radius = Math.hypot(radius, v);
  if (!(radius > 0)) return null;
  const candidate = model.proposeProjectedGradient(radius, system.stepConstraints(state));
  if (!candidate.direction) return null;
  const product = sparseProduct(modelMatrix, candidate.direction);
  const slope = modelResidual.reduce((sum, v, i) => sum + v * product[i], 0);
  return { direction: candidate.direction, diagnostics: {
    method: 'constrained-residual-gradient', predictedReduction: candidate.predictedReduction,
    linearizedMeritSlope: slope, beforeSquaredNorm: radius * radius,
    scaledStepNorm: candidate.scaledStepNorm, projection: candidate.projection, equationsChanged: false,
  } };
}

export function coupledNewtonRayStalled(history, proposal, event, decrease) {
  const recent = [...history.slice(-2), { ...proposal, activeChange: event.changed, residualDecrease: decrease }];
  return recent.length === 3 && recent.every(h => h.step > 0 && h.step < 1e-3
    && !h.activeChange && !h.directionRecovery && h.residualDecrease
    && h.residualDecrease.afterSquaredNorm > .99 * h.residualDecrease.beforeSquaredNorm);
}

// A bounded quasi-Newton search direction, never a change to the residual.
// Dampen local BL feedback while retaining all Euler/BL cross derivatives.
// Only directions descending the unregularized linear model can enter the ordinary
// physical, geometry-maintenance and nonlinear line search.
export function coupledRegularizedDirection(matrix, residual, ne, strength, { weights, preferredOrdering, pivotTolerance } = {}) {
  if (!Number.isInteger(ne) || ne < 0 || ne >= matrix.n || (matrix.n - ne) % 4
    || !Number.isFinite(strength) || strength <= 0) throw new Error('Invalid coupled direction regularization.');
  const shifted = { ...matrix, values: matrix.values.slice() };
  let changed = 0;
  for (let row = ne; row < matrix.n; row++) if ((row - ne) % 4 < 3) {
    for (let p = matrix.rowPtr[row]; p < matrix.rowPtr[row + 1]; p++) if (matrix.colIndex[p] === row) {
      shifted.values[p] *= 1 + strength;
      if (matrix.values[p] !== 0) changed++;
    }
  }
  if (!changed) return null;
  const linear = solveSparseDirect(shifted, residual.map(v => -v), { preferredOrdering, pivotTolerance });
  const jd = sparseProduct(matrix, linear.x);
  let slope = 0, squaredNorm = 0;
  for (let i = 0; i < matrix.n; i++) {
    const w = weights?.[i] ?? 1;
    slope += w * w * residual[i] * jd[i];
    squaredNorm += (w * residual[i]) ** 2;
  }
  if (!Number.isFinite(slope) || slope >= -1e-4 * squaredNorm) return null;
  return { linear, diagnostics: { method: 'regularized-bl-direction', strength, changedDiagonals: changed,
    linearizedMeritSlope: slope, beforeSquaredNorm: squaredNorm,
    relativeLinearResidual: linear.relativeResidual, equationsChanged: false } };
}

// After guarded directions have been tried, roundoff-sized accepted steps
// cannot justify spending the remaining startup budget on the same state.
export function coupledStepStagnation(history, tolerance) {
  const rows = history.slice(-3);
  if (rows.length !== 3 || rows.some(h => !(h.step > 0 && h.step < 1e-10)
    || h.activeChange || h.accepted === false || !(h.residual > tolerance)
    || !h.residualDecrease || !(h.residualDecrease.afterSquaredNorm >= .999 * h.residualDecrease.beforeSquaredNorm))) return null;
  return { iterations: rows.map(h => h.iteration), steps: rows.map(h => h.step),
    residual: rows.at(-1).residual, equationsChanged: false };
}
