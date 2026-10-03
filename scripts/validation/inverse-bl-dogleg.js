// SPDX-License-Identifier: GPL-2.0-or-later
// Research trust model on the conditional inverse-BL manifold. All Euler
// rows and all surface/wake edge-matching rows remain. No normal equations.
import { createDoglegModel } from '../../src/numerics/dogleg.js';
import { sparseProduct } from '../../src/numerics/sparse.js';
import { reduceInverseBL } from './inverse-bl-schur.js';

export function createInverseBLDogleg(system, matrix, residual, newton) {
  const reduced = reduceInverseBL(system, matrix, residual.map(v => -v));
  const model = createDoglegModel(reduced.matrix, reduced.rhs.map(v => -v),
    Float64Array.from(reduced.columns, col => newton[col]));
  const merit = .5 * residual.reduce((sum, v) => sum + v * v, 0);
  const resolvedDecrease = 64 * Number.EPSILON * merit;
  const propose = radius => {
    const p = model.propose(radius);
    if (!p.direction) throw new Error(p.reason);
    return { ...p, reducedDirection: p.direction, direction: reduced.lift(p.direction) };
  };
  const assess = (retainedDirection, nextResidual, radius) => {
    if (nextResidual.length !== residual.length || !nextResidual.every(Number.isFinite)
      || !(radius > 0) || !Number.isFinite(radius)) throw new Error('Invalid inverse dogleg assessment.');
    // BL substitution includes its tiny remaining inhomogeneous defect.
    // Score against EVERY original equation, independently of Schur scoring.
    const jd = sparseProduct(matrix, reduced.lift(retainedDirection));
    const predictedReduction = -residual.reduce((sum, v, i) => sum + v * jd[i], 0)
      - .5 * jd.reduce((sum, v) => sum + v * v, 0);
    const actualReduction = .5 * residual.reduce((sum, v, i) => sum + (v - nextResidual[i]) * (v + nextResidual[i]), 0);
    const scaledStepNorm = Math.sqrt(retainedDirection.reduce((sum, v, i) => sum + (v * model.scales[i]) ** 2, 0));
    const ratio = predictedReduction > resolvedDecrease ? actualReduction / predictedReduction : null;
    const accepted = Number.isFinite(ratio) && ratio > 1e-4 && actualReduction > resolvedDecrease;
    let nextRadius = radius;
    if (!accepted || ratio < .25) nextRadius = .25 * radius;
    else if (ratio > .75 && scaledStepNorm >= .9 * radius) nextRadius = 2 * radius;
    return { accepted, radius, nextRadius, predictedReduction, actualReduction, ratio,
      scaledStepNorm, scaledNewtonNorm: model.newtonNorm, cauchyNorm: model.cauchyNorm };
  };
  return { propose, assess, initialRadius: Math.min(model.newtonNorm, model.cauchyNorm),
    reduced, model, merit, resolvedDecrease };
}
