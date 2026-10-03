// SPDX-License-Identifier: GPL-2.0-or-later
// Scaled dogleg trust-region globalization of a square residual system.
// The Newton direction is supplied by the caller's certified linear solve.
// This module combines it with the residual's scaled steepest-descent
// direction; it does not change the equations or certify a nonlinear root.
import { projectHalfspaces } from './halfspace-projection.js';
import { correctSecondOrderStep } from './second-order-step.js';
const dot = (a, b) => a.reduce((sum, v, i) => sum + v * b[i], 0);
const norm = a => { let value = 0; for (const v of a) value = Math.hypot(value, v); return value; };

export function createDoglegModel(matrix, residual, newtonDirection) {
  const n = residual.length, sparse = matrix.rowPtr !== undefined;
  if (!n || newtonDirection.length !== n || !residual.every(Number.isFinite) || !newtonDirection.every(Number.isFinite)
    || (sparse ? matrix.n !== n || matrix.rowPtr.length !== n + 1 || matrix.rowPtr[0] !== 0
      || matrix.rowPtr[n] !== matrix.values.length || matrix.colIndex.length !== matrix.values.length
      || !matrix.values.every(Number.isFinite) : matrix.length !== n * n || !matrix.every(Number.isFinite)))
    throw new Error('Invalid dogleg linear model.');
  const visit = callback => {
    for (let row = 0; row < n; row++) {
      if (sparse) {
        if (!Number.isInteger(matrix.rowPtr[row + 1]) || matrix.rowPtr[row + 1] < matrix.rowPtr[row]) throw new Error('Invalid dogleg sparse rows.');
        for (let k = matrix.rowPtr[row]; k < matrix.rowPtr[row + 1]; k++) {
          const col = matrix.colIndex[k]; if (!Number.isInteger(col) || col < 0 || col >= n) throw new Error('Invalid dogleg sparse column.');
          callback(row, col, matrix.values[k]);
        }
      } else for (let col = 0; col < n; col++) callback(row, col, matrix[row * n + col]);
    }
  };
  const scales = new Float64Array(n), gradient = new Float64Array(n);
  visit((row, col, value) => { scales[col] = Math.hypot(scales[col], value); gradient[col] += value * residual[row]; });
  if (!scales.every(v => v > 0 && Number.isFinite(v))) throw new Error('Dogleg model has an unresolved Jacobian column.');
  const product = step => { const value = new Float64Array(n); visit((i, j, a) => { value[i] += a * step[j]; }); return value; };
  // z = D*dx, A = J/D. This metric is invariant to changes of units of
  // individual state variables. No J^T*J matrix is formed or factorized.
  const scaledGradient = gradient.map((v, i) => v / scales[i]);
  const gradientNorm = norm(scaledGradient), ag = product(scaledGradient.map((v, i) => v / scales[i])), agNorm = norm(ag);
  const scaledNewton = Float64Array.from(newtonDirection, (v, i) => v * scales[i]), newtonNorm = norm(scaledNewton);
  const alpha = gradientNorm > 0 && agNorm > 0 ? (gradientNorm / agNorm) ** 2 : 0;
  const cauchy = scaledGradient.map(v => -alpha * v), cauchyNorm = norm(cauchy);
  if (![gradientNorm, agNorm, alpha, newtonNorm, cauchyNorm].every(Number.isFinite)) throw new Error('Nonfinite dogleg model.');
  const score = (z, kind) => {
    const direction = z.map((v, i) => v / scales[i]), jd = product(direction);
    const predictedReduction = -dot(residual, jd) - .5 * dot(jd, jd), scaledStepNorm = norm(z);
    if (!(predictedReduction > 0) || !Number.isFinite(predictedReduction) || !direction.every(Number.isFinite))
      return { direction: null, reason: 'No resolved dogleg model decrease.' };
    return { direction, kind, predictedReduction, scaledStepNorm,
      linearizedResidualNorm: norm(residual.map((v, i) => v + jd[i])) };
  };
  const validRadius = radius => {
    if (!(radius > 0) || !Number.isFinite(radius)) throw new Error('Invalid dogleg trust radius.');
  };
  const propose = radius => {
    validRadius(radius);
    if (newtonNorm <= radius) return score(scaledNewton.slice(), 'newton');
    if (!(cauchyNorm > 0)) return { direction: null, reason: 'Stationary linearized merit; residual convergence unproven.' };
    if (cauchyNorm >= radius) return score(cauchy.map(v => radius * v / cauchyNorm), 'gradient');
    const change = scaledNewton.map((v, i) => v - cauchy[i]);
    const a = dot(change, change), b = dot(cauchy, change), c = (cauchyNorm - radius) * (cauchyNorm + radius);
    const root = Math.sqrt(b * b - a * c);
    const t = b >= 0 ? -c / (b + root) : (-b + root) / a;
    return score(cauchy.map((v, i) => v + t * change[i]), 'dogleg');
  };
  // Near an admissibility boundary, the Cauchy direction can point out of
  // the domain even though Newton points inward. Its truncated ray is an
  // alternative model-descent direction inside the SAME trust region.
  // A good measured reduction is still required; no geometry is projected.
  const proposeNewtonRay = radius => {
    validRadius(radius);
    return score(scaledNewton.map(v => Math.min(1, radius / newtonNorm) * v), 'newton-ray');
  };
  // An optional feasible linearized Cauchy direction for domain boundaries.
  // Project -g in the same D metric, then minimize the residual model along
  // that segment. Zero belongs to every halfspace, so scaling toward zero
  // preserves linearized feasibility and the trust-radius bound. Nonlinear
  // geometry/admissibility and actual merit must still be checked by callers.
  const projectAndMinimize = (target, constraints, kind) => {
    const projection = projectHalfspaces(target, constraints.map(({ gradient, lower }) => ({ lower,
      gradient: new Map([...gradient].map(([j, v]) => [j, v / scales[j]])) })));
    const { point: z, ...diagnostics } = projection;
    if (!projection.converged) return { direction: null, reason: 'Linearized constraint projection did not converge.', projection: diagnostics };
    const jd = product(z.map((v, i) => v / scales[i])), slope = dot(scaledGradient, z), curvature = dot(jd, jd);
    if (!(slope < 0) || !(curvature > 0)) return { direction: null, reason: 'No projected model-descent direction.', projection: diagnostics };
    const fraction = Math.min(1, -slope / curvature);
    return { ...score(z.map(v => fraction * v), kind), projection: diagnostics };
  };
  const proposeProjectedGradient = (radius, constraints) => {
    validRadius(radius);
    if (!(gradientNorm > 0)) return { direction: null, reason: 'Stationary linearized merit; residual convergence unproven.' };
    return projectAndMinimize(scaledGradient.map(v => -radius * v / gradientNorm), constraints, 'projected-gradient');
  };
  // A metric projection of the truncated Newton direction is an additional
  // proposal, not the solution of a constrained least-squares subproblem.
  // Projection can destroy descent; the full residual model must score it.
  const proposeProjectedNewton = (radius, constraints) => {
    validRadius(radius);
    const fraction = Math.min(1, radius / newtonNorm);
    return projectAndMinimize(scaledNewton.map(v => fraction * v), constraints, 'projected-newton');
  };
  const scoreDirection = (direction, kind) => {
    if (direction.length !== n || !direction.every(Number.isFinite)) throw new Error('Invalid corrected model direction.');
    return score(Float64Array.from(direction, (v, i) => v * scales[i]), kind);
  };
  const proposeConstrainedDogleg = (radius, constraints, cauchy) => {
    validRadius(radius);
    if (!cauchy?.direction) return { direction: null };
    const end = propose(radius);
    if (!end.direction) return { direction: null };
    // Both endpoints are in the scaled trust ball. Starting at the feasible
    // projected Cauchy point, stop its segment toward ordinary dogleg at the
    // first linearized face. Convexity preserves every halfspace and the
    // trust bound. Minimize the full quadratic model on this bounded segment.
    const change = end.direction.map((v, i) => v - cauchy.direction[i]);
    let maximumFraction = 1, limitingConstraint = null;
    constraints.forEach(({ gradient, lower }, id) => {
      let slack = -lower, slope = 0;
      for (const [j, a] of gradient) { slack += a * cauchy.direction[j]; slope += a * change[j]; }
      if (slope < 0) {
        const bound = Math.max(0, slack) / -slope;
        if (bound < maximumFraction) { maximumFraction = bound; limitingConstraint = id; }
      }
    });
    if (!(maximumFraction > 0)) return { direction: null };
    const jc = product(cauchy.direction), jd = product(change), curvature = dot(jd, jd);
    const slope = dot(residual.map((v, i) => v + jc[i]), jd);
    if (!(curvature > 0) || !(slope < 0)) return { direction: null };
    const fraction = Math.min(maximumFraction, -slope / curvature);
    const direction = cauchy.direction.map((v, i) => v + fraction * change[i]);
    const proposal = scoreDirection(direction, 'constrained-dogleg');
    // Do not add a duplicate proposal at a roundoff-sized segment. The
    // original Cauchy step remains available if the extension is rejected.
    if (!proposal.direction || !(proposal.predictedReduction > cauchy.predictedReduction * (1 + 64 * Number.EPSILON)))
      return { direction: null };
    return { ...proposal, projection: cauchy.projection, doglegSegment: { fraction, maximumFraction,
      limitingConstraint, cauchyPredictedReduction: cauchy.predictedReduction } };
  };
  return { scales, gradientNorm, newtonNorm, cauchyNorm, propose, proposeNewtonRay, proposeProjectedGradient,
    proposeProjectedNewton, proposeConstrainedDogleg, scoreDirection };
}

export function takeDoglegStep({ initial, currentResidual, matrix, newtonDirection, residual, admissible = () => true,
  radius = 1, maximumRadius = 1e6, maxTrials = 30, activeSet, linearizedConstraints, constraintValues }) {
  if (initial.length !== currentResidual.length || !initial.every(Number.isFinite) || !(radius > 0) || !Number.isFinite(radius)
    || !Number.isFinite(maximumRadius) || maximumRadius < radius || !Number.isInteger(maxTrials) || maxTrials < 1)
    throw new Error('Invalid dogleg step controls.');
  if (activeSet && !['snapshot', 'restore', 'prepare'].every(k => typeof activeSet[k] === 'function'))
    throw new Error('Invalid dogleg active-set callbacks.');
  if (linearizedConstraints !== undefined && typeof linearizedConstraints !== 'function')
    throw new Error('Invalid dogleg constraint callback.');
  if (constraintValues !== undefined && (typeof constraintValues !== 'function' || !linearizedConstraints))
    throw new Error('Nonlinear constraint values require linearized constraint callbacks.');
  const model = createDoglegModel(matrix, currentResidual, newtonDirection), trials = [];
  const merit = .5 * dot(currentResidual, currentResidual), resolvedDecrease = 64 * Number.EPSILON * merit;
  const saved = activeSet?.snapshot(); let keepActive = false, constraints;
  try {
  while (trials.length < maxTrials) {
    if (!(radius > 0)) return { accepted: false, reason: 'Unresolved trust radius; residual convergence unproven.', radius, trials };
    let proposal = model.propose(radius), domainRejected = false, projectedProposals, gradientProposal;
    const trialRadius = radius;
    for (let candidate = 0; candidate < (constraintValues ? 6 : linearizedConstraints ? 5 : 2) && trials.length < maxTrials; candidate++) {
      if (activeSet) activeSet.restore(saved);
      if (candidate === 1) {
        if (proposal.kind === 'newton') { if (linearizedConstraints) continue; else break; }
        proposal = model.proposeNewtonRay(trialRadius);
      }
      if (candidate >= 2 && candidate <= 4) {
        // Projection repairs a domain-obstructed direction. If the ordinary
        // proposals were admissible but their nonlinear model was poor,
        // contract the trust radius instead. An active constraint on an
        // unrelated gradient proposal must not preempt that contraction.
        if (!domainRejected) break;
        try {
          constraints ??= linearizedConstraints(initial);
          if (!projectedProposals) {
            const newton = { kind: 'projected-newton', ...model.proposeProjectedNewton(trialRadius, constraints) };
            // Use the residual model's Cauchy length, bounded by the trust
            // radius, before projection and segment minimization. Projecting
            // a radius-sized gradient can activate distant bounds and fail
            // even when the much shorter Cauchy point is already feasible.
            const cauchyRadius = Math.min(trialRadius, model.cauchyNorm);
            gradientProposal = { kind: 'projected-gradient', ...model.proposeProjectedGradient(cauchyRadius, constraints) };
            const newtonActive = newton.projection?.active > 0, gradientActive = gradientProposal.projection?.active > 0;
            // Empty/inactive constraints preserve ordinary globalization.
            // If Newton's projection IS active, a feasible unprojected Cauchy
            // point is a valid competing proposal, labeled accordingly.
            if (!newtonActive && !gradientActive) break;
            if (!newtonActive && newton.direction) newton.direction = null;
            if (!gradientActive && gradientProposal.direction) gradientProposal.kind = 'cauchy';
            const dogleg = model.proposeConstrainedDogleg(trialRadius, constraints, gradientProposal);
            projectedProposals = [newton, ...(dogleg.direction ? [dogleg] : []), gradientProposal];
            // A projected Newton direction has no Cauchy-decrease guarantee.
            // Try the better model prediction first, preserving Newton order
            // for roundoff ties. Nonlinear failures still try the alternative
            // and then contract; model dominance is not actual acceptance.
            projectedProposals.sort((a, b) => {
              if (!a.direction) return b.direction ? -1 : 0;
              if (!b.direction) return 1;
              if (Math.abs(a.predictedReduction - b.predictedReduction)
                <= 64 * Number.EPSILON * Math.max(a.predictedReduction, b.predictedReduction)) return 0;
              return b.predictedReduction - a.predictedReduction;
            });
          }
          proposal = projectedProposals[candidate - 2];
        } catch (error) {
          return { accepted: false, reason: `Linearized constraint assembly failed: ${error.message}`, radius, trials };
        }
        if (!proposal) continue;
        if (!proposal.direction) {
          if (proposal.reason) trials.push({ radius: trialRadius, kind: proposal.kind, reason: proposal.reason, projection: proposal.projection });
          continue;
        }
      }
      if (candidate === 5) {
        // Sorting the proposals must not change which direction receives the
        // existing second-order correction: it is the Cauchy proposal.
        proposal = gradientProposal;
        if (!proposal?.direction) continue;
        const correctedKind = `${proposal.kind}-soc`;
        try {
          // Discard any auxiliary remapping made by the rejected trial.
          // The correction starts from the original continuous proposal;
          // active-event preparation is rerun on its corrected candidate.
          const corrected = correctSecondOrderStep({ initial, direction: proposal.direction, scales: model.scales,
            radius: trialRadius, constraints, values: constraintValues });
          const { direction, ...correction } = corrected;
          if (!corrected.corrected) { trials.push({ radius: trialRadius, kind: correctedKind, reason: corrected.reason, correction }); continue; }
          proposal = { ...model.scoreDirection(direction, correctedKind), projection: proposal.projection, correction };
          // A feasible curvature correction can lose model descent. Reject
          // this candidate and retry a smaller region, not the whole solve.
          if (!proposal.direction) {
            trials.push({ radius: trialRadius, kind: correctedKind, reason: proposal.reason,
              projection: proposal.projection, correction }); continue;
          }
        } catch (error) { trials.push({ radius: trialRadius, kind: correctedKind, reason: error.message }); continue; }
      }
      if (!proposal.direction) return { accepted: false, reason: proposal.reason, radius, trials };
      if (proposal.predictedReduction <= resolvedDecrease) {
        if (candidate >= 2) {
          trials.push({ radius: trialRadius, kind: proposal.kind, reason: 'Projected model decrease is at roundoff.' }); continue;
        }
        return { accepted: false, reason: 'Trust-region decrease is at roundoff; residual convergence unproven.', radius, trials };
      }
      const x = Float64Array.from(initial, (v, i) => v + proposal.direction[i]);
      let next, failure, event;
      try {
        event = activeSet?.prepare(x, initial);
        if (!x.every(Number.isFinite) || !admissible(x)) throw new Error('Inadmissible trust-region trial.');
        next = Float64Array.from(residual(x));
        if (next.length !== initial.length || !next.every(Number.isFinite)) throw new Error('Invalid trust-region residual.');
      } catch (error) { failure = error.message; }
      if (failure) { domainRejected = true; trials.push({ radius: trialRadius, kind: proposal.kind, reason: failure }); continue; }
      if (event?.changed) {
        // An active event can change an auxiliary variable's meaning (N vs
        // turbulent shear). The old linear model cannot score that remapping.
        // The caller's prepare callback must bound physical event travel;
        // accept an admissible event, keep its new active set and rebuild the
        // model on the next iteration. Never report a fictitious merit ratio.
        trials.push({ radius: trialRadius, kind: proposal.kind, activeChange: true, changes: event.changes });
        keepActive = true;
        return { accepted: true, x, residual: next, radius: trialRadius, trialRadius, kind: proposal.kind,
          activeChange: true, changes: event.changes, meritComparable: false,
          ratio: null, actualReduction: null, predictedReduction: null, linearizedResidualNorm: null,
          scaledStepNorm: proposal.scaledStepNorm, scaledNewtonNorm: model.newtonNorm,
          ...(proposal.projection ? { projection: proposal.projection } : {}),
          ...(proposal.doglegSegment ? { doglegSegment: proposal.doglegSegment } : {}),
          ...(proposal.correction ? { correction: proposal.correction } : {}), trials };
      }
      // Difference of squared norms without subtracting two large sums.
      const actualReduction = .5 * currentResidual.reduce((sum, v, i) => sum + (v - next[i]) * (v + next[i]), 0);
      const ratio = actualReduction / proposal.predictedReduction;
      trials.push({ radius, kind: proposal.kind, ratio, actualReduction, predictedReduction: proposal.predictedReduction });
      if (Number.isFinite(ratio) && ratio > 1e-4 && actualReduction > resolvedDecrease) {
        if (ratio < .25) radius *= .25;
        else if (ratio > .75 && proposal.scaledStepNorm >= .9 * radius) radius = Math.min(maximumRadius, 2 * radius);
        keepActive = true;
        return { accepted: true, x, residual: next,
        radius, trialRadius, kind: proposal.kind, ratio, actualReduction, predictedReduction: proposal.predictedReduction,
        scaledStepNorm: proposal.scaledStepNorm, scaledNewtonNorm: model.newtonNorm,
        linearizedResidualNorm: proposal.linearizedResidualNorm,
        ...(proposal.projection ? { projection: proposal.projection } : {}),
        ...(proposal.doglegSegment ? { doglegSegment: proposal.doglegSegment } : {}),
        ...(proposal.correction ? { correction: proposal.correction } : {}), trials };
      }
    }
    radius *= .25;
  }
  return { accepted: false, reason: 'Trust-region trial limit; residual convergence unproven.', radius, trials };
  } finally { if (activeSet && !keepActive) activeSet.restore(saved); }
}
