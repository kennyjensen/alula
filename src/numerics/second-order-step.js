// SPDX-License-Identifier: GPL-2.0-or-later
// Bounded minimum-norm correction for nonlinear constraint curvature.
// Frozen constraint normals are pulled into the caller's D metric. This
// only constructs a proposal; actual residual/admissibility acceptance is
// still required. No accepted state, active regime or mesh is modified.
import { projectHalfspaces } from './halfspace-projection.js';

const norm = x => { let n = 0; for (const v of x) n = Math.hypot(n, v); return n; };
export function correctSecondOrderStep({ initial, direction, scales, radius, constraints, values,
  maxCorrections = 4, maximumCorrectionFraction = .25 }) {
  const n = initial.length;
  if (!n || direction.length !== n || scales.length !== n || !initial.every(Number.isFinite)
    || !direction.every(Number.isFinite) || !scales.every(v => Number.isFinite(v) && v > 0)
    || !(radius > 0) || !Number.isFinite(radius) || typeof values !== 'function'
    || !Number.isInteger(maxCorrections) || maxCorrections < 1
    || !(maximumCorrectionFraction > 0 && maximumCorrectionFraction < 1)) throw new Error('Invalid second-order step controls.');
  const original = Float64Array.from(direction, (v, i) => v * scales[i]), originalNorm = norm(original);
  if (!(originalNorm > 0) || originalNorm > radius * (1 + 1e-11)) throw new Error('Second-order proposal exceeds its trust radius.');
  const rows = constraints.map(({ gradient, value, lower }) => {
    if (!(gradient instanceof Map) || !(value > 0) || !Number.isFinite(value) || !Number.isFinite(lower)
      || !(value + lower > 0)) throw new Error('Second-order correction requires resolved positive initial margins.');
    let prediction = value, arithmeticScale = value;
    for (const [j, d] of gradient) { prediction += d * direction[j]; arithmeticScale += Math.abs(d * direction[j]); }
    if (prediction < value + lower - 64 * Number.EPSILON * arithmeticScale)
      throw new Error('Second-order correction requires a linearly feasible proposal.');
    return { target: value + lower, gradient: new Map([...gradient].map(([j, d]) => [j, d / scales[j]])) };
  });
  let step = original.slice(), correctionNorm = 0, evaluations = 0;
  const history = [], maximumCorrection = maximumCorrectionFraction * originalNorm;
  for (let iteration = 0; iteration <= maxCorrections; iteration++) {
    let current;
    try { current = values(initial.map((v, i) => v + step[i] / scales[i])); evaluations++; }
    catch (error) { return { corrected: false, reason: error.message, history, evaluations }; }
    if (current.length !== rows.length || !current.every(Number.isFinite))
      return { corrected: false, reason: 'Invalid second-order constraint values.', history, evaluations };
    const violated = current.reduce((sum, v) => sum + (v <= 0 ? 1 : 0), 0);
    if (!violated) return iteration ? { corrected: true, direction: step.map((v, i) => v / scales[i]),
      correctionNorm, scaledStepNorm: norm(step), originalNorm, history, evaluations }
      : { corrected: false, reason: 'The rejected trial has no supported nonlinear constraint violation.', history, evaluations };
    if (iteration === maxCorrections) break;
    let projection;
    try { projection = projectHalfspaces(new Float64Array(n), rows.map((r, i) => ({ gradient: r.gradient, lower: r.target - current[i] })),
      { maximumNorm: maximumCorrection }); }
    catch (error) { return { corrected: false, reason: error.message, history, evaluations }; }
    const { point: change, ...diagnostics } = projection;
    if (!projection.converged) return { corrected: false, reason: 'Second-order correction projection failed.', history, evaluations, projection: diagnostics };
    const proposed = step.map((v, i) => v + change[i]), proposedNorm = norm(proposed), fraction = Math.min(1, radius / proposedNorm);
    step = proposed.map(v => fraction * v);
    correctionNorm = norm(step.map((v, i) => v - original[i]));
    history.push({ violated, correctionNorm, radiusFraction: fraction, projection: diagnostics });
    if (correctionNorm > maximumCorrection * (1 + 1e-11))
      return { corrected: false, reason: 'Second-order correction exceeds its step bound.', history, evaluations };
  }
  return { corrected: false, reason: 'Second-order correction did not recover positive nonlinear margins.', history, evaluations };
}
