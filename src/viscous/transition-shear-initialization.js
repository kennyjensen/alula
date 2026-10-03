// SPDX-License-Identifier: GPL-2.0-or-later
// Initial shear for a newly selected mixed transition interval. Fixed
// theta/delta/Ue/upstream N fix XT and the native shear row coefficients:
// R(S) = A - B*S - C*log(S/ST), B >= 0, C > 0. The positive root is unique.
// This local event preparation does not solve the coupled flow or alter its
// equations. In particular, mixed energy can change with the initialized S.
import { evaluateTransitionInterval } from './transition-interval.js';

export function initializeMixedTransitionShear(kernel, { upstream, downstream, tripS, expected }) {
  const tolerance = 1e-14, maximumExpansions = 160, maximumBisections = 80;
  let evaluations = 0, reference;
  const sample = logShear => {
    const aux = Math.exp(logShear);
    if (!(aux > 0) || !Number.isFinite(aux)) throw new Error('Mixed-transition shear bracket left the finite positive domain.');
    const value = evaluateTransitionInterval(kernel, { upstream, downstream: { ...downstream, aux }, tripS });
    evaluations++;
    const transition = value.transition;
    if (!transition || !Number.isFinite(value.residual[0])) throw new Error('Mixed-transition shear initialization has no finite selected row.');
    if (expected && (transition.s !== expected.s || transition.forced !== (expected.kind === 'forced')))
      throw new Error('Mixed-transition shear initialization changed the selected transition.');
    reference ??= { s: transition.s, forced: transition.forced };
    if (transition.s !== reference.s || transition.forced !== reference.forced)
      throw new Error('Mixed-transition shear bracket changed the transition location or branch.');
    return { logShear, aux, residual: value.residual[0] };
  };
  const finish = (point, method, iterations) => {
    if (!Number.isFinite(point.residual) || Math.abs(point.residual) > tolerance)
      throw new Error('Mixed-transition shear initialization did not close its native row.');
    return { aux: point.aux, diagnostics: { method, residual: point.residual, tolerance, iterations, evaluations,
      transition: reference, equationsChanged: false, initialGuessOnly: true } };
  };
  if (expected?.kind === 'forced' && expected.s === downstream.s) {
    const aux = kernel.station({ ...downstream, aux: .03 }, 'turbulent').transitionShear;
    if (!(aux > 0) || !Number.isFinite(aux)) throw new Error('Invalid native terminal-transition shear.');
    // Preserve the native expression's value rather than a log/exp roundtrip.
    const value = evaluateTransitionInterval(kernel, { upstream, downstream: { ...downstream, aux }, tripS });
    evaluations++;
    if (value.transition?.s !== downstream.s || value.transition.forced !== true)
      throw new Error('Terminal-transition shear initialization changed the selected transition.');
    reference = { s: value.transition.s, forced: true };
    return finish({ aux, residual: value.residual[0] }, 'native-terminal-shear', 0);
  }
  const seed = Number.isFinite(downstream.aux) && downstream.aux > 0 ? downstream.aux : .03;
  let low = sample(Math.log(seed)), high = low;
  if (Math.abs(low.residual) <= tolerance) return finish(low, 'native-mixed-shear-log-bisection', 0);
  for (let k = 0; low.residual < 0 && k < maximumExpansions; k++) low = sample(low.logShear - 2);
  for (let k = 0; high.residual > 0 && k < maximumExpansions; k++) high = sample(high.logShear + 2);
  if (!(low.residual >= 0 && high.residual <= 0)) throw new Error('Could not bracket the native mixed-transition shear root.');
  let best = Math.abs(low.residual) < Math.abs(high.residual) ? low : high;
  for (let k = 0; k < maximumBisections; k++) {
    const middle = (low.logShear + high.logShear) / 2;
    if (middle === low.logShear || middle === high.logShear) return finish(best, 'native-mixed-shear-log-bisection', k);
    const point = sample(middle);
    if (Math.abs(point.residual) < Math.abs(best.residual)) best = point;
    if (Math.abs(point.residual) <= tolerance) return finish(point, 'native-mixed-shear-log-bisection', k + 1);
    if (point.residual > 0) low = point; else high = point;
  }
  return finish(best, 'native-mixed-shear-log-bisection', maximumBisections);
}
