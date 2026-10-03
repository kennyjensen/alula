// SPDX-License-Identifier: GPL-2.0-or-later
// Native BLDIF shear initialization for an event's affected turbulent tail/wake.
// With upstream Ctau and both stations' theta/delta/Ue fixed,
// R(S) = A - B*S - C*log(S/S1), B >= 0, C > 0. Its positive root is unique.
// This prepares a starting value; all three BL equations remain in Newton.
export function initializeTransportShear(kernel, { upstream, downstream, regime }) {
  if (!['turbulent', 'wake'].includes(regime)) throw new Error('Shear transport initialization requires a turbulent or wake interval.');
  if (![upstream?.aux, downstream?.aux].every(s => Number.isFinite(s) && s > 0))
    throw new Error('Shear transport initialization requires positive finite upstream and downstream shear.');
  const tolerance = 1e-14, maximumExpansions = 160, maximumBisections = 80;
  let evaluations = 0;
  const sample = (logShear, aux = Math.exp(logShear)) => {
    if (!(aux > 0) || !Number.isFinite(aux)) throw new Error('Native shear bracket left the positive finite domain.');
    const value = kernel.interval({ upstream, downstream: { ...downstream, aux }, regime });
    evaluations++;
    const residual = value.residual[0], logDerivative = value.downstream[0][0] * aux;
    if (!Number.isFinite(residual) || !Number.isFinite(logDerivative) || !(logDerivative < 0))
      throw new Error('Native shear transport row must be finite and strictly decreasing.');
    return { logShear, aux, residual };
  };
  const finish = (point, iterations) => {
    if (Math.abs(point.residual) > tolerance) throw new Error('Native shear initialization did not close its original row.');
    return { aux: point.aux, diagnostics: { method: 'native-transport-shear-log-bisection', regime,
      residual: point.residual, tolerance, iterations, evaluations, equationsChanged: false, initialGuessOnly: true } };
  };
  // Preserve an already closed seed exactly, including its floating-point bits.
  let low = sample(Math.log(downstream.aux), downstream.aux), high = low;
  if (Math.abs(low.residual) <= tolerance) return finish(low, 0);
  for (let k = 0; low.residual < 0 && k < maximumExpansions; k++) low = sample(low.logShear - 2);
  for (let k = 0; high.residual > 0 && k < maximumExpansions; k++) high = sample(high.logShear + 2);
  if (!(low.residual >= 0 && high.residual <= 0)) throw new Error('Could not bracket the original native transport shear root.');
  let best = Math.abs(low.residual) < Math.abs(high.residual) ? low : high;
  for (let k = 0; k < maximumBisections; k++) {
    const middle = (low.logShear + high.logShear) / 2;
    if (middle === low.logShear || middle === high.logShear) return finish(best, k);
    const point = sample(middle);
    if (Math.abs(point.residual) < Math.abs(best.residual)) best = point;
    if (Math.abs(point.residual) <= tolerance) return finish(point, k + 1);
    if (point.residual > 0) low = point; else high = point;
  }
  return finish(best, maximumBisections);
}
